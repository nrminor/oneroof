#!/usr/bin/env python3

import argparse
import shlex
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).parent.parent

# ANSI escape codes for terminal formatting
BOLD = "\033[1m"
UNDERLINE = "\033[4m"
RESET = "\033[0m"


class BoldGroupHelpFormatter(argparse.RawDescriptionHelpFormatter):
    """Custom formatter that renders argument group headings in bold and underlined."""

    def start_section(self, heading: str | None) -> None:
        """Override to wrap section headings in bold + underline ANSI codes."""
        if heading:
            heading = f"{BOLD}{UNDERLINE}{heading}{RESET}"
        super().start_section(heading)


def parse_command_line_args() -> argparse.Namespace:  # noqa: PLR0915
    """Parse command line arguments for the OneRoof pipeline.

    This function sets up an argument parser with various subcommands and their
    respective arguments. The subcommands include 'env', 'validate', 'resume',
    and 'run', each with its own set of options.

    Returns:
        argparse.Namespace: An object containing all the parsed arguments.

    Subcommands:
        env: Check that all dependencies are available in the environment.
        validate: Validate provided inputs.
        resume: Resume the previous run.
        run: Run the full pipeline.

    The 'run' subcommand has numerous options for configuring the pipeline,
    including input files, reference sequences, basecalling parameters,
    and various thresholds for filtering and analysis.
    """
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(title="Subcommands", dest="subcommands")

    # create the env subcommand
    _env = subparsers.add_parser(
        "env",
        help="Check that all dependencies are available in the environment",
    )

    # create the validate subcommand
    _validate = subparsers.add_parser("validate", help="Validate provided inputs.")

    # create the resume subcommand
    resume = subparsers.add_parser("resume", help="Resume the previous run.")
    resume.add_argument(
        "-v",
        "--verbose",
        action="count",
        default=0,
        help="Increase verbosity level (-v for WARNING, -vv for INFO, -vvv for DEBUG)",
        required=False,
    )

    # create the run subcommand
    run = subparsers.add_parser(
        "run",
        help="Run the full pipeline.",
        formatter_class=BoldGroupHelpFormatter,
    )

    # -------------------------------------------------------------------------
    # Reference & Primers
    # -------------------------------------------------------------------------
    ref_group = run.add_argument_group(
        "Reference & Primers",
        "Reference sequences and primer definitions",
    )
    ref_group.add_argument(
        "--refseq",
        type=str,
        required=True,
        help="The reference sequence to be used for mapping in FASTA format.",
    )
    ref_group.add_argument(
        "--ref_gbk",
        type=str,
        default=None,
        help="The reference sequence to be used for variant annotation in Genbank format.",
    )
    ref_group.add_argument(
        "--primer_bed",
        type=str,
        default=None,
        help="A BED file of primer coordinates relative to the reference.",
    )
    ref_group.add_argument(
        "--primer_tsv",
        type=str,
        default=None,
        help="A TSV file of primer sequences (alternative to --primer_bed).",
    )
    ref_group.add_argument(
        "--fwd_suffix",
        type=str,
        default=None,
        help="Suffix in the primer BED file denoting forward primers [default: '_LEFT'].",
    )
    ref_group.add_argument(
        "--rev_suffix",
        type=str,
        default=None,
        help="Suffix in the primer BED file denoting reverse primers [default: '_RIGHT'].",
    )

    # -------------------------------------------------------------------------
    # Input Data
    # -------------------------------------------------------------------------
    input_group = run.add_argument_group(
        "Input Data",
        "Specify input data sources (Nanopore or Illumina)",
    )
    input_group.add_argument(
        "--pod5_dir",
        type=str,
        default=None,
        help="Local directory containing POD5 files.",
    )
    input_group.add_argument(
        "--prepped_data",
        type=str,
        default=None,
        help="Directory with already basecalled and demultiplexed Nanopore data.",
    )
    input_group.add_argument(
        "--precalled_staging",
        type=str,
        default=None,
        help="Directory to watch for Nanopore FASTQs/BAMs as they become available.",
    )
    input_group.add_argument(
        "--illumina_fastq_dir",
        type=str,
        default=None,
        help="Directory containing Illumina paired-end FASTQ files.",
    )
    input_group.add_argument(
        "--remote_pod5_location",
        type=str,
        default=None,
        help="Remote SSH location to watch for POD5 files in real-time.",
    )
    input_group.add_argument(
        "--file_watcher_config",
        type=str,
        default=None,
        help="Configuration file for remote file monitoring.",
    )
    input_group.add_argument(
        "--pod5_staging",
        type=str,
        default=None,
        help="Local cache directory for POD5 files arriving from remote location.",
    )

    # -------------------------------------------------------------------------
    # Nanopore Basecalling
    # -------------------------------------------------------------------------
    basecall_group = run.add_argument_group(
        "Basecalling",
        "Nanopore basecalling options (only applies to POD5 input)",
    )
    basecall_group.add_argument(
        "--model",
        type=str,
        default=None,
        help="Nanopore basecalling model [default: 'sup@latest'].",
    )
    basecall_group.add_argument(
        "--model_cache",
        type=str,
        default=None,
        help="Directory to cache basecalling models.",
    )
    basecall_group.add_argument(
        "--kit",
        type=str,
        default=None,
        help="Nanopore barcoding kit used for library preparation.",
    )
    basecall_group.add_argument(
        "--pod5_batch_size",
        type=int,
        default=None,
        help="Number of POD5 files to basecall at once.",
    )
    basecall_group.add_argument(
        "--basecall_max",
        type=int,
        default=None,
        help="Maximum parallel basecaller instances [default: 1].",
    )

    # -------------------------------------------------------------------------
    # Read Filtering
    # -------------------------------------------------------------------------
    filter_group = run.add_argument_group(
        "Read Filtering",
        "Quality and length filters for input reads",
    )
    filter_group.add_argument(
        "--min_len",
        type=int,
        default=None,
        help="Minimum acceptable read length [default: 1].",
    )
    filter_group.add_argument(
        "--max_len",
        type=int,
        default=None,
        help="Maximum acceptable read length [default: unlimited].",
    )
    filter_group.add_argument(
        "--min_qual",
        type=int,
        default=None,
        help="Minimum acceptable average read quality [default: 20].",
    )

    # -------------------------------------------------------------------------
    # Primer Finding & Trimming
    # -------------------------------------------------------------------------
    primer_group = run.add_argument_group(
        "Primer Finding & Trimming",
        "Options for primer detection and removal",
    )
    primer_group.add_argument(
        "--max_mismatch",
        type=int,
        default=None,
        help="Maximum mismatches allowed when finding primers [default: 0].",
    )
    primer_group.add_argument(
        "--forward_window",
        type=int,
        default=None,
        help="Search for forward primers in first N bases only (0 = entire read) [default: 0].",
    )
    primer_group.add_argument(
        "--reverse_window",
        type=int,
        default=None,
        help="Search for reverse primers in last N bases only (0 = entire read) [default: 0].",
    )

    # -------------------------------------------------------------------------
    # Alignment & Coverage
    # -------------------------------------------------------------------------
    align_group = run.add_argument_group(
        "Alignment & Coverage",
        "Alignment behavior and coverage options",
    )
    align_group.add_argument(
        "--secondary",
        action="store_true",
        help="Enable secondary alignments for each amplicon.",
    )
    align_group.add_argument(
        "--dedup",
        action="store_true",
        help="Enable read deduplication.",
    )
    align_group.add_argument(
        "--early_downsample_to",
        type=int,
        default=None,
        help="Early downsampling coverage target (0 = no early downsampling) [default: 0].",
    )
    align_group.add_argument(
        "--downsample_to",
        type=int,
        default=None,
        help="Target coverage for downsampling (0 = no downsampling) [default: 0].",
    )
    align_group.add_argument(
        "--min_depth_coverage",
        type=int,
        default=None,
        help="Minimum depth of coverage required [default: 10].",
    )

    # -------------------------------------------------------------------------
    # Consensus & Variants
    # -------------------------------------------------------------------------
    consensus_group = run.add_argument_group(
        "Consensus & Variants",
        "Consensus calling and variant annotation options",
    )
    consensus_group.add_argument(
        "--min_consensus_freq",
        type=float,
        default=None,
        help="Minimum variant frequency for consensus inclusion [default: 0.5].",
    )
    consensus_group.add_argument(
        "--snpeff_cache",
        type=str,
        default=None,
        help="Directory to cache custom snpEff database.",
    )
    consensus_group.add_argument(
        "--snpEff_config",
        type=str,
        default=None,
        help="Custom snpEff configuration file.",
    )

    # -------------------------------------------------------------------------
    # Haplotyping
    # -------------------------------------------------------------------------
    haplo_group = run.add_argument_group(
        "Haplotyping",
        "Amplicon haplotype phasing options",
    )
    haplo_group.add_argument(
        "--min_haplo_reads",
        type=int,
        default=None,
        help="Minimum read support to report an amplicon-haplotype [default: 2].",
    )
    haplo_group.add_argument(
        "--devider_preset",
        type=str,
        default=None,
        choices=["old-long-reads", "nanopore-r9", "nanopore-r10", "hi-fi"],
        help="Haplotype phasing preset for devider [default: 'nanopore-r10'].",
    )

    # -------------------------------------------------------------------------
    # Metagenomics
    # -------------------------------------------------------------------------
    meta_group = run.add_argument_group(
        "Metagenomics",
        "Metagenomic profiling options",
    )
    meta_group.add_argument(
        "--meta_ref",
        type=str,
        default=None,
        help="Reference dataset (local FASTA or pre-built Sylph database) for profiling.",
    )
    meta_group.add_argument(
        "--sylph_tax_db",
        type=str,
        default=None,
        help="Taxonomic annotation database for the Sylph dataset.",
    )
    meta_group.add_argument(
        "--sylph_db_link",
        type=str,
        default=None,
        help="URL to download the Sylph database.",
    )
    meta_group.add_argument(
        "--k",
        type=int,
        default=None,
        help="K-mer size for metagenomic sketching [default: 31].",
    )

    # -------------------------------------------------------------------------
    # Decontamination
    # -------------------------------------------------------------------------
    decon_group = run.add_argument_group(
        "Decontamination",
        "Host/contaminant removal options",
    )
    decon_group.add_argument(
        "--contam_fasta",
        type=str,
        default=None,
        help="Contamination FASTA dataset to scrub from reads.",
    )
    decon_group.add_argument(
        "--contam_link",
        type=str,
        default=None,
        help="URL to download the contamination FASTA dataset.",
    )

    # -------------------------------------------------------------------------
    # Phylogenetics
    # -------------------------------------------------------------------------
    phylo_group = run.add_argument_group(
        "Phylogenetics",
        "Nextclade phylogenetic placement options",
    )
    phylo_group.add_argument(
        "--nextclade_dataset",
        type=str,
        default=None,
        help="Nextclade dataset name or path.",
    )
    phylo_group.add_argument(
        "--nextclade_cache",
        type=str,
        default=None,
        help="Directory to cache Nextclade datasets.",
    )

    # -------------------------------------------------------------------------
    # Output & Execution
    # -------------------------------------------------------------------------
    output_group = run.add_argument_group(
        "Output & Execution",
        "Output location and execution options",
    )
    output_group.add_argument(
        "--results",
        type=str,
        default=None,
        help="Directory to place results [default: './results'].",
    )
    output_group.add_argument(
        "--log",
        action="store_true",
        help="Use log scale for multisample coverage plots.",
    )
    output_group.add_argument(
        "--email",
        type=str,
        default=None,
        help="Email address(es) for completion notification (comma-delimited).",
    )
    output_group.add_argument(
        "--low_memory",
        action="store_true",
        help="Run in low-memory mode, limiting parallel high-memory processes.",
    )
    output_group.add_argument(
        "--cleanup",
        action="store_true",
        help="Clean up work directory after successful run.",
    )
    output_group.add_argument(
        "--resume",
        action="store_true",
        help="Resume from a previous run.",
    )
    output_group.add_argument(
        "-profile",
        type=str,
        nargs="+",
        choices=["standard", "docker", "singularity", "apptainer", "containerless"],
        default=None,
        help="Nextflow execution profile(s) to use.",
    )

    return parser.parse_args()


def generate_nextflow_command(args: argparse.Namespace) -> str:
    """Generate a Nextflow command based on the provided arguments.

    This function takes the parsed command-line arguments and converts them into
    a Nextflow command string. It filters out None values and handles boolean flags
    appropriately.

    Args:
        args (argparse.Namespace): The parsed command-line arguments.

    Returns:
        str: A string representing the complete Nextflow command to be executed.

    The function performs the following steps:
    1. Converts the argparse.Namespace to a dictionary.
    2. Filters out None values and creates Nextflow parameter strings.
    3. Handles boolean flags by including only the parameter name if True.
    4. Joins all parameters into a single string.
    5. Constructs and returns the full Nextflow command.
    """
    # Convert the namespace to a dictionary
    args_dict = {k: v for k, v in vars(args).items() if k != "subcommands"}

    # Filter out None values and create Nextflow parameter strings
    nextflow_params = []
    for arg, value in args_dict.items():
        if value is not None:
            if isinstance(value, bool):
                # For boolean flags, just include the parameter name if True
                if value:
                    nextflow_params.append(f"--{arg}")
            else:
                # For other types, include the value
                nextflow_params.append(f"--{arg} {shlex.quote(str(value))}")

    # Join all parameters into a single string
    params_str = " ".join(nextflow_params)

    # Construct and return the full Nextflow command
    return f"nextflow run . {params_str}"


def run_nextflow(run_command: str) -> None:
    """Run the Nextflow command.

    This function executes the provided Nextflow command, handles resuming a previous run,
    and writes the command to a file for future resumption.

    Args:
        run_command (str): The Nextflow command to be executed.

    Returns:
        None

    The function performs the following steps:
    1. Writes the run command to a '.nfresume' file, appending '-resume' if not already present.
    2. If the command ends with '-resume', it only writes to the file and returns.
    3. Splits the command into a list and executes it using subprocess.run().
    """

    with Path(".nfresume").open("w", encoding="utf8") as resume_handle:
        if run_command.endswith("-resume"):
            resume_handle.write(run_command)
            return
        resume_handle.write(f"{run_command} -resume")

    split_command = run_command.split("")
    subprocess.run(split_command, check=True)  # noqa: S603
    return


def resume_nextflow() -> None:
    """Resume a previous Nextflow run.

    This function resumes a previously executed Nextflow run by reading the command
    from a '.nfresume' file and executing it. It checks for the existence of the
    '.nfresume' file before proceeding.

    Raises:
        AssertionError: If the '.nfresume' file does not exist, indicating that no
        previous run was detected.

    Returns:
        None

    The function performs the following steps:
    1. Checks for the existence of the '.nfresume' file.
    2. Reads the Nextflow command from the '.nfresume' file.
    3. Calls the run_nextflow() function with the read command to resume the run.
    """

    assert Path(
        ".nfresume",
    ).exists(), "Previous run not detected. Make sure you start with `oneroof run`"
    "before switching to `oneroof resume`."

    with Path(".nfresume").open("r", encoding="utf8") as resume_handle:
        run_command = resume_handle.readline().strip()
    run_nextflow(run_command)


def main() -> None:
    """Main function serving as the entry point for the OneRoof pipeline.

    This function orchestrates the execution of the pipeline based on the chosen subcommand.
    It parses command-line arguments, validates the subcommand, and invokes the appropriate
    function for running or resuming the Nextflow pipeline.

    The function supports the following subcommands:
    - 'run': Executes a new pipeline run with the provided arguments.
    - 'resume': Resumes a previously executed pipeline run.
    - Other subcommands (e.g., 'env', 'validate') are recognized but not yet implemented.

    Raises:
        AssertionError: If an invalid number of subcommands is detected.
        SystemExit: If an unsupported subcommand is provided.

    Returns:
        None

    The function performs the following steps:
    1. Parses command-line arguments.
    2. Validates that exactly one subcommand is chosen.
    3. Executes the appropriate action based on the chosen subcommand:
       - For 'run', generates and executes a Nextflow command.
       - For 'resume', resumes a previous Nextflow run.
       - For other subcommands, exits with an informative message.
    """

    args = parse_command_line_args()
    chosen_subcommand = [v for k, v in vars(args).items() if k == "subcommands"]
    assert len(chosen_subcommand) == 1, (
        "Invalid state registered in the form of requesting more than one or zero sub"
    )
    f"commands: {chosen_subcommand}"

    subcommand = chosen_subcommand[0]
    match subcommand:
        case "run":
            run_command = generate_nextflow_command(args)
            run_nextflow(run_command)
        case "resume":
            resume_nextflow()
        case _:
            sys.exit(
                f"The subcommand {subcommand} is not yet supported but will be soon!",
            )


if __name__ == "__main__":
    main()
