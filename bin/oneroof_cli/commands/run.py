# ruff: noqa: PLR0913, FBT002, UP045
"""
The 'run' command for the OneRoof CLI.

This is the main command that executes the Nextflow pipeline with all
configurable parameters organized into logical groups.

Note: PLR0913 (too many arguments) is disabled because CLI commands legitimately
need many parameters. FBT002 (boolean default in function) is disabled because
Typer uses boolean defaults for flag options.

Note: We intentionally do NOT use `from __future__ import annotations` here
because Typer needs to introspect the type annotations at runtime, and PEP 563
deferred evaluation breaks this.
"""

from enum import Enum
from pathlib import Path
from typing import Annotated, Optional

import typer
from oneroof_cli.app import app
from oneroof_cli.utils import generate_nextflow_command, info, run_nextflow

# =============================================================================
# Enums for Constrained Choices
# =============================================================================


class DeviderPreset(str, Enum):
    """Haplotype phasing presets for devider."""

    old_long_reads = "old-long-reads"
    nanopore_r9 = "nanopore-r9"
    nanopore_r10 = "nanopore-r10"
    hi_fi = "hi-fi"


# =============================================================================
# Available Profiles (for help text)
# =============================================================================

# Core execution profiles
CORE_PROFILES = [
    "standard",
    "docker",
    "singularity",
    "apptainer",
    "containerless",
]

# Test profiles (Illumina)
ILLUMINA_TEST_PROFILES = [
    "illumina_test_with_primers",
    "illumina_test_without_primers",
    "illumina_test_with_metagenomics",
    "illumina_test_with_phylo",
    "illumina_test_missing_fastq",
    "illumina_test_bad_primers",
    "illumina_test_decon",
]

# Test profiles (Nanopore)
NANOPORE_TEST_PROFILES = [
    "nanopore_test_with_primers",
    "nanopore_test_without_primers",
    "nanopore_test_with_haplo",
    "nanopore_test_with_metagenomics",
    "nanopore_test_with_phylo",
    "nanopore_test_missing_input",
    "nanopore_test_bad_primers",
]

ALL_PROFILES = CORE_PROFILES + ILLUMINA_TEST_PROFILES + NANOPORE_TEST_PROFILES


# =============================================================================
# Help Panel Names (for organizing --help output)
# =============================================================================

PANEL_REFERENCE = "Reference & Primers"
PANEL_INPUT = "Input Data"
PANEL_BASECALLING = "Basecalling"
PANEL_FILTERING = "Read Filtering"
PANEL_PRIMER_TRIM = "Primer Finding & Trimming"
PANEL_ALIGNMENT = "Alignment & Coverage"
PANEL_CONSENSUS = "Consensus & Variants"
PANEL_HAPLOTYPING = "Haplotyping"
PANEL_METAGENOMICS = "Metagenomics"
PANEL_DECONTAMINATION = "Decontamination"
PANEL_PHYLOGENETICS = "Phylogenetics"
PANEL_OUTPUT = "Output & Execution"


# =============================================================================
# The Run Command
# =============================================================================


@app.command("run")
def run_pipeline(
    # -------------------------------------------------------------------------
    # Reference & Primers
    # -------------------------------------------------------------------------
    refseq: Annotated[
        Optional[Path],
        typer.Option(
            "--refseq",
            help="Reference sequence for mapping (FASTA format). Required unless using a test profile.",
            exists=True,
            file_okay=True,
            dir_okay=False,
            readable=True,
            resolve_path=True,
            rich_help_panel=PANEL_REFERENCE,
        ),
    ] = None,
    ref_gbk: Annotated[
        Optional[Path],
        typer.Option(
            "--ref-gbk",
            help="Reference sequence for variant annotation (GenBank format).",
            exists=True,
            file_okay=True,
            dir_okay=False,
            readable=True,
            resolve_path=True,
            rich_help_panel=PANEL_REFERENCE,
        ),
    ] = None,
    primer_bed: Annotated[
        Optional[Path],
        typer.Option(
            "--primer-bed",
            help="BED file of primer coordinates relative to the reference.",
            exists=True,
            file_okay=True,
            dir_okay=False,
            readable=True,
            resolve_path=True,
            rich_help_panel=PANEL_REFERENCE,
        ),
    ] = None,
    primer_tsv: Annotated[
        Optional[Path],
        typer.Option(
            "--primer-tsv",
            help="TSV file of primer sequences (alternative to --primer-bed).",
            exists=True,
            file_okay=True,
            dir_okay=False,
            readable=True,
            resolve_path=True,
            rich_help_panel=PANEL_REFERENCE,
        ),
    ] = None,
    fwd_suffix: Annotated[
        Optional[str],
        typer.Option(
            "--fwd-suffix",
            help="Suffix denoting forward primers in BED file.",
            show_default="'_LEFT'",
            rich_help_panel=PANEL_REFERENCE,
        ),
    ] = None,
    rev_suffix: Annotated[
        Optional[str],
        typer.Option(
            "--rev-suffix",
            help="Suffix denoting reverse primers in BED file.",
            show_default="'_RIGHT'",
            rich_help_panel=PANEL_REFERENCE,
        ),
    ] = None,
    # -------------------------------------------------------------------------
    # Input Data
    # -------------------------------------------------------------------------
    pod5_dir: Annotated[
        Optional[Path],
        typer.Option(
            "--pod5-dir",
            help="Local directory containing POD5 files.",
            exists=True,
            file_okay=False,
            dir_okay=True,
            readable=True,
            resolve_path=True,
            rich_help_panel=PANEL_INPUT,
        ),
    ] = None,
    prepped_data: Annotated[
        Optional[Path],
        typer.Option(
            "--prepped-data",
            help="Directory with already basecalled and demultiplexed Nanopore data.",
            exists=True,
            file_okay=False,
            dir_okay=True,
            readable=True,
            resolve_path=True,
            rich_help_panel=PANEL_INPUT,
        ),
    ] = None,
    precalled_staging: Annotated[
        Optional[Path],
        typer.Option(
            "--precalled-staging",
            help="Directory to watch for Nanopore FASTQs/BAMs as they become available.",
            exists=True,
            file_okay=False,
            dir_okay=True,
            readable=True,
            resolve_path=True,
            rich_help_panel=PANEL_INPUT,
        ),
    ] = None,
    illumina_fastq_dir: Annotated[
        Optional[Path],
        typer.Option(
            "--illumina-fastq-dir",
            help="Directory containing Illumina paired-end FASTQ files.",
            exists=True,
            file_okay=False,
            dir_okay=True,
            readable=True,
            resolve_path=True,
            rich_help_panel=PANEL_INPUT,
        ),
    ] = None,
    remote_pod5_location: Annotated[
        Optional[str],
        typer.Option(
            "--remote-pod5-location",
            help="Remote SSH location to watch for POD5 files in real-time.",
            rich_help_panel=PANEL_INPUT,
        ),
    ] = None,
    file_watcher_config: Annotated[
        Optional[Path],
        typer.Option(
            "--file-watcher-config",
            help="Configuration file for remote file monitoring.",
            exists=True,
            file_okay=True,
            dir_okay=False,
            readable=True,
            resolve_path=True,
            rich_help_panel=PANEL_INPUT,
        ),
    ] = None,
    pod5_staging: Annotated[
        Optional[Path],
        typer.Option(
            "--pod5-staging",
            help="Local cache directory for POD5 files arriving from remote location.",
            file_okay=False,
            dir_okay=True,
            resolve_path=True,
            rich_help_panel=PANEL_INPUT,
        ),
    ] = None,
    # -------------------------------------------------------------------------
    # Basecalling
    # -------------------------------------------------------------------------
    model: Annotated[
        Optional[str],
        typer.Option(
            "--model",
            help="Nanopore basecalling model.",
            show_default="'sup@latest'",
            rich_help_panel=PANEL_BASECALLING,
        ),
    ] = None,
    model_cache: Annotated[
        Optional[Path],
        typer.Option(
            "--model-cache",
            help="Directory to cache basecalling models.",
            file_okay=False,
            dir_okay=True,
            resolve_path=True,
            rich_help_panel=PANEL_BASECALLING,
        ),
    ] = None,
    kit: Annotated[
        Optional[str],
        typer.Option(
            "--kit",
            help="Nanopore barcoding kit used for library preparation.",
            rich_help_panel=PANEL_BASECALLING,
        ),
    ] = None,
    pod5_batch_size: Annotated[
        Optional[int],
        typer.Option(
            "--pod5-batch-size",
            help="Number of POD5 files to basecall at once.",
            min=1,
            rich_help_panel=PANEL_BASECALLING,
        ),
    ] = None,
    basecall_max: Annotated[
        Optional[int],
        typer.Option(
            "--basecall-max",
            help="Maximum parallel basecaller instances.",
            show_default="1",
            min=1,
            rich_help_panel=PANEL_BASECALLING,
        ),
    ] = None,
    # -------------------------------------------------------------------------
    # Read Filtering
    # -------------------------------------------------------------------------
    min_len: Annotated[
        Optional[int],
        typer.Option(
            "--min-len",
            help="Minimum acceptable read length.",
            show_default="1",
            min=1,
            rich_help_panel=PANEL_FILTERING,
        ),
    ] = None,
    max_len: Annotated[
        Optional[int],
        typer.Option(
            "--max-len",
            help="Maximum acceptable read length.",
            show_default="unlimited",
            min=1,
            rich_help_panel=PANEL_FILTERING,
        ),
    ] = None,
    min_qual: Annotated[
        Optional[int],
        typer.Option(
            "--min-qual",
            help="Minimum acceptable average read quality.",
            show_default="20",
            min=0,
            rich_help_panel=PANEL_FILTERING,
        ),
    ] = None,
    # -------------------------------------------------------------------------
    # Primer Finding & Trimming
    # -------------------------------------------------------------------------
    max_mismatch: Annotated[
        Optional[int],
        typer.Option(
            "--max-mismatch",
            help="Maximum mismatches allowed when finding primers.",
            show_default="0",
            min=0,
            rich_help_panel=PANEL_PRIMER_TRIM,
        ),
    ] = None,
    forward_window: Annotated[
        Optional[int],
        typer.Option(
            "--forward-window",
            help="Search for forward primers in first N bases only (0 = entire read).",
            show_default="0",
            min=0,
            rich_help_panel=PANEL_PRIMER_TRIM,
        ),
    ] = None,
    reverse_window: Annotated[
        Optional[int],
        typer.Option(
            "--reverse-window",
            help="Search for reverse primers in last N bases only (0 = entire read).",
            show_default="0",
            min=0,
            rich_help_panel=PANEL_PRIMER_TRIM,
        ),
    ] = None,
    # -------------------------------------------------------------------------
    # Alignment & Coverage
    # -------------------------------------------------------------------------
    secondary: Annotated[
        bool,
        typer.Option(
            "--secondary/--no-secondary",
            help="Enable secondary alignments for each amplicon.",
            rich_help_panel=PANEL_ALIGNMENT,
        ),
    ] = False,
    dedup: Annotated[
        bool,
        typer.Option(
            "--dedup/--no-dedup",
            help="Enable read deduplication.",
            rich_help_panel=PANEL_ALIGNMENT,
        ),
    ] = False,
    early_downsample_to: Annotated[
        Optional[int],
        typer.Option(
            "--early-downsample-to",
            help="Early downsampling coverage target (0 = no early downsampling).",
            show_default="0",
            min=0,
            rich_help_panel=PANEL_ALIGNMENT,
        ),
    ] = None,
    downsample_to: Annotated[
        Optional[int],
        typer.Option(
            "--downsample-to",
            help="Target coverage for downsampling (0 = no downsampling).",
            show_default="0",
            min=0,
            rich_help_panel=PANEL_ALIGNMENT,
        ),
    ] = None,
    min_depth_coverage: Annotated[
        Optional[int],
        typer.Option(
            "--min-depth-coverage",
            help="Minimum depth of coverage required.",
            show_default="10",
            min=1,
            rich_help_panel=PANEL_ALIGNMENT,
        ),
    ] = None,
    # -------------------------------------------------------------------------
    # Consensus & Variants
    # -------------------------------------------------------------------------
    min_consensus_freq: Annotated[
        Optional[float],
        typer.Option(
            "--min-consensus-freq",
            help="Minimum variant frequency for consensus inclusion.",
            show_default="0.5",
            min=0.0,
            max=1.0,
            rich_help_panel=PANEL_CONSENSUS,
        ),
    ] = None,
    snpeff_cache: Annotated[
        Optional[Path],
        typer.Option(
            "--snpeff-cache",
            help="Directory to cache custom snpEff database.",
            file_okay=False,
            dir_okay=True,
            resolve_path=True,
            rich_help_panel=PANEL_CONSENSUS,
        ),
    ] = None,
    snpeff_config: Annotated[
        Optional[Path],
        typer.Option(
            "--snpeff-config",
            help="Custom snpEff configuration file.",
            exists=True,
            file_okay=True,
            dir_okay=False,
            readable=True,
            resolve_path=True,
            rich_help_panel=PANEL_CONSENSUS,
        ),
    ] = None,
    # -------------------------------------------------------------------------
    # Haplotyping
    # -------------------------------------------------------------------------
    min_haplo_reads: Annotated[
        Optional[int],
        typer.Option(
            "--min-haplo-reads",
            help="Minimum read support to report an amplicon-haplotype.",
            show_default="2",
            min=1,
            rich_help_panel=PANEL_HAPLOTYPING,
        ),
    ] = None,
    devider_preset: Annotated[
        Optional[DeviderPreset],
        typer.Option(
            "--devider-preset",
            help="Haplotype phasing preset for devider.",
            show_default="'nanopore-r10'",
            case_sensitive=False,
            rich_help_panel=PANEL_HAPLOTYPING,
        ),
    ] = None,
    # -------------------------------------------------------------------------
    # Metagenomics
    # -------------------------------------------------------------------------
    meta_ref: Annotated[
        Optional[Path],
        typer.Option(
            "--meta-ref",
            help="Reference dataset (local FASTA or pre-built Sylph database) for profiling.",
            exists=True,
            readable=True,
            resolve_path=True,
            rich_help_panel=PANEL_METAGENOMICS,
        ),
    ] = None,
    sylph_tax_db: Annotated[
        Optional[Path],
        typer.Option(
            "--sylph-tax-db",
            help="Taxonomic annotation database for the Sylph dataset.",
            exists=True,
            readable=True,
            resolve_path=True,
            rich_help_panel=PANEL_METAGENOMICS,
        ),
    ] = None,
    sylph_db_link: Annotated[
        Optional[str],
        typer.Option(
            "--sylph-db-link",
            help="URL to download the Sylph database.",
            rich_help_panel=PANEL_METAGENOMICS,
        ),
    ] = None,
    k: Annotated[
        Optional[int],
        typer.Option(
            "--k",
            "-k",
            help="K-mer size for metagenomic sketching.",
            show_default="31",
            min=1,
            rich_help_panel=PANEL_METAGENOMICS,
        ),
    ] = None,
    # -------------------------------------------------------------------------
    # Decontamination
    # -------------------------------------------------------------------------
    contam_fasta: Annotated[
        Optional[Path],
        typer.Option(
            "--contam-fasta",
            help="Contamination FASTA dataset to scrub from reads.",
            exists=True,
            file_okay=True,
            dir_okay=False,
            readable=True,
            resolve_path=True,
            rich_help_panel=PANEL_DECONTAMINATION,
        ),
    ] = None,
    contam_link: Annotated[
        Optional[str],
        typer.Option(
            "--contam-link",
            help="URL to download the contamination FASTA dataset.",
            rich_help_panel=PANEL_DECONTAMINATION,
        ),
    ] = None,
    # -------------------------------------------------------------------------
    # Phylogenetics
    # -------------------------------------------------------------------------
    nextclade_dataset: Annotated[
        Optional[str],
        typer.Option(
            "--nextclade-dataset",
            help="Nextclade dataset name or path.",
            rich_help_panel=PANEL_PHYLOGENETICS,
        ),
    ] = None,
    nextclade_cache: Annotated[
        Optional[Path],
        typer.Option(
            "--nextclade-cache",
            help="Directory to cache Nextclade datasets.",
            file_okay=False,
            dir_okay=True,
            resolve_path=True,
            rich_help_panel=PANEL_PHYLOGENETICS,
        ),
    ] = None,
    # -------------------------------------------------------------------------
    # Output & Execution
    # -------------------------------------------------------------------------
    results: Annotated[
        Optional[Path],
        typer.Option(
            "--results",
            help="Directory to place results.",
            show_default="'./results'",
            file_okay=False,
            dir_okay=True,
            resolve_path=True,
            rich_help_panel=PANEL_OUTPUT,
        ),
    ] = None,
    log: Annotated[
        bool,
        typer.Option(
            "--log/--no-log",
            help="Use log scale for multisample coverage plots.",
            rich_help_panel=PANEL_OUTPUT,
        ),
    ] = False,
    email: Annotated[
        Optional[str],
        typer.Option(
            "--email",
            help="Email address(es) for completion notification (comma-delimited).",
            rich_help_panel=PANEL_OUTPUT,
        ),
    ] = None,
    low_memory: Annotated[
        bool,
        typer.Option(
            "--low-memory/--no-low-memory",
            help="Run in low-memory mode, limiting parallel high-memory processes.",
            rich_help_panel=PANEL_OUTPUT,
        ),
    ] = False,
    cleanup: Annotated[
        bool,
        typer.Option(
            "--cleanup/--no-cleanup",
            help="Clean up work directory after successful run.",
            rich_help_panel=PANEL_OUTPUT,
        ),
    ] = False,
    resume: Annotated[
        bool,
        typer.Option(
            "--resume/--no-resume",
            help="Resume from a previous run (uses Nextflow's -resume flag).",
            rich_help_panel=PANEL_OUTPUT,
        ),
    ] = False,
    profile: Annotated[
        Optional[str],
        typer.Option(
            "--profile",
            "-p",
            help=(
                "Nextflow execution profile(s), comma-separated. "
                "Core: standard, docker, singularity, apptainer, containerless. "
                "Test profiles also available (e.g., illumina_test_with_primers)."
            ),
            rich_help_panel=PANEL_OUTPUT,
        ),
    ] = None,
    # -------------------------------------------------------------------------
    # Dry Run (CLI-only, not passed to Nextflow)
    # -------------------------------------------------------------------------
    dry_run: Annotated[
        bool,
        typer.Option(
            "--dry-run",
            "-n",
            help="Show the Nextflow command without executing it.",
            rich_help_panel=PANEL_OUTPUT,
        ),
    ] = False,
) -> None:
    """
    [bold green]Run[/bold green] the OneRoof pipeline.

    Execute the complete amplicon sequencing analysis workflow, from raw reads
    through consensus calling and variant annotation.

    [dim]Examples:[/dim]

        # Basic Illumina run
        oneroof run --refseq ref.fasta --illumina-fastq-dir ./reads/

        # Nanopore with primers
        oneroof run --refseq ref.fasta --prepped-data ./basecalled/ --primer-bed primers.bed

        # With Docker profile and resume
        oneroof run --refseq ref.fasta --illumina-fastq-dir ./reads/ --profile docker --resume
    """
    # Build the arguments dictionary for Nextflow
    # We convert Python naming (snake_case with hyphens in CLI) to Nextflow naming (snake_case)
    args: dict[str, object] = {
        # Reference & Primers
        "refseq": refseq,
        "ref_gbk": ref_gbk,
        "primer_bed": primer_bed,
        "primer_tsv": primer_tsv,
        "fwd_suffix": fwd_suffix,
        "rev_suffix": rev_suffix,
        # Input Data
        "pod5_dir": pod5_dir,
        "prepped_data": prepped_data,
        "precalled_staging": precalled_staging,
        "illumina_fastq_dir": illumina_fastq_dir,
        "remote_pod5_location": remote_pod5_location,
        "file_watcher_config": file_watcher_config,
        "pod5_staging": pod5_staging,
        # Basecalling
        "model": model,
        "model_cache": model_cache,
        "kit": kit,
        "pod5_batch_size": pod5_batch_size,
        "basecall_max": basecall_max,
        # Read Filtering
        "min_len": min_len,
        "max_len": max_len,
        "min_qual": min_qual,
        # Primer Finding & Trimming
        "max_mismatch": max_mismatch,
        "forward_window": forward_window,
        "reverse_window": reverse_window,
        # Alignment & Coverage
        "secondary": secondary if secondary else None,
        "dedup": dedup if dedup else None,
        "early_downsample_to": early_downsample_to,
        "downsample_to": downsample_to,
        "min_depth_coverage": min_depth_coverage,
        # Consensus & Variants
        "min_consensus_freq": min_consensus_freq,
        "snpeff_cache": snpeff_cache,
        "snpEff_config": snpeff_config,  # Note: Nextflow uses snpEff_config (camelCase)
        # Haplotyping
        "min_haplo_reads": min_haplo_reads,
        "devider_preset": devider_preset.value if devider_preset else None,
        # Metagenomics
        "meta_ref": meta_ref,
        "sylph_tax_db": sylph_tax_db,
        "sylph_db_link": sylph_db_link,
        "k": k,
        # Decontamination
        "contam_fasta": contam_fasta,
        "contam_link": contam_link,
        # Phylogenetics
        "nextclade_dataset": nextclade_dataset,
        "nextclade_cache": nextclade_cache,
        # Output & Execution
        "results": results,
        "log": log if log else None,
        "email": email,
        "low_memory": low_memory if low_memory else None,
        "cleanup": cleanup if cleanup else None,
        # Profile is handled specially in generate_nextflow_command
        # It's passed as a comma-separated string directly to Nextflow
        "profile": profile,
    }

    # Handle resume flag - this goes to Nextflow as -resume, not --resume
    command = generate_nextflow_command(args)
    if resume:
        command = f"{command} -resume"

    if dry_run:
        info("Dry run mode - command not executed")
        from oneroof_cli.utils import console  # noqa: PLC0415

        console.print(f"\n[bold]Command:[/bold]\n{command}\n")
        return

    run_nextflow(command)
