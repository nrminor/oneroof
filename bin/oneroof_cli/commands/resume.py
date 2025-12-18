"""
The 'resume' command for the OneRoof CLI.

Resumes a previous pipeline run using the cached .nfresume file.
"""

from oneroof_cli.app import app
from oneroof_cli.utils import resume_nextflow


@app.command("resume")
def resume_pipeline() -> None:
    """
    [bold blue]Resume[/bold blue] a previous pipeline run.

    This command reads the cached Nextflow command from the .nfresume file
    (created during a previous 'oneroof run') and re-executes it with the
    -resume flag to continue from where it left off.

    [dim]Note:[/dim] This is different from 'oneroof run --resume', which
    requires you to re-specify all parameters. The 'resume' subcommand
    uses the exact same parameters from your previous run.

    [dim]Example:[/dim]

        # First, run the pipeline
        oneroof run --refseq ref.fasta --illumina-fastq-dir ./reads/

        # If it fails or is interrupted, resume with:
        oneroof resume
    """
    resume_nextflow()
