"""
OneRoof CLI - A Typer-based command-line interface for the OneRoof pipeline.

This package provides a user-friendly wrapper around the Nextflow-based
OneRoof pipeline for amplicon sequencing analysis.

Usage:
    oneroof run --refseq ref.fasta --illumina-fastq-dir ./reads/
    oneroof resume
    oneroof --help
"""

import sys

import typer
from rich.console import Console

from oneroof_cli.app import app

# Import commands to register them with the app
from oneroof_cli.commands import env, resume, run, validate  # noqa: F401

__all__ = ["app", "main"]

console = Console()


def main() -> None:
    """Main entry point for the OneRoof CLI."""
    try:
        app()
    except KeyboardInterrupt:
        console.print("\n[yellow]Interrupted by user[/yellow]")
        sys.exit(130)
    except typer.Exit:
        # Normal exit from Typer - re-raise to preserve exit code
        raise
    except typer.Abort:
        # User aborted (e.g., Ctrl+C during prompt)
        sys.exit(1)
