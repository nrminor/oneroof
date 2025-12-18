#!/usr/bin/env python3
"""
Extract metrics from OneRoof pipeline outputs.

This CLI provides subcommands for extracting metrics from various pipeline
stages into validated JSON files. These per-sample JSON files are later
assembled into the final OneRoof report.

Usage:
    extract_metrics.py coverage --sample-id SAMPLE --bed coverage.bed --output metrics.json
    extract_metrics.py --help

Future subcommands (Phase 2):
    extract_metrics.py alignment --sample-id SAMPLE --bam aligned.bam --output metrics.json
    extract_metrics.py variants --sample-id SAMPLE --effects-tsv variants.tsv --output metrics.json
    extract_metrics.py consensus --sample-id SAMPLE --fasta consensus.fa --output metrics.json
"""

from pathlib import Path
from typing import Annotated

import typer
from rich.console import Console

from reporting.extractors.coverage import extract as extract_coverage

app = typer.Typer(
    name="extract_metrics",
    help="Extract metrics from OneRoof pipeline outputs into validated JSON.",
    add_completion=False,
    rich_markup_mode="rich",
    no_args_is_help=True,
    context_settings={"help_option_names": ["-h", "--help"]},
)
console = Console()


@app.callback()
def main() -> None:
    """
    Extract metrics from OneRoof pipeline outputs into validated JSON.

    Each subcommand extracts metrics from a specific pipeline stage.
    """
    pass


@app.command("coverage")
def coverage(
    sample_id: Annotated[
        str,
        typer.Option("--sample-id", "-s", help="Sample identifier"),
    ],
    bed: Annotated[
        Path,
        typer.Option(
            "--bed",
            "-b",
            help="Path to bedtools genomecov BED file",
            exists=True,
            file_okay=True,
            dir_okay=False,
            readable=True,
            resolve_path=True,
        ),
    ],
    output: Annotated[
        Path,
        typer.Option(
            "--output",
            "-o",
            help="Output JSON file path",
            file_okay=True,
            dir_okay=False,
            resolve_path=True,
        ),
    ],
) -> None:
    """
    Extract coverage metrics from bedtools genomecov output.

    Parses a per-base BED file (from `bedtools genomecov -bga`) and computes
    coverage statistics including mean/median depth and genome coverage at
    various thresholds.
    """
    metrics = extract_coverage(sample_id, bed)
    output.write_text(metrics.model_dump_json(indent=2))
    console.print(f"[green]Wrote coverage metrics to {output}[/green]")


if __name__ == "__main__":
    app()
