#!/usr/bin/env python3
"""
Extract metrics from OneRoof pipeline outputs.

This CLI provides subcommands for extracting metrics from various pipeline
stages into validated JSON files. These per-sample JSON files are later
assembled into the final OneRoof report.

Usage:
    extract_metrics.py coverage --sample-id SAMPLE --bed coverage.bed --output metrics.json
    extract_metrics.py alignment --sample-id SAMPLE --bam aligned.bam --output metrics.json
    extract_metrics.py variants --sample-id SAMPLE --effects-tsv variants.tsv --output metrics.json
    extract_metrics.py consensus --sample-id SAMPLE --fasta consensus.fa --output metrics.json
    extract_metrics.py metagenomics --sample-id SAMPLE --profile-tsv profile.tsv --output metrics.json
    extract_metrics.py haplotyping --sample-id SAMPLE --devider-dir sample_devider/ --output metrics.json
    extract_metrics.py --help
"""

from pathlib import Path
from typing import Annotated

import typer
from reporting.extractors.alignment import extract as extract_alignment
from reporting.extractors.consensus import extract as extract_consensus
from reporting.extractors.coverage import extract as extract_coverage
from reporting.extractors.haplotyping import extract as extract_haplotyping
from reporting.extractors.metagenomics import extract as extract_metagenomics
from reporting.extractors.variants import extract as extract_variants
from rich.console import Console

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


@app.command("alignment")
def alignment(
    sample_id: Annotated[
        str,
        typer.Option("--sample-id", "-s", help="Sample identifier"),
    ],
    bam: Annotated[
        Path,
        typer.Option(
            "--bam",
            "-b",
            help="Path to indexed BAM file",
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
    Extract alignment metrics from a BAM file.

    Parses an indexed BAM file and extracts read counts and mapping statistics
    including total reads, mapped reads, unmapped reads, and mapping rate.
    """
    metrics = extract_alignment(sample_id, bam)
    output.write_text(metrics.model_dump_json(indent=2))
    console.print(f"[green]Wrote alignment metrics to {output}[/green]")


@app.command("variants")
def variants(
    sample_id: Annotated[
        str,
        typer.Option("--sample-id", "-s", help="Sample identifier"),
    ],
    effects_tsv: Annotated[
        Path,
        typer.Option(
            "--effects-tsv",
            "-e",
            help="Path to SnpSift extractFields TSV file",
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
    consensus_threshold: Annotated[
        float,
        typer.Option(
            "--consensus-threshold",
            "-t",
            help="Allele frequency threshold for consensus calls",
        ),
    ] = 0.5,
) -> None:
    """
    Extract variant metrics from a SnpSift effects TSV file.

    Parses the variant effects TSV and extracts counts by mutation type
    (SNP, insertion, deletion, MNP), effect type, and consensus vs subclonal
    classification based on allele frequency.
    """
    metrics = extract_variants(sample_id, effects_tsv, consensus_threshold)
    output.write_text(metrics.model_dump_json(indent=2))
    console.print(f"[green]Wrote variant metrics to {output}[/green]")


@app.command("consensus")
def consensus(
    sample_id: Annotated[
        str,
        typer.Option("--sample-id", "-s", help="Sample identifier"),
    ],
    fasta: Annotated[
        Path,
        typer.Option(
            "--fasta",
            "-f",
            help="Path to consensus FASTA file",
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
    Extract consensus metrics from a FASTA file.

    Parses the consensus FASTA and extracts sequence quality statistics
    including N content, completeness, ambiguous bases, and GC content.
    """
    metrics = extract_consensus(sample_id, fasta)
    output.write_text(metrics.model_dump_json(indent=2))
    console.print(f"[green]Wrote consensus metrics to {output}[/green]")


@app.command("metagenomics")
def metagenomics(
    sample_id: Annotated[
        str,
        typer.Option("--sample-id", "-s", help="Sample identifier"),
    ],
    profile_tsv: Annotated[
        Path,
        typer.Option(
            "--profile-tsv",
            "-p",
            help="Path to Sylph profile TSV file",
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
    top_n: Annotated[
        int,
        typer.Option(
            "--top-n",
            "-n",
            help="Number of top hits to include",
        ),
    ] = 5,
) -> None:
    """
    Extract metagenomics metrics from a Sylph profile TSV file.

    Parses the Sylph profile output and extracts top taxonomic hits,
    total abundance, and unknown fraction.
    """
    metrics = extract_metagenomics(sample_id, profile_tsv, top_n)
    output.write_text(metrics.model_dump_json(indent=2))
    console.print(f"[green]Wrote metagenomics metrics to {output}[/green]")


@app.command("haplotyping")
def haplotyping(
    sample_id: Annotated[
        str,
        typer.Option("--sample-id", "-s", help="Sample identifier"),
    ],
    devider_dir: Annotated[
        Path,
        typer.Option(
            "--devider-dir",
            "-d",
            help="Path to Devider output directory",
            exists=True,
            file_okay=False,
            dir_okay=True,
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
    Extract haplotyping metrics from a Devider output directory.

    Parses the Devider output files (snp_haplotypes.fasta, ids.txt, hap_info.txt)
    and extracts haplotype counts, read assignments, and SNP statistics.
    Nanopore only.
    """
    metrics = extract_haplotyping(sample_id, devider_dir)
    output.write_text(metrics.model_dump_json(indent=2))
    console.print(f"[green]Wrote haplotyping metrics to {output}[/green]")


if __name__ == "__main__":
    app()
