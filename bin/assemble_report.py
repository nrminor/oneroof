#!/usr/bin/env python3
"""
Assemble OneRoof report from collected metrics.

This CLI aggregates per-sample metrics JSON files into a unified report,
generates MultiQC custom content files, and produces visualizations.

Usage:
    assemble_report.py assemble \\
        --metrics-dir /path/to/metrics \\
        --output-dir /path/to/output \\
        --platform ont \\
        --reference-name "SARS-CoV-2"

    assemble_report.py --help
"""

import json
from datetime import datetime
from pathlib import Path
from typing import Annotated

import typer
from rich.console import Console

from reporting import (
    generate_coverage_table_tsv,
    generate_general_stats_tsv,
    generate_multiqc_config,
)
from reporting.schema import (
    AlignmentMetrics,
    OneRoofReport,
    Platform,
    PrimerInfo,
    QCStatus,
    ReferenceInfo,
    RunMetadata,
    RunParameters,
    SampleMetrics,
    Summary,
)
from reporting.visualizations import coverage_bar_chart, coverage_summary_heatmap

app = typer.Typer(
    name="assemble_report",
    help="Assemble OneRoof report from collected metrics.",
    add_completion=False,
    rich_markup_mode="rich",
    no_args_is_help=True,
    context_settings={"help_option_names": ["-h", "--help"]},
)
console = Console()


# Default QC thresholds
DEFAULT_QC_THRESHOLDS = {
    "coverage_pass": 0.95,  # ≥95% genome at 10x for pass
    "coverage_warn": 0.80,  # ≥80% genome at 10x for warn
    "completeness_pass": 0.98,  # ≥98% non-N bases for pass
    "completeness_warn": 0.90,  # ≥90% non-N bases for warn
}


@app.callback()
def main() -> None:
    """
    Assemble OneRoof report from collected metrics.

    Aggregates per-sample metrics, determines QC status, generates JSON reports,
    MultiQC files, and visualizations.
    """
    pass


def load_sample_metrics(metrics_dir: Path) -> dict[str, dict]:
    """
    Load all per-sample metrics JSON files from a directory.

    Globs for `*_metrics.json` files and merges metrics by sample_id.
    Different metric types (coverage, alignment, variants, etc.) are
    stored under their respective keys.

    Args:
        metrics_dir: Directory containing `*_metrics.json` files

    Returns:
        Dict mapping sample_id to merged metrics dict
    """
    samples: dict[str, dict] = {}

    for metrics_file in metrics_dir.glob("*_metrics.json"):
        with open(metrics_file) as f:
            data = json.load(f)

        sample_id = data.get("sample_id")
        if not sample_id:
            console.print(
                f"[yellow]Warning: No sample_id in {metrics_file}, skipping[/yellow]"
            )
            continue

        if sample_id not in samples:
            samples[sample_id] = {"sample_id": sample_id}

        # Determine metric type from filename
        # e.g., "sample1_coverage_metrics.json" -> "coverage"
        filename = metrics_file.stem  # e.g., "sample1_coverage_metrics"
        if "_coverage_metrics" in filename:
            samples[sample_id]["coverage"] = data
        elif "_alignment_metrics" in filename:
            samples[sample_id]["alignment"] = data
        elif "_variant_metrics" in filename:
            samples[sample_id]["variants"] = data
        elif "_consensus_metrics" in filename:
            samples[sample_id]["consensus"] = data
        elif "_metagenomics_metrics" in filename:
            samples[sample_id]["metagenomics"] = data
        elif "_haplotyping_metrics" in filename:
            samples[sample_id]["haplotyping"] = data
        else:
            # Unknown metric type - store under generic key
            samples[sample_id]["other"] = data

    return samples


def determine_qc_status(sample: dict, thresholds: dict) -> tuple[QCStatus, list[str]]:
    """
    Determine QC status and notes for a sample based on thresholds.

    Args:
        sample: Sample metrics dict with coverage, consensus, etc.
        thresholds: Dict with coverage_pass, coverage_warn, completeness_pass,
                   completeness_warn thresholds

    Returns:
        Tuple of (QCStatus, list of human-readable notes)
    """
    notes: list[str] = []
    status = QCStatus.PASS

    # Check coverage threshold (genome coverage at 10x)
    coverage = sample.get("coverage", {})
    genome_cov_10x = coverage.get("genome_coverage_at_10x", 0)

    coverage_pass = thresholds.get(
        "coverage_pass", DEFAULT_QC_THRESHOLDS["coverage_pass"]
    )
    coverage_warn = thresholds.get(
        "coverage_warn", DEFAULT_QC_THRESHOLDS["coverage_warn"]
    )

    if genome_cov_10x < coverage_warn:
        status = QCStatus.FAIL
        notes.append(
            f"Low coverage: {genome_cov_10x:.1%} genome at ≥10x (threshold: {coverage_warn:.0%})"
        )
    elif genome_cov_10x < coverage_pass:
        if status != QCStatus.FAIL:
            status = QCStatus.WARN
        notes.append(
            f"Marginal coverage: {genome_cov_10x:.1%} genome at ≥10x (threshold: {coverage_pass:.0%})"
        )

    # Check completeness threshold (if consensus metrics available)
    consensus = sample.get("consensus", {})
    if consensus:
        completeness = consensus.get("completeness", 0)
        completeness_pass = thresholds.get(
            "completeness_pass", DEFAULT_QC_THRESHOLDS["completeness_pass"]
        )
        completeness_warn = thresholds.get(
            "completeness_warn", DEFAULT_QC_THRESHOLDS["completeness_warn"]
        )

        if completeness < completeness_warn:
            status = QCStatus.FAIL
            notes.append(
                f"Low completeness: {completeness:.1%} (threshold: {completeness_warn:.0%})"
            )
        elif completeness < completeness_pass:
            if status != QCStatus.FAIL:
                status = QCStatus.WARN
            notes.append(
                f"Marginal completeness: {completeness:.1%} (threshold: {completeness_pass:.0%})"
            )

    return status, notes


def compute_summary(samples: dict[str, dict]) -> Summary:
    """
    Compute aggregate summary statistics from sample metrics.

    Args:
        samples: Dict mapping sample_id to metrics dict

    Returns:
        Summary model with aggregate statistics
    """
    sample_count = len(samples)
    samples_pass = 0
    samples_warn = 0
    samples_fail = 0

    coverage_depths: list[float] = []
    genome_coverages: list[float] = []
    total_variants = 0

    for sample_data in samples.values():
        # Count by QC status
        qc_status = sample_data.get("qc_status")
        if qc_status == QCStatus.PASS or qc_status == "pass":
            samples_pass += 1
        elif qc_status == QCStatus.WARN or qc_status == "warn":
            samples_warn += 1
        else:
            samples_fail += 1

        # Collect coverage metrics
        coverage = sample_data.get("coverage", {})
        if coverage:
            mean_cov = coverage.get("mean_coverage", 0)
            genome_cov = coverage.get("genome_coverage_at_10x", 0)
            coverage_depths.append(mean_cov)
            genome_coverages.append(genome_cov)

        # Collect variant metrics (if available)
        variants = sample_data.get("variants", {})
        if variants:
            total_variants += variants.get("total_called", 0)

    # Calculate means
    mean_coverage_depth = (
        sum(coverage_depths) / len(coverage_depths) if coverage_depths else 0.0
    )
    mean_genome_coverage = (
        sum(genome_coverages) / len(genome_coverages) if genome_coverages else 0.0
    )

    return Summary(
        sample_count=sample_count,
        samples_pass=samples_pass,
        samples_warn=samples_warn,
        samples_fail=samples_fail,
        mean_coverage_depth=mean_coverage_depth,
        mean_genome_coverage=mean_genome_coverage,
        total_variants_called=total_variants,
        unique_variants=0,  # Placeholder - requires variant deduplication logic
    )


def build_sample_metrics(sample_id: str, sample_data: dict) -> SampleMetrics:
    """
    Build a SampleMetrics model from raw sample data.

    Args:
        sample_id: Sample identifier
        sample_data: Dict with coverage, alignment, variants, etc.

    Returns:
        Validated SampleMetrics model
    """
    coverage = sample_data.get("coverage", {})

    # Build AlignmentMetrics from coverage data
    # In Phase 1, we only have coverage metrics, so read counts are placeholders
    alignment = AlignmentMetrics(
        total_reads=sample_data.get("alignment", {}).get("total_reads", 0),
        mapped_reads=sample_data.get("alignment", {}).get("mapped_reads", 0),
        mapping_rate=sample_data.get("alignment", {}).get("mapping_rate", 0.0),
        mean_coverage=coverage.get("mean_coverage", 0.0),
        median_coverage=coverage.get("median_coverage"),
        genome_coverage_at_1x=coverage.get("genome_coverage_at_1x", 0.0),
        genome_coverage_at_10x=coverage.get("genome_coverage_at_10x", 0.0),
        genome_coverage_at_100x=coverage.get("genome_coverage_at_100x", 0.0),
    )

    return SampleMetrics(
        sample_id=sample_id,
        qc_status=sample_data.get("qc_status", QCStatus.FAIL),
        qc_notes=sample_data.get("qc_notes", []),
        alignment=alignment,
        variants=sample_data.get("variants"),
        consensus=sample_data.get("consensus"),
        metagenomics=sample_data.get("metagenomics"),
        haplotyping=sample_data.get("haplotyping"),
    )


@app.command("assemble")
def assemble(
    metrics_dir: Annotated[
        Path,
        typer.Option(
            "--metrics-dir",
            "-m",
            help="Directory containing *_metrics.json files",
            exists=True,
            file_okay=False,
            dir_okay=True,
            readable=True,
            resolve_path=True,
        ),
    ],
    output_dir: Annotated[
        Path,
        typer.Option(
            "--output-dir",
            "-o",
            help="Output directory for reports",
            file_okay=False,
            dir_okay=True,
            resolve_path=True,
        ),
    ],
    platform: Annotated[
        str,
        typer.Option(
            "--platform",
            "-p",
            help="Sequencing platform (ont or illumina)",
        ),
    ],
    reference_name: Annotated[
        str,
        typer.Option(
            "--reference-name",
            "-r",
            help="Reference sequence name for report title",
        ),
    ],
    multiqc_dir: Annotated[
        Path | None,
        typer.Option(
            "--multiqc-dir",
            help="Output directory for MultiQC files (defaults to output-dir)",
            file_okay=False,
            dir_okay=True,
            resolve_path=True,
        ),
    ] = None,
    viz_dir: Annotated[
        Path | None,
        typer.Option(
            "--viz-dir",
            help="Output directory for visualizations (defaults to output-dir/visualizations)",
            file_okay=False,
            dir_okay=True,
            resolve_path=True,
        ),
    ] = None,
    config_template: Annotated[
        Path | None,
        typer.Option(
            "--config-template",
            help="Path to MultiQC config template YAML",
            exists=True,
            file_okay=True,
            dir_okay=False,
            readable=True,
            resolve_path=True,
        ),
    ] = None,
    qc_thresholds: Annotated[
        str | None,
        typer.Option(
            "--qc-thresholds",
            help="QC thresholds as JSON string",
        ),
    ] = None,
) -> None:
    """
    Assemble OneRoof report from collected metrics.

    Loads per-sample metrics JSON files, determines QC status, builds the
    unified report, and generates MultiQC files and visualizations.
    """
    # Parse QC thresholds
    thresholds = DEFAULT_QC_THRESHOLDS.copy()
    if qc_thresholds:
        try:
            user_thresholds = json.loads(qc_thresholds)
            thresholds.update(user_thresholds)
        except json.JSONDecodeError as e:
            console.print(f"[red]Error parsing --qc-thresholds: {e}[/red]")
            raise typer.Exit(1) from e

    # Set default directories
    if multiqc_dir is None:
        multiqc_dir = output_dir
    if viz_dir is None:
        viz_dir = output_dir / "visualizations"

    # Create output directories
    output_dir.mkdir(parents=True, exist_ok=True)
    multiqc_dir.mkdir(parents=True, exist_ok=True)
    viz_dir.mkdir(parents=True, exist_ok=True)

    # Load sample metrics
    console.print(f"[blue]Loading metrics from {metrics_dir}[/blue]")
    samples = load_sample_metrics(metrics_dir)

    if not samples:
        console.print("[red]No sample metrics found![/red]")
        raise typer.Exit(1)

    console.print(f"[green]Loaded metrics for {len(samples)} samples[/green]")

    # Determine QC status for each sample
    for sample_id, sample_data in samples.items():
        status, notes = determine_qc_status(sample_data, thresholds)
        sample_data["qc_status"] = status
        sample_data["qc_notes"] = notes

    # Build run metadata
    # Note: In Phase 1, we use minimal metadata. Full metadata comes from
    # Nextflow params in later phases.
    run_metadata = RunMetadata(
        platform=Platform(platform.lower()),
        reference=ReferenceInfo(
            name=reference_name,
            path="",  # Not available in Phase 1
            length=0,  # Not available in Phase 1
        ),
        primers=PrimerInfo(provided=False),
        parameters=RunParameters(
            min_depth_coverage=10,
            min_consensus_freq=0.5,
            min_variant_frequency=0.2,
        ),
    )

    # Build sample metrics models
    sample_metrics = {
        sample_id: build_sample_metrics(sample_id, sample_data)
        for sample_id, sample_data in samples.items()
    }

    # Compute summary
    summary = compute_summary(samples)

    # Build the report
    report = OneRoofReport(
        generated_at=datetime.now(),
        run_metadata=run_metadata,
        summary=summary,
        samples=sample_metrics,
    )

    # Write full JSON report
    full_report_path = output_dir / "oneroof_report.json"
    full_report_path.write_text(report.model_dump_json(indent=2))
    console.print(f"[green]Wrote full report to {full_report_path}[/green]")

    # Write summary JSON report (no per-sample details)
    summary_report = report.model_copy(update={"samples": {}})
    summary_report_path = output_dir / "oneroof_report_summary.json"
    summary_report_path.write_text(summary_report.model_dump_json(indent=2))
    console.print(f"[green]Wrote summary report to {summary_report_path}[/green]")

    # Generate MultiQC files
    console.print("[blue]Generating MultiQC files...[/blue]")

    # General stats TSV
    general_stats_path = multiqc_dir / "oneroof_general_stats_mqc.tsv"
    generate_general_stats_tsv(samples, general_stats_path)
    console.print(f"[green]Wrote {general_stats_path}[/green]")

    # Coverage table TSV
    coverage_table_path = multiqc_dir / "oneroof_coverage_table_mqc.tsv"
    generate_coverage_table_tsv(samples, coverage_table_path)
    console.print(f"[green]Wrote {coverage_table_path}[/green]")

    # MultiQC config
    multiqc_config_path = multiqc_dir / "multiqc_config.yaml"
    generate_multiqc_config(
        run_metadata.model_dump(),
        multiqc_config_path,
        template_path=config_template,
    )
    console.print(f"[green]Wrote {multiqc_config_path}[/green]")

    # Generate visualizations
    console.print("[blue]Generating visualizations...[/blue]")

    # Prepare coverage metrics list for visualizations
    coverage_metrics_list = [
        sample_data.get("coverage", {})
        for sample_data in samples.values()
        if sample_data.get("coverage")
    ]

    if coverage_metrics_list:
        # Coverage heatmap
        heatmap_paths = coverage_summary_heatmap(
            coverage_metrics_list,
            viz_dir / "coverage_heatmap",
            formats=["html"],
        )
        for path in heatmap_paths:
            console.print(f"[green]Wrote {path}[/green]")

        # Coverage bar chart
        bar_paths = coverage_bar_chart(
            coverage_metrics_list,
            viz_dir / "coverage_bar",
            formats=["html"],
        )
        for path in bar_paths:
            console.print(f"[green]Wrote {path}[/green]")

    console.print("[bold green]Report assembly complete![/bold green]")


if __name__ == "__main__":
    app()
