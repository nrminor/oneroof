"""
MultiQC custom content file generators for OneRoof reporting.

MultiQC supports custom content via specially formatted files with embedded
YAML configuration in comments. This module generates TSV files that MultiQC
will automatically detect and include in reports.

Reference: https://multiqc.info/docs/development/modules/#custom-content

File format:
    # id: 'oneroof_coverage'
    # section_name: 'Coverage Summary'
    # plot_type: 'table'
    # pconfig:
    #     namespace: 'OneRoof'
    # headers:
    #     col1:
    #         title: 'Column 1'
    #         format: '{:,.0f}'
    Sample	col1	col2
    sample1	100	200
"""

from pathlib import Path
from typing import Any

import polars as pl
import yaml


def _write_tsv_with_header(
    output_path: Path,
    yaml_config: dict[str, Any],
    headers: list[str],
    rows: list[list[Any]],
) -> None:
    """
    Write a TSV file with embedded YAML header for MultiQC.

    Args:
        output_path: Path to write the TSV file
        yaml_config: YAML configuration dict to embed in comments
        headers: Column headers
        rows: Data rows (list of lists)
    """
    lines = []

    # Serialize YAML config and prefix each line with #
    yaml_str = yaml.dump(yaml_config, default_flow_style=False, sort_keys=False)
    for yaml_line in yaml_str.splitlines():
        lines.append(f"# {yaml_line}")

    # Write header row
    lines.append("\t".join(headers))

    # Write data rows
    for row in rows:
        lines.append("\t".join(str(v) for v in row))

    output_path.write_text("\n".join(lines) + "\n")


def generate_general_stats_tsv(
    samples: dict[str, dict],
    output_path: Path,
) -> Path:
    """
    Generate a TSV file for MultiQC General Statistics table.

    This adds OneRoof metrics to the main General Statistics table that
    appears at the top of every MultiQC report.

    Args:
        samples: Dict mapping sample_id to metrics dict
        output_path: Path to write the TSV file

    Returns:
        Path to the written file
    """
    yaml_config = {
        "id": "oneroof_general_stats",
        "plot_type": "generalstats",
        "pconfig": {
            "namespace": "OneRoof",
        },
        "headers": {
            "mean_coverage": {
                "title": "Mean Cov",
                "description": "Mean coverage depth",
                "format": "{:,.1f}",
                "scale": "Blues",
            },
            "genome_coverage_pct": {
                "title": "Genome %",
                "description": "Percentage of genome with ≥10x coverage",
                "format": "{:,.1f}",
                "suffix": "%",
                "scale": "RdYlGn",
                "min": 0,
                "max": 100,
                "cond_formatting_rules": {
                    "pass": [{"gt": 95}],
                    "warn": [{"gt": 80}],
                    "fail": [{"lt": 80}],
                },
            },
            "completeness_pct": {
                "title": "Complete %",
                "description": "Consensus completeness (non-N bases)",
                "format": "{:,.1f}",
                "suffix": "%",
                "scale": "RdYlGn",
                "min": 0,
                "max": 100,
                "cond_formatting_rules": {
                    "pass": [{"gt": 98}],
                    "warn": [{"gt": 90}],
                    "fail": [{"lt": 90}],
                },
            },
            "n_pct": {
                "title": "N %",
                "description": "Percentage of N bases in consensus",
                "format": "{:,.1f}",
                "suffix": "%",
                "scale": "OrRd",
                "min": 0,
                "max": 100,
                "cond_formatting_rules": {
                    "pass": [{"lt": 1}],
                    "warn": [{"lt": 5}],
                    "fail": [{"gt": 5}],
                },
            },
        },
    }

    headers = [
        "Sample",
        "mean_coverage",
        "genome_coverage_pct",
        "completeness_pct",
        "n_pct",
    ]
    rows = []

    for sample_id, metrics in samples.items():
        # Extract coverage metrics
        coverage = metrics.get("coverage", {})
        consensus = metrics.get("consensus", {})

        mean_cov = coverage.get("mean_coverage", 0)
        genome_pct = coverage.get("genome_coverage_at_10x", 0) * 100
        completeness = consensus.get("completeness", 0) * 100 if consensus else 0
        n_pct = consensus.get("n_percentage", 0) if consensus else 0

        rows.append([sample_id, mean_cov, genome_pct, completeness, n_pct])

    _write_tsv_with_header(output_path, yaml_config, headers, rows)
    return output_path


def generate_coverage_table_tsv(
    samples: dict[str, dict],
    output_path: Path,
) -> Path:
    """
    Generate a TSV file for a detailed coverage table section.

    Args:
        samples: Dict mapping sample_id to metrics dict
        output_path: Path to write the TSV file

    Returns:
        Path to the written file
    """
    yaml_config = {
        "id": "oneroof_coverage_table",
        "section_name": "Coverage Summary",
        "description": "Detailed coverage statistics from OneRoof analysis",
        "plot_type": "table",
        "pconfig": {
            "namespace": "OneRoof",
            "id": "oneroof_coverage_table",
            "title": "Coverage Statistics",
        },
        "headers": {
            "mean_coverage": {
                "title": "Mean Coverage",
                "format": "{:,.1f}",
                "scale": "Blues",
            },
            "median_coverage": {
                "title": "Median Coverage",
                "format": "{:,.1f}",
                "scale": "Blues",
            },
            "genome_1x": {
                "title": "≥1x %",
                "format": "{:,.1f}",
                "suffix": "%",
                "scale": "RdYlGn",
            },
            "genome_10x": {
                "title": "≥10x %",
                "format": "{:,.1f}",
                "suffix": "%",
                "scale": "RdYlGn",
            },
            "genome_100x": {
                "title": "≥100x %",
                "format": "{:,.1f}",
                "suffix": "%",
                "scale": "RdYlGn",
            },
            "min_coverage": {
                "title": "Min",
                "format": "{:,.0f}",
            },
            "max_coverage": {
                "title": "Max",
                "format": "{:,.0f}",
            },
        },
    }

    headers = [
        "Sample",
        "mean_coverage",
        "median_coverage",
        "genome_1x",
        "genome_10x",
        "genome_100x",
        "min_coverage",
        "max_coverage",
    ]
    rows = []

    for sample_id, metrics in samples.items():
        coverage = metrics.get("coverage", {})

        rows.append(
            [
                sample_id,
                coverage.get("mean_coverage", 0),
                coverage.get("median_coverage", 0),
                coverage.get("genome_coverage_at_1x", 0) * 100,
                coverage.get("genome_coverage_at_10x", 0) * 100,
                coverage.get("genome_coverage_at_100x", 0) * 100,
                coverage.get("min_coverage", 0),
                coverage.get("max_coverage", 0),
            ],
        )

    _write_tsv_with_header(output_path, yaml_config, headers, rows)
    return output_path


def generate_variant_bargraph_tsv(
    samples: dict[str, dict],
    output_path: Path,
) -> Path | None:
    """
    Generate a TSV file for a MultiQC stacked bar graph of variant types.

    Creates a horizontal stacked bar chart showing SNPs, insertions, deletions,
    and MNPs per sample. Samples with zero total variants are excluded.

    Args:
        samples: Dict mapping sample_id to metrics dict (must have "variants" key)
        output_path: Path to write the TSV file

    Returns:
        Path to the written file, or None if no variant data available
    """
    # Collect samples with variant data
    rows = []
    for sample_id, metrics in samples.items():
        variants = metrics.get("variants", {})
        if not variants:
            continue

        snps = variants.get("snps", 0)
        insertions = variants.get("insertions", 0)
        deletions = variants.get("deletions", 0)
        mnps = variants.get("mnps", 0)

        # Skip samples with no variants
        if snps + insertions + deletions + mnps == 0:
            continue

        rows.append([sample_id, snps, insertions, deletions, mnps])

    if not rows:
        return None

    yaml_config = {
        "id": "oneroof_variants",
        "section_name": "Variant Summary",
        "description": "Variant types called per sample",
        "plot_type": "bargraph",
        "pconfig": {
            "id": "oneroof_variant_bargraph",
            "title": "OneRoof: Variant Types",
            "ylab": "Count",
        },
    }

    headers = ["Sample", "SNPs", "Insertions", "Deletions", "MNPs"]

    _write_tsv_with_header(output_path, yaml_config, headers, rows)
    return output_path


# Performance tier thresholds (as fraction of median across all amplicons)
# Matches thresholds in lib/reporting/visualizations/amplicon_efficiency.py
_TIER_THRESHOLDS = {
    "good": 0.5,  # >= 50% of overall median
    "moderate": 0.1,  # >= 10% of overall median
    # below 10% = poor
}


def generate_amplicon_efficiency_tsv(
    summary_path: Path,
    output_path: Path,
) -> Path | None:
    """
    Generate a TSV file for a MultiQC table of amplicon efficiency metrics.

    Creates a sortable table showing per-amplicon performance with columns for
    median reads, dropout rate, sample count, and performance tier. Tiers are
    color-coded via conditional formatting (good=green, moderate=yellow, poor=red).

    Args:
        summary_path: Path to amplicon_summary.tsv
        output_path: Path to write the TSV file

    Returns:
        Path to the written file, or None if no data available
    """
    if not summary_path.exists():
        return None

    data = pl.read_csv(summary_path, separator="\t")

    if len(data) == 0:
        return None

    # Compute per-amplicon stats (same logic as amplicon_efficiency.py)
    stats = data.group_by("amplicon_name").agg(
        pl.col("reads").median().alias("median_reads"),
        (pl.col("reads") == 0).mean().alias("dropout_rate"),
        pl.col("reads").count().alias("sample_count"),
    )

    if len(stats) == 0:
        return None

    # Compute overall median for tier classification
    raw_median = stats["median_reads"].median()
    if (
        raw_median is None
        or not isinstance(raw_median, (int, float))
        or raw_median == 0
    ):
        overall_median = 1.0
    else:
        overall_median = float(raw_median)

    # Assign performance tiers
    good_threshold = _TIER_THRESHOLDS["good"]
    moderate_threshold = _TIER_THRESHOLDS["moderate"]

    stats = stats.with_columns(
        pl.when(pl.col("median_reads") >= overall_median * good_threshold)
        .then(pl.lit("good"))
        .when(pl.col("median_reads") >= overall_median * moderate_threshold)
        .then(pl.lit("moderate"))
        .otherwise(pl.lit("poor"))
        .alias("performance_tier"),
    )

    # Sort by median reads descending
    stats = stats.sort("median_reads", descending=True)

    yaml_config = {
        "id": "oneroof_amplicon_efficiency",
        "section_name": "Amplicon Efficiency",
        "description": "Per-amplicon performance metrics across all samples",
        "plot_type": "table",
        "pconfig": {
            "id": "oneroof_amplicon_efficiency_table",
            "title": "OneRoof: Amplicon Efficiency",
            "sortRows": "true",
        },
        "headers": {
            "median_reads": {
                "title": "Median Reads",
                "format": "{:,.0f}",
                "scale": "Blues",
            },
            "dropout_pct": {
                "title": "Dropout %",
                "format": "{:,.1f}",
                "suffix": "%",
                "scale": "OrRd",
            },
            "sample_count": {
                "title": "Samples",
                "format": "{:,.0f}",
            },
            "performance_tier": {
                "title": "Tier",
                "cond_formatting_rules": {
                    "pass": [{"eq": "good"}],
                    "warn": [{"eq": "moderate"}],
                    "fail": [{"eq": "poor"}],
                },
            },
        },
    }

    headers = [
        "Amplicon",
        "median_reads",
        "dropout_pct",
        "sample_count",
        "performance_tier",
    ]

    rows = []
    for row in stats.iter_rows(named=True):
        rows.append(
            [
                row["amplicon_name"],
                row["median_reads"] or 0,
                (row["dropout_rate"] or 0) * 100,  # Convert to percentage
                row["sample_count"],
                row["performance_tier"],
            ]
        )

    _write_tsv_with_header(output_path, yaml_config, headers, rows)
    return output_path


def generate_amplicon_heatmap_tsv(
    summary_path: Path,
    output_path: Path,
) -> Path | None:
    """
    Generate a TSV file for a MultiQC heatmap of amplicon read counts.

    Creates a heatmap with samples on the Y-axis (alphabetically sorted) and
    amplicons on the X-axis (sorted by genomic position). Values are read counts.

    Args:
        summary_path: Path to amplicon_summary.tsv
        output_path: Path to write the TSV file

    Returns:
        Path to the written file, or None if no data available
    """
    if not summary_path.exists():
        return None

    data = pl.read_csv(summary_path, separator="\t")

    if len(data) == 0:
        return None

    # Get amplicon order by genomic position
    amplicon_order = (
        data.lazy()
        .select("amplicon_name", "start_pos")
        .unique()
        .sort("start_pos")
        .collect()["amplicon_name"]
        .to_list()
    )

    # Get sample order alphabetically
    sample_order = sorted(data["sample_name"].unique().to_list())

    if not amplicon_order or not sample_order:
        return None

    # Build the heatmap data as a matrix
    # For MultiQC heatmap TSV: first column is row labels (samples),
    # subsequent columns are amplicons (in position order)
    # Header row has amplicon names

    # Pivot the data to get sample × amplicon matrix
    pivot = data.pivot(
        on="amplicon_name",
        index="sample_name",
        values="reads",
    )

    # Reorder columns to match amplicon genomic order
    # First column should be sample_name, then amplicons in order
    ordered_cols = ["sample_name"] + [a for a in amplicon_order if a in pivot.columns]
    pivot = pivot.select(ordered_cols)

    # Sort rows alphabetically by sample name
    pivot = pivot.sort("sample_name")

    # Write as TSV with MultiQC header
    yaml_config = {
        "id": "oneroof_amplicon_heatmap",
        "section_name": "Amplicon Coverage Heatmap",
        "description": "Read counts per sample and amplicon (amplicons sorted by genomic position)",
        "plot_type": "heatmap",
        "pconfig": {
            "id": "oneroof_amplicon_heatmap_plot",
            "title": "OneRoof: Amplicon Coverage",
            "xlab": "Amplicon",
            "ylab": "Sample",
            "square": "false",
            "xcats_samples": "false",
            "min": "0",
        },
    }

    # Build header and rows
    headers = ["Sample"] + [a for a in amplicon_order if a in pivot.columns]
    rows = []
    for row in pivot.iter_rows(named=True):
        row_data = [row["sample_name"]]
        for amp in amplicon_order:
            if amp in row:
                row_data.append(row[amp] if row[amp] is not None else 0)
        rows.append(row_data)

    _write_tsv_with_header(output_path, yaml_config, headers, rows)
    return output_path


# Default relative path to the MultiQC config template from project root
DEFAULT_MULTIQC_TEMPLATE = Path("conf/multiqc_config.yaml")


def generate_multiqc_config(
    run_metadata: dict,
    output_path: Path,
    template_path: Path | None = None,
) -> Path:
    """
    Generate a MultiQC configuration file from template.

    Loads the base template and substitutes dynamic values like platform
    and reference name.

    Args:
        run_metadata: Run metadata dict with platform, reference info, etc.
        output_path: Path to write the YAML config file
        template_path: Path to template file. If None, uses DEFAULT_MULTIQC_TEMPLATE
                      relative to current working directory.

    Returns:
        Path to the written file

    Raises:
        FileNotFoundError: If template_path does not exist
    """
    if template_path is None:
        template_path = DEFAULT_MULTIQC_TEMPLATE

    if not template_path.exists():
        msg = f"MultiQC config template not found: {template_path}"
        raise FileNotFoundError(msg)

    platform = run_metadata.get("platform", "unknown").upper()
    reference_name = run_metadata.get("reference", {}).get("name", "unknown")

    # Load template and substitute placeholders
    template_content = template_path.read_text()
    config_content = template_content.format(
        platform=platform,
        reference=reference_name,
    )

    output_path.write_text(config_content)
    return output_path
