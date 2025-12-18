"""
Coverage heatmap visualization for OneRoof reporting.

Generates a multi-sample coverage heatmap showing coverage depth across
samples. Useful for identifying samples with low coverage or systematic
dropout patterns.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import altair as alt
import polars as pl

from .utils import register_oneroof_theme, save_chart

if TYPE_CHECKING:
    from pathlib import Path


def prepare_heatmap_data(coverage_metrics: list[dict]) -> pl.DataFrame:
    """
    Prepare coverage metrics for heatmap visualization.

    For the initial implementation, this creates a simple summary table
    with one row per sample showing key coverage statistics. Future versions
    may support per-position or per-amplicon breakdowns.

    Args:
        coverage_metrics: List of coverage metrics dicts (one per sample),
                         each containing sample_id, mean_coverage,
                         genome_coverage_at_1x, etc.

    Returns:
        DataFrame with columns suitable for heatmap visualization
    """
    if not coverage_metrics:
        return pl.DataFrame(
            schema={
                "sample_id": pl.Utf8,
                "metric": pl.Utf8,
                "value": pl.Float64,
            },
        )

    # Convert to long format for heatmap (sample x metric)
    rows = []
    for metrics in coverage_metrics:
        sample_id = metrics["sample_id"]
        rows.extend(
            [
                {
                    "sample_id": sample_id,
                    "metric": "Mean Coverage",
                    "value": float(metrics.get("mean_coverage", 0)),
                },
                {
                    "sample_id": sample_id,
                    "metric": "Genome % ≥1x",
                    "value": float(metrics.get("genome_coverage_at_1x", 0)) * 100,
                },
                {
                    "sample_id": sample_id,
                    "metric": "Genome % ≥10x",
                    "value": float(metrics.get("genome_coverage_at_10x", 0)) * 100,
                },
                {
                    "sample_id": sample_id,
                    "metric": "Genome % ≥100x",
                    "value": float(metrics.get("genome_coverage_at_100x", 0)) * 100,
                },
            ],
        )

    return pl.DataFrame(rows)


def coverage_summary_heatmap(
    coverage_metrics: list[dict],
    output_path: Path,
    formats: list[str] | None = None,
    title: str = "Coverage Summary",
) -> list[Path]:
    """
    Generate a coverage summary heatmap showing key metrics across samples.

    Args:
        coverage_metrics: List of coverage metrics dicts (one per sample)
        output_path: Base output path (without extension)
        formats: Output formats (default: ["html"])
        title: Chart title

    Returns:
        List of paths to saved files
    """
    if formats is None:
        formats = ["html"]

    register_oneroof_theme()

    data = prepare_heatmap_data(coverage_metrics)

    if len(data) == 0:
        # Return empty list if no data
        return []

    # Define metric order for consistent display
    metric_order = [
        "Mean Coverage",
        "Genome % ≥1x",
        "Genome % ≥10x",
        "Genome % ≥100x",
    ]

    # Create heatmap
    # Use symlog scale to handle zero values gracefully
    chart = (
        alt.Chart(data)
        .mark_rect()
        .encode(
            alt.X("metric:N")
            .axis(labelAngle=-45, labelOverlap=False)
            .sort(metric_order)
            .title(None),
            alt.Y("sample_id:N").title("Sample"),
            alt.Color("value:Q").scale(type="symlog", scheme="viridis", constant=1).title("Value"),
            tooltip=[
                alt.Tooltip("sample_id:N", title="Sample"),
                alt.Tooltip("metric:N", title="Metric"),
                alt.Tooltip("value:Q", title="Value", format=".2f"),
            ],
        )
        .properties(
            width=300,
            height=max(100, len(coverage_metrics) * 25),
            title=title,
        )
        .configure_view(strokeWidth=0)
        .configure_axis(domain=False)
    )

    return save_chart(chart, output_path, formats)


def coverage_bar_chart(
    coverage_metrics: list[dict],
    output_path: Path,
    formats: list[str] | None = None,
    title: str = "Mean Coverage by Sample",
) -> list[Path]:
    """
    Generate a bar chart showing mean coverage per sample.

    This is often more readable than a heatmap for single-metric comparisons.

    Args:
        coverage_metrics: List of coverage metrics dicts (one per sample)
        output_path: Base output path (without extension)
        formats: Output formats (default: ["html"])
        title: Chart title

    Returns:
        List of paths to saved files
    """
    if formats is None:
        formats = ["html"]

    register_oneroof_theme()

    if not coverage_metrics:
        return []

    data = pl.DataFrame(
        [
            {
                "sample_id": m["sample_id"],
                "mean_coverage": float(m.get("mean_coverage", 0)),
            }
            for m in coverage_metrics
        ],
    )

    chart = (
        alt.Chart(data)
        .mark_bar()
        .encode(
            alt.X("mean_coverage:Q").title("Mean Coverage"),
            alt.Y("sample_id:N").sort("-x").title("Sample"),
            alt.Color("mean_coverage:Q")
            .scale(type="symlog", scheme="viridis", constant=1)
            .legend(None),
            tooltip=[
                alt.Tooltip("sample_id:N", title="Sample"),
                alt.Tooltip("mean_coverage:Q", title="Mean Coverage", format=".1f"),
            ],
        )
        .properties(
            width=400,
            height=max(100, len(coverage_metrics) * 25),
            title=title,
        )
    )

    return save_chart(chart, output_path, formats)
