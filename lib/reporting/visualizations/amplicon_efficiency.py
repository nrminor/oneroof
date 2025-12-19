"""
Amplicon efficiency visualizations for OneRoof reporting.

Generates visualizations for amplicon performance analysis:
- Amplicon ranking bar chart (sorted by median reads, colored by performance tier)
- Amplicon dropout scatter plot (median reads vs dropout rate)
- Amplicon heatmap (samples x position-sorted amplicons, colored by read count)

These visualizations help identify problematic amplicons that may need
primer redesign or protocol optimization.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import altair as alt
import polars as pl

from .utils import COLORS, register_oneroof_theme, save_chart

if TYPE_CHECKING:
    from pathlib import Path

# Performance tier thresholds (as fraction of median across all amplicons)
TIER_THRESHOLDS = {
    "good": 0.5,  # >= 50% of overall median
    "moderate": 0.1,  # >= 10% of overall median
    # below 10% = poor
}


def load_amplicon_summary(summary_path: Path) -> pl.DataFrame:
    """
    Load amplicon summary TSV file.

    Args:
        summary_path: Path to amplicon_summary.tsv

    Returns:
        DataFrame with columns: sample_name, amplicon_name, start_pos, end_pos, reads
    """
    return pl.read_csv(summary_path, separator="\t")


def prepare_amplicon_stats(data: pl.DataFrame) -> pl.DataFrame:
    """
    Compute per-amplicon statistics across all samples.

    Args:
        data: DataFrame with columns: sample_name, amplicon_name, reads

    Returns:
        DataFrame with columns: amplicon_name, median_reads, dropout_rate,
                                sample_count, performance_tier
    """
    if len(data) == 0:
        return pl.DataFrame(
            schema={
                "amplicon_name": pl.Utf8,
                "median_reads": pl.Float64,
                "dropout_rate": pl.Float64,
                "sample_count": pl.Int64,
                "performance_tier": pl.Utf8,
            },
        )

    # Compute per-amplicon stats
    stats = data.group_by("amplicon_name").agg(
        pl.col("reads").median().alias("median_reads"),
        (pl.col("reads") == 0).mean().alias("dropout_rate"),
        pl.col("reads").count().alias("sample_count"),
    )

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

    # Thresholds as floats for multiplication
    good_threshold = TIER_THRESHOLDS["good"]
    moderate_threshold = TIER_THRESHOLDS["moderate"]

    # Assign performance tiers
    return stats.with_columns(
        pl.when(pl.col("median_reads") >= overall_median * good_threshold)
        .then(pl.lit("good"))
        .when(pl.col("median_reads") >= overall_median * moderate_threshold)
        .then(pl.lit("moderate"))
        .otherwise(pl.lit("poor"))
        .alias("performance_tier"),
    )


def amplicon_ranking_bar(
    summary_path: Path,
    output_path: Path,
    formats: list[str] | None = None,
    title: str = "Amplicon Efficiency Ranking",
) -> list[Path]:
    """
    Generate a bar chart of amplicons sorted by median read count.

    Amplicons are colored by performance tier (good/moderate/poor) based
    on their median reads relative to the overall median.

    Args:
        summary_path: Path to amplicon_summary.tsv
        output_path: Base output path (without extension)
        formats: Output formats (default: ["html"])
        title: Chart title

    Returns:
        List of paths to saved files
    """
    if formats is None:
        formats = ["html"]

    register_oneroof_theme()

    data = load_amplicon_summary(summary_path)
    stats = prepare_amplicon_stats(data)

    if len(stats) == 0:
        return []

    # Sort by median reads descending
    stats = stats.sort("median_reads", descending=True)

    # Define tier colors
    tier_colors = {
        "good": COLORS["pass"],
        "moderate": COLORS["warn"],
        "poor": COLORS["fail"],
    }

    chart = (
        alt.Chart(stats)
        .mark_bar()
        .encode(
            alt.X("amplicon_name:N")
            .sort(field="median_reads", order="descending")
            .title("Amplicon")
            .axis(labelAngle=-45),
            alt.Y("median_reads:Q").title("Median Reads"),
            alt.Color("performance_tier:N")
            .scale(
                domain=list(tier_colors.keys()),
                range=list(tier_colors.values()),
            )
            .title("Performance"),
            tooltip=[
                alt.Tooltip("amplicon_name:N", title="Amplicon"),
                alt.Tooltip("median_reads:Q", title="Median Reads", format=".0f"),
                alt.Tooltip("dropout_rate:Q", title="Dropout Rate", format=".1%"),
                alt.Tooltip("sample_count:Q", title="Samples"),
                alt.Tooltip("performance_tier:N", title="Tier"),
            ],
        )
        .properties(
            width=max(400, len(stats) * 20),  # Scale width with amplicon count
            height=300,
            title=title,
        )
    )

    return save_chart(chart, output_path, formats)


def amplicon_dropout_scatter(
    summary_path: Path,
    output_path: Path,
    formats: list[str] | None = None,
    title: str = "Amplicon Dropout Analysis",
) -> list[Path]:
    """
    Generate a scatter plot of median reads vs dropout rate.

    Helps identify problematic amplicons that have both low coverage
    and high dropout rates across samples.

    Args:
        summary_path: Path to amplicon_summary.tsv
        output_path: Base output path (without extension)
        formats: Output formats (default: ["html"])
        title: Chart title

    Returns:
        List of paths to saved files
    """
    if formats is None:
        formats = ["html"]

    register_oneroof_theme()

    data = load_amplicon_summary(summary_path)
    stats = prepare_amplicon_stats(data)

    if len(stats) == 0:
        return []

    # Define tier colors
    tier_colors = {
        "good": COLORS["pass"],
        "moderate": COLORS["warn"],
        "poor": COLORS["fail"],
    }

    chart = (
        alt.Chart(stats)
        .mark_circle(size=100)
        .encode(
            alt.X("median_reads:Q")
            .title("Median Reads")
            .scale(type="symlog"),  # Handle wide range of values
            alt.Y("dropout_rate:Q").title("Dropout Rate").scale(domain=[0, 1]),
            alt.Color("performance_tier:N")
            .scale(
                domain=list(tier_colors.keys()),
                range=list(tier_colors.values()),
            )
            .title("Performance"),
            tooltip=[
                alt.Tooltip("amplicon_name:N", title="Amplicon"),
                alt.Tooltip("median_reads:Q", title="Median Reads", format=".0f"),
                alt.Tooltip("dropout_rate:Q", title="Dropout Rate", format=".1%"),
                alt.Tooltip("sample_count:Q", title="Samples"),
                alt.Tooltip("performance_tier:N", title="Tier"),
            ],
        )
        .properties(width=400, height=300, title=title)
        .interactive()
    )

    return save_chart(chart, output_path, formats)


def prepare_amplicon_heatmap_data(data: pl.LazyFrame) -> pl.DataFrame:
    """
    Prepare data for heatmap visualization.

    Adds log-transformed reads column for color scale while preserving
    raw reads for tooltips.

    Args:
        data: LazyFrame with columns: sample_name, amplicon_name, start_pos, reads

    Returns:
        DataFrame with additional log_reads column (log1p transformed)
    """
    return data.with_columns(
        pl.col("reads").cast(pl.Float64).log1p().alias("log_reads"),
    ).collect()


def amplicon_heatmap(
    summary_path: Path,
    output_path: Path,
    formats: list[str] | None = None,
    title: str = "Amplicon Coverage Heatmap",
) -> list[Path]:
    """
    Generate a heatmap of read counts across samples and amplicons.

    Y-axis shows samples (alphabetically sorted), X-axis shows amplicons
    (sorted by genomic position). Color intensity represents read count
    on a log scale, with tooltips showing the raw values.

    Args:
        summary_path: Path to amplicon_summary.tsv
        output_path: Base output path (without extension)
        formats: Output formats (default: ["html"])
        title: Chart title

    Returns:
        List of paths to saved files
    """
    if formats is None:
        formats = ["html"]

    register_oneroof_theme()

    # Load lazily and prepare data
    data_lf = pl.scan_csv(summary_path, separator="\t")
    data = prepare_amplicon_heatmap_data(data_lf)

    if len(data) == 0:
        return []

    # Get amplicon order by start_pos
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

    chart = (
        alt.Chart(data)
        .mark_rect()
        .encode(
            alt.X("amplicon_name:N")
            .sort(amplicon_order)
            .title("Amplicon (by genomic position)")
            .axis(labelAngle=-45),
            alt.Y("sample_name:N").sort(sample_order).title("Sample"),
            alt.Color("log_reads:Q").scale(scheme="blues").title("Reads (log)"),
            tooltip=[
                alt.Tooltip("sample_name:N", title="Sample"),
                alt.Tooltip("amplicon_name:N", title="Amplicon"),
                alt.Tooltip("reads:Q", title="Reads", format=","),
                alt.Tooltip("start_pos:Q", title="Position", format=","),
            ],
        )
        .properties(
            width=max(400, len(amplicon_order) * 15),
            height=max(300, len(sample_order) * 20),
            title=title,
        )
    )

    return save_chart(chart, output_path, formats)
