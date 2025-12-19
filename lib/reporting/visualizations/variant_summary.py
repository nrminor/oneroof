"""
Variant summary visualizations for OneRoof reporting.

Generates stacked bar charts showing variant counts by mutation type
(SNP, insertion, deletion, MNP) and by predicted effect (missense,
synonymous, etc.) across samples.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import altair as alt
import polars as pl

from .utils import register_oneroof_theme, save_chart

if TYPE_CHECKING:
    from pathlib import Path


# Color scheme for mutation types
MUTATION_TYPE_COLORS = {
    "SNP": "#3b82f6",  # Blue
    "insertion": "#22c55e",  # Green
    "deletion": "#ef4444",  # Red
    "MNP": "#a855f7",  # Purple
}

# Color scheme for effect types (using a categorical palette)
EFFECT_TYPE_SCHEME = "tableau20"


def prepare_variant_type_data(samples: dict[str, dict]) -> pl.DataFrame:
    """
    Prepare variant metrics for mutation type bar chart.

    Transforms per-sample variant metrics into long format suitable
    for a stacked bar chart.

    Args:
        samples: Dict mapping sample_id to metrics dict, where each
                 metrics dict has a "variants" key containing
                 snps, insertions, deletions, mnps counts.

    Returns:
        DataFrame with columns: sample_id, mutation_type, count
    """
    rows = []

    for sample_id, metrics in samples.items():
        variants = metrics.get("variants", {})
        if not variants:
            continue

        rows.extend(
            [
                {
                    "sample_id": sample_id,
                    "mutation_type": "SNP",
                    "count": variants.get("snps", 0),
                },
                {
                    "sample_id": sample_id,
                    "mutation_type": "insertion",
                    "count": variants.get("insertions", 0),
                },
                {
                    "sample_id": sample_id,
                    "mutation_type": "deletion",
                    "count": variants.get("deletions", 0),
                },
                {
                    "sample_id": sample_id,
                    "mutation_type": "MNP",
                    "count": variants.get("mnps", 0),
                },
            ],
        )

    if not rows:
        return pl.DataFrame(
            schema={
                "sample_id": pl.Utf8,
                "mutation_type": pl.Utf8,
                "count": pl.Int64,
            },
        )

    return pl.DataFrame(rows)


def prepare_variant_effect_data(samples: dict[str, dict]) -> pl.DataFrame:
    """
    Prepare variant metrics for effect type bar chart.

    Transforms per-sample variant effect counts into long format suitable
    for a stacked bar chart.

    Args:
        samples: Dict mapping sample_id to metrics dict, where each
                 metrics dict has a "variants" key containing
                 a "by_effect" dict mapping effect names to counts.

    Returns:
        DataFrame with columns: sample_id, effect_type, count
    """
    rows = []

    for sample_id, metrics in samples.items():
        variants = metrics.get("variants", {})
        by_effect = variants.get("by_effect", {})
        if not by_effect:
            continue

        for effect_type, count in by_effect.items():
            rows.append(
                {
                    "sample_id": sample_id,
                    "effect_type": effect_type,
                    "count": count,
                },
            )

    if not rows:
        return pl.DataFrame(
            schema={
                "sample_id": pl.Utf8,
                "effect_type": pl.Utf8,
                "count": pl.Int64,
            },
        )

    return pl.DataFrame(rows)


def variant_type_bar(
    samples: dict[str, dict],
    output_path: Path,
    formats: list[str] | None = None,
    title: str = "Variants by Mutation Type",
) -> list[Path]:
    """
    Generate a stacked bar chart showing variant counts by mutation type.

    Args:
        samples: Dict mapping sample_id to metrics dict
        output_path: Base output path (without extension)
        formats: Output formats (default: ["html"])
        title: Chart title

    Returns:
        List of paths to saved files
    """
    if formats is None:
        formats = ["html"]

    register_oneroof_theme()

    data = prepare_variant_type_data(samples)

    if len(data) == 0:
        return []

    # Define mutation type order for consistent stacking
    type_order = ["SNP", "insertion", "deletion", "MNP"]

    # Calculate total variants per sample for sorting
    sample_totals = data.group_by("sample_id").agg(pl.col("count").sum().alias("total"))
    sample_order = sample_totals.sort("total", descending=True)["sample_id"].to_list()

    # Build color scale from our defined colors
    domain = list(MUTATION_TYPE_COLORS.keys())
    range_ = list(MUTATION_TYPE_COLORS.values())

    chart = (
        alt.Chart(data)
        .mark_bar()
        .encode(
            alt.X("count:Q").stack("zero").title("Variant Count"),
            alt.Y("sample_id:N").sort(sample_order).title("Sample"),
            alt.Color("mutation_type:N")
            .scale(domain=domain, range=range_)
            .sort(type_order)
            .title("Mutation Type"),
            alt.Order("mutation_type:N"),
            tooltip=[
                alt.Tooltip("sample_id:N", title="Sample"),
                alt.Tooltip("mutation_type:N", title="Type"),
                alt.Tooltip("count:Q", title="Count", format=",d"),
            ],
        )
        .properties(
            width=400,
            height=max(100, len(sample_order) * 25),
            title=title,
        )
    )

    return save_chart(chart, output_path, formats)


def variant_effect_bar(
    samples: dict[str, dict],
    output_path: Path,
    formats: list[str] | None = None,
    title: str = "Variants by Effect Type",
) -> list[Path]:
    """
    Generate a stacked bar chart showing variant counts by effect type.

    Args:
        samples: Dict mapping sample_id to metrics dict
        output_path: Base output path (without extension)
        formats: Output formats (default: ["html"])
        title: Chart title

    Returns:
        List of paths to saved files
    """
    if formats is None:
        formats = ["html"]

    register_oneroof_theme()

    data = prepare_variant_effect_data(samples)

    if len(data) == 0:
        return []

    # Calculate total variants per sample for sorting
    sample_totals = data.group_by("sample_id").agg(pl.col("count").sum().alias("total"))
    sample_order = sample_totals.sort("total", descending=True)["sample_id"].to_list()

    # Sort effect types by total count across all samples for legend order
    effect_totals = data.group_by("effect_type").agg(
        pl.col("count").sum().alias("total"),
    )
    effect_order = effect_totals.sort("total", descending=True)["effect_type"].to_list()

    chart = (
        alt.Chart(data)
        .mark_bar()
        .encode(
            alt.X("count:Q").stack("zero").title("Variant Count"),
            alt.Y("sample_id:N").sort(sample_order).title("Sample"),
            alt.Color("effect_type:N")
            .scale(scheme=EFFECT_TYPE_SCHEME)
            .sort(effect_order)
            .title("Effect Type"),
            alt.Order("effect_type:N"),
            tooltip=[
                alt.Tooltip("sample_id:N", title="Sample"),
                alt.Tooltip("effect_type:N", title="Effect"),
                alt.Tooltip("count:Q", title="Count", format=",d"),
            ],
        )
        .properties(
            width=400,
            height=max(100, len(sample_order) * 25),
            title=title,
        )
    )

    return save_chart(chart, output_path, formats)
