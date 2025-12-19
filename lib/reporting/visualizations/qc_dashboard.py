"""
QC dashboard visualizations for OneRoof reporting.

Generates summary visualizations for quality control status across samples:
- QC status pie/donut chart (pass/warn/fail distribution)
- Coverage distribution histogram
- Completeness distribution histogram
- QC scatter plot (coverage vs completeness, colored by status)
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import altair as alt
import polars as pl

from .utils import COLORS, register_oneroof_theme, save_chart

if TYPE_CHECKING:
    from pathlib import Path


def prepare_qc_status_data(samples: dict[str, dict]) -> pl.DataFrame:
    """
    Prepare QC status counts for pie chart.

    Args:
        samples: Dict mapping sample_id to metrics dict, where each
                 metrics dict has a "qc_status" key.

    Returns:
        DataFrame with columns: status, count
    """
    status_counts = {"pass": 0, "warn": 0, "fail": 0}

    for metrics in samples.values():
        status = metrics.get("qc_status", "fail")
        # Handle both enum and string values
        if hasattr(status, "value"):
            status = status.value
        status = str(status).lower()
        if status in status_counts:
            status_counts[status] += 1

    rows = [
        {"status": status, "count": count}
        for status, count in status_counts.items()
        if count > 0
    ]

    if not rows:
        return pl.DataFrame(schema={"status": pl.Utf8, "count": pl.Int64})

    return pl.DataFrame(rows)


def prepare_distribution_data(samples: dict[str, dict]) -> pl.DataFrame:
    """
    Prepare coverage and completeness data for distribution plots.

    Args:
        samples: Dict mapping sample_id to metrics dict

    Returns:
        DataFrame with columns: sample_id, mean_coverage, completeness, qc_status
    """
    rows = []

    for sample_id, metrics in samples.items():
        coverage = metrics.get("coverage", {})
        consensus = metrics.get("consensus", {})
        qc_status = metrics.get("qc_status", "fail")

        # Handle both enum and string values
        if hasattr(qc_status, "value"):
            qc_status = qc_status.value
        qc_status = str(qc_status).lower()

        mean_cov = coverage.get("mean_coverage", 0)
        completeness = consensus.get("completeness", 0) if consensus else 0

        rows.append(
            {
                "sample_id": sample_id,
                "mean_coverage": float(mean_cov),
                "completeness": float(completeness) * 100,  # Convert to percentage
                "qc_status": qc_status,
            },
        )

    if not rows:
        return pl.DataFrame(
            schema={
                "sample_id": pl.Utf8,
                "mean_coverage": pl.Float64,
                "completeness": pl.Float64,
                "qc_status": pl.Utf8,
            },
        )

    return pl.DataFrame(rows)


def qc_status_summary(
    samples: dict[str, dict],
    output_path: Path,
    formats: list[str] | None = None,
    title: str = "QC Status Summary",
) -> list[Path]:
    """
    Generate a donut chart showing QC status distribution.

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

    data = prepare_qc_status_data(samples)

    if len(data) == 0:
        return []

    # Define colors for each status
    status_colors = {
        "pass": COLORS["pass"],
        "warn": COLORS["warn"],
        "fail": COLORS["fail"],
    }

    domain = list(status_colors.keys())
    range_ = list(status_colors.values())

    # Calculate total for percentage labels
    total = data["count"].sum()

    # Add percentage column
    data = data.with_columns((pl.col("count") / total * 100).alias("percentage"))

    chart = (
        alt.Chart(data)
        .mark_arc(innerRadius=50)
        .encode(
            alt.Theta("count:Q"),
            alt.Color("status:N").scale(domain=domain, range=range_).title("Status"),
            tooltip=[
                alt.Tooltip("status:N", title="Status"),
                alt.Tooltip("count:Q", title="Count"),
                alt.Tooltip("percentage:Q", title="Percentage", format=".1f"),
            ],
        )
        .properties(width=250, height=250, title=title)
    )

    return save_chart(chart, output_path, formats)


def coverage_distribution(
    samples: dict[str, dict],
    output_path: Path,
    formats: list[str] | None = None,
    title: str = "Coverage Distribution",
) -> list[Path]:
    """
    Generate a histogram of mean coverage values across samples.

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

    data = prepare_distribution_data(samples)

    if len(data) == 0:
        return []

    # For single sample, show a simple bar instead of histogram
    if len(data) == 1:
        chart = (
            alt.Chart(data)
            .mark_bar(width=50)
            .encode(
                alt.X("sample_id:N").title("Sample"),
                alt.Y("mean_coverage:Q").title("Mean Coverage"),
                alt.Color("qc_status:N")
                .scale(
                    domain=["pass", "warn", "fail"],
                    range=[COLORS["pass"], COLORS["warn"], COLORS["fail"]],
                )
                .title("QC Status"),
                tooltip=[
                    alt.Tooltip("sample_id:N", title="Sample"),
                    alt.Tooltip("mean_coverage:Q", title="Mean Coverage", format=".1f"),
                ],
            )
            .properties(width=200, height=300, title=title)
        )
    else:
        chart = (
            alt.Chart(data)
            .mark_bar()
            .encode(
                alt.X("mean_coverage:Q").bin(maxbins=20).title("Mean Coverage"),
                alt.Y("count()").title("Number of Samples"),
                tooltip=[
                    alt.Tooltip("mean_coverage:Q", bin=True, title="Coverage Range"),
                    alt.Tooltip("count()", title="Count"),
                ],
            )
            .properties(width=400, height=300, title=title)
        )

    return save_chart(chart, output_path, formats)


def completeness_distribution(
    samples: dict[str, dict],
    output_path: Path,
    formats: list[str] | None = None,
    title: str = "Completeness Distribution",
) -> list[Path]:
    """
    Generate a histogram of completeness values across samples.

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

    data = prepare_distribution_data(samples)

    if len(data) == 0:
        return []

    # For single sample, show a simple bar instead of histogram
    if len(data) == 1:
        chart = (
            alt.Chart(data)
            .mark_bar(width=50)
            .encode(
                alt.X("sample_id:N").title("Sample"),
                alt.Y("completeness:Q")
                .title("Completeness (%)")
                .scale(domain=[0, 100]),
                alt.Color("qc_status:N")
                .scale(
                    domain=["pass", "warn", "fail"],
                    range=[COLORS["pass"], COLORS["warn"], COLORS["fail"]],
                )
                .title("QC Status"),
                tooltip=[
                    alt.Tooltip("sample_id:N", title="Sample"),
                    alt.Tooltip(
                        "completeness:Q",
                        title="Completeness (%)",
                        format=".1f",
                    ),
                ],
            )
            .properties(width=200, height=300, title=title)
        )
    else:
        chart = (
            alt.Chart(data)
            .mark_bar()
            .encode(
                alt.X("completeness:Q").bin(maxbins=20).title("Completeness (%)"),
                alt.Y("count()").title("Number of Samples"),
                tooltip=[
                    alt.Tooltip("completeness:Q", bin=True, title="Completeness Range"),
                    alt.Tooltip("count()", title="Count"),
                ],
            )
            .properties(width=400, height=300, title=title)
        )

    return save_chart(chart, output_path, formats)


def qc_scatter(
    samples: dict[str, dict],
    output_path: Path,
    formats: list[str] | None = None,
    title: str = "Coverage vs Completeness",
) -> list[Path]:
    """
    Generate a scatter plot of coverage vs completeness, colored by QC status.

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

    data = prepare_distribution_data(samples)

    if len(data) == 0:
        return []

    chart = (
        alt.Chart(data)
        .mark_circle(size=100)
        .encode(
            alt.X("mean_coverage:Q").title("Mean Coverage").scale(type="symlog"),
            alt.Y("completeness:Q").title("Completeness (%)").scale(domain=[0, 100]),
            alt.Color("qc_status:N")
            .scale(
                domain=["pass", "warn", "fail"],
                range=[COLORS["pass"], COLORS["warn"], COLORS["fail"]],
            )
            .title("QC Status"),
            tooltip=[
                alt.Tooltip("sample_id:N", title="Sample"),
                alt.Tooltip("mean_coverage:Q", title="Mean Coverage", format=".1f"),
                alt.Tooltip("completeness:Q", title="Completeness (%)", format=".1f"),
                alt.Tooltip("qc_status:N", title="QC Status"),
            ],
        )
        .properties(width=400, height=300, title=title)
        .interactive()
    )

    return save_chart(chart, output_path, formats)
