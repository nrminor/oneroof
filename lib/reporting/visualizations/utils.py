"""
Altair theme and chart utilities for OneRoof visualizations.

Provides a consistent visual theme and helper functions for saving charts
in multiple formats (HTML, SVG, PNG).
"""

from __future__ import annotations

from pathlib import Path

import altair as alt

# OneRoof color palette
COLORS = {
    "primary": "#2563eb",  # Blue
    "secondary": "#64748b",  # Slate
    "success": "#22c55e",  # Green
    "warning": "#f59e0b",  # Amber
    "danger": "#ef4444",  # Red
    "pass": "#22c55e",
    "warn": "#f59e0b",
    "fail": "#ef4444",
}

# Sequential color scheme for heatmaps
HEATMAP_SCHEME = "viridis"


@alt.theme.register("oneroof", enable=True)
def _oneroof_theme() -> alt.theme.ThemeConfig:
    """Return the OneRoof Altair theme configuration."""
    return alt.theme.ThemeConfig(
        {
            "background": "#ffffff",
            "config": {
                "title": {
                    "fontSize": 16,
                    "fontWeight": "bold",
                    "anchor": "start",
                    "color": "#1e293b",
                },
                "axis": {
                    "labelFontSize": 11,
                    "titleFontSize": 12,
                    "titleColor": "#475569",
                    "labelColor": "#64748b",
                    "gridColor": "#e2e8f0",
                    "domainColor": "#cbd5e1",
                },
                "legend": {
                    "labelFontSize": 11,
                    "titleFontSize": 12,
                    "titleColor": "#475569",
                    "labelColor": "#64748b",
                },
                "view": {
                    "strokeWidth": 0,
                },
                "range": {
                    "category": [
                        COLORS["primary"],
                        COLORS["secondary"],
                        COLORS["success"],
                        COLORS["warning"],
                        COLORS["danger"],
                    ],
                },
            },
        }
    )


def register_oneroof_theme() -> None:
    """
    Register and enable the OneRoof Altair theme.

    Call this once at the start of visualization generation to ensure
    consistent styling across all charts.

    Note: The theme is auto-registered via decorator, but this function
    ensures it's enabled when called explicitly.
    """
    alt.theme.enable("oneroof")


def save_chart(
    chart: alt.Chart,
    output_path: Path,
    formats: list[str] | None = None,
) -> list[Path]:
    """
    Save an Altair chart in one or more formats.

    Args:
        chart: The Altair chart to save
        output_path: Base output path (without extension)
        formats: List of formats to save. Defaults to ["html"].
                 Supported: "html", "svg", "png"

    Returns:
        List of paths to saved files
    """
    if formats is None:
        formats = ["html"]

    saved_paths = []
    output_path = Path(output_path)

    for fmt in formats:
        if fmt == "html":
            path = output_path.with_suffix(".html")
            chart.save(path, embed_options={"renderer": "svg"})
            saved_paths.append(path)
        elif fmt == "svg":
            path = output_path.with_suffix(".svg")
            chart.save(path)
            saved_paths.append(path)
        elif fmt == "png":
            path = output_path.with_suffix(".png")
            chart.save(path, scale_factor=2)  # 2x for retina
            saved_paths.append(path)
        else:
            msg = f"Unsupported format: {fmt}. Use 'html', 'svg', or 'png'."
            raise ValueError(msg)

    return saved_paths
