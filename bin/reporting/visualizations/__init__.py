"""
Visualization generators for OneRoof Reporting.

This subpackage provides Altair-based chart generators for the OneRoof report.
Charts are output as self-contained HTML files (with embedded Vega-Lite spec)
or as static SVG/PNG images.

Modules:
    utils: Theme registration and chart saving utilities
    coverage_heatmap: Multi-sample coverage heatmap (Phase 1)
    variant_summary: Variant type and effect bar charts (Phase 3)
    qc_dashboard: QC status and distribution visualizations (Phase 3)
    amplicon_efficiency: Amplicon performance charts (Phase 3)
"""

from .coverage_heatmap import (
    coverage_bar_chart,
    coverage_summary_heatmap,
    prepare_heatmap_data,
)
from .utils import COLORS, register_oneroof_theme, save_chart

__all__ = [
    "COLORS",
    "coverage_bar_chart",
    "coverage_summary_heatmap",
    "prepare_heatmap_data",
    "register_oneroof_theme",
    "save_chart",
]
