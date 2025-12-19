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

from .amplicon_efficiency import (
    amplicon_dropout_scatter,
    amplicon_heatmap,
    amplicon_ranking_bar,
    load_amplicon_summary,
    prepare_amplicon_stats,
    prepare_heatmap_data,
)
from .coverage_heatmap import (
    coverage_bar_chart,
    coverage_summary_heatmap,
    prepare_heatmap_data,
)
from .qc_dashboard import (
    completeness_distribution,
    coverage_distribution,
    qc_scatter,
    qc_status_summary,
)
from .utils import COLORS, register_oneroof_theme, save_chart
from .variant_summary import (
    prepare_variant_effect_data,
    prepare_variant_type_data,
    variant_effect_bar,
    variant_type_bar,
)

__all__ = [
    "COLORS",
    "amplicon_dropout_scatter",
    "amplicon_heatmap",
    "amplicon_ranking_bar",
    "completeness_distribution",
    "coverage_bar_chart",
    "coverage_distribution",
    "coverage_summary_heatmap",
    "load_amplicon_summary",
    "prepare_amplicon_stats",
    "prepare_heatmap_data",
    "prepare_variant_effect_data",
    "prepare_variant_type_data",
    "qc_scatter",
    "qc_status_summary",
    "register_oneroof_theme",
    "save_chart",
    "variant_effect_bar",
    "variant_type_bar",
]
