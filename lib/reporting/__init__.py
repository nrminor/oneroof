"""
OneRoof Reporting Package.

This package provides metrics extraction, report assembly, and visualization
generation for the OneRoof amplicon sequencing pipeline.

Subpackages:
    extractors: Per-sample metrics extraction from pipeline outputs
    visualizations: Altair-based chart generation (Phase 1+)

Modules:
    schema: Pydantic models defining the canonical report JSON structure
    multiqc: MultiQC custom content file generation (Phase 1+)
"""

from .multiqc import (
    DEFAULT_MULTIQC_TEMPLATE,
    generate_amplicon_efficiency_tsv,
    generate_amplicon_heatmap_tsv,
    generate_coverage_table_tsv,
    generate_general_stats_tsv,
    generate_multiqc_config,
    generate_variant_bargraph_tsv,
)

__all__ = [
    "DEFAULT_MULTIQC_TEMPLATE",
    "generate_amplicon_efficiency_tsv",
    "generate_amplicon_heatmap_tsv",
    "generate_coverage_table_tsv",
    "generate_general_stats_tsv",
    "generate_multiqc_config",
    "generate_variant_bargraph_tsv",
]
