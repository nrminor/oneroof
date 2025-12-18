"""
Metrics Extractors for OneRoof Reporting.

Each extractor module parses output from a specific pipeline stage and produces
a validated Pydantic model containing the extracted metrics. These per-sample
JSON files are later assembled into the final OneRoof report.

Extractors:
    coverage: Parse bedtools genomecov output for coverage statistics
    alignment: Parse BAM files for read mapping statistics (Phase 2)
    variants: Parse SnpSift variant effects TSV (Phase 2)
    consensus: Parse consensus FASTA for sequence statistics (Phase 2)
    metagenomics: Parse Sylph profile output (Phase 2)
    haplotyping: Parse Devider output for haplotype statistics (Phase 2, ONT only)
"""

from .coverage import extract as extract_coverage

__all__ = ["extract_coverage"]
