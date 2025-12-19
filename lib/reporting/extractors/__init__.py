"""
Metrics Extractors for OneRoof Reporting.

Each extractor module parses output from a specific pipeline stage and produces
a validated Pydantic model containing the extracted metrics. These per-sample
JSON files are later assembled into the final OneRoof report.

Extractors:
    coverage: Parse bedtools genomecov output for coverage statistics
    alignment: Parse BAM files for read mapping statistics
    variants: Parse SnpSift variant effects TSV
    consensus: Parse consensus FASTA for sequence statistics
    metagenomics: Parse Sylph profile output
    haplotyping: Parse Devider output for haplotype statistics (Nanopore only)
"""

from .alignment import extract as extract_alignment
from .consensus import extract as extract_consensus
from .coverage import extract as extract_coverage
from .haplotyping import extract as extract_haplotyping
from .metagenomics import extract as extract_metagenomics
from .variants import extract as extract_variants

__all__ = [
    "extract_alignment",
    "extract_consensus",
    "extract_coverage",
    "extract_haplotyping",
    "extract_metagenomics",
    "extract_variants",
]
