"""
Alignment metrics extractor for OneRoof reporting.

Parses BAM files using pysam to extract read mapping statistics.
This extractor focuses on read-level metrics (counts, mapping rate),
while coverage metrics are handled by the coverage extractor.

Note: The BAM file must be indexed (.bai file present) for efficient
index statistics retrieval. If index statistics are unavailable,
the extractor falls back to iterating through all reads.
"""

from pathlib import Path

import pysam
from pydantic import BaseModel, Field


class AlignmentExtractedMetrics(BaseModel):
    """Metrics extracted from aligned BAM files."""

    sample_id: str
    total_reads: int = Field(ge=0, description="Total reads in BAM file")
    mapped_reads: int = Field(ge=0, description="Reads aligned to reference")
    unmapped_reads: int = Field(ge=0, description="Reads not aligned")
    mapping_rate: float = Field(
        ge=0,
        le=1,
        description="Fraction of reads successfully mapped",
    )


def extract(sample_id: str, bam_path: Path) -> AlignmentExtractedMetrics:
    """
    Extract alignment metrics from a BAM file.

    Args:
        sample_id: Sample identifier
        bam_path: Path to the indexed BAM file

    Returns:
        Validated AlignmentExtractedMetrics model
    """
    mapped_reads = 0
    unmapped_reads = 0

    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        # Try to get counts from BAM index (more efficient for large files)
        # The .mapped and .unmapped attributes use the index when available
        try:
            mapped_reads = bam.mapped
            unmapped_reads = bam.unmapped
        except (ValueError, AttributeError):
            # Index not available, count directly
            mapped_reads, unmapped_reads = _count_reads_directly(bam)

    total_reads = mapped_reads + unmapped_reads
    mapping_rate = mapped_reads / total_reads if total_reads > 0 else 0.0

    return AlignmentExtractedMetrics(
        sample_id=sample_id,
        total_reads=total_reads,
        mapped_reads=mapped_reads,
        unmapped_reads=unmapped_reads,
        mapping_rate=mapping_rate,
    )


def _count_reads_directly(bam: pysam.AlignmentFile) -> tuple[int, int]:
    """
    Count mapped and unmapped reads by iterating through the BAM.

    Args:
        bam: Open pysam AlignmentFile

    Returns:
        Tuple of (mapped_reads, unmapped_reads)
    """
    mapped = 0
    unmapped = 0

    bam.reset()
    for read in bam.fetch(until_eof=True):
        if read.is_unmapped:
            unmapped += 1
        else:
            mapped += 1

    return mapped, unmapped
