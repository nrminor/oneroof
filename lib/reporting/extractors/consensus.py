"""
Consensus metrics extractor for OneRoof reporting.

Parses consensus FASTA files from ivar consensus to compute sequence
quality statistics including N content, completeness, and GC content.

Example input (from ivar consensus):
    >Consensus_sample1.consensus_threshold_0.5_quality_0
    ACGTNNNNACGT...
"""

from pathlib import Path

from Bio import SeqIO
from pydantic import BaseModel, Field


class ConsensusExtractedMetrics(BaseModel):
    """Metrics extracted from consensus FASTA files."""

    sample_id: str
    length: int = Field(ge=0, description="Total consensus sequence length")
    n_count: int = Field(ge=0, description="Number of N bases")
    n_percentage: float = Field(ge=0, le=100, description="Percentage of N bases")
    completeness: float = Field(
        ge=0,
        le=1,
        description="Fraction of non-N bases (1 - n_count/length)",
    )
    ambiguous_bases: int = Field(
        ge=0,
        description="Number of ambiguous bases (non-ACGTN IUPAC codes)",
    )
    gc_content: float = Field(ge=0, le=100, description="GC percentage of non-N bases")


# Standard nucleotides
STANDARD_BASES = set("ACGTN")


def extract(sample_id: str, fasta_path: Path) -> ConsensusExtractedMetrics:
    """
    Extract consensus metrics from a FASTA file.

    If the FASTA contains multiple sequences, metrics are computed across
    all sequences combined (summed lengths, counts, etc.).

    Args:
        sample_id: Sample identifier
        fasta_path: Path to the consensus FASTA file

    Returns:
        Validated ConsensusExtractedMetrics model
    """
    total_length = 0
    n_count = 0
    g_count = 0
    c_count = 0
    a_count = 0
    t_count = 0
    ambiguous_count = 0

    for record in SeqIO.parse(fasta_path, "fasta"):
        seq = str(record.seq).upper()
        total_length += len(seq)

        for base in seq:
            if base == "N":
                n_count += 1
            elif base == "G":
                g_count += 1
            elif base == "C":
                c_count += 1
            elif base == "A":
                a_count += 1
            elif base == "T":
                t_count += 1
            elif base not in STANDARD_BASES:
                # IUPAC ambiguity codes: R, Y, S, W, K, M, B, D, H, V
                ambiguous_count += 1

    # Handle empty file
    if total_length == 0:
        return ConsensusExtractedMetrics(
            sample_id=sample_id,
            length=0,
            n_count=0,
            n_percentage=0.0,
            completeness=0.0,
            ambiguous_bases=0,
            gc_content=0.0,
        )

    # Calculate metrics
    n_percentage = (n_count / total_length) * 100
    completeness = 1 - (n_count / total_length)

    # GC content is calculated from non-N bases only
    acgt_count = a_count + c_count + g_count + t_count
    if acgt_count > 0:
        gc_content = ((g_count + c_count) / acgt_count) * 100
    else:
        # All Ns or empty - GC content is undefined, report as 0
        gc_content = 0.0

    return ConsensusExtractedMetrics(
        sample_id=sample_id,
        length=total_length,
        n_count=n_count,
        n_percentage=n_percentage,
        completeness=completeness,
        ambiguous_bases=ambiguous_count,
        gc_content=gc_content,
    )
