"""
Variant metrics extractor for OneRoof reporting.

Parses SnpSift extractFields TSV output to compute variant statistics.
The TSV contains annotated variant calls with columns for position, allele
frequencies, and predicted effects.

Example input columns (from SnpSift extractFields):
    CHROM, REF, POS, ALT, AF, AC, DP, GEN[0].REF_DP, GEN[0].ALT_DP,
    GEN[0].ALT_FREQ, MQ, ANN[0].GENE, ANN[0].EFFECT, ANN[0].HGVS_P,
    ANN[0].CDS_POS, ANN[0].AA_POS
"""

from pathlib import Path

import polars as pl
from pydantic import BaseModel, Field


class VariantExtractedMetrics(BaseModel):
    """Metrics extracted from SnpSift variant effects TSV."""

    sample_id: str
    total_called: int = Field(ge=0, description="Total variants called")
    consensus_variants: int = Field(
        ge=0,
        description="Variants at or above consensus threshold",
    )
    subclonal_variants: int = Field(
        ge=0,
        description="Variants below consensus threshold",
    )
    snps: int = Field(ge=0, description="Single nucleotide polymorphisms")
    insertions: int = Field(ge=0, description="Insertion variants")
    deletions: int = Field(ge=0, description="Deletion variants")
    mnps: int = Field(ge=0, description="Multi-nucleotide polymorphisms")
    by_effect: dict[str, int] = Field(
        default_factory=dict,
        description="Variant counts by SnpEff effect type",
    )


def load_variant_effects_tsv(tsv_path: Path) -> pl.DataFrame:
    """
    Load a SnpSift extractFields TSV file.

    Handles the column name transformations from SnpSift's output format
    (e.g., 'GEN[0].REF_DP' -> 'gen_0_ref_dp').

    Args:
        tsv_path: Path to the variant effects TSV file

    Returns:
        DataFrame with standardized column names, or empty DataFrame if no variants
    """
    # Read with polars, handling empty files
    try:
        df = pl.read_csv(
            tsv_path,
            separator="\t",
            has_header=True,
            null_values=["", "."],
            infer_schema_length=1000,
        )
    except pl.exceptions.NoDataError:
        # Empty file (no header even)
        return pl.DataFrame()

    if len(df) == 0:
        # Header only, no data rows
        return pl.DataFrame()

    # Standardize column names: lowercase, replace brackets and dots
    return df.select(
        pl.all().name.map(
            lambda c: c.lower().replace("[", "_").replace("]", "").replace(".", "_"),
        ),
    )


def classify_mutation_type(ref: str, alt: str) -> str:
    """
    Classify mutation type based on reference and alternate allele lengths.

    Args:
        ref: Reference allele sequence
        alt: Alternate allele sequence

    Returns:
        One of: "SNP", "insertion", "deletion", "MNP"
    """
    ref_len = len(ref)
    alt_len = len(alt)

    if ref_len == alt_len:
        return "SNP" if ref_len == 1 else "MNP"
    if ref_len > alt_len:
        return "deletion"
    return "insertion"


def extract(
    sample_id: str,
    effects_tsv: Path,
    consensus_threshold: float = 0.5,
) -> VariantExtractedMetrics:
    """
    Extract variant metrics from a SnpSift effects TSV file.

    Args:
        sample_id: Sample identifier
        effects_tsv: Path to the SnpSift extractFields TSV file
        consensus_threshold: Allele frequency threshold for consensus calls (default 0.5)

    Returns:
        Validated VariantExtractedMetrics model
    """
    df = load_variant_effects_tsv(effects_tsv)

    # Handle empty files
    if len(df) == 0:
        return VariantExtractedMetrics(
            sample_id=sample_id,
            total_called=0,
            consensus_variants=0,
            subclonal_variants=0,
            snps=0,
            insertions=0,
            deletions=0,
            mnps=0,
            by_effect={},
        )

    total_called = len(df)

    # Count consensus vs subclonal variants
    # Use 'af' column (allele frequency from VCF)
    if "af" in df.columns:
        af_col = df["af"].cast(pl.Float64, strict=False).fill_null(0.0)
        consensus_variants = int((af_col >= consensus_threshold).sum())
        subclonal_variants = total_called - consensus_variants
    else:
        # If no AF column, assume all are consensus
        consensus_variants = total_called
        subclonal_variants = 0

    # Classify mutation types
    snps = 0
    insertions = 0
    deletions = 0
    mnps = 0

    if "ref" in df.columns and "alt" in df.columns:
        for ref, alt in zip(df["ref"].to_list(), df["alt"].to_list()):
            if ref is None or alt is None:
                continue
            mut_type = classify_mutation_type(str(ref), str(alt))
            if mut_type == "SNP":
                snps += 1
            elif mut_type == "insertion":
                insertions += 1
            elif mut_type == "deletion":
                deletions += 1
            elif mut_type == "MNP":
                mnps += 1

    # Count by effect type
    by_effect: dict[str, int] = {}
    effect_col_name = "ann_0_effect"
    if effect_col_name in df.columns:
        effect_counts = (
            df.group_by(effect_col_name)
            .agg(pl.len().alias("count"))
            .filter(pl.col(effect_col_name).is_not_null())
        )
        for row in effect_counts.iter_rows(named=True):
            effect = row[effect_col_name]
            if effect:
                by_effect[str(effect)] = int(row["count"])

    return VariantExtractedMetrics(
        sample_id=sample_id,
        total_called=total_called,
        consensus_variants=consensus_variants,
        subclonal_variants=subclonal_variants,
        snps=snps,
        insertions=insertions,
        deletions=deletions,
        mnps=mnps,
        by_effect=by_effect,
    )
