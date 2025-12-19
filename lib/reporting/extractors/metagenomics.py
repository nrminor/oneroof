"""
Metagenomics metrics extractor for OneRoof reporting.

Parses Sylph profile TSV output to extract metagenomic profiling results.
Sylph performs species-level taxonomic profiling with ANI estimates.

Example input columns (from `sylph profile --estimate-unknown`):
    Sample_file, Genome_file, Taxonomic_abundance, Sequence_abundance,
    Adjusted_ANI, Eff_cov, ANI_5-95_percentile, Eff_lambda, ...

The --estimate-unknown flag scales abundances to account for unknown sequences,
so if Taxonomic_abundance values sum to less than 100%, the remainder represents
the unknown/unclassified fraction.
"""

from pathlib import Path

import polars as pl
from pydantic import BaseModel, Field


class MetagenomicsHit(BaseModel):
    """A single hit from metagenomic profiling."""

    taxon: str = Field(description="Genome/taxon identifier")
    ani: float = Field(ge=0, le=100, description="Adjusted ANI estimate")
    relative_abundance: float = Field(
        ge=0,
        le=100,
        description="Taxonomic abundance percentage",
    )


class MetagenomicsExtractedMetrics(BaseModel):
    """Metrics extracted from Sylph profile output."""

    sample_id: str
    top_hits: list[MetagenomicsHit] = Field(
        default_factory=list,
        description="Top hits by taxonomic abundance",
    )
    total_hits: int = Field(ge=0, description="Total number of hits in profile")
    total_abundance: float = Field(
        ge=0,
        le=100,
        description="Sum of all taxonomic abundances",
    )
    unknown_fraction: float = Field(
        ge=0,
        le=1,
        description="Fraction of reads not classified (1 - total_abundance/100)",
    )


def load_sylph_profile(tsv_path: Path) -> pl.DataFrame:
    """
    Load a Sylph profile TSV file.

    Args:
        tsv_path: Path to the Sylph profile TSV file

    Returns:
        DataFrame with profile data, or empty DataFrame if no hits
    """
    try:
        df = pl.read_csv(
            tsv_path,
            separator="\t",
            has_header=True,
            null_values=["", ".", "NA", "NA-NA"],
            infer_schema_length=1000,
        )
    except pl.exceptions.NoDataError:
        return pl.DataFrame()

    if len(df) == 0:
        return pl.DataFrame()

    # Standardize column names to lowercase with underscores
    return df.select(
        pl.all().name.map(lambda c: c.lower().replace("-", "_").replace(" ", "_")),
    )


def extract(
    sample_id: str,
    profile_tsv: Path,
    top_n: int = 5,
) -> MetagenomicsExtractedMetrics:
    """
    Extract metagenomics metrics from a Sylph profile TSV file.

    Args:
        sample_id: Sample identifier
        profile_tsv: Path to the Sylph profile TSV file
        top_n: Number of top hits to include (default 5)

    Returns:
        Validated MetagenomicsExtractedMetrics model
    """
    df = load_sylph_profile(profile_tsv)

    # Handle empty profiles
    if len(df) == 0:
        return MetagenomicsExtractedMetrics(
            sample_id=sample_id,
            top_hits=[],
            total_hits=0,
            total_abundance=0.0,
            unknown_fraction=1.0,
        )

    total_hits = len(df)

    # Get taxonomic abundance column
    abundance_col = "taxonomic_abundance"
    if abundance_col not in df.columns:
        # Fallback if column name differs
        for col in df.columns:
            if "taxonomic" in col.lower() and "abundance" in col.lower():
                abundance_col = col
                break
        else:
            # No abundance column found, return empty metrics
            return MetagenomicsExtractedMetrics(
                sample_id=sample_id,
                top_hits=[],
                total_hits=total_hits,
                total_abundance=0.0,
                unknown_fraction=1.0,
            )

    # Calculate total abundance
    total_abundance = float(df[abundance_col].sum())

    # Unknown fraction is what's left after accounting for classified reads
    # Abundances are percentages (0-100), so divide by 100 for fraction
    unknown_fraction = max(0.0, 1.0 - (total_abundance / 100.0))

    # Sort by abundance and get top N hits
    df_sorted = df.sort(abundance_col, descending=True).head(top_n)

    # Extract top hits
    top_hits = []
    genome_col = _find_column(df_sorted, ["genome_file", "genome", "contig_name"])
    ani_col = _find_column(df_sorted, ["adjusted_ani", "ani", "naive_ani"])

    if genome_col and ani_col:
        for row in df_sorted.iter_rows(named=True):
            taxon = str(row[genome_col]) if row[genome_col] else "unknown"
            ani = float(row[ani_col]) if row[ani_col] is not None else 0.0
            abundance = (
                float(row[abundance_col]) if row[abundance_col] is not None else 0.0
            )

            top_hits.append(
                MetagenomicsHit(
                    taxon=taxon,
                    ani=ani,
                    relative_abundance=abundance,
                ),
            )

    return MetagenomicsExtractedMetrics(
        sample_id=sample_id,
        top_hits=top_hits,
        total_hits=total_hits,
        total_abundance=total_abundance,
        unknown_fraction=unknown_fraction,
    )


def _find_column(df: pl.DataFrame, candidates: list[str]) -> str | None:
    """Find the first matching column name from a list of candidates."""
    for candidate in candidates:
        if candidate in df.columns:
            return candidate
    return None
