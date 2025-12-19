"""
Coverage metrics extractor for OneRoof reporting.

Parses bedtools genomecov output (BED format with -bga flag) to compute
coverage statistics for a sample. The output is a per-base BED file with
columns: chrom, start, end, depth.

Example input (from `bedtools genomecov -bga`):
    MN908947.3    0       54      0
    MN908947.3    54      100     15
    MN908947.3    100     200     150
    ...
"""

from pathlib import Path

import polars as pl
from pydantic import BaseModel, Field


class CoverageMetrics(BaseModel):
    """Intermediate coverage metrics extracted from genomecov BED output."""

    sample_id: str
    total_bases: int = Field(ge=0, description="Total bases in reference")
    mean_coverage: float = Field(ge=0, description="Weighted mean depth")
    median_coverage: float = Field(ge=0, description="Weighted median depth")
    genome_coverage_at_1x: float = Field(
        ge=0, le=1, description="Fraction of genome with ≥1x coverage"
    )
    genome_coverage_at_10x: float = Field(
        ge=0, le=1, description="Fraction of genome with ≥10x coverage"
    )
    genome_coverage_at_100x: float = Field(
        ge=0, le=1, description="Fraction of genome with ≥100x coverage"
    )
    min_coverage: int = Field(ge=0, description="Minimum depth observed")
    max_coverage: int = Field(ge=0, description="Maximum depth observed")


def load_coverage_bed(bed_path: Path) -> pl.DataFrame:
    """
    Load a bedtools genomecov BED file.

    Args:
        bed_path: Path to the per-base BED file from `bedtools genomecov -bga`

    Returns:
        DataFrame with columns: chrom, start, end, depth
    """
    return pl.read_csv(
        bed_path,
        separator="\t",
        has_header=False,
        new_columns=["chrom", "start", "end", "depth"],
        schema={
            "chrom": pl.Utf8,
            "start": pl.Int64,
            "end": pl.Int64,
            "depth": pl.Int64,
        },
    )


def compute_coverage_stats(df: pl.DataFrame) -> dict:
    """
    Compute coverage statistics from a genomecov DataFrame.

    Args:
        df: DataFrame with columns: chrom, start, end, depth

    Returns:
        Dictionary with computed coverage statistics
    """
    # Calculate region lengths
    df = df.with_columns((pl.col("end") - pl.col("start")).alias("length"))

    total_bases = df["length"].sum()

    if total_bases == 0:
        return {
            "total_bases": 0,
            "mean_coverage": 0.0,
            "median_coverage": 0.0,
            "genome_coverage_at_1x": 0.0,
            "genome_coverage_at_10x": 0.0,
            "genome_coverage_at_100x": 0.0,
            "min_coverage": 0,
            "max_coverage": 0,
        }

    # Weighted mean: sum(depth * length) / total_bases
    weighted_sum = (df["depth"] * df["length"]).sum()
    mean_coverage = weighted_sum / total_bases

    # Weighted median: expand depths by length and find median
    # For efficiency with large genomes, we compute this from the cumulative distribution
    sorted_df = df.sort("depth")
    sorted_df = sorted_df.with_columns(
        (pl.col("length").cum_sum() / total_bases).alias("cumulative_frac")
    )
    # Find the first row where cumulative fraction >= 0.5
    median_row = sorted_df.filter(pl.col("cumulative_frac") >= 0.5).head(1)
    median_coverage = float(median_row["depth"][0]) if len(median_row) > 0 else 0.0

    # Coverage at thresholds
    bases_at_1x = df.filter(pl.col("depth") >= 1)["length"].sum()
    bases_at_10x = df.filter(pl.col("depth") >= 10)["length"].sum()
    bases_at_100x = df.filter(pl.col("depth") >= 100)["length"].sum()

    # Min and max coverage
    # Note: Polars min/max return PythonLiteral which includes int at runtime
    min_cov_value = df["depth"].min()
    max_cov_value = df["depth"].max()
    min_coverage = 0 if min_cov_value is None else int(min_cov_value)  # type: ignore[arg-type]
    max_coverage = 0 if max_cov_value is None else int(max_cov_value)  # type: ignore[arg-type]

    return {
        "total_bases": int(total_bases),
        "mean_coverage": float(mean_coverage),
        "median_coverage": float(median_coverage),
        "genome_coverage_at_1x": float(bases_at_1x / total_bases),
        "genome_coverage_at_10x": float(bases_at_10x / total_bases),
        "genome_coverage_at_100x": float(bases_at_100x / total_bases),
        "min_coverage": min_coverage,
        "max_coverage": max_coverage,
    }


def extract(sample_id: str, bed_path: Path) -> CoverageMetrics:
    """
    Extract coverage metrics from a bedtools genomecov BED file.

    Args:
        sample_id: Sample identifier
        bed_path: Path to the per-base BED file

    Returns:
        Validated CoverageMetrics model
    """
    df = load_coverage_bed(bed_path)
    stats = compute_coverage_stats(df)
    return CoverageMetrics(sample_id=sample_id, **stats)
