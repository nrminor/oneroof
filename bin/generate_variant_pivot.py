#!/usr/bin/env python3

# /// script
# requires-python = ">= 3.10"
# dependencies = [
#     "patito",
#     "polars-lts-cpu",
#     "loguru",
# ]
# ///

from __future__ import annotations

import argparse
from pathlib import Path

import patito as pt
import polars as pl
from loguru import logger


def parse_command_line_args() -> argparse.Namespace:
    """
    TODO
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input_table",
        "-i",
        type=Path,
        required=True,
        help="Input table of VCF fields extracted from an annotated VCF with SnpSift",
    )

    return parser.parse_args()


class VariantPivotSchema(pt.Model):
    """Schema for validating variant pivot table input data."""

    contig: str = pt.Field(description="Contig/chromosome name")
    ref: str = pt.Field(min_length=1, description="Reference allele")
    pos: int = pt.Field(gt=0, description="Genomic position")
    alt: str = pt.Field(min_length=1, description="Alternative allele")
    af: float | None = pt.Field(ge=0.0, le=1.0, description="Allele frequency")
    ac: int | None = pt.Field(ge=0, description="Allele count")
    dp: int | None = pt.Field(ge=0, description="Read depth")
    mq: float | None = pt.Field(ge=0, description="Mapping quality")
    gene: str | None = pt.Field(description="Gene name")
    aa_effect: str | None = pt.Field(description="Amino acid effect")
    ref_codon_alt: str | None = pt.Field(description="Reference codon/alternative codon")
    cds_pos: int | None = pt.Field(gt=0, description="CDS position")
    aa_pos: int | None = pt.Field(gt=0, description="Amino acid position")


def validate_variant_data(df: pl.DataFrame) -> pl.DataFrame:
    """Validate variant pivot data using patito schema.

    Args:
        df: DataFrame with variant data

    Returns:
        Validated DataFrame

    Raises:
        ValueError: If validation fails
    """
    # Validate using patito
    try:
        VariantPivotSchema.validate(df)
    except Exception as e:
        msg = f"Variant data validation failed: {e!s}"
        raise ValueError(msg) from e

    return df


def main() -> None:
    """
    Generate variant pivot table from VCF data extracted with SnpSift.
    """
    args = parse_command_line_args()
    input_table = args.input_table

    # Read the data
    pivot_df = pl.read_csv(
        input_table,
        separator="\t",
        has_header=False,
        skip_rows=1,
        columns=[
            "contig",
            "ref",
            "pos",
            "alt",
            "af",
            "ac",
            "dp",
            "mq",
            "gene",
            "aa_effect",
            "ref_codon_alt",
            "cds_pos",
            "aa_pos",
        ],
    ).with_columns(pl.col("aa_effect").str.replace(".p", "").alias("aa_effect"))

    # Validate the data
    logger.info(f"Validating data from {input_table}")
    pivot_df = validate_variant_data(pivot_df)
    logger.success(f"Validation passed for {len(pivot_df)} variants")

    # TODO: Add actual pivot table generation logic here
    logger.info("Variant pivot table generation complete")


if __name__ == "__main__":
    main()
