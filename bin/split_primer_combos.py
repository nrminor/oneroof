#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "patito",
#     "polars",
# ]
# ///

"""
Short Python script for splitting a BED file with many possible primer
combinations into one BED per-combination.

```
usage: split_primer_combos.py [-h] --input_bed INPUT_BED [--forward_suffix FORWARD_SUFFIX] [--reverse_suffix REVERSE_SUFFIX]

options:
  -h, --help            show this help message and exit
  --input_bed INPUT_BED, -i INPUT_BED
                        BED file to be split per primer-combination.
  --forward_suffix FORWARD_SUFFIX, -f FORWARD_SUFFIX
                        Suffix in amplicon bed for forward primers.
  --reverse_suffix REVERSE_SUFFIX, -r REVERSE_SUFFIX
                        Suffix in amplicon bed for reverse primers.
```
"""

from __future__ import annotations

import argparse
from pathlib import Path

import patito as pt
import polars as pl


def parse_command_line_args() -> tuple[Path, str, str]:
    """
        Parse command line arguments while passing errors onto main.

    Args:
        `None`

    Returns:
        `Tuple[Path, str]`: a tuple containing the path to the input BED
        file, the forward suffix, and the reverse suffix for each primer name.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input_bed",
        "-i",
        type=Path,
        required=True,
        help="BED file to be split per primer-combination.",
    )
    parser.add_argument(
        "--forward_suffix",
        "-f",
        type=str,
        required=False,
        default="_LEFT",
        help="Suffix in amplicon bed for forward primers.",
    )
    parser.add_argument(
        "--reverse_suffix",
        "-r",
        type=str,
        required=False,
        default="_RIGHT",
        help="Suffix in amplicon bed for reverse primers.",
    )
    args = parser.parse_args()

    return args.input_bed, args.forward_suffix, args.reverse_suffix


class PrimerBedSchema(pt.Model):
    """Schema for validating primer BED format data."""

    Ref: str = pt.Field(description="Reference sequence name")
    Start_Position: int = pt.Field(ge=0, description="Start position (0-based)")
    Stop_Position: int = pt.Field(gt=0, description="Stop position (exclusive)")
    NAME: str = pt.Field(description="Primer name")
    INDEX: int = pt.Field(ge=0, description="Primer index")
    SENSE: str = pt.Field(description="Primer sense/strand")


def validate_primer_bed(df: pl.DataFrame) -> pl.DataFrame:
    """Validate primer BED format data using patito schema.

    Args:
        df: DataFrame with primer BED data

    Returns:
        Validated DataFrame

    Raises:
        ValueError: If validation fails
    """
    # Rename columns temporarily for validation
    df_for_validation = df.rename(
        {
            "Start Position": "Start_Position",
            "Stop Position": "Stop_Position",
        },
    )

    # Validate using patito
    try:
        PrimerBedSchema.validate(df_for_validation)
    except pt.exceptions.DataFrameValidationError as e:
        msg = f"Primer BED data validation failed: {e}"
        raise ValueError(msg) from e

    # Additional validation: ensure start < stop
    invalid_rows = df.filter(pl.col("Start Position") >= pl.col("Stop Position"))
    if len(invalid_rows) > 0:
        msg = f"Found {len(invalid_rows)} rows where start >= stop"
        raise ValueError(msg)

    # Validate SENSE values
    valid_senses = {"+", "-", "PLUS", "MINUS", "FORWARD", "REVERSE"}
    invalid_sense = df.filter(~pl.col("SENSE").is_in(valid_senses))
    if len(invalid_sense) > 0:
        unique_invalid = invalid_sense.select("SENSE").unique().to_series().to_list()
        msg = f"Invalid SENSE values found: {unique_invalid}. Expected one of: {valid_senses}"
        raise ValueError(msg)

    return df


def main() -> None:
    """
    Main handles data I/O and ultimately splits the BED file.
    """

    bed_to_split, fwd_suff, rev_suff = parse_command_line_args()

    # Read the BED file
    bed_df = pl.read_csv(
        bed_to_split,
        separator="\t",
        has_header=False,
        new_columns=[
            "Ref",
            "Start Position",
            "Stop Position",
            "NAME",
            "INDEX",
            "SENSE",
        ],
    )

    # Validate the data
    bed_df = validate_primer_bed(bed_df)

    bed_dfs = bed_df.with_columns(
        pl.col("NAME").str.replace(fwd_suff, "").str.replace(rev_suff, "").alias("NAME"),
    ).partition_by("NAME")

    for df in bed_dfs:
        splicing = df.select("NAME").unique().item()
        assert (
            len(df) == 2  # noqa: PLR2004
        ), f"Problematic splicing occurred with {splicing}"
        df.write_csv(
            file=f"{splicing}.bed",
            separator="\t",
            include_header=False,
        )


if __name__ == "__main__":
    main()
