#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = [
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


def main() -> None:
    """
    Main handles data I/O and ultimately splits the BED file.
    """

    bed_to_split, fwd_suff, rev_suff = parse_command_line_args()

    bed_dfs = (
        pl.read_csv(
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
        .drop_nulls()
        .with_columns(
            pl.col("NAME")
            .str.replace(fwd_suff, "")
            .str.replace(rev_suff, "")
            .alias("NAME"),
        )
        .partition_by("NAME")
    )

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
