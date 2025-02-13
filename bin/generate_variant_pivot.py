#!/usr/bin/env python3

# /// script
# requires-python = ">= 3.10"
# dependencies = [
#     "polars-lts-cpu",
#     "loguru",
# ]
# ///


import argparse
from pathlib import Path

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


def main() -> None:
    """
    TODO
    """
    args = parse_command_line_args()
    input_table = args.input_table

    _ = pl.read_csv(
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
    logger.info("Hi mom!")


if __name__ == "__main__":
    main()
