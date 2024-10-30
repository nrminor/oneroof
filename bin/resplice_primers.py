#!/usr/bin/env python3

"""
Module that finds all possible combinations of spike-in primers in an amplicon scheme.

Example usage:
```
usage: resplice_primers.py [-h] --input_bed INPUT_BED [--output_prefix OUTPUT_PREFIX] [--config CONFIG]

options:
    -h, --help            show this help message and exit
    --input_bed INPUT_BED, -i INPUT_BED
                        BED file with one-off spike-in primers to be respliced into possible amplicons.
    --output_prefix OUTPUT_PREFIX, -o OUTPUT_PREFIX
                        Output prefix for final respliced amplicon BED file.
    --config CONFIG, -c CONFIG
                        YAML file used to configure module such that it avoids hardcoding
```
"""

from __future__ import annotations

import argparse
import shutil
from itertools import product
from pathlib import Path

import polars as pl
from loguru import logger


def parse_command_line_args() -> argparse.Namespace:
    """
        Parse command line arguments while passing errors onto main.

    Args:
        `None`

    Returns:
        `tuple[Path, str]`: a tuple containing the path to the input BED
        file and a string representing the desired output name.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input_bed",
        "-i",
        type=Path,
        required=True,
        help="BED file with one-off spike-in primers to be respliced into possible amplicons.",
    )
    parser.add_argument(
        "--output_prefix",
        "-o",
        type=str,
        required=False,
        default="respliced",
        help="Output prefix for final respliced amplicon BED file.",
    )
    parser.add_argument(
        "--fwd_suffix",
        "-f",
        type=str,
        required=False,
        default="_LEFT",
        help="The suffix to be expected in the names for forward primers",
    )
    parser.add_argument(
        "--rev_suffix",
        "-r",
        type=str,
        required=False,
        default="_RIGHT",
        help="The suffix to be expected in the names for reverse primers",
    )
    return parser.parse_args()


def dedup_primers(partitioned_bed: list[pl.DataFrame]) -> list[pl.DataFrame]:
    """
        `dedup_primers()` primers finds repeated instances of the same primer
        name and adds a unique identifier for each. This ensures that joins
        downstream are one-to-many as opposed to many-to-many.

    Args:
        `partitioned_bed: list[pl.DataFrame]`: A list of Polars DataFrames that
        have been partitioned by amplicon.

    Returns:
        `list[pl.DataFrame]`: A partitioned list of Polars dataframes with no
        repeat primer names.
    """

    for i, primer_df in enumerate(partitioned_bed):
        if True in primer_df.select("NAME").is_duplicated().to_list():
            new_dfs = primer_df.with_columns(
                primer_df.select("NAME").is_duplicated().alias("duped"),
            ).partition_by("duped")

            for j, dup_frame in enumerate(new_dfs):
                if True in dup_frame.select("duped").to_series().to_list():
                    renamed = (
                        dup_frame.with_row_count(offset=1)
                        .cast({"row_nr": pl.Utf8})
                        .with_columns(
                            pl.concat_str(
                                [pl.col("NAME"), pl.col("row_nr")],
                                separator="-",
                            ).alias("NAME"),
                        )
                        .select(
                            "Ref",
                            "Start Position",
                            "Stop Position",
                            "ORIG_NAME",
                            "NAME",
                            "INDEX",
                            "SENSE",
                            "Amplicon",
                            "duped",
                        )
                    )
                    new_dfs[j] = renamed
            partitioned_bed[i] = pl.concat(new_dfs)

    return pl.concat(partitioned_bed).drop("duped").partition_by("Amplicon")


def resolve_primer_names(
    to_combine: list[str],
    combine_to: list[str],
) -> tuple[list[str], list[str]]:
    """
        `resolve_primer_names()` names each possible pairing of primers in
        amplicons where singletons, forward or reverse, have been added to
        increase template coverage.

    Args:
        `to_combine: list[str]`: A list of forward primers to resolve.
        `combine_to: list[str]`: A list of reverse primers to resolve.

    Returns:
        `tuple[list[str], list[str]]`: A tuple containing two lists, the first
        being a list of primer names to use with joining, and the second being
        a list of new primer names to use once left-joining is complete.
    """

    primer_pairs = list(product(to_combine, combine_to))

    new_primer_pairs = []
    for fwd_primer, rev_primer in primer_pairs:
        fwd_suffix = fwd_primer.split("_")[-1]
        rev_suffix = rev_primer.split("_")[-1]
        primer_label = fwd_primer.replace(f"_{fwd_suffix}", "").split("-")[-1]
        try:
            int(primer_label)
        except TypeError:
            primer_label = "1"
        amplicon = "_".join(
            fwd_primer.replace(f"_{fwd_suffix}", "")
            .replace(f"-{primer_label}", "")
            .split("_")[0:2],
        )
        new_fwd_primer = f"{amplicon}_splice{primer_label}_{fwd_suffix}"
        new_rev_primer = f"{amplicon}_splice{primer_label}_{rev_suffix}"
        new_primer_pairs.append([new_fwd_primer, new_rev_primer])

    primers_to_join = [item[0] for item in primer_pairs] + [
        item[1] for item in primer_pairs
    ]

    new_primer_names = [item[0] for item in new_primer_pairs] + [
        item[1] for item in new_primer_pairs
    ]

    return (primers_to_join, new_primer_names)


@logger.catch
def resplice_primers(dedup_partitioned: list[pl.DataFrame]) -> list[pl.DataFrame]:
    """
        `resplice_primers()` determines whether spike-ins are forward or reverse
        primers (or both) and uses that information to handle resplicing
        possible combinations.

    Args:
        `dedup_partitioned: list[pl.DataFrame]`: A Polars dataframe with no
        duplicate primer names.

    Returns:
        `list[pl.DataFrame]`: A list of Polars DataFrames where each dataframe
        is represents all possible pairings of primers within a single amplicon.
    """

    mutated_frames: list[pl.DataFrame] = []
    for i, dedup_df in enumerate(dedup_partitioned):
        primer_pair_number = 2
        if dedup_df.shape[0] != primer_pair_number:
            primers = dedup_df["NAME"]

            fwd_primers = [primer for primer in primers if "LEFT" in primer]
            rev_primers = [primer for primer in primers if "RIGHT" in primer]

            if len(fwd_primers) == 0 and len(rev_primers) == 0:
                return mutated_frames

            if len(fwd_primers) > len(rev_primers):
                to_combine = rev_primers
                combine_to = fwd_primers
            elif len(fwd_primers) < len(rev_primers):
                to_combine = fwd_primers
                combine_to = rev_primers
            else:
                break

            primers_to_join, new_primer_names = resolve_primer_names(
                to_combine,
                combine_to,
            )

            assert len(primers_to_join) == len(
                new_primer_names,
            ), f"Insufficient number of replacement names generated for partition {i}"
            new_df = (
                pl.DataFrame({"NAME": primers_to_join})
                .cast(pl.Utf8)
                .join(
                    dedup_df.with_columns(pl.col("NAME").cast(pl.Utf8)),
                    how="left",
                    on="NAME",
                    validate="m:1",
                    coalesce=False,
                )
                .with_columns(pl.Series(new_primer_names).alias("NAME"))
                .select(
                    "Ref",
                    "Start Position",
                    "Stop Position",
                    "ORIG_NAME",
                    "NAME",
                    "INDEX",
                    "SENSE",
                    "Amplicon",
                )
            )
            mutated_frames.append(new_df)
        elif dedup_df.shape[0] == 1:
            logger.error("There is a single primer without an amplicon running around!")
            raise ValueError
        else:
            mutated_frames.append(dedup_df)

    return mutated_frames


def finalize_primer_pairings(
    mutated_frames: list[pl.DataFrame],
    fwd_suffix: str,
    rev_suffix: str,
) -> pl.DataFrame:
    """
        `finalize_primer_pairings()` removes any spikeins with possible pairings
        that could not be determined.

    Args:
        `mutated_frames: list[pl.DataFrame]`: A list of Polars DataFrames, each
        representing a respliced amplicon.

    Returns:
        `pl.DataFrame`: A concatenated Polars dataframe that will be written
        out to a new BED file.
    """

    final_frames: list[pl.DataFrame] = []
    for df in mutated_frames:
        fwd_keepers = [
            primer
            for primer in df.select("NAME").to_series().to_list()
            if fwd_suffix in primer
        ]
        rev_keepers = [
            primer
            for primer in df.select("NAME").to_series().to_list()
            if rev_suffix in primer
        ]
        if len(fwd_keepers) > 0 and len(rev_keepers) > 0:
            final_frames.append(df)

    return pl.concat(final_frames)


def main() -> None:
    """
    `main()` coordinates the flow of data through the module's functions.
    """

    args = parse_command_line_args()

    partitioned_bed = (
        pl.read_csv(
            args.bed_file,
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
        .with_columns(pl.col("NAME").alias("ORIG_NAME"))
        .with_columns(
            pl.col("NAME")
            .str.replace_all(args.fwd_prefix, "")
            .str.replace_all(args.rev_suffix, "")
            .str.replace_all(r"-\d+", "")
            .alias("Amplicon"),
        )
        .select(
            "Ref",
            "Start Position",
            "Stop Position",
            "ORIG_NAME",
            "NAME",
            "INDEX",
            "SENSE",
            "Amplicon",
        )
        .with_columns(pl.col("NAME").is_duplicated().alias("duped"))
        .partition_by("Amplicon")
    )

    dedup_partitioned = dedup_primers(partitioned_bed)

    mutated_frames = resplice_primers(dedup_partitioned)

    if len(mutated_frames) == 0:
        shutil.copy(args.bed_file, f"{args.output_prefix}.bed")
        return

    final_df = finalize_primer_pairings(mutated_frames)

    final_df.drop("Amplicon").drop("NAME").sort(
        "Start Position",
        "Stop Position",
    ).write_csv(
        f"{args.output_prefix}.bed",
        separator="\t",
        include_header=False,
    )


if __name__ == "__main__":
    main()
