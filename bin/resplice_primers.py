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
from itertools import product
from pathlib import Path
from typing import Any

import polars as pl
from loguru import logger


def parse_command_line_args() -> argparse.Namespace:
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


def partition_by_amplicon(parsed_bed: pl.DataFrame) -> list[pl.DataFrame]:
    return (
        parsed_bed.with_columns(
            # make sure to save the original name for debugging purposes
            pl.col("NAME").alias("ORIG_NAME"),
            # this gets rid of any pre-existing spike-in primer indices, delimited by a hyphen,
            # that might need to be corrected, .e.g, having "amplicon01_left" and
            # "amplicon02-2_LEFT", but not "amplicon02-1_LEFT"
            pl.col("NAME").str.replace_all(r"-\d+", "").alias("NAME"),
        )
        .with_columns(
            pl.col("NAME")
            .str.replace_all("_LEFT", "")
            .str.replace_all("_RIGHT", "")
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
        .partition_by("Amplicon")
    )


def normalize_indices(
    partitioned_bed: list[pl.DataFrame],
    fwd_suffix: str,
    rev_suffix: str,
) -> dict[Any, pl.DataFrame]:
    normalized_dfs: list[pl.DataFrame] = []
    for primer_df in partitioned_bed:
        new_dfs = primer_df.partition_by("NAME")

        # TODO: Should there be a step here that appends and continues if the name
        # dataframe has only one row? Also, in the current setup where `partition_by_amplicon()`
        # does what it does before this function is run, can't we assume that any dataframes
        # with more than one row *are guaranteed* to have duplicates, in which case the dupe
        # check here is unnecessary (and verbose)?

        for j, name_df in enumerate(new_dfs):
            # run a check for whether there's >= 1 spikein. If there are none, there will just
            # be one primer per name partition, in which case we can skip to the next primer
            # set without renaming anything
            dupe_check = (
                name_df.with_columns(pl.col("NAME").is_duplicated().alias("duped"))
                .select("duped")
                .to_series()
                .to_list()
            )
            if True not in dupe_check:
                normalized_dfs.append(name_df)
                continue

            # otherwise, we'll need to rename the primers to account for spike-ins explicitly.
            corrected_indices = (
                # add a row that is the 1-based index of the primer and cast it as a string
                # to make concatenating columns more type-safe
                name_df.with_row_index("index", offset=1)
                .cast({"index": pl.Utf8})
                # overwrite the old NAME column with a new column that takes the core
                # amplicon label, the index of the primer, and the forward or reverse suffix
                .with_columns(
                    # use the polars.concat_str argument to concatenate the amplicon value,
                    # a hyphen, the row index column, and the forward or reverse suffix
                    # depending on what the original primer was.
                    pl.concat_str(
                        [
                            pl.col("Amplicon"),
                            pl.lit("-"),
                            pl.col("index"),
                            pl.when(pl.col("NAME").str.contains(fwd_suffix))
                            .then(pl.lit(fwd_suffix))
                            .otherwise(pl.lit(rev_suffix)),
                        ],
                        separator="",
                    ).alias(
                        "NAME",
                    ),
                )
                # reorder the relevant columns with a select expression.
                .select(
                    "Ref",
                    "Start Position",
                    "Stop Position",
                    "ORIG_NAME",
                    "NAME",
                    "INDEX",
                    "SENSE",
                    "Amplicon",
                ),
            )

            # overwrite the previous entry in this position of the dataframe list with the
            # updated entry, which has renamed primers to account for spike-ins.
            new_dfs[j] = corrected_indices[0]

        # after going through the forward and reverse primers for this amplicon in the above
        # control flow, append the updated dataframes to the accumulating normalized dataframes.
        normalized_dfs.append(pl.concat(new_dfs))

    return pl.concat(normalized_dfs).partition_by("Amplicon", as_dict=True)


def convertable_to_int(s: str) -> bool:
    try:
        int(s)
    except ValueError:
        return False
    except TypeError:
        return False
    else:
        return True


def resolve_primer_names(
    deficit_primers: list[str],
    excess_primers: list[str],
    fwd_suffix: str = "_LEFT",
    rev_suffix: str = "_RIGHT",
) -> tuple[list[str], list[str]]:
    # TODO: Is it possible that we no longer need to keep track of whether a primer is
    # in the deficit or the excess groups? This would spare us having to check yet again
    # for which primers are forward and reverse below.

    # find all possible combinations of the provided primers
    all_possible_pairs = list(product(deficit_primers, excess_primers))

    # initialize some mutable state to keep track of the new primer labels and their
    # pairing, which pairs have previously been handled, which resplice combo the current
    # iteration is handling, and the which old primer names should be in which order for
    # a join downstream
    new_primer_pairs: list[tuple[str, str]] = []
    handled_pairs: list[list[int]] = []
    old_primer_pairs: list[tuple[str, str]] = []
    combo_ticker = 0

    # loop through both primers, where primer1 is from the "deficit primer" list, and
    # primer2 is from the "excess primer" list.
    for primer1, primer2 in all_possible_pairs:
        # pull of the last element, delimited by hyphen, on both primer names
        primer1_final_element = primer1.split("-")[-1]
        primer2_final_element = primer2.split("-")[-1]

        # run checks to make sure indices that could be used for tracking each pair are
        # parseable from the primer name
        assert convertable_to_int(
            primer1_final_element,
        ), f"The primer {primer1} does not end with a hyphen-delimited integer, e.g. \
        '-1', which is required for properly handling different possible primer combinations. Aborting."
        assert convertable_to_int(
            primer2_final_element,
        ), f"The primer {primer2} does not end with a hyphen-delimited integer, e.g. \
        '-1', which is required for properly handling different possible primer combinations. Aborting."

        # figure out which combination of primers is being handled and check
        # whether it has already been handled
        primer1_index = int(primer1_final_element)
        primer2_index = int(primer2_final_element)
        current_pair = sorted((primer1_index, primer2_index))
        if current_pair in handled_pairs:
            continue

        handled_pairs.append(current_pair)

        # now that we know we're working with a previously unhandled pairing, incrememt
        # the combo ticker by one
        combo_ticker += 1

        # determine which of the primers are forward and which are reverse
        if fwd_suffix in primer1:
            old_fwd_primer = primer1

            # crash if the other primer doesn't contain a reverse suffix
            assert (
                rev_suffix in primer2
            ), f"Could not find the expected reverse suffix {rev_suffix} in the primer {primer2}. Aborting."
            old_rev_primer = primer2

        elif fwd_suffix in primer2:
            old_fwd_primer = primer2

            # crash if the other primer doesn't contain a reverse suffix
            assert (
                rev_suffix in primer1
            ), f"Could not find the expected reverse suffix {rev_suffix} in the primer {primer1}. Aborting."
            old_rev_primer = primer1

        else:
            # if the suffixes aren't found in either of the primers, something has gone very wrong.
            message = f"Neither {primer1} nor {primer2} contained the required \
            forward suffix {fwd_suffix} or the required reverse suffix {rev_suffix}. \
            Aborting."
            raise ValueError(message)

        # use f-strings to construct new names that make the combinations explicit
        new_fwd_primer = f"{old_fwd_primer}_splice{combo_ticker}"
        new_rev_primer = f"{old_rev_primer}_splice{combo_ticker}"

        # continue accumulating old and new primer pair lists
        old_primer_pairs.append((old_fwd_primer, old_rev_primer))
        new_primer_pairs.append((new_fwd_primer, new_rev_primer))

    # flatten the tuples at each position of the pair lists with two comprehensions
    # to make it explicit to the reader that forward primers come before reverse primers
    # in the flattened list. These comprehensions handle the old primer names.
    old_primer_names = [old_fwd for old_fwd, _ in old_primer_pairs] + [
        old_rev for _, old_rev in old_primer_pairs
    ]

    # do the same thing for the new primer names
    new_primer_names = [new_fwd for new_fwd, _ in new_primer_pairs] + [
        new_rev for _, new_rev in new_primer_pairs
    ]

    return (old_primer_names, new_primer_names)


def resplice_primers(
    amplicon_dfs: dict[Any, pl.DataFrame],
    fwd_suffix: str,
    rev_suffix: str,
) -> list[pl.DataFrame]:
    mutated_frames: list[pl.DataFrame] = []
    for amplicon, primer_df in amplicon_dfs.items():
        # our usual expectation is that there are two primers per amplicon. If this is
        # the case, append it to our output list and skip to the next amplicon's primers
        primer_pair_number = 2
        if primer_df.shape[0] == primer_pair_number:
            mutated_frames.append(primer_df)
            continue

        # If an amplicon only has one primer associated with it, something has gone wrong.
        if primer_df.shape[0] == 1:
            logger.error(
                f"There is a single primer without an amplicon running around! Here's \
                the parsed BED file for this amplicon:\n\n{primer_df}",
            )
            raise ValueError

        # pull out the primer labels and separate out forward and reverse primers
        primers = primer_df["NAME"]
        fwd_primers = [primer for primer in primers if fwd_suffix in primer]
        rev_primers = [primer for primer in primers if rev_suffix in primer]

        # determine which of the primer groupings, between the forward and reverse
        # primers, is the larger group.
        # TODO: Restating the above, this may no longer be necessary because of the
        # use of the product function. Instead, we could just replace it for an equality
        # check for the lengths of the two primer lists. If the check is true, append
        # the current primer_df, which won't need to go through the combinatorics, and
        # continue, i.e.,:
        # ```
        # if len(fwd_primers) == len(rev_primers):
        #   mutated_frames.append(primer_df)  # noqa: ERA001
        #   continue  # noqa: ERA001
        # ````
        if len(fwd_primers) > len(rev_primers):
            deficit_primers = rev_primers
            excess_primers = fwd_primers
        elif len(fwd_primers) < len(rev_primers):
            deficit_primers = fwd_primers
            excess_primers = rev_primers
        else:
            assert (
                (len(fwd_primers) + len(rev_primers)) % 2 == 0
            ), f"Invalid primer partitioning encountered with the forward primers \
            {fwd_primers} and the reverse primers {rev_primers}. Aborting."
            mutated_frames.append(primer_df)
            continue

        # compute the new names for the primers that will be used to expand shared primers
        # across all spike-ins
        old_primer_names, new_primer_names = resolve_primer_names(
            deficit_primers,
            excess_primers,
            fwd_suffix,
            rev_suffix,
        )

        # double check that we have the correct number of primer labels, though with the
        # current design this should be impossible
        assert len(old_primer_names) == len(
            new_primer_names,
        ), f"Insufficient number of replacement names ({new_primer_names}) \
        generated for partition for amplicon {amplicon}: {primer_df}"

        # run a join on the old primer names to bring in the new primer names in their
        # proper locations
        new_df = (
            pl.DataFrame({"NAME": old_primer_names})
            .cast(pl.Utf8)
            .join(
                primer_df.with_columns(pl.col("NAME").cast(pl.Utf8)),
                how="left",
                on="NAME",
                validate="1:1",
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

    bed_df = pl.read_csv(
        args.input_bed,
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

    # identify amplicons/primer pairings and create a separate dataframe for the primers
    # associated with each
    partitioned_bed = partition_by_amplicon(bed_df)

    # make sure that all primers have informative names in amplicons where there are
    # spiked-in primers
    indexed_primer_dfs = normalize_indices(
        partitioned_bed,
        args.fwd_suffix,
        args.rev_suffix,
    )

    respliced_dfs = resplice_primers(
        indexed_primer_dfs,
        args.fwd_suffix,
        args.rev_suffix,
    )

    if len(respliced_dfs) == len(indexed_primer_dfs):
        logger.warning(
            f"The number of amplicon primer sets that made it through resplicing,\
            {len(respliced_dfs)}, does not match the number of input amplicons, \
            {len(indexed_primer_dfs)}. Data loss may have occurred for this primer set.",
        )

    # TODO: Is this function necessary? It seems to be controlling for something that
    # should no longer be possible, right?
    final_df = finalize_primer_pairings(
        respliced_dfs,
        args.fwd_suffix,
        args.rev_suffix,
    )

    # write the now unnecessary columns, sort by start and stop position, and write
    # out with the user-provided output prefix
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
