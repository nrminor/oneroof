#!/usr/bin/env python3
# /// script
# requires-python = ">= 3.10"

# dependencies = [
#     "polars-lts-cpu",
#     "loguru",
# ]
# ///

"""
Module that finds all possible combinations of spike-in primers in an amplicon scheme.

Example usage:
```
usage: resplice_primers.py [-h] --input_bed INPUT_BED [--output_prefix OUTPUT_PREFIX]
                           [--fwd_suffix FWD_SUFFIX] [--rev_suffix REV_SUFFIX]
                           [--idx_delim IDX_DELIM] [--idx_position IDX_POSITION] [-v]

options:
  -h, --help            show this help message and exit
  --input_bed INPUT_BED, -i INPUT_BED
                        BED file with one-off spike-in primers to be respliced into
                        possible amplicons.
  --output_prefix OUTPUT_PREFIX, -o OUTPUT_PREFIX
                        Output prefix for final respliced amplicon BED file.
  --fwd_suffix FWD_SUFFIX, -f FWD_SUFFIX
                        The suffix to be expected in the names for forward primers
  --rev_suffix REV_SUFFIX, -r REV_SUFFIX
                        The suffix to be expected in the names for reverse primers
  --idx_delim IDX_DELIM, -d IDX_DELIM
                        The symbol used to delimit the index of a spike-in primer,
                        which differentiates it from the primer original around the
                        same position. Defaults to a hyphen/dash: '-'.
  --idx_position IDX_POSITION, -p IDX_POSITION
                        The position where the primer spike-in index should be
                        expected, defaulting to the final position after splitting by
                        the specified index delimiter symbol.
  -v, --verbose         Increase verbosity level (-v for WARNING, -vv for INFO, -vvv
                        for DEBUG)
```
"""

from __future__ import annotations

import argparse
import sys
from itertools import product
from pathlib import Path
from pprint import pformat
from typing import Any

import polars as pl
from loguru import logger


def parse_command_line_args() -> argparse.Namespace:
    """
    Creates an argument parser object to handle command line arguments.

    Returns:
        argparse.Namespace: The parsed command line arguments.

    The available arguments are:
        --input_bed/-i: Path to the BED file containing spike-in primers.
        --output_prefix/-o: Output prefix for the respliced amplicon BED file.
        --fwd_suffix/-f: The suffix expected in forward primer names.
        --rev_suffix/-r: The suffix expected in reverse primer names.
        --idx_delim/-d: The symbol delimiting spike-in primer indices.
        --idx_position/-p: Position to expect primer spike-in index after splitting by delimiter.
        --verbose/-v: Increase verbosity level (-v for WARNING, -vv for INFO, -vvv for DEBUG)
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
    parser.add_argument(
        "--idx_delim",
        "-d",
        type=str,
        required=False,
        default="-",
        help="The symbol used to delimit the index of a spike-in primer, which differentiates it from the primer original around the same position. Defaults to a hyphen/dash: '-'.",
    )
    parser.add_argument(
        "--idx_position",
        "-p",
        type=int,
        required=False,
        default=-1,
        help="The position where the primer spike-in index should be expected, defaulting to the final position after splitting by the specified index delimiter symbol.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="count",
        default=0,
        help="Increase verbosity level (-v for WARNING, -vv for INFO, -vvv for DEBUG)",
    )

    return parser.parse_args()


def check_bed_existence(provided_bed_path: str | Path) -> None:
    """
    Checks for the existence of a BED file at the provided path.

    Args:
        provided_bed_path (str | Path): Path to the BED file to verify

    Returns:
        None

    Raises:
        SystemExit: If BED file does not exist or is not a file. Lists available .bed
        files in the same directory if any are found.
    """
    logger.info(f"Parsing input BED file at {provided_bed_path}...")
    if not Path(provided_bed_path).is_file():
        available_bed_files = [
            str(file) for file in Path(provided_bed_path).parent.glob("*.bed")
        ]
        if len(available_bed_files) == 0:
            logger.error(
                f"The provided input BED file, {provided_bed_path}, does not exist or is not a file. Aborting.",
            )
            sys.exit(1)
        logger.error(
            f"The provided input BED file, {provided_bed_path}, does not exist or is not a file. Perhaps you meant one of the following files?\n\n{pformat(available_bed_files, indent=1)}",
        )
        sys.exit(1)


def check_idx_delims(
    parsed_bed: pl.DataFrame,
    idx_delim: str = "-",
    _idx_position: int = -1,
) -> None:
    """
    Checks primer names for multiple occurrences of the index delimiter symbol.

    This function examines all primer names in the BED file to identify any that contain
    multiple instances of the specified index delimiter symbol. Having multiple
    delimiters could cause issues with proper primer indexing and matching.

    Args:
        parsed_bed (pl.DataFrame): A polars DataFrame containing the parsed BED file
        data
        idx_delim (str, optional): The symbol used to delimit spike-in primer indices.
            Defaults to "-".
        _idx_position (int, optional): Position to expect the index after splitting by
            delimiter. Defaults to -1.

    Returns:
        None

    If any primers are found with multiple delimiter symbols, a warning message is
    logged listing the problematic primer names.
    """
    primer_names = parsed_bed.select("NAME").to_series().to_list()

    warning_list = []
    for name in primer_names:
        delim_count = name.count(idx_delim)

        if delim_count <= 1:
            continue

        warning_list.append(name)

    if len(warning_list) == 0:
        return

    logger.warning(
        f"{len(warning_list)} primer names contained more than one of the symbol used to delimit spike-in primer index, '{idx_delim}', which is a special value in this program. Unexpected behavior, such as merging multiple amplicons into one or incorrect resplicing, is likely to occur.\n\nTo fix it, either use a different delimiter symbol than '{idx_delim}' or modify the primer names so that the symbol '{idx_delim}' is only used to denote a spike-in index.\n\nBelow is the list of {len(warning_list)} primers that generated this warning:\n\n{pformat(warning_list, indent=1)}",
    )


def partition_by_amplicon(
    parsed_bed: pl.DataFrame,
    fwd_suffix: str,
    rev_suffix: str,
) -> list[pl.DataFrame]:
    """
    Partitions a BED file of PCR primer coordinates by amplicon name.

    Args:
        parsed_bed (pl.DataFrame): A polars DataFrame containing parsed primer coordinates
            from a BED file. Expected columns are Ref, Start Position, Stop Position, NAME,
            INDEX and SENSE.
        fwd_suffix (str): The suffix used to identify forward PCR primers.
        rev_suffix (str): The suffix used to identify reverse PCR primers.

    Returns:
        list[pl.DataFrame]: A list of Polars DataFrames, each containing primers from a
            shared amplicon, identified by the amplicon name derived from stripping
            suffixes from the primer NAME column.
    """
    return (
        parsed_bed.with_columns(
            # make sure to save the original name for debugging purposes
            pl.col("NAME").alias("ORIG_NAME"),
            # this gets rid of any pre-existing spike-in primer indices, delimited by a
            # hyphen, that might need to be corrected, .e.g, having "amplicon01_left"
            # and "amplicon02-2_LEFT", but not "amplicon02-1_LEFT"
            pl.col("NAME").str.replace_all(r"-\d+", "").alias("NAME"),
        )
        .with_columns(
            # take the name column, which now has no hyphen-delimited indices, and also
            # remove the forward and reverse suffixes, calling the resulting column
            # 'Amplicon'. This should result in an amplicon column that will group
            # together associated forward and reverse primers, as well as primer name
            # groupings with >= 1 forward or reverse primer. If the forward and reverse
            # groupings have more than one entry, they should all now have the same name
            # and thus be ready for indexing downstream.
            pl.col("NAME")
            .str.replace_all(fwd_suffix, "")
            .str.replace_all(rev_suffix, "")
            .alias("Amplicon"),
        )
        # use a select expression to reorder columns
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


def assign_new_indices(
    primer_name_df: pl.DataFrame,
    fwd_suffix: str,
    rev_suffix: str,
    idx_delim: str = "-",
    _idx_position: int = -1,
) -> pl.DataFrame:
    """
    Assigns new indices to primers within an amplicon scheme when spike-in primers are present.

    Given a DataFrame of primers within a single amplicon that may contain spike-ins,
    assigns a new index to each primer that could potentially pair with another primer
    from the same amplicon. This indexing scheme ensures that spike-in primers can be
    differentiated from their original primers by adding a numerical suffix.

    Args:
        primer_name_df (pl.DataFrame): Dataframe containing all primers belonging to a
            single amplicon.
        fwd_suffix (str): The suffix used to identify forward PCR primers.
        rev_suffix (str): The suffix used to identify reverse PCR primers.
        idx_delim  (str): The character used to delimit the index, defaulting to "-".
        _idx_position (int): The as-yet unused expected position within a split primer
            name where the index is expected. Defaults to the final position.s

    Returns:
        pl.DataFrame: DataFrame with a new NAME column that uniquely identifies each
        primer with a numerical suffix based on its position in the amplicon scheme.

    The function assigns indices to primers by:
    - Adding a 1-based row index
    - Creating a new NAME column by concatenating:
        - Original amplicon name
        - Forward/reverse suffix
        - Hyphen delimiter by default, but the user can provide their own delimiter
        - Row index
    """
    return (
        # add a row that is the 1-based index of the primer and cast it as a string
        # to make concatenating columns more type-safe
        primer_name_df.with_row_index("index", offset=1)
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
                    pl.when(pl.col("NAME").str.contains(fwd_suffix))
                    .then(pl.lit(fwd_suffix))
                    .otherwise(pl.lit(rev_suffix)),
                    pl.lit(idx_delim),
                    pl.col("index"),
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
        )
    )


def normalize_indices(
    partitioned_bed: list[pl.DataFrame],
    fwd_suffix: str,
    rev_suffix: str,
) -> dict[Any, pl.DataFrame]:
    """
    Normalizes indices within a group of primers, assigning new indices to distinguish spike-in primers.

    This function performs several key steps:
    1. Partitions primers by amplicon and checks for presence of spike-ins
    2. For each amplicon group with spike-ins, assigns sequential indices to differentiate primers
    3. Handles special cases like single primers or duplicates
    4. Returns a dictionary mapping amplicon names to their primer DataFrames

    The control flow involves nested loops:
    - Outer loop processes each amplicon's primers
    - Inner loop handles forward/reverse primer groups within each amplicon
    - For each group, checks are performed to ensure valid primer naming
    - New indices are assigned only when spike-ins are detected

    Args:
        partitioned_bed (list[pl.DataFrame]): List of DataFrames, each containing primers from
            a single amplicon
        fwd_suffix (str): Suffix identifying forward primers
        rev_suffix (str): Suffix identifying reverse primers

    Returns:
        dict[Any, pl.DataFrame]: Dictionary mapping amplicon names to DataFrames containing
        their normalized primers. Each primer now has a unique index if it is part of a
        spike-in group.

    The function maintains separation between original and spike-in primers through careful
    index management, ensuring that downstream analyses can properly pair complementary
    primers while preserving information about which primers are spike-ins versus originals.
    """
    # initialize a new dataframe to be mutated in the following for loops
    normalized_dfs: list[pl.DataFrame] = []

    # The outer loop goes through each amplicon, which should have at least two primers,
    # a forward primer and a reverse primer.
    for primer_df in partitioned_bed:
        # partition the dataframe for the current amplicon by name, which should
        # effectively mean partitioning by forward or reverse.
        new_dfs: list[pl.DataFrame] = primer_df.partition_by("NAME")

        # the inner loop iterates through each group of potentially redundant primer
        # names, assigning new indices to differentiate them if they are. For example,
        # if there are three forward primers (i.e., one original primer and two spike-in
        # primers), they should all have the same name at this stage, and will thus have
        # -1, -2, and -3 appended to them respectively.
        for j, primer_name_df in enumerate(new_dfs):

            # having passed the above checks, we can now proceed to assigning new names
            # to the primers to account for spike-ins explicitly.
            corrected_indices = assign_new_indices(
                primer_name_df,
                fwd_suffix,
                rev_suffix,
            )

            # overwrite the previous entry in this position of the dataframe list with
            # the updated entry, which has renamed primers to account for spike-ins.
            new_dfs[j] = corrected_indices

            # end of inner NAME loop

        # after going through the forward and reverse primers for this amplicon in the
        # above control flow, append the updated dataframes to the accumulating
        # normalized dataframes.
        normalized_dfs.append(pl.concat(new_dfs))

        # end of outer Amplicon loop

    return pl.concat(normalized_dfs).partition_by("Amplicon", as_dict=True)


def _convertable_to_int(s: str) -> bool:
    try:
        int(s)
    except ValueError:
        return False
    except TypeError:
        return False
    else:
        return True


def resolve_primer_names(
    old_fwd_primers: list[str],
    old_rev_primers: list[str],
    idx_delim: str = "-",
    idx_position: int = -1,
) -> tuple[list[str], list[str]]:
    """
    Resolves primer names by generating new labels for all possible combinations of forward and reverse primers.

    This function takes lists of forward and reverse primer names and generates unique identifiers for each
    possible pairing combination. It handles primers that have been marked with indices to denote spike-ins
    or variants.

    Args:
        old_fwd_primers (list[str]): List of forward primer names to be paired
        old_rev_primers (list[str]): List of reverse primer names to be paired
        idx_delim (str, optional): Delimiter used to separate the index in primer names. Defaults to "-"
        idx_position (int, optional): Position of the index when splitting by delimiter. Defaults to -1

    Returns:
        tuple[list[str], list[str]]: Two lists - one containing the original primer names in order, and
            one containing the newly generated primer names with unique combination identifiers

    The function:
    1. Finds all possible forward/reverse primer combinations
    2. Assigns a unique splice index to each valid pairing
    3. Handles primers with existing hyphen-delimited indices
    4. Returns both old and new primer names for mapping the changes
    """
    # find all possible combinations of the provided primers
    all_possible_pairs = list(product(old_fwd_primers, old_rev_primers))

    # initialize some mutable state to keep track of the new primer labels and their
    # pairing, which pairs have previously been handled, which resplice combo the
    # current iteration is handling, and which old primer names should be in which
    # order for a join downstream.
    new_primer_pairs: list[tuple[str, str]] = []
    handled_pairs: list[list[int]] = []
    old_primer_pairs: list[tuple[str, str]] = []
    combo_ticker = 0

    # loop through both primers, where primer1 is from the "deficit primer" list, and
    # primer2 is from the "excess primer" list.
    for old_fwd_primer, old_rev_primer in all_possible_pairs:
        # pull of the last element, delimited by hyphen, on both primer names
        fwd_final_element = old_fwd_primer.split(idx_delim)[idx_position]
        rev_final_element = old_rev_primer.split(idx_delim)[idx_position]

        # run checks to make sure indices that could be used for tracking each pair are
        # parseable from the primer name
        if not _convertable_to_int(fwd_final_element):
            logger.error(
                f"The primer {old_fwd_primer}, which has been paired with {old_rev_primer} does not end with a hyphen-delimited integer, e.g. '-1', which is required for properly handling different possible primer combinations. Printing the set of combinations being handled and aborting:\n\n{all_possible_pairs}",
            )
            sys.exit(1)
        if not _convertable_to_int(rev_final_element):
            logger.error(
                f"The primer {old_rev_primer}, which has been paired with {old_fwd_primer} does not end with a hyphen-delimited integer, e.g. '-1', which is required for properly handling different possible primer combinations. Printing the set of combinations being handled and aborting:\n\n{all_possible_pairs}",
            )
            sys.exit(1)

        # figure out which combination of primers is being handled and check
        # whether it has already been handled
        primer1_index = int(fwd_final_element)
        primer2_index = int(rev_final_element)
        current_pair = sorted((primer1_index, primer2_index))
        if current_pair in handled_pairs:
            continue

        # now that we know we're working with a previously unhandled pairing, incrememt
        # the combo ticker by one
        combo_ticker += 1

        # use f-strings to construct new names that make the combinations explicit
        new_fwd_primer = f"{old_fwd_primer}_splice{combo_ticker}"
        new_rev_primer = f"{old_rev_primer}_splice{combo_ticker}"

        # continue accumulating old and new primer pair lists
        old_primer_pairs.append((old_fwd_primer, old_rev_primer))
        new_primer_pairs.append((new_fwd_primer, new_rev_primer))

        # now that we know nothing has gone awry, at this pair to the accumulating list
        # of handled primer pairs
        handled_pairs.append(current_pair)

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
    idx_delim: str = "-",
    idx_position: int = -1,
) -> list[pl.DataFrame]:
    """
    Generates all possible combinations of forward and reverse primers within each amplicon.

    This function takes amplicon DataFrames, identifies spike-in primers, and creates new
    primer combinations by pairing forward and reverse primers. It handles special cases like
    standard primer pairs and single primers.

    Args:
        amplicon_dfs (dict[Any, pl.DataFrame]): Dictionary mapping amplicon names to their
            primer DataFrames
        fwd_suffix (str): Suffix identifying forward primers
        rev_suffix (str): Suffix identifying reverse primers
        idx_delim (str, optional): Delimiter used for spike-in indices. Defaults to "-"
        idx_position (int, optional): Position of spike-in index after splitting. Defaults to -1

    Returns:
        list[pl.DataFrame]: List of DataFrames containing valid primer combinations for each
            amplicon, with updated names reflecting the possible pairings.

    The function:
    1. Checks for standard two-primer amplicons
    2. Identifies forward and reverse primers in amplicons with spike-ins
    3. Generates valid primer combinations and assigns new names
    4. Skips invalid or incomplete primer sets
    5. Returns DataFrames ready for final processing
    """
    mutated_frames: list[pl.DataFrame] = []
    for amplicon, primer_df in amplicon_dfs.items():
        # our usual expectation is that there are two primers per amplicon. If this is
        # the case, append it to our output list and skip to the next amplicon's primers
        if primer_df.shape[0] == 2:  # noqa: PLR2004
            # double check that there's an identifiable forward and reverse primer
            # within the pair
            pair_labels = primer_df.select("NAME").to_series().to_list()
            logger.debug(
                f"Pair of primers within the amplicon {amplicon} detected: {pair_labels}. No resplicing will be needed here, though double check that the necessary forward and reverse suffixes, {fwd_suffix} and {rev_suffix}, are present.",
            )
            assert any(
                fwd_suffix in primer for primer in pair_labels
            ), f"The forward suffix {fwd_suffix} is missing in the provided primer pairs: {pair_labels}. Aborting."
            assert any(
                rev_suffix in primer for primer in pair_labels
            ), f"The reverse suffix {rev_suffix} is missing in the provided primer pairs: {pair_labels}. Aborting."

            # if so, all is well. Append the primers into the growing list of correct
            # dataframes and move onto the next amplicon
            mutated_frames.append(primer_df)
            continue

        # If an amplicon only has one primer associated with it, something has gone
        # wrong.
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

        # compute the new names for the primers that will be used to expand shared
        # primers across all spike-ins
        logger.debug(
            f"Resplicing will be performed for the following primers:\n\n{pformat(fwd_primers + rev_primers, indent=1)}",
        )
        old_primer_names, new_primer_names = resolve_primer_names(
            fwd_primers,
            rev_primers,
            idx_delim,
            idx_position,
        )

        # double check that we have the correct number of primer labels, though with the
        # current design this should be impossible
        assert len(old_primer_names) == len(
            new_primer_names,
        ), f"Insufficient number of replacement names ({new_primer_names}) generated \
        for partition for amplicon {amplicon}: {primer_df}"

        # run a join on the old primer names to bring in the new primer names in their
        # proper locations
        new_df = (
            pl.DataFrame({"NAME": old_primer_names})
            .cast(pl.Utf8)
            .join(
                primer_df.with_columns(pl.col("NAME").cast(pl.Utf8)),
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
    This script generates all possible combinations of amplicons based on an input BED
    file of spike-in primers and original primers.

    The workflow consists of the following steps:
    1. Parse command line arguments specifying input BED, output prefix, and primer
       naming conventions
    2. Read in BED file of primers as a polars DataFrame
    3. Group primers by amplicon and partition them into separate dataframes
    4. Normalize primer indices within each amplicon group to handle spike-ins
    5. Generate all valid primer pair combinations for each amplicon
    6. Filter out invalid pairings and finalize primer names
    7. Output a new BED file containing all valid primer pair combinations

    Each primer is required to have forward/reverse suffixes and a spike-in index
    delimited by a specified character. The script checks that all primers follow this
    naming convention to ensure reliable pairing.

    The output maintains the original BED file format while expanding out all possible
    amplicon combinations from the given spike-in primers.
    """
    args = parse_command_line_args()

    # set up logger based on the requested verbosity
    match args.verbose:
        case 0:
            level = "WARNING"
        case 1:
            level = "SUCCESS"
        case 2:
            level = "INFO"
        case 3:
            level = "DEBUG"
        case _:
            level = "WARNING"
    logger.remove()
    logger.add(sys.stderr, colorize=True, level=level)
    logger.info(f"Logging set to {level}")

    # read in the BED file
    check_bed_existence(args.input_bed)
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

    # look through primer names before proceeding to make sure primers aren't named
    # inconsistently with respect the symbols they contain
    check_idx_delims(
        bed_df,
        args.idx_delim,
        args.idx_position,
    )

    # identify amplicons/primer pairings and create a separate dataframe for the primers
    # associated with each
    partitioned_bed = partition_by_amplicon(
        bed_df,
        args.fwd_suffix,
        args.rev_suffix,
    )
    logger.info(
        f"{len(bed_df)} input primers split into {len(partitioned_bed)} discrete amplicons.",
    )
    if len(bed_df) % len(partitioned_bed) != 0:
        logger.info(
            f"There appear to be {len(bed_df) % len(partitioned_bed)} spike-in primers.",
        )
    else:
        logger.info(
            "There don't appear to be any spike-in primers. Assuming that amplicons can be parsed properly, this means that no resplicing will occur, though primer indices will still be updated.",
        )

    # make sure that all primers have informative names in amplicons where there are
    # spiked-in primers
    logger.info(
        f"Re-assigning 1-based indices to account for current or future spike-ins amongst the {len(bed_df)} input primers.",
    )
    indexed_primer_dfs = normalize_indices(
        partitioned_bed,
        args.fwd_suffix,
        args.rev_suffix,
    )

    logger.info(
        "Re-indexing successful. Proceeding to resplicing any implicit combinations that are not accounted for with a row in the input BED file.",
    )
    respliced_dfs = resplice_primers(
        indexed_primer_dfs,
        args.fwd_suffix,
        args.rev_suffix,
        args.idx_delim,
        args.idx_position,
    )

    if len(respliced_dfs) != len(indexed_primer_dfs):
        logger.warning(
            f"The number of amplicon primer sets that made it through resplicing, {len(respliced_dfs)}, does not match the number of input amplicons, {len(indexed_primer_dfs)}. Data loss may have occurred for this primer set.",
        )

    # run one last check that each group has at least one identifiable forward primer
    # and one identifiable reverse primer, if not multiple.
    logger.info(
        f"Finalizing primer pairings amongst {len(respliced_dfs)} amplicons by verifying that every forward primer comes with a reverse primer.",
    )
    final_df = finalize_primer_pairings(
        respliced_dfs,
        args.fwd_suffix,
        args.rev_suffix,
    )

    # drop the now unnecessary columns, sort by start and stop position, and write
    # out with the user-provided output prefix
    logger.info(
        f"The input BED file contained {bed_df.shape[0]} primers. After searching for any spike-ins to resplice, the final BED file contains {final_df.shape[0]}",
    )
    final_df.drop("Amplicon").drop("ORIG_NAME").sort(
        "Ref",
        "Start Position",
        "Stop Position",
    ).write_csv(
        f"{args.output_prefix}.bed",
        separator="\t",
        include_header=False,
    )

    logger.success(
        f"Reading and processing of primers in '{args.input_bed}' was successful. The output BED file was written to '{Path.cwd()!s}/{args.output_prefix}.bed'. Goodbye!",
    )


if __name__ == "__main__":
    main()
