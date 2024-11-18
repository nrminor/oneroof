#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "loguru",
# ]
# ///

"""
This module provides functionality to validate and normalize BED files containing primer information.

It includes functions to:
1. Parse command-line arguments
2. Check for correct primer suffixes
3. Orient primer coordinates
4. Validate and normalize the input BED file

The main function processes the input BED file, checks for valid primer labels,
and reorients the coordinates if necessary. The resulting validated and normalized
BED file is then written to the specified output file.

Usage:
    python validate_primer_bed.py --input_bed <input_bed_file> [--output_prefix <output_prefix>]
                          [--fwd_suffix <forward_suffix>] [--rev_suffix <reverse_suffix>]

For more information on usage and options, use the --help flag.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from loguru import logger


def parse_command_line_args() -> argparse.Namespace:
    """
    Parse command line arguments for the BED file validation script.

    This function sets up the argument parser and defines the following arguments:
    --input_bed (-i): Input BED file to be normalized and validated (required)
    --output_prefix (-o): Prefix for the output validated BED file (optional, default: "validated")
    --fwd_suffix (-f): Suffix expected in the names for forward primers (optional, default: "_LEFT")
    --rev_suffix (-r): Suffix expected in the names for reverse primers (optional, default: "_RIGHT")

    Returns:
        argparse.Namespace: An object containing the parsed command-line arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input_bed",
        "-i",
        type=Path,
        required=True,
        help="Input BED file to be normalized and validated before any further usage.",
    )
    parser.add_argument(
        "--output_prefix",
        "-o",
        type=str,
        required=False,
        default="validated",
        help="The prefix to use when naming the output validated BED file.",
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


def check_for_suffixes(
    row: list[str],
    fwd_suffix: str = "_LEFT",
    rev_suffix: str = "_RIGHT",
) -> None | str:
    """
    Check if the primer label in the given row contains the expected forward or reverse suffix.

    Args:
        row (list[str]): A list representing a row from the BED file.
        fwd_suffix (str, optional): The expected suffix for forward primers. Defaults to "_LEFT".
        rev_suffix (str, optional): The expected suffix for reverse primers. Defaults to "_RIGHT".

    Returns:
        None | str: None if the primer label contains the expected suffix, otherwise returns the invalid primer label.

    Raises:
        AssertionError: If the row has fewer than the required 6 columns.
    """

    required_col_count = 6
    assert (
        len(row) >= required_col_count
    ), f"BED file row has fewer than the required 6 columns: \n\n {row}"

    primer_label = row[3]
    if fwd_suffix in primer_label or rev_suffix in primer_label:
        return None

    logger.error(
        f"Primer label without the expected forward suffix, {fwd_suffix}, or the expected reverse suffix, {rev_suffix}, encountered: {primer_label}. Primer labels must be standardized such that all primers have either the forward or reverse suffixes before proceeding.",
    )
    return primer_label


def check_for_pairs(
    rows: list[list[str]],
    fwd_suffix: str = "_LEFT",
    rev_suffix: str = "_RIGHT",
) -> list[str]:
    """
    Check for primer pairs in the given rows of a BED file.

    This function identifies primer pairs based on the provided forward and reverse suffixes.
    It counts the occurrences of each amplicon label and identifies singletons (unpaired primers).

    Args:
        rows (list[list[str]]): A list of rows from the BED file, where each row is a list of strings.
        fwd_suffix (str, optional): The suffix for forward primers. Defaults to "_LEFT".
        rev_suffix (str, optional): The suffix for reverse primers. Defaults to "_RIGHT".

    Returns:
        list[str]: A list of amplicon labels that appear only once (singletons).

    Note:
        This function assumes that primer labels follow a specific naming convention,
        where the amplicon label is separated from the primer direction suffix by either
        the provided suffixes or a hyphen.
    """
    primer_matching = []
    amplicon_tally = {}
    singletons = []
    pair_limit = 2
    for row in rows:
        primer_label = row[3]
        if primer_label.endswith((fwd_suffix, rev_suffix)):
            amplicon_label = primer_label.split(
                fwd_suffix if primer_label.endswith(fwd_suffix) else rev_suffix,
            )
            hyphen_detector = str(amplicon_label[0])
            if "-" in hyphen_detector:
                amplicon_label = hyphen_detector.split("-")
            primer_matching.append(amplicon_label[0])

    for primer_label in primer_matching:
        if primer_label not in amplicon_tally:
            amplicon_tally[primer_label] = 0
        amplicon_tally[primer_label] += 1

    for amplicon, tally in amplicon_tally.items():
        if tally < pair_limit:
            singletons.append(amplicon)

    return singletons


def orient_primer_coords(row: list[str], row_index: int) -> list[str]:
    """
    Orient primer coordinates to ensure start position precedes stop position.

    This function checks if the start position is greater than the stop position
    and swaps them if necessary. It also adjusts the strand information accordingly.

    Args:
        row (list[str]): A list representing a row from the BED file.
        row_index (int): The index of the current row in the BED file.

    Returns:
        list[str]: The row with correctly oriented coordinates and strand information.

    Raises:
        ValueError: If an unsupported value is encountered in the strand column.
    """

    start = int(row[1])
    stop = int(row[2])
    sign = row[5]
    if start < stop:
        return row

    row[2] = str(start)
    row[1] = str(stop)

    if sign == "+":
        row[5] = "-"
    elif sign == "-":
        row[5] = "+"
    else:
        message = f"unsupported value encountered in the strand column of BED file row {row_index}: {row}"
        raise ValueError(message)

    return row


def normalize_bed_lines(
    input_bed: Path,
    output_prefix: str = "validated",
    fwd_suffix: str = "_LEFT",
    rev_suffix: str = "_RIGHT",
) -> None:
    """
    Validate and normalize the input BED file containing primer information.

    This function reads the input BED file, checks for valid primer labels,
    orients the coordinates correctly, and writes the normalized data to a new file.

    Args:
        input_bed (Path): Path to the input BED file.
        output_prefix (str, optional): Prefix for the output file name. Defaults to "validated".
        fwd_suffix (str, optional): Suffix for forward primers. Defaults to "_LEFT".
        rev_suffix (str, optional): Suffix for reverse primers. Defaults to "_RIGHT".

    Raises:
        AssertionError: If any invalid primer labels are found or if a row has fewer than 6 columns.
        ValueError: If an unsupported value is encountered in the strand column.

    Returns:
        None: The function writes the normalized data to a new file and doesn't return anything.
    """
    with open(input_bed, encoding="utf8") as file, open(
        f"{output_prefix}.bed",
        "w",
        encoding="utf8",
    ) as output:
        # collect all the lines from the bed file
        lines = [row.strip().split("\t") for row in file if row != ""]

        # check that all lines contain primer labels with either the expected forward
        # suffix or the expected reverse suffix, which is necessary to be able to
        # properly handle cases where one primer may end up paired with multiple other
        # primers
        label_test = [
            check_for_suffixes(line, fwd_suffix, rev_suffix) for line in lines
        ]
        invalid_labels = [label for label in label_test if label is not None]

        # crash if any invalid primer labels were provided
        assert (
            len(invalid_labels) == 0
        ), f"Invalid primer label(s) without the expected forward suffix, {fwd_suffix}, and the expected reverse suffix, {rev_suffix}detected: {invalid_labels}."

        # make sure that each primer label can be mapped back to at least two primers
        pair_test = check_for_pairs(lines, fwd_suffix, rev_suffix)
        assert (
            len(pair_test) == 0
        ), f"These amplicon labels do not appear in more than one primer: {pair_test}."

        # for all lines, make sure there are at least 6 columns, make sure all primers
        # are oriented correctly such that all start positions precede stop positions,
        # an append to a list of fixed lines
        fixed_lines = []
        for i, row in enumerate(lines):
            required_col_count = 6
            assert (
                len(row) >= required_col_count
            ), f"BED file row {i} has fewer than the required 6 columns: \n\n {row}"
            new_row = orient_primer_coords(row, i)
            fixed_lines.append("\t".join(new_row))

        output.write("\n".join(fixed_lines))


def main() -> None:
    # parse command line arguments
    args = parse_command_line_args()

    # make sure that the provided BED file exists and is a file
    assert args.input_bed.is_file(), f"Inputed path is not file: {args.input_bed}"

    # run normalizations
    normalize_bed_lines(args.input_bed, args.output_prefix)


if __name__ == "__main__":
    main()
