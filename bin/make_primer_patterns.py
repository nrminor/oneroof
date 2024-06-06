#!/usr/bin/env python3

"""
Usage:
```
python3 make_primer_patterns.py <FASTA> <OUTPUT_PREFIX>
```
"""

import os
import sys
from typing import List


def generate_regex_patterns(primer_fasta: str, label: str) -> None:
    """
    Generates regular expression patterns based on primer sequences from a FASTA file.

    This function reads a FASTA file containing two primer sequences, extracts the sequences,
    and generates two regular expression patterns:
    1. A pattern to match all characters after the forward primer sequence (excluding newlines and whitespace).
    2. A pattern to match all characters before the reverse primer sequence (excluding newlines and whitespace).

    The generated patterns are written to an output file named based on the provided label.

    Parameters:
    primer_fasta (str): The path to the FASTA file containing the primer sequences.
    label (str): The label used to name the output file where the patterns will be saved.

    Raises:
    AssertionError: If the FASTA file does not contain exactly two sequences.
    AssertionError: If the FASTA file is not formatted like an output from `bedtools getfasta`.

    The format expected for the FASTA file is similar to the following:
    ```
    >identifier:0-16
    ATCG...
    >identifier:17-32
    GCTA...
    ```

    Example:
    ```
    generate_regex_patterns("primers.fasta", "output_patterns")
    ```
    This would read the primers from "primers.fasta" and write the patterns to "output_patterns.txt".
    """
    # initialize a list of strings and parse the lines from the primer FASTA\
    # into it
    lines: List[str]
    with open(primer_fasta, "r", encoding="utf8") as primer_handle:
        lines = [line.strip() for line in primer_handle]

    # pull out the sequences
    seqs = [line for line in lines if not line.startswith(">")]

    # make sure there are only two sequences present, and thus that we are only
    # considering a single amplicon
    assert (
        len(seqs) == 2
    ), f"The provided FASTA does not contain exactly two sequences:\n{seqs}"

    # pull out the start coordinates, assuming `bedtools getfasta` formatting,
    # and run an assert that checks this assumption, a corollary of which is
    # that the forward primer starts before the reverse primer starts
    start_coords = [
        int(line.split(":")[-1].split("-")[0]) for line in lines if line.startswith(">")
    ]
    assert (
        start_coords[0] < start_coords[1]
    ), "Please double check that the provided FASTA is formatted like an output from `bedtools getfasta`, e.g.\n\n'>PP599462.1:0-16'"

    # unpack the sequences
    forward_primer, reverse_primer = seqs

    # use raw strings and '+' operator overloading to construct the patterns
    forward_pattern = r"(?<=" + forward_primer + r")[^\s]*"
    reverse_pattern = r"[^\s]*(?=" + reverse_primer + r")"

    # iterate through the patterns to write the lines of the pattern file
    with open(f"{label}.txt", "w", encoding="utf8") as output_handle:
        output_handle.write(forward_pattern)
        output_handle.write("\n")
        output_handle.write(reverse_pattern)


def main() -> None:
    """
    Main function to generate regular expression patterns from a primer FASTA file.

    This function reads the primer FASTA file path and output label from the command line arguments,
    ensures the file exists, and calls `generate_regex_patterns` to create the regex patterns and
    save them to a file.

    Command Line Arguments:
    primer_fasta_path (str): The path to the primer FASTA file.
    label (str): The label for naming the output file.

    Raises:
    AssertionError: If the provided primer FASTA file does not exist.
    """
    # pull  the primer FASTA from the command line
    primer_fasta_path = sys.argv[1]
    label = sys.argv[2]

    # make sure the file provided by the user exists
    assert os.path.isfile(
        primer_fasta_path
    ), f"The provided primer FASTA, {primer_fasta_path}, does not exist"

    # generate the new regex patterns
    generate_regex_patterns(primer_fasta_path, label)


if __name__ == "__main__":
    main()
