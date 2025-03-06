#! /usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "plotnine",
#     "polars",
#     "pyarrow",
# ]
# ///

"""
This module provides functionality for processing and visualizing coverage data from BED files.

It includes functions for parsing command-line arguments, reading sample lookup data,
accumulating coverage data from multiple files, fixing data frames, and plotting coverage
across different samples and chromosomes.

The main components of this module are:
1. Command-line argument parsing
2. Sample lookup data reading
3. Coverage data accumulation and processing
4. Data frame fixing for analysis
5. Coverage plotting
6. Main execution flow

This module requires external libraries such as argparse, json, os, pathlib, polars, and plotnine.

Usage:
    Run this script from the command line with appropriate arguments to process
    coverage data and generate visualizations.

Example:
    python multisample_plot.py --input_dir /path/to/bed_files --sample_lookup /path/to/lookup.json --min_coverage 20
"""

from __future__ import annotations

import argparse
from math import log10
from pathlib import Path

import polars as pl
from plotnine import (
    aes,
    facet_wrap,
    geom_hline,
    geom_line,
    ggplot,
    ggsave,
    guide_legend,
    guides,
    labs,
    theme_minimal,
)


def parse_command_line_args() -> argparse.Namespace:
    """
    Parse command line arguments for the script.

    This function sets up the argument parser and defines the required and optional
    arguments for the script. It handles the following arguments:
    - input_dir: Directory to scan for BED files.
    - sample_lookup: JSON-formatted sample lookup file.
    - min_coverage: Minimum coverage threshold (optional, default: 20).

    Returns:
        argparse.Namespace: An object containing the parsed command line arguments.

    Example:
        args = parse_command_line_args()
        input_directory = args.input_dir
        sample_lookup_file = args.sample_lookup
        min_coverage_threshold = args.min_coverage
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input_dir",
        "-i",
        type=Path,
        required=True,
        help="Directory to scan for BED files.",
    )
    parser.add_argument(
        "--min_coverage",
        "-m",
        type=int,
        required=False,
        default=20,
        help="The minimum expected coverage.",
    )
    parser.add_argument(
        "--log",
        action="store_true",
        help="Whether to plot the log of coverage.",
    )
    return parser.parse_args()


def accumulate_cov_dfs(directory: str, sample_list: list[str]) -> pl.DataFrame:
    """
    Accumulate and concatenate multiple CSV files from a specified directory into a single Polars DataFrame.

    This function reads all CSV files in the given directory, each file should contain columns named
    'chromosome', 'start', 'stop', and 'coverage'. It adds an additional column called 'sample',
    with values taken from the provided 'sample_lookup' list. The length of the 'sample_lookup' list must match
    the number of CSV files in the directory.

    The function performs the following steps:
    1. Lists all files in the specified directory.
    2. Reads each CSV file into a Polars DataFrame.
    3. Adds a 'sample' column to each DataFrame using values from the 'sample_lookup' list.
    4. Concatenates all DataFrames into a single DataFrame.

    Parameters:
    - directory (str): The path to the directory containing the CSV files. Each CSV file should have columns
                       named 'chromosome', 'start', 'stop', and 'coverage'.
    - sample_lookup (list): A list of values to be added as the 'sample' column in the resulting DataFrame.
                      The length of this list should match the number of CSV files in the directory.

    Returns:
    - pl.DataFrame: A Polars DataFrame containing data from all files, with an additional 'sample' column.

    Raises:
    - ValueError: If the number of files in the directory does not match the length of the 'sample_lookup' list.

    Example:
    >>> directory = 'path/to/csv_files'
    >>> sample_lookup = ['Sample1', 'Sample2', 'Sample3']
    >>> result_df = accumulate_cov_dfs(directory, sample_lookup)
    >>> print(result_df)
    shape: (number_of_rows, 5)
    ┌───────────┬───────┬───────┬──────────┬────────┐
    │ chromosome │ start │ stop  │ coverage │ sample │
    │ ---       │ ---   │ ---   │ ---      │ ---    │
    │ str       │ i64   │ i64   │ i64      │ str    │
    ├───────────┼───────┼───────┼──────────┼────────┤
    │ ...       │ ...   │ ...   │ ...      │ ...    │
    └───────────┴───────┴───────┴──────────┴────────┘
    """
    assert Path(
        directory,
    ).is_dir(), f"The provided input directory {directory} does not exist."

    sample_lookup = {
        sample_id: Path(directory) / Path(f"{sample_id}.per-base.bed")
        for sample_id in sample_list
    }

    df_list = []
    for sample_id, bed_file in sample_lookup.items():
        bc_df = (
            pl.read_csv(
                bed_file,
                separator="\t",
                has_header=False,
                new_columns=["chromosome", "start", "stop", "coverage"],
            )
            .with_columns(sample=pl.lit(sample_id))
            .with_columns(
                pl.int_ranges(start=pl.col("start"), end=pl.col("stop")).alias(
                    "position",
                ),
            )
            .drop("start", "stop")
            .explode("position")
        )

        df_list.append(bc_df)

    bc_stacked = df_list[0]
    for df in df_list[1:]:
        bc_stacked = bc_stacked.vstack(df)
    return bc_stacked


def plot_log_coverages(
    all_barcodes: pl.DataFrame,
    min_desired_depth: int = 20,
) -> ggplot:
    """
    Create a line plot of log-transformed coverage depth across different samples and chromosomes.

    This function generates a plot showing the log-transformed coverage depth for different
    samples across chromosomes or segments. The coverage values are log10-transformed before
    plotting to better visualize large variations in coverage.

    Parameters:
        all_barcodes (pl.DataFrame): A Polars DataFrame containing the coverage data. Must include
                                    columns 'coverage', 'sample', 'position', and 'chromosome'.

    Returns:
        ggplot: A `ggplot` object representing the log-transformed coverage plot.

    Example:
        >>> sample_df = pl.DataFrame({
        >>>     'position': [1, 2, 3, 4],
        >>>     'coverage': [10, 100, 1000, 10000],
        >>>     'sample': ['A', 'A', 'B', 'B'],
        >>>     'chromosome': ['chr1', 'chr1', 'chr2', 'chr2']
        >>> })
        >>> plot = plot_log_coverages(sample_df)
    """
    barcodes_log = all_barcodes.with_columns(
        pl.col("coverage").log(base=10).alias("log_coverage"),
    )
    log_desired_depth = log10(min_desired_depth)
    return (
        ggplot(
            barcodes_log.to_pandas(),
            aes(
                x="position",
                y="log_coverage",
                color="sample",
            ),
        )
        + geom_line()
        + geom_hline(yintercept=log_desired_depth, linetype="dashed")
        + labs(
            title="Log Coverage across Samples",
            x="Position on Chromosome/Segment",
            y="Log10 Depth of Coverage",
        )
        + facet_wrap("~chromosome", scales="free_x")
        + theme_minimal()
        + guides(color=guide_legend(ncol=3))
    )


def plot_coverages(all_barcodes: pl.DataFrame, min_desired_depth: int = 20) -> ggplot:
    """
    Create a line plot of coverage depth across different barcodes and chromosomes.

    This function generates a plot showing the coverage depth for different samples (barcodes) across
    chromosomes or segments. The coverage depth is plotted as a line, and a horizontal dashed line
    indicates the minimum desired depth.

    The input DataFrame should have the following columns:
    - 'start': The start position on the chromosome/segment.
    - 'stop': The stop position on the chromosome/segment.
    - 'coverage': The depth of coverage at each position.
    - 'sample': The sample or barcode identifier.
    - 'chromosome': The chromosome or segment name.

    Parameters:
    - all_barcodes (pl.DataFrame): A Polars DataFrame containing the coverage data. Must include columns
      'start', 'stop', 'coverage', 'sample', and 'chromosome'.
    - min_desired_depth (int): The minimum desired depth of coverage to be indicated by a horizontal dashed line
      on the plot. Default is 20.

    Returns:
    - ggplot: A `ggplot` object representing the coverage plot.

    Example:
    >>> import polars as pl
    >>> from plotnine import ggplot
    >>> # Sample DataFrame
    >>> df = pl.DataFrame({
    >>>     'start': [1, 10, 20, 30],
    >>>     'stop': [5, 15, 25, 35],
    >>>     'coverage': [15, 25, 5, 30],
    >>>     'sample': ['A', 'A', 'B', 'B'],
    >>>     'chromosome': ['chr1', 'chr1', 'chr2', 'chr2']
    >>> })
    >>> plot = plot_coverages(df)
    >>> print(plot)
    """
    return (
        ggplot(
            all_barcodes.to_pandas(),
            aes(
                x="position",
                y="coverage",
                color="sample",
            ),
        )
        + geom_line()
        + geom_hline(yintercept=min_desired_depth, linetype="dashed")
        + labs(
            title="Coverage across Samples",
            x="Position on Chromosome/Segment",
            y="Depth of Coverage (read count)",
        )
        + facet_wrap("~chromosome", scales="free_x")
        + theme_minimal()
        + guides(color=guide_legend(ncol=3))
    )


def main() -> None:
    """
    Main function to execute the sequence of data processing and plotting tasks.

    This function orchestrates the process of loading sample data, accumulating coverage data frames,
    fixing the data frame, plotting coverage, and saving the plot to a PDF file.

    The function performs the following steps:
    1. Reads a JSON file containing a list of samples.
    2. Accumulates coverage data frames based on the loaded sample list.
    3. Fixes the data frame for further analysis.
    4. Plots the coverage data with a specified minimum depth.
    5. Saves the generated plot to a PDF file.

    The function does not return any value.

    Steps:
    1. **Read Samples**: Reads sample data from a JSON file using `read_sample_lookup()`.
    2. **Accumulate Coverage Data Frames**: Processes the data frames using `accumulate_cov_dfs()`.
    3. **Plot Coveraged**: Generates a plot using `plot_coverages()` or `plot_log_coverages()`.
    4. **Save Plot**: Saves the plot to a PDF file using `ggsave()`.

    Example:
    >>> main()

    Note:
    Ensure that the necessary functions (`read_sample_lookup`, `accumulate_cov_dfs`, `fix_dataframe`, `plot_coverage`, and `ggsave`) are correctly implemented and imported for this function to work properly.
    """
    args = parse_command_line_args()
    min_desired_depth = args.min_coverage
    sample_list = [
        path.name.replace(".per-base.bed", "")
        for path in Path(args.input_dir).glob("*.per-base.bed")
    ]
    sample_dataframe = accumulate_cov_dfs(args.input_dir, sample_list)

    if args.log:
        plot_instance = plot_log_coverages(
            sample_dataframe,
            min_desired_depth,
        )
    else:
        plot_instance = plot_coverages(
            sample_dataframe,
            min_desired_depth,
        )

    ggsave(
        plot_instance,
        "multisample.fixed_y.coverage.pdf",
        format="pdf",
        height=6,
        width=11,
    )


if __name__ == "__main__":
    main()
