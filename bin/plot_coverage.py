#!/usr/bin/env python3

"""
usage: plot_coverage.py [-h] --input INPUT [--label LABEL]

options:
  -h, --help            show this help message and exit
  --input INPUT, -i INPUT
                        Mosdepth data file to use in plotting.
  --label LABEL, -l LABEL
                        Label to use as a prefix for the output plot
"""

import argparse
import os
from pathlib import Path
from typing import Tuple

import polars as pl
from plotnine import (
    aes,
    facet_wrap,
    geom_hline,
    geom_rect,
    ggplot,
    ggsave,
    labs,
    theme_minimal,
)


def parse_command_line_args() -> argparse.Namespace:
    """
    Parse command line arguments for the plotting script.

    This function sets up and parses command line arguments needed for
    the plotting script, including the input Mosdepth data file and an
    optional label for the output plot.

    Returns:
        argparse.Namespace: An object containing the parsed command line arguments.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        "-i",
        type=Path,
        required=True,
        help="Mosdepth data file to use in plotting.",
    )
    parser.add_argument(
        "--label",
        "-l",
        type=str,
        required=False,
        default="sample",
        help="Label to use as a prefix for the output plot",
    )
    parser.add_argument(
        "--depth",
        "-d",
        type=int,
        required=False,
        default=20,
        help="Minimum depth of coverage to visualize on the plots.",
    )
    args = parser.parse_args()

    return args


def construct_plot(coverage_lf: pl.LazyFrame, label: str, depth: int) -> ggplot:
    """
    Render a coverage plot using the grammar of graphics implementation in plotnine.

    This function takes a LazyFrame of coverage data and a label, converts the
    LazyFrame to a pandas DataFrame, and creates a plot of the coverage data
    that can be faceted by chromosome or segment. It will also highlight areas
    of low coverage.

    Args:
        coverage_lf (pl.LazyFrame): A Polars LazyFrame containing the coverage data.
        label (str): A label to use in the plot title to identify the sample.

    Returns:
        ggplot: A ggplot object representing the coverage plot.
    """
    low_cov_df = coverage_lf.filter(pl.col("coverage") < 0).collect().to_pandas()
    base_plot = (
        ggplot(
            coverage_lf.collect().to_pandas(),
            aes(xmin="start", xmax="stop", ymin=0, ymax="coverage"),
        )
        + geom_rect(fill="black", color="black")
        + geom_hline(yintercept=depth, linetype="dashed")
        + labs(
            title=f"Coverage for Sample ID {label}",
            x="Position on Chromosome/Segment",
            y="Depth of Coverage (read count)",
        )
        + theme_minimal()
    )
    if low_cov_df.shape[0] == 0:
        return base_plot

    return base_plot + geom_rect(
        data=low_cov_df,
        mapping=aes(ymin=0, ymax=float("inf")),
        fill="gray",
        alpha=0.3,
    )


def finish_plot(core_plot: ggplot, contig_count: int) -> Tuple[ggplot]:
    """
    Finalize the plot by adding faceting based on the number of contigs.

    This function takes a core plot and the number of contigs, then applies appropriate
    faceting using facet_wrap. For a single contig, it returns one plot with free x-axis
    scaling. For multiple contigs, it returns two plots: one with free x-axis scaling and
    another with both free x and y-axis scaling.

    Parameters:
    core_plot (ggplot): The base plot to be faceted.
    contig_count (int): The number of contigs in the dataset.

    Returns:
    Tuple[ggplot]: A tuple containing either one or two ggplot objects:
        - For single contig: (plot_with_free_x,)
        - For multiple contigs: (plot_with_free_x, plot_with_free_xy)

    Raises:
    AssertionError: If the provided contig_count is less than 1.

    Note:
    This function uses the facet_wrap function from the plotnine library to create
    separate panels for each chromosome/contig.
    """
    assert contig_count >= 1, f"Invalid contig count provided: {contig_count}"

    # make two plots with free and fixed Y-axis for multisegment references. Otherwise,
    # just make 1
    match contig_count:
        case 1:
            return core_plot + facet_wrap("~chromosome", scales="free_x")
        case _:
            return (
                core_plot + facet_wrap("~chromosome", scales="free_x"),
                core_plot + facet_wrap("~chromosome", scales="free"),
            )


def compute_perc_cov(
    coverage_lf: pl.LazyFrame, label: str, contig_count: int, depth: int
) -> pl.LazyFrame:
    """
    Compute the percentage of coverage above a specified depth for each chromosome/contig.

    This function processes a LazyFrame containing coverage data to calculate the proportion
    of each chromosome/contig that has coverage above a specified depth threshold.

    Parameters:
    coverage_lf (pl.LazyFrame): A Polars LazyFrame containing coverage data with columns
                                'chromosome', 'start', 'stop', and 'coverage'.
    contig_count (int): The expected number of unique contigs/chromosomes in the data.
    depth (int): The coverage depth threshold.

    Returns:
    pl.LazyFrame: A LazyFrame with columns 'chromosome' and 'proportion_above_cutoff',
                    where 'proportion_above_cutoff' represents the fraction of the
                    chromosome/contig with coverage above the specified depth.

    Raises:
    AssertionError: If the resulting LazyFrame does not contain the expected number of rows
                    (one per contig/chromosome).

    Notes:
    - The function assumes the input LazyFrame has columns 'chromosome', 'start', 'stop',
        and 'coverage'.
    - The computation is done lazily and the result is not materialized until collected.
    """

    percent_passing_lf = (
        coverage_lf.with_columns(
            pl.lit(label).alias("sample id"),
            pl.col("stop").max().over("chromosome").name.suffix("_max"),
            pl.col("start").min().over("chromosome").name.suffix("_min"),
            (pl.col("stop") - pl.col("start")).alias("block_length"),
        )
        .with_columns(
            (pl.col("stop_max") - pl.col("start_min")).alias("reference length")
        )
        .drop("stop_max", "start_min")
        .filter(pl.col("coverage") >= depth)
        .drop("coverage")
        .unique()
        .with_columns(pl.col("block_length").sum().over("chromosome").alias("sum"))
        .drop("start", "stop", "block_length")
        .with_columns(
            (pl.col("sum") / pl.col("reference length"))
            .over("chromosome")
            .alias(f"proportion ≥ {depth}X coverage")
        )
        .drop("sum")
        .select(
            "sample id",
            "chromosome",
            "reference length",
            f"proportion ≥ {depth}X coverage",
        )
        .unique()
    )

    if percent_passing_lf.collect().shape[0] == contig_count:
        return percent_passing_lf

    return (
        coverage_lf.with_columns(
            pl.lit(label).alias("sample id"),
            pl.col("stop").max().over("chromosome").name.suffix("_max"),
            pl.col("start").min().over("chromosome").name.suffix("_min"),
            (pl.col("stop") - pl.col("start")).alias("block_length"),
        )
        .with_columns(
            (pl.col("stop_max") - pl.col("start_min")).alias("reference length")
        )
        .select("chromosome", "reference length")
        .unique()
        .join(
            percent_passing_lf,
            on=["chromosome", "reference length"],
            how="left",
            coalesce=True,
        )
        .with_columns(
            pl.when(pl.col(f"proportion ≥ {depth}X coverage").is_null())
            .then(0)
            .otherwise(pl.col(f"proportion ≥ {depth}X coverage"))
            .alias(f"proportion ≥ {depth}X coverage"),
            pl.when(pl.col("sample id").is_null())
            .then(pl.lit(label))
            .otherwise(pl.col("sample id"))
            .alias("sample id"),
        )
    )


def main() -> None:
    """
    Script entrypoint. Main parses args, checks that the provided input file exists,
    lazily scans the coverage BED file, renders a facet-wrapped visualization of
    coverage, and saves it to a PDF.
    """
    # collect parsed command line arguments and make sure the input exists
    args = parse_command_line_args()
    assert args.input and os.path.isfile(
        args.input
    ), f"The provided input file {args.input} does not exist."

    # open a Polars query to lazily read the file with explicit column names
    coverage_lf = pl.scan_csv(
        args.input,
        separator="\t",
        new_columns=["chromosome", "start", "stop", "coverage"],
    )

    # determine how many contigs are present
    contig_count = coverage_lf.select("chromosome").unique().collect().shape[0]

    # construct the core plot
    core_plot = construct_plot(coverage_lf, args.label, args.depth)

    # finish plot by handling faceting differently depending on whether there are multiple contigs
    # in the reference
    rendered = finish_plot(core_plot, contig_count)

    # if the rendered can be unpacked into two plots, write both separately. Otherwise,
    # just output the one plot with a fixed Y-axis
    try:
        fixed_plot, free_plot = rendered
    except TypeError:
        ggsave(rendered, f"{args.label}.coverage.pdf", format="pdf", height=6, width=11)
    else:
        ggsave(
            fixed_plot,
            f"{args.label}.fixed_y.coverage.pdf",
            format="pdf",
            height=6,
            width=11,
        )
        ggsave(
            free_plot,
            f"{args.label}.free_scales.coverage.pdf",
            format="pdf",
            height=6,
            width=11,
        )

    # write out a table of the percentage of positions that are greater than the cutoff
    (
        compute_perc_cov(coverage_lf, args.label, contig_count, args.depth)
        .collect()
        .write_csv(f"{args.label}.passing_cov.tsv", separator="\t")
    )


if __name__ == "__main__":
    main()
