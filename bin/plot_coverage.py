#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "patito",
#     "plotnine",
#     "polars",
#     "pyarrow",
# ]
# ///

"""
usage: plot_coverage.py [-h] --input INPUT [--label LABEL]

options:
  -h, --help            show this help message and exit
  --input INPUT, -i INPUT
                        Mosdepth data file to use in plotting.
  --label LABEL, -l LABEL
                        Label to use as a prefix for the output plot
"""

from __future__ import annotations

import argparse
from pathlib import Path

import patito as pt
import polars as pl
from plotnine import (
    aes,
    facet_wrap,
    geom_hline,
    geom_line,
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

    return parser.parse_args()


class MosdepthBedSchema(pt.Model):
    """Schema for validating Mosdepth BED format coverage data."""

    chromosome: str = pt.Field(description="Chromosome or contig name")
    start: int = pt.Field(ge=0, description="Start position (0-based)")
    stop: int = pt.Field(gt=0, description="Stop position (exclusive)")
    coverage: float = pt.Field(ge=0, description="Coverage depth")


def validate_bed_data(bed_df: pl.LazyFrame) -> pl.LazyFrame:
    """Validate BED format coverage data using patito schema.

    Args:
        bed_df: LazyFrame with BED coverage data

    Returns:
        Validated LazyFrame

    Raises:
        ValueError: If validation fails
    """
    # Collect a small sample to validate schema
    sample_df = bed_df.head(1000).collect()

    # Validate using patito
    try:
        MosdepthBedSchema.validate(sample_df)
    except pt.exceptions.DataFrameValidationError as e:
        msg = f"BED coverage data validation failed:\n{e}"
        raise ValueError(msg) from e

    # Additional validation: ensure start < stop
    return bed_df.filter(pl.col("start") < pl.col("stop"))


def adapt_bed_to_df(input_bed: str | Path, depth: int) -> pl.LazyFrame:
    """
    Transform Mosdepth BED output into a Polars LazyFrame for visualization.

    This function takes a BED-format file output by Mosdepth and converts it into a
    LazyFrame suitable for coverage visualization. It performs initial data loading
    and transformations including position exploding and coverage aggregation.

    Args:
        input_bed (str | Path): Path to the BED file containing coverage data.
        depth (int): Minimum depth threshold for coverage analysis.

    Returns:
        pl.LazyFrame: A LazyFrame containing processed coverage data with columns for
                     chromosome, position, coverage, and depth assessment.
    """
    # Read the BED file
    raw_bed_lf = pl.scan_csv(
        input_bed,
        separator="\t",
        has_header=False,
        new_columns=["chromosome", "start", "stop", "coverage"],
    )

    # Validate the data
    bed_lf = validate_bed_data(raw_bed_lf)

    exploded_df = (
        bed_lf.with_columns(
            pl.int_ranges(start=pl.col("start"), end=pl.col("stop")).alias(
                "position",
            ),
        )
        .drop("start", "stop")
        .explode("position")
    )

    return exploded_df.join(
        exploded_df.group_by("chromosome")
        .agg(
            pl.col("position").count().alias("reference length"),
            pl.col("coverage").max().alias("max_coverage"),
        )
        .with_columns(
            pl.when(pl.col("max_coverage") < depth)
            .then(False)  # noqa: FBT003
            .otherwise(True)  # noqa: FBT003
            .alias("passes_depth_cutoff"),
        )
        .drop("max_coverage"),
        on="chromosome",
        how="left",
    )


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
    low_cov_df = (
        coverage_lf.filter(
            ~pl.col("passes_depth_cutoff") | (pl.col("coverage") < depth),
        )
        .with_columns(
            ((pl.col("position") - pl.col("position").shift(1)) != 1)
            .fill_null(True)  # noqa: FBT003
            .alias("new_block_flag"),
        )
        .with_columns(pl.col("new_block_flag").cum_sum().alias("block"))
        .drop("new_block_flag")
        .group_by(["chromosome", "block"])
        .agg(
            pl.col("position").min().alias("start"),
            pl.col("position").max().alias("stop"),
        )
        .collect()
    )

    assert low_cov_df.shape[0] != 0

    return (
        ggplot(
            coverage_lf.collect().to_pandas(),
        )
        + geom_rect(
            data=low_cov_df,
            mapping=aes(xmin="start", xmax="stop", ymin=0, ymax=float("inf")),
            fill="gray",
            alpha=0.3,
        )
        + geom_line(
            aes(
                x="position",
                y="coverage",
            ),
        )
        + geom_hline(yintercept=depth, linetype="dashed")
        + labs(
            title=f"Coverage for Sample ID {label}",
            x="Position on Chromosome/Segment",
            y="Depth of Coverage (read count)",
        )
        + facet_wrap("~chromosome", scales="free_x")
        + theme_minimal()
    )


def finish_plot(core_plot: ggplot, contig_count: int) -> ggplot | tuple[ggplot, ggplot]:
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
    coverage_lf: pl.LazyFrame,
    label: str,
    contig_count: int,
    depth: int,
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
        coverage_lf.join(
            (
                coverage_lf.filter(
                    ~pl.col("passes_depth_cutoff") | pl.col("coverage").ge(depth),
                )
                .group_by("chromosome")
                .agg(pl.col("position").count().alias("count"))
            ),
            on="chromosome",
            how="left",
        )
        .with_columns(
            pl.when(pl.col("passes_depth_cutoff"))
            .then(pl.col("count") / pl.col("reference length"))
            .otherwise(pl.lit(0.0))
            .alias(f"proportion ≥ {depth}X coverage"),
        )
        .drop("count", "position", "passes_depth_cutoff")
        .unique()
        .with_columns(
            pl.lit(label).alias("sample id"),
        )
        .select(
            "sample id",
            "chromosome",
            "reference length",
            f"proportion ≥ {depth}X coverage",
        )
        .unique()
    )

    assert percent_passing_lf.collect().shape[0] == contig_count, (
        f"The number of rows does not match the expected number of contigs, {contig_count}."
    )

    return percent_passing_lf


def main() -> None:
    """
    Script entrypoint. Main parses args, checks that the provided input file exists,
    lazily scans the coverage BED file, renders a facet-wrapped visualization of
    coverage, and saves it to a PDF.
    """
    # collect parsed command line arguments and make sure the input exists
    args = parse_command_line_args()
    assert args.input, f"The provided input file {args.input} does not exist."
    assert Path(
        args.input,
    ).is_file(), f"The provided input file {args.input} does not exist."

    # open a Polars query to lazily read the file with explicit column names
    coverage_lf = adapt_bed_to_df(args.input, args.depth)

    # determine how many contigs are present
    contig_count = coverage_lf.select("chromosome").unique().collect().shape[0]

    # construct the core plot
    core_plot = construct_plot(coverage_lf, args.label, args.depth)

    # finish plot by handling faceting differently depending on whether there are
    # multiple contigs in the reference
    rendered = finish_plot(core_plot, contig_count)

    # if the rendered can be unpacked into two plots, write both separately. Otherwise,
    # just output the one plot with a fixed Y-axis
    multi_segment_plot_count = 2
    if isinstance(rendered, tuple) and len(rendered) == multi_segment_plot_count:
        fixed_plot, free_plot = rendered
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
    elif isinstance(rendered, ggplot):
        ggsave(
            rendered,
            f"{args.label}.coverage.pdf",
            format="pdf",
            height=6,
            width=11,
        )
    else:
        assert isinstance(rendered, ggplot) or len(rendered) >= 1, (
            f"Unexpected behavior for plot rendering: {rendered}"
        )

    # write out a table of the percentage of positions that are greater than the cutoff
    (
        compute_perc_cov(coverage_lf, args.label, contig_count, args.depth)
        .collect()
        .write_csv(f"{args.label}.passing_cov.tsv", separator="\t")
    )


if __name__ == "__main__":
    main()
