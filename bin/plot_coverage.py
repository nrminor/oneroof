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
from pathlib import Path

import polars as pl
from plotnine import (
    aes,
    facet_wrap,
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
    args = parser.parse_args()

    return args


def render_plot(coverage_lf: pl.LazyFrame, label: str) -> ggplot:
    """
    Render a coverage plot using the grammar of graphics implementation in plotnine.

    This function takes a LazyFrame of coverage data and a label, converts the
    LazyFrame to a pandas DataFrame, and creates a plot of the coverage data
    faceted by chromosome or segment.

    Args:
        coverage_lf (pl.LazyFrame): A Polars LazyFrame containing the coverage data.
        label (str): A label to use in the plot title to identify the sample.

    Returns:
        ggplot: A ggplot object representing the coverage plot.
    """
    chroms = coverage_lf.select("chromosome").unique().collect()
    coverage_pd = coverage_lf.collect().to_pandas()

    if chroms.shape[0] > 1:
        return (
            ggplot(coverage_pd, aes(xmin="start", xmax="stop", ymin=0, ymax="coverage"))
            + geom_rect(fill="black", color="black")
            + labs(
                title=f"Coverage for Sample ID {label}",
                x="Position on Chromosome/Segment",
                y="Depth of Coverage (read count)",
            )
            + theme_minimal()
            + facet_wrap("~chromosome", scales="free_x")
        )

    return (
        ggplot(coverage_pd, aes(xmin="start", xmax="stop", ymin=0, ymax="coverage"))
        + geom_rect(fill="black", color="black")
        + labs(
            title=f"Coverage for Sample ID {label}",
            x="Position on Chromosome/Segment",
            y="Depth of Coverage (read count)",
        )
        + theme_minimal()
    )


def main() -> None:
    """
    Script entrypoint. Main parses args, checks that the provided input file exists,
    lazily scans the coverage BED file, renders a facet-wrapped visualization of
    coverage, and saves it to a PDF.
    """
    args = parse_command_line_args()

    assert args.input, f"The provided input file {args.input} does not exist."

    coverage_lf = pl.scan_csv(
        args.input,
        separator="\t",
        new_columns=["chromosome", "start", "stop", "coverage"],
    )

    plot = render_plot(coverage_lf, args.label)

    ggsave(plot, f"{args.label}.coverage.pdf", format="pdf", height=6, width=11)


if __name__ == "__main__":
    main()
