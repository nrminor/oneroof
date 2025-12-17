#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "polars",
# ]
# ///

"""Summarize amplicon coverage from stats files and BED file."""

from __future__ import annotations

import argparse

import polars as pl


def main(bed_file: str, output_file: str, stats_pattern: str) -> None:
    """
    Create amplicon summary TSV from stats files and BED file.

    Reads all stats files matching the glob pattern, extracts sample and amplicon
    names from the 'file' column, joins with BED file to get amplicon positions,
    and writes a summary TSV.

    Parameters:
        bed_file: Path to BED file with primer positions
        output_file: Path to output TSV file
        stats_pattern: Glob pattern for stats files (e.g., "stats_*.tsv")
    """
    # Read all stats files with glob, extract info from the 'file' column
    # (which contains filenames like "SAMPLE.QIAseq_X-Y.no_downsampling.fasta.gz")
    stats = (
        pl.scan_csv(
            stats_pattern,
            separator="\t",
            glob=True,
        )
        .with_columns(
            # Normalize .QIAseq_ to _QIAseq_ for consistent parsing
            pl.col("file")
            .str.replace(r"\.QIAseq_", "_QIAseq_")
            .alias("normalized_file")
        )
        .with_columns(
            # Extract sample name: everything before _QIAseq_, or first dot-segment
            pl.when(pl.col("normalized_file").str.contains("_QIAseq_"))
            .then(
                pl.col("normalized_file").str.extract(
                    r"^([^_]+(?:_[^_]+)*?)_QIAseq_", group_index=1
                )
            )
            .otherwise(
                pl.col("normalized_file").str.extract(r"^([^.]+)", group_index=1)
            )
            .alias("sample_name"),
            # Extract amplicon name: QIAseq_XXX (including any suffix like -1)
            pl.col("file")
            .str.extract(r"(QIAseq_[^.]+)", group_index=1)
            .alias("amplicon_name"),
            # Extract base amplicon for joining: QIAseq_N (just the number, no suffix)
            pl.col("file")
            .str.extract(r"(QIAseq_\d+)", group_index=1)
            .alias("base_amplicon"),
        )
        .select(
            "sample_name",
            "amplicon_name",
            "base_amplicon",
            pl.col("num_seqs").alias("reads"),
        )
    )

    # Read BED file, compute amplicon start/end positions
    bed = (
        pl.scan_csv(
            bed_file,
            separator="\t",
            has_header=False,
            new_columns=["chrom", "start", "end", "name", "score", "strand"],
        )
        .with_columns(
            # Extract base amplicon from primer name (e.g., QIAseq_2_LEFT -> QIAseq_2)
            pl.col("name")
            .str.extract(r"(QIAseq_\d+)", group_index=1)
            .alias("base_amplicon"),
            # Flag LEFT vs RIGHT primers
            pl.col("name").str.contains("_LEFT").alias("is_left"),
            pl.col("name").str.contains("_RIGHT").alias("is_right"),
        )
        .filter(pl.col("base_amplicon").is_not_null())
    )

    # Aggregate to get min(start) for LEFT primers, max(end) for RIGHT primers
    amplicon_positions = bed.group_by("base_amplicon").agg(
        pl.col("start").filter(pl.col("is_left")).min().alias("start_pos"),
        pl.col("end").filter(pl.col("is_right")).max().alias("end_pos"),
    )

    # Join stats with positions and format output
    result = (
        stats.join(amplicon_positions, on="base_amplicon", how="left")
        .select(
            "sample_name",
            "amplicon_name",
            pl.col("start_pos").cast(pl.String).fill_null("NA"),
            pl.col("end_pos").cast(pl.String).fill_null("NA"),
            "reads",
        )
        .collect()
    )

    result.write_csv(output_file, separator="\t")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Summarize amplicon coverage from stats files and BED file.",
    )
    parser.add_argument("--bed", required=True, help="Path to BED file")
    parser.add_argument(
        "--output",
        default="amplicon_summary.tsv",
        help="Output TSV file",
    )
    parser.add_argument(
        "--pattern",
        default="stats_*.tsv",
        help="Glob pattern for stats files (default: stats_*.tsv)",
    )

    args = parser.parse_args()
    main(args.bed, args.output, args.pattern)
