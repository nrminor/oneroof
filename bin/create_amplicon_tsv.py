#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "polars",
# ]
# ///

"""
Summarize amplicon coverage from stats files and primer position data.

Supports two input modes for position data:
  1. --primer-tsv: Use primer_pairs.tsv from prepare_primers.py (preferred)
  2. --bed: Use BED file with primer coordinates (legacy)

The primer TSV approach is cleaner as it already contains pre-computed
amplicon positions and works regardless of primer naming conventions.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import polars as pl


def _load_positions_from_tsv(tsv_path: str) -> pl.LazyFrame:
    """
    Load amplicon positions from primer_pairs.tsv.

    The TSV has columns: amplicon_name, fwd_sequence, rev_sequence, chrom,
    amplicon_start, amplicon_end

    Returns a LazyFrame with: amplicon_name, start_pos, end_pos
    """
    return pl.scan_csv(tsv_path, separator="\t").select(
        pl.col("amplicon_name"),
        pl.col("amplicon_start").alias("start_pos"),
        pl.col("amplicon_end").alias("end_pos"),
    )


def _load_positions_from_bed(
    bed_path: str, fwd_suffix: str, rev_suffix: str
) -> pl.LazyFrame:
    """
    Load amplicon positions from BED file (legacy approach).

    Extracts base amplicon name from primer names and computes positions
    from LEFT/RIGHT primer coordinates.

    Returns a LazyFrame with: amplicon_name, start_pos, end_pos
    """
    bed = (
        pl.scan_csv(
            bed_path,
            separator="\t",
            has_header=False,
            new_columns=["chrom", "start", "end", "name", "score", "strand"],
        )
        .with_columns(
            # Extract base amplicon name by removing suffixes and indices
            pl.col("name")
            .str.replace(fwd_suffix, "")
            .str.replace(rev_suffix, "")
            .str.replace(r"_splice\d+$", "")
            .str.replace(r"-\d+$", "")
            .alias("amplicon_name"),
            # Flag primer direction
            pl.col("name").str.contains(fwd_suffix).alias("is_fwd"),
            pl.col("name").str.contains(rev_suffix).alias("is_rev"),
        )
        .filter(pl.col("amplicon_name").is_not_null())
    )

    # Aggregate to get amplicon span
    return bed.group_by("amplicon_name").agg(
        pl.col("start").filter(pl.col("is_fwd")).min().alias("start_pos"),
        pl.col("end").filter(pl.col("is_rev")).max().alias("end_pos"),
    )


def _parse_stats_files(stats_pattern: str) -> pl.LazyFrame:
    """
    Parse amplicon stats files and extract sample/amplicon info.

    Stats files have columns including 'file' and 'num_seqs'.
    The 'file' column contains filenames like:
      "SAMPLE.amplicon_name.no_downsampling.fasta.gz"

    Returns a LazyFrame with: sample_name, amplicon_name, reads
    """
    return (
        pl.scan_csv(stats_pattern, separator="\t", glob=True)
        .with_columns(
            # Extract sample name: first segment before '.'
            pl.col("file").str.extract(r"^([^.]+)", group_index=1).alias("sample_name"),
            # Extract amplicon name: second segment (between first and second '.')
            pl.col("file")
            .str.extract(r"^[^.]+\.([^.]+)", group_index=1)
            .alias("amplicon_name"),
        )
        .select(
            "sample_name",
            "amplicon_name",
            pl.col("num_seqs").alias("reads"),
        )
    )


def _create_summary(
    stats_lf: pl.LazyFrame,
    positions_lf: pl.LazyFrame,
    output_path: str,
) -> None:
    """
    Join stats with positions and write summary TSV.
    """
    result = (
        stats_lf.join(positions_lf, on="amplicon_name", how="left")
        .select(
            "sample_name",
            "amplicon_name",
            pl.col("start_pos").cast(pl.String).fill_null("NA"),
            pl.col("end_pos").cast(pl.String).fill_null("NA"),
            "reads",
        )
        .sort("sample_name", "amplicon_name")
        .collect()
    )

    result.write_csv(output_path, separator="\t")
    print(f"Wrote {len(result)} rows to {output_path}", file=sys.stderr)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # Position data input (mutually exclusive)
    pos_group = parser.add_mutually_exclusive_group(required=True)
    pos_group.add_argument(
        "--primer-tsv",
        type=Path,
        help="Primer pairs TSV from prepare_primers.py (preferred)",
    )
    pos_group.add_argument(
        "--bed",
        type=Path,
        help="BED file with primer coordinates (legacy)",
    )

    # Stats input
    parser.add_argument(
        "--pattern",
        default="stats_*.tsv",
        help="Glob pattern for stats files (default: stats_*.tsv)",
    )

    # Output
    parser.add_argument(
        "--output",
        default="amplicon_summary.tsv",
        help="Output TSV file (default: amplicon_summary.tsv)",
    )

    # BED-specific options
    parser.add_argument(
        "--fwd-suffix",
        default="_LEFT",
        help="Forward primer suffix for BED parsing (default: _LEFT)",
    )
    parser.add_argument(
        "--rev-suffix",
        default="_RIGHT",
        help="Reverse primer suffix for BED parsing (default: _RIGHT)",
    )

    args = parser.parse_args()

    # Load position data
    if args.primer_tsv:
        assert args.primer_tsv.is_file(), f"Primer TSV not found: {args.primer_tsv}"
        positions_lf = _load_positions_from_tsv(str(args.primer_tsv))
    else:
        assert args.bed.is_file(), f"BED file not found: {args.bed}"
        positions_lf = _load_positions_from_bed(
            str(args.bed),
            args.fwd_suffix,
            args.rev_suffix,
        )

    # Parse stats files
    stats_lf = _parse_stats_files(args.pattern)

    # Create and write summary
    _create_summary(stats_lf, positions_lf, args.output)


if __name__ == "__main__":
    main()
