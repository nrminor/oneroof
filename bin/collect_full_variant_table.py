#!/usr/bin/env python3

# /// script
# requires-python = ">= 3.10"
# dependencies = [
#     "polars>=1.22.0,<2",
#     "loguru>=0.7.0,<1",
#     "biopython>=1.80,<2",
# ]
# ///

"""
Collect and enrich variant data from all samples into a final, queryable table.

Combines per-sample variant effect files (from SnpSift) into a single table with
derived columns for easier exploration in Jupyter, DuckDB, or similar tools.
"""

import argparse
import re
from pathlib import Path

import polars as pl
from Bio.Data.IUPACData import protein_letters_3to1_extended
from loguru import logger

DEFAULT_CONSENSUS_THRESHOLD = 0.8
VARIANT_EFFECTS_SUFFIX = "_variant_effects.tsv"

# Single source of truth for SnpSift column definitions.
# Each tuple: (original_header_name, final_column_name, polars_dtype)
# This ensures consistent type inference across all input files,
# avoiding errors when columns are empty in some files but numeric in others.
SNPSIFT_COLUMN_DEFS = [
    ("CHROM", "chrom", pl.String),
    ("REF", "ref", pl.String),
    ("POS", "pos", pl.Int64),
    ("ALT", "alt", pl.String),
    ("AF", "af", pl.Float64),
    ("AC", "ac", pl.Int64),
    ("DP", "dp", pl.Int64),
    ("GEN[0].REF_DP", "ref_dp", pl.Int64),
    ("GEN[0].ALT_DP", "alt_dp", pl.Int64),
    ("GEN[0].ALT_FREQ", "alt_freq", pl.Float64),
    ("MQ", "mq", pl.Float64),
    ("ANN[0].GENE", "gene", pl.String),
    ("ANN[0].EFFECT", "effect", pl.String),
    ("ANN[0].HGVS_P", "hgvs_p", pl.String),
    ("ANN[0].CDS_POS", "cds_pos", pl.Int64),
    ("ANN[0].AA_POS", "aa_pos", pl.Int64),
]

# Derived constants from the single source of truth
SNPSIFT_SCHEMA = {orig: dtype for orig, _, dtype in SNPSIFT_COLUMN_DEFS}
SNPSIFT_RENAME_MAP = {orig: final for orig, final, _ in SNPSIFT_COLUMN_DEFS}
SNPSIFT_COLUMNS = [final for _, final, _ in SNPSIFT_COLUMN_DEFS]

FINAL_COLUMNS = [
    "sample_id",
    "chrom",
    "pos",
    "ref",
    "alt",
    "af",
    "dp",
    "ref_dp",
    "alt_dp",
    "mq",
    "gene",
    "effect",
    "hgvs_p",
    "cds_pos",
    "aa_pos",
    "variant_id",
    "aa_change",
    "mutation_type",
    "is_consensus",
    "sample_count",
    "is_shared",
]

AA_THREE_TO_ONE = {**protein_letters_3to1_extended, "Ter": "*"}

HGVS_PATTERN = re.compile(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|=)")


def parse_command_line_args() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Collect variant data from all samples into a final enriched table.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--input-dir",
        "-i",
        type=Path,
        required=True,
        help="Directory containing *_variant_effects.tsv files",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=Path,
        required=True,
        help="Output basename (produces .tsv and .parquet files)",
    )
    parser.add_argument(
        "--consensus-threshold",
        "-t",
        type=float,
        default=DEFAULT_CONSENSUS_THRESHOLD,
        help="Allele frequency threshold for consensus calls",
    )
    return parser.parse_args()


def _three_letter_to_one(aa: str) -> str:
    """Convert three-letter amino acid code to one-letter."""
    return AA_THREE_TO_ONE.get(aa.capitalize(), "X")


def _parse_hgvs_to_short(hgvs: str | None) -> str | None:
    """
    Parse HGVS protein notation to short form.

    Examples:
        p.Asp614Gly -> D614G
        p.Phe924= -> F924F (synonymous)
        None -> None
    """
    if hgvs is None:
        return None

    match = HGVS_PATTERN.match(hgvs)
    if match is None:
        return None

    ref_aa = _three_letter_to_one(match.group(1))
    position = match.group(2)
    alt_part = match.group(3)

    alt_aa = ref_aa if alt_part == "=" else _three_letter_to_one(alt_part)

    return f"{ref_aa}{position}{alt_aa}"


def discover_input_files(input_dir: Path) -> list[tuple[str, Path]]:
    """
    Discover variant effect TSV files and extract sample IDs.

    Returns list of (sample_id, file_path) tuples.
    """
    files = sorted(input_dir.glob(f"*{VARIANT_EFFECTS_SUFFIX}"))

    assert len(files) > 0, f"No *{VARIANT_EFFECTS_SUFFIX} files found in {input_dir}"

    result = []
    for f in files:
        sample_id = f.name.removesuffix(VARIANT_EFFECTS_SUFFIX)
        result.append((sample_id, f))

    logger.info(f"Discovered {len(result)} sample files in {input_dir}")
    return result


def load_single_file(sample_id: str, file_path: Path) -> pl.LazyFrame:
    """Load a single variant effects TSV file as a LazyFrame."""
    return (
        pl.scan_csv(
            file_path,
            separator="\t",
            has_header=False,
            skip_rows=1,
            schema=SNPSIFT_SCHEMA,
            null_values=["", "."],
        )
        .rename(SNPSIFT_RENAME_MAP)
        .select(SNPSIFT_COLUMNS)
        .with_columns(pl.lit(sample_id).alias("sample_id"))
    )


def load_and_concat(files: list[tuple[str, Path]]) -> pl.LazyFrame:
    """Load all variant effect files and concatenate into single LazyFrame."""
    frames = [load_single_file(sample_id, path) for sample_id, path in files]
    return pl.concat(frames, how="vertical_relaxed")


def variant_id_expr() -> pl.Expr:
    """Expression to create variant_id column."""
    return pl.concat_str(
        [
            pl.col("chrom"),
            pl.lit(":"),
            pl.col("pos").cast(pl.String),
            pl.lit(":"),
            pl.col("ref"),
            pl.lit(">"),
            pl.col("alt"),
        ],
    ).alias("variant_id")


def mutation_type_expr() -> pl.Expr:
    """Expression to classify mutation type based on ref/alt lengths."""
    ref_len = pl.col("ref").str.len_chars()
    alt_len = pl.col("alt").str.len_chars()

    return (
        pl.when(ref_len == alt_len)
        .then(pl.when(ref_len == 1).then(pl.lit("SNP")).otherwise(pl.lit("MNP")))
        .when(ref_len > alt_len)
        .then(pl.lit("deletion"))
        .otherwise(pl.lit("insertion"))
        .alias("mutation_type")
    )


def is_consensus_expr(threshold: float) -> pl.Expr:
    """Expression to determine if variant is at consensus level."""
    return (pl.col("af") >= threshold).alias("is_consensus")


def aa_change_expr() -> pl.Expr:
    """
    Expression to create human-readable amino acid change.

    Uses map_elements for HGVS parsing since regex isn't available in expressions.
    """
    return (
        pl.when(pl.col("hgvs_p").is_not_null() & (pl.col("hgvs_p") != ""))
        .then(
            pl.concat_str(
                [
                    pl.col("gene"),
                    pl.lit(":"),
                    pl.col("hgvs_p").map_elements(
                        _parse_hgvs_to_short,
                        return_dtype=pl.Utf8,
                    ),
                ],
            ),
        )
        .otherwise(pl.lit(None))
        .alias("aa_change")
    )


def add_derived_columns(lf: pl.LazyFrame, threshold: float) -> pl.LazyFrame:
    """Add derived columns to the LazyFrame."""
    return lf.with_columns(
        variant_id_expr(),
        mutation_type_expr(),
        is_consensus_expr(threshold),
        aa_change_expr(),
    )


def add_cross_sample_metrics(lf: pl.LazyFrame) -> pl.LazyFrame:
    """Add columns that depend on cross-sample aggregation."""
    return lf.with_columns(
        pl.col("variant_id").count().over("variant_id").alias("sample_count"),
    ).with_columns(
        (pl.col("sample_count") > 1).alias("is_shared"),
    )


def select_and_sort(lf: pl.LazyFrame) -> pl.LazyFrame:
    """Select final columns and sort output."""
    return lf.select(FINAL_COLUMNS).sort(["chrom", "pos", "sample_id"])


def write_outputs(df: pl.DataFrame, output_basename: Path) -> None:
    """Write DataFrame to TSV and Parquet files.

    Note: We collect once before this function rather than using sink_csv/sink_parquet
    twice, which would execute the lazy query twice.
    """
    tsv_path = output_basename.with_suffix(".tsv")
    parquet_path = output_basename.with_suffix(".parquet")

    df.write_csv(tsv_path, separator="\t")
    logger.info(f"Wrote TSV output to {tsv_path}")

    df.write_parquet(parquet_path, compression="zstd")
    logger.info(f"Wrote Parquet output to {parquet_path}")


def main() -> None:
    """Collect variant data from all samples into a final enriched table."""
    args = parse_command_line_args()

    input_files = discover_input_files(args.input_dir)

    result: pl.DataFrame = (
        load_and_concat(input_files)
        .pipe(add_derived_columns, args.consensus_threshold)
        .pipe(add_cross_sample_metrics)
        .pipe(select_and_sort)
        .collect()
    )

    logger.info(
        f"Generated table with {len(result)} rows across {result['sample_id'].n_unique()} samples",
    )

    write_outputs(result, args.output)


if __name__ == "__main__":
    main()
