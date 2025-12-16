#!/usr/bin/env python

import argparse
import glob
import re

import pandas as pd


def extract_sample_name(filename):
    base_name = filename.replace(".QIAseq_", "_QIAseq_")
    if "_QIAseq_" in base_name:
        return base_name.split("_QIAseq_")[0]
    return base_name.split(".")[0]


def extract_amplicon_name(filename):
    match = re.search(r"(QIAseq_[^._\s]+)", filename)
    if match:
        return match.group(1)
    return "unknown"


def get_amplicon_positions(bed_df, amplicon_name):
    # Match base form (e.g., QIAseq_92 from QIAseq_92_splice9)
    base_amplicon = re.sub(r"[_-](splice\d+|\d+)$", "", amplicon_name)

    # Find matching LEFT and RIGHT primers (include splices)
    left_primers = bed_df[
        bed_df["name"].str.contains(f"{base_amplicon}_LEFT", na=False)
    ]
    right_primers = bed_df[
        bed_df["name"].str.contains(f"{base_amplicon}_RIGHT", na=False)
    ]

    if not left_primers.empty and not right_primers.empty:
        start_pos = left_primers["start"].min()
        end_pos = right_primers["end"].max()
        return start_pos, end_pos

    return None, None


def main(bed_file, output_file, stats_pattern):
    # Read stats files
    stats_files = glob.glob(stats_pattern)
    if not stats_files:
        raise FileNotFoundError(
            f"No stats files found matching pattern: {stats_pattern}",
        )

    all_stats_data = [pd.read_csv(f, sep="\t") for f in stats_files]
    combined_stats_df = pd.concat(all_stats_data, ignore_index=True)

    # Read BED file
    bed_columns = ["chrom", "start", "end", "name", "score", "strand"]
    bed_df = pd.read_csv(bed_file, sep="\t", names=bed_columns, header=None)

    output_data = []
    for _, row in combined_stats_df.iterrows():
        filename = row["file"]
        sample_name = extract_sample_name(filename)
        amplicon_name = extract_amplicon_name(filename)
        reads = row["num_seqs"]
        start_pos, end_pos = get_amplicon_positions(bed_df, amplicon_name)

        output_data.append(
            {
                "sample_name": sample_name,
                "amplicon_name": amplicon_name,
                "start_pos": start_pos if start_pos is not None else "NA",
                "end_pos": end_pos if end_pos is not None else "NA",
                "reads": reads,
            },
        )

    output_df = pd.DataFrame(output_data)
    output_df.to_csv(output_file, sep="\t", index=False)


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
