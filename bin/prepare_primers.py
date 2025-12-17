#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = ["polars", "loguru"]
# ///

"""
Prepare primer information for amplicon sequencing analysis.

Consolidates the functionality of:
  - validate_primer_bed.py (BED validation and coordinate orientation)
  - resplice_primers.py (spike-in primer combination handling)
  - split_primer_combos.py (amplicon partitioning)
  - make_primer_patterns.py (pattern generation)
  - bedtools getfasta (sequence extraction)

Supports two input modes:
  1. BED + Reference FASTA: Full validation, resplicing, sequence extraction
  2. Primer TSV: Direct passthrough with validation

Outputs:
  - primer_pairs.tsv: amplicon_name, fwd_sequence, rev_sequence, chrom, start, end
  - respliced.bed (optional): The respliced BED file for reference/visualization

Example usage:

    # From BED + FASTA (full pipeline)
    prepare_primers.py \\
        --input-bed primers.bed \\
        --reference reference.fasta \\
        --output-tsv primer_pairs.tsv \\
        --output-bed respliced.bed

    # From existing primer TSV (passthrough with validation)
    prepare_primers.py \\
        --input-tsv primers.tsv \\
        --output-tsv primer_pairs.tsv
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from itertools import product
from pathlib import Path
from typing import TYPE_CHECKING

import polars as pl
from loguru import logger

if TYPE_CHECKING:
    from collections.abc import Sequence


@dataclass(frozen=True, slots=True)
class PrimerPair:
    """A validated forward/reverse primer pair with sequences and coordinates."""

    amplicon_name: str
    fwd_sequence: str
    rev_sequence: str
    chrom: str
    amplicon_start: int
    amplicon_end: int

    def __post_init__(self) -> None:
        """Validate primer pair on construction."""
        assert self.amplicon_name, "amplicon_name cannot be empty"
        assert self.fwd_sequence, "fwd_sequence cannot be empty"
        assert self.rev_sequence, "rev_sequence cannot be empty"
        assert self.chrom, "chrom cannot be empty"
        assert self.amplicon_start >= 0, f"Invalid start: {self.amplicon_start}"
        assert self.amplicon_end > self.amplicon_start, (
            f"End ({self.amplicon_end}) must be > start ({self.amplicon_start})"
        )


@dataclass(frozen=True, slots=True)
class BedRecord:
    """A single BED file record."""

    chrom: str
    start: int
    end: int
    name: str
    score: int | str
    strand: str

    def __post_init__(self) -> None:
        """Validate BED record on construction."""
        assert self.chrom, "chrom cannot be empty"
        assert self.start >= 0, f"Invalid start position: {self.start}"
        assert self.end > self.start, f"End ({self.end}) must be > start ({self.start})"
        assert self.name, "name cannot be empty"
        assert self.strand in ("+", "-"), f"Invalid strand: {self.strand}"


MIN_READS_THRESHOLD = 10
MIN_BED_FIELDS = 6
PRIMER_PAIR_SIZE = 2


def parse_fasta(fasta_path: Path) -> dict[str, str]:
    """
    Parse a FASTA file into a dictionary of {sequence_id: sequence}.

    Args:
        fasta_path: Path to the FASTA file

    Returns:
        Dictionary mapping sequence IDs to their sequences (uppercase)
    """
    assert fasta_path.is_file(), f"FASTA file not found: {fasta_path}"

    sequences: dict[str, str] = {}
    current_id: str | None = None
    current_seq: list[str] = []

    def _save_current() -> None:
        nonlocal current_id, current_seq
        if current_id is None:
            return
        seq = "".join(current_seq).upper()
        assert seq, f"Empty sequence for ID: {current_id}"
        sequences[current_id] = seq
        current_id = None
        current_seq = []

    with fasta_path.open(encoding="utf8") as handle:
        for line_num, line in enumerate(handle, 1):
            line = line.strip()  # noqa: PLW2901
            if not line:
                continue

            if line.startswith(">"):
                _save_current()
                current_id = line[1:].split()[0]
                assert current_id, f"Empty sequence ID at line {line_num}"
                continue

            assert current_id is not None, f"Sequence data before header at line {line_num}"
            current_seq.append(line)

        _save_current()

    assert sequences, f"No sequences found in {fasta_path}"
    return sequences


def extract_sequence(
    fasta_dict: dict[str, str],
    chrom: str,
    start: int,
    end: int,
) -> str:
    """
    Extract a subsequence from the parsed FASTA dictionary.

    Args:
        fasta_dict: Dictionary from parse_fasta()
        chrom: Chromosome/sequence ID
        start: 0-based start position
        end: 0-based end position (exclusive, BED-style)

    Returns:
        The extracted sequence (uppercase)
    """
    assert chrom in fasta_dict, (
        f"Chromosome '{chrom}' not found. Available: {list(fasta_dict.keys())[:5]}..."
    )

    seq = fasta_dict[chrom]
    assert start >= 0, f"Start position cannot be negative: {start}"
    assert end <= len(seq), f"End ({end}) exceeds sequence length ({len(seq)})"
    assert start < end, f"Start ({start}) must be < end ({end})"

    return seq[start:end]


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    assert seq, "Cannot reverse complement empty sequence"
    complement = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")
    return seq.translate(complement)[::-1]


def parse_bed_file(bed_path: Path) -> list[BedRecord]:
    """
    Parse a BED file into a list of BedRecord objects.

    Args:
        bed_path: Path to the BED file

    Returns:
        List of validated BedRecord objects
    """
    assert bed_path.is_file(), f"BED file not found: {bed_path}"

    records = []
    with bed_path.open(encoding="utf8") as handle:
        for line_num, line in enumerate(handle, 1):
            line = line.strip()  # noqa: PLW2901
            if not line or line.startswith("#"):
                continue

            fields = line.split("\t")
            assert len(fields) >= MIN_BED_FIELDS, (
                f"BED line {line_num} has {len(fields)} fields, expected >= {MIN_BED_FIELDS}: {line}"
            )

            chrom, start_str, end_str, name, score, strand = fields[:MIN_BED_FIELDS]

            start = int(start_str)
            end = int(end_str)

            # Orient coordinates if needed (start should be < end)
            if start > end:
                start, end = end, start
                strand = "-" if strand == "+" else "+"

            records.append(
                BedRecord(
                    chrom=chrom,
                    start=start,
                    end=end,
                    name=name,
                    score=score,
                    strand=strand,
                ),
            )

    assert records, f"No valid records found in {bed_path}"
    return records


def has_valid_suffix(name: str, fwd_suffix: str, rev_suffix: str) -> bool:
    """Check if a primer name contains a valid forward or reverse suffix."""
    return fwd_suffix in name or rev_suffix in name


def is_forward_primer(name: str, fwd_suffix: str) -> bool:
    """Check if a primer name indicates a forward primer."""
    return fwd_suffix in name


def is_reverse_primer(name: str, rev_suffix: str) -> bool:
    """Check if a primer name indicates a reverse primer."""
    return rev_suffix in name


def extract_amplicon_base(name: str, fwd_suffix: str, rev_suffix: str) -> str:
    """
    Extract the base amplicon name from a primer name.

    Removes the forward/reverse suffix and any spike-in index.

    Examples:
        "QIAseq_21_LEFT" -> "QIAseq_21"
        "QIAseq_171-2_LEFT" -> "QIAseq_171"
        "QIAseq_171_LEFT_splice1" -> "QIAseq_171"
    """
    # First remove any splice suffix
    base = name.split("_splice")[0]

    # Remove forward/reverse suffix
    if fwd_suffix in base:
        base = base.split(fwd_suffix)[0]
    elif rev_suffix in base:
        base = base.split(rev_suffix)[0]

    # Remove spike-in index (trailing -N)
    if "-" in base:
        parts = base.rsplit("-", 1)
        if parts[-1].isdigit():
            base = parts[0]

    assert base, f"Could not extract amplicon base from: {name}"
    return base


def validate_primer_names(
    records: Sequence[BedRecord],
    fwd_suffix: str,
    rev_suffix: str,
) -> None:
    """
    Validate that all primer names have required suffixes.

    Raises:
        ValueError: If any primer lacks a valid suffix
    """
    invalid = [r.name for r in records if not has_valid_suffix(r.name, fwd_suffix, rev_suffix)]

    if not invalid:
        return

    msg = (
        f"Found {len(invalid)} primer(s) without required suffix "
        f"({fwd_suffix} or {rev_suffix}): {invalid}"
    )
    raise ValueError(msg)


def validate_primer_pairs(
    records: Sequence[BedRecord],
    fwd_suffix: str,
    rev_suffix: str,
) -> None:
    """
    Validate that each amplicon has at least one forward and one reverse primer.

    Raises:
        ValueError: If any amplicon is missing a forward or reverse primer
    """
    amplicons_with_fwd: set[str] = set()
    amplicons_with_rev: set[str] = set()

    for record in records:
        base = extract_amplicon_base(record.name, fwd_suffix, rev_suffix)
        if is_forward_primer(record.name, fwd_suffix):
            amplicons_with_fwd.add(base)
        if is_reverse_primer(record.name, rev_suffix):
            amplicons_with_rev.add(base)

    all_amplicons = amplicons_with_fwd | amplicons_with_rev
    missing_fwd = all_amplicons - amplicons_with_fwd
    missing_rev = all_amplicons - amplicons_with_rev

    errors = []
    if missing_fwd:
        errors.append(f"Amplicons missing forward primer: {sorted(missing_fwd)}")
    if missing_rev:
        errors.append(f"Amplicons missing reverse primer: {sorted(missing_rev)}")

    if errors:
        raise ValueError("\n".join(errors))


def validate_bed_records(
    records: Sequence[BedRecord],
    fwd_suffix: str,
    rev_suffix: str,
) -> list[BedRecord]:
    """
    Validate BED records and return them unchanged if valid.

    Args:
        records: List of BedRecord objects
        fwd_suffix: Forward primer suffix
        rev_suffix: Reverse primer suffix

    Returns:
        The same records if validation passes

    Raises:
        ValueError: If validation fails
    """
    validate_primer_names(records, fwd_suffix, rev_suffix)
    validate_primer_pairs(records, fwd_suffix, rev_suffix)
    return list(records)


def group_by_amplicon(
    records: Sequence[BedRecord],
    fwd_suffix: str,
    rev_suffix: str,
) -> dict[str, list[BedRecord]]:
    """
    Group BED records by their base amplicon name.

    Args:
        records: List of BedRecord objects
        fwd_suffix: Forward primer suffix
        rev_suffix: Reverse primer suffix

    Returns:
        Dictionary mapping amplicon names to their primer records
    """
    groups: dict[str, list[BedRecord]] = {}

    for record in records:
        base = extract_amplicon_base(record.name, fwd_suffix, rev_suffix)
        groups.setdefault(base, []).append(record)

    # Validate each group has at least 2 primers
    for amplicon, group_records in groups.items():
        assert len(group_records) >= PRIMER_PAIR_SIZE, (
            f"Amplicon {amplicon} has only {len(group_records)} primer(s), need >= {PRIMER_PAIR_SIZE}"
        )

    return groups


def strip_existing_indices(name: str) -> str:
    """
    Remove any existing spike-in indices from a primer name.

    Examples:
        "QIAseq_171-2_LEFT" -> "QIAseq_171_LEFT"
        "QIAseq_21_LEFT" -> "QIAseq_21_LEFT"
    """
    # Find the suffix position
    for suffix in ("_LEFT", "_RIGHT"):
        if suffix not in name:
            continue

        prefix, rest = name.split(suffix, 1)

        # Remove trailing -N from prefix
        if "-" in prefix:
            parts = prefix.rsplit("-", 1)
            if parts[-1].isdigit():
                prefix = parts[0]

        return prefix + suffix + rest

    return name


def assign_indices_to_group(
    records: Sequence[BedRecord],
    fwd_suffix: str,
    rev_suffix: str,
) -> list[tuple[BedRecord, str]]:
    """
    Assign sequential indices to primers within an amplicon group.

    Forward and reverse primers are indexed separately.

    Args:
        records: Primer records for a single amplicon
        fwd_suffix: Forward primer suffix
        rev_suffix: Reverse primer suffix

    Returns:
        List of (record, indexed_name) tuples
    """
    fwd_records = [r for r in records if is_forward_primer(r.name, fwd_suffix)]
    rev_records = [r for r in records if is_reverse_primer(r.name, rev_suffix)]

    assert fwd_records, "No forward primers in group"
    assert rev_records, "No reverse primers in group"

    result: list[tuple[BedRecord, str]] = []

    for idx, record in enumerate(fwd_records, 1):
        clean_name = strip_existing_indices(record.name)
        indexed_name = f"{clean_name}-{idx}"
        result.append((record, indexed_name))

    for idx, record in enumerate(rev_records, 1):
        clean_name = strip_existing_indices(record.name)
        indexed_name = f"{clean_name}-{idx}"
        result.append((record, indexed_name))

    return result


def generate_splice_combinations(
    indexed_primers: Sequence[tuple[BedRecord, str]],
    fwd_suffix: str,
    rev_suffix: str,
) -> list[tuple[BedRecord, str]]:
    """
    Generate all valid forward/reverse primer combinations.

    For amplicons with multiple forward or reverse primers (spike-ins),
    creates all possible pairings with unique splice identifiers.

    Args:
        indexed_primers: List of (record, indexed_name) tuples
        fwd_suffix: Forward primer suffix
        rev_suffix: Reverse primer suffix

    Returns:
        List of (record, spliced_name) tuples for all combinations
    """
    fwd_primers = [(r, n) for r, n in indexed_primers if is_forward_primer(n, fwd_suffix)]
    rev_primers = [(r, n) for r, n in indexed_primers if is_reverse_primer(n, rev_suffix)]

    assert fwd_primers, "No forward primers found"
    assert rev_primers, "No reverse primers found"

    # Simple case: exactly one of each - no splicing needed
    if len(fwd_primers) == 1 and len(rev_primers) == 1:
        return list(indexed_primers)

    # Generate all combinations
    all_pairs = list(product(fwd_primers, rev_primers))
    assert all_pairs, "No valid primer pairs generated"

    result: list[tuple[BedRecord, str]] = []

    for splice_idx, ((fwd_rec, fwd_name), (rev_rec, rev_name)) in enumerate(all_pairs, 1):
        # Strip the index suffix and add splice identifier
        fwd_base = fwd_name.rsplit("-", 1)[0]
        rev_base = rev_name.rsplit("-", 1)[0]

        fwd_spliced = f"{fwd_base}_splice{splice_idx}"
        rev_spliced = f"{rev_base}_splice{splice_idx}"

        result.append((fwd_rec, fwd_spliced))
        result.append((rev_rec, rev_spliced))

    return result


def resplice_all_amplicons(
    records: Sequence[BedRecord],
    fwd_suffix: str,
    rev_suffix: str,
) -> list[tuple[BedRecord, str]]:
    """
    Process all amplicons through the resplicing pipeline.

    Args:
        records: All validated BED records
        fwd_suffix: Forward primer suffix
        rev_suffix: Reverse primer suffix

    Returns:
        List of (record, final_name) tuples for all respliced primers
    """
    groups = group_by_amplicon(records, fwd_suffix, rev_suffix)

    all_respliced: list[tuple[BedRecord, str]] = []

    for _amplicon, group_records in sorted(groups.items()):
        indexed = assign_indices_to_group(group_records, fwd_suffix, rev_suffix)
        spliced = generate_splice_combinations(indexed, fwd_suffix, rev_suffix)
        all_respliced.extend(spliced)

    assert all_respliced, "No primers after resplicing"
    return all_respliced


def _group_respliced_by_pair(
    respliced: Sequence[tuple[BedRecord, str]],
    fwd_suffix: str,
    rev_suffix: str,
) -> dict[str, list[tuple[BedRecord, str]]]:
    """
    Group respliced primers by their pair identifier.

    Primers with matching splice indices belong to the same pair.

    Args:
        respliced: List of (record, final_name) tuples
        fwd_suffix: Forward primer suffix
        rev_suffix: Reverse primer suffix

    Returns:
        Dictionary mapping pair identifiers to their primer tuples
    """

    def get_pair_key(name: str) -> str:
        """Extract the pair grouping key from a primer name."""
        # Extract splice suffix if present
        if "_splice" not in name:
            return extract_amplicon_base(name, fwd_suffix, rev_suffix)

        base = name.split(fwd_suffix)[0] if fwd_suffix in name else name.split(rev_suffix)[0]
        splice = name.split("_splice")[1]
        return f"{base}_splice{splice}"

    groups: dict[str, list[tuple[BedRecord, str]]] = {}

    for record, name in respliced:
        key = get_pair_key(name)
        groups.setdefault(key, []).append((record, name))

    # Validate each group has exactly 2 primers (one fwd, one rev)
    for key, group in groups.items():
        assert len(group) == PRIMER_PAIR_SIZE, (
            f"Pair {key} has {len(group)} primers, expected {PRIMER_PAIR_SIZE}"
        )

        names = [n for _, n in group]
        has_fwd = any(is_forward_primer(n, fwd_suffix) for n in names)
        has_rev = any(is_reverse_primer(n, rev_suffix) for n in names)
        assert has_fwd, f"Pair {key} missing forward primer"
        assert has_rev, f"Pair {key} missing reverse primer"

    return groups


def extract_primer_pair(
    pair_key: str,
    primers: Sequence[tuple[BedRecord, str]],
    fasta_dict: dict[str, str],
    fwd_suffix: str,
    rev_suffix: str,
) -> PrimerPair:
    """
    Extract sequences and create a PrimerPair from a pair of primers.

    Args:
        pair_key: The pair identifier (used as amplicon name)
        primers: Exactly 2 (record, name) tuples
        fasta_dict: Parsed reference FASTA
        fwd_suffix: Forward primer suffix
        rev_suffix: Reverse primer suffix

    Returns:
        A validated PrimerPair object
    """
    assert len(primers) == PRIMER_PAIR_SIZE, (
        f"Expected {PRIMER_PAIR_SIZE} primers, got {len(primers)}"
    )

    # Separate forward and reverse
    fwd_primer = next((r, n) for r, n in primers if is_forward_primer(n, fwd_suffix))
    rev_primer = next((r, n) for r, n in primers if is_reverse_primer(n, rev_suffix))

    fwd_record, _ = fwd_primer
    rev_record, _ = rev_primer

    # Validate same chromosome
    assert fwd_record.chrom == rev_record.chrom, (
        f"Primers on different chromosomes: {fwd_record.chrom} vs {rev_record.chrom}"
    )

    # Extract sequences
    fwd_seq = extract_sequence(
        fasta_dict,
        fwd_record.chrom,
        fwd_record.start,
        fwd_record.end,
    )
    rev_seq = extract_sequence(
        fasta_dict,
        rev_record.chrom,
        rev_record.start,
        rev_record.end,
    )

    # Reverse complement the reverse primer
    rev_seq_rc = reverse_complement(rev_seq)

    # Compute amplicon span
    amplicon_start = min(fwd_record.start, rev_record.start)
    amplicon_end = max(fwd_record.end, rev_record.end)

    return PrimerPair(
        amplicon_name=pair_key,
        fwd_sequence=fwd_seq,
        rev_sequence=rev_seq_rc,
        chrom=fwd_record.chrom,
        amplicon_start=amplicon_start,
        amplicon_end=amplicon_end,
    )


def generate_primer_pairs(
    respliced: Sequence[tuple[BedRecord, str]],
    fasta_dict: dict[str, str],
    fwd_suffix: str,
    rev_suffix: str,
) -> list[PrimerPair]:
    """
    Generate PrimerPair objects from respliced primers.

    Args:
        respliced: List of (record, final_name) tuples
        fasta_dict: Parsed reference FASTA
        fwd_suffix: Forward primer suffix
        rev_suffix: Reverse primer suffix

    Returns:
        List of PrimerPair objects
    """
    groups = _group_respliced_by_pair(respliced, fwd_suffix, rev_suffix)

    pairs = [
        extract_primer_pair(key, primers, fasta_dict, fwd_suffix, rev_suffix)
        for key, primers in sorted(groups.items())
    ]

    assert pairs, "No primer pairs generated"
    return pairs


def _primer_pairs_to_dataframe(pairs: Sequence[PrimerPair]) -> pl.LazyFrame:
    """Convert a list of PrimerPair objects to a Polars LazyFrame."""
    assert pairs, "Cannot create DataFrame from empty primer pairs"

    return pl.LazyFrame(
        [
            {
                "amplicon_name": p.amplicon_name,
                "fwd_sequence": p.fwd_sequence,
                "rev_sequence": p.rev_sequence,
                "chrom": p.chrom,
                "amplicon_start": p.amplicon_start,
                "amplicon_end": p.amplicon_end,
            }
            for p in pairs
        ],
    )


def _respliced_to_bed_dataframe(
    respliced: Sequence[tuple[BedRecord, str]],
) -> pl.LazyFrame:
    """
    Convert respliced primers to a BED-format LazyFrame.

    Args:
        respliced: List of (record, final_name) tuples

    Returns:
        LazyFrame with BED columns, sorted by position
    """
    assert respliced, "Cannot create BED from empty respliced primers"

    rows = [
        {
            "chrom": record.chrom,
            "start": record.start,
            "end": record.end,
            "name": name,
            "score": record.score,
            "strand": record.strand,
        }
        for record, name in respliced
    ]

    return pl.LazyFrame(rows).sort("chrom", "start", "end")


VALID_BASES = frozenset("ACGTRYSWKMBDHVN")


def _validate_sequence(seq: str, context: str) -> None:
    """Validate that a sequence contains only valid nucleotide characters."""
    invalid = set(seq.upper()) - VALID_BASES
    assert not invalid, f"Invalid characters in {context}: {invalid}"


def validate_primer_tsv(tsv_path: Path) -> pl.LazyFrame:
    """
    Validate and load a primer TSV file.

    Expected columns: amplicon_name, fwd_sequence, rev_sequence
    Optional columns: chrom, amplicon_start, amplicon_end

    Args:
        tsv_path: Path to the TSV file

    Returns:
        Validated LazyFrame with standardized columns
    """
    assert tsv_path.is_file(), f"TSV file not found: {tsv_path}"

    # Read eagerly first for validation, then convert to lazy
    df = pl.read_csv(tsv_path, separator="\t")
    assert not df.is_empty(), f"TSV file is empty: {tsv_path}"

    # Column name normalization mapping
    column_aliases = {
        "amplicon_name": ["name", "amplicon"],
        "fwd_sequence": ["fwd_seq", "forward_sequence", "forward"],
        "rev_sequence": ["rev_seq", "reverse_sequence", "reverse"],
    }

    # Normalize column names
    for canonical, aliases in column_aliases.items():
        if canonical in df.columns:
            continue
        for alias in aliases:
            if alias not in df.columns:
                continue
            df = df.rename({alias: canonical})
            break

    # Check required columns
    required = {"amplicon_name", "fwd_sequence", "rev_sequence"}
    missing = required - set(df.columns)
    assert not missing, f"TSV missing required columns: {missing}"

    # Validate sequences
    for row_idx, row in enumerate(df.iter_rows(named=True)):
        _validate_sequence(row["fwd_sequence"], f"row {row_idx} fwd_sequence")
        _validate_sequence(row["rev_sequence"], f"row {row_idx} rev_sequence")

    # Convert to lazy and add optional columns with defaults
    lf = df.lazy()

    if "chrom" not in df.columns:
        lf = lf.with_columns(pl.lit("unknown").alias("chrom"))
    if "amplicon_start" not in df.columns:
        lf = lf.with_columns(pl.lit(0).alias("amplicon_start"))
    if "amplicon_end" not in df.columns:
        lf = lf.with_columns(pl.lit(0).alias("amplicon_end"))

    return lf.select(
        "amplicon_name",
        "fwd_sequence",
        "rev_sequence",
        "chrom",
        "amplicon_start",
        "amplicon_end",
    )


def process_bed_input(
    bed_path: Path,
    reference_path: Path,
    fwd_suffix: str,
    rev_suffix: str,
) -> tuple[pl.LazyFrame, pl.LazyFrame]:
    """
    Process BED + FASTA input through the full pipeline.

    Pipeline stages:
      1. Parse BED file
      2. Validate primer names and pairs
      3. Parse reference FASTA
      4. Resplice primer combinations
      5. Extract sequences and generate primer pairs
      6. Generate output LazyFrames

    Args:
        bed_path: Path to primer BED file
        reference_path: Path to reference FASTA
        fwd_suffix: Forward primer suffix
        rev_suffix: Reverse primer suffix

    Returns:
        Tuple of (primer_pairs_lf, respliced_bed_lf)
    """
    # Stage 1: Parse BED
    logger.info(f"Parsing BED file: {bed_path}")
    records = parse_bed_file(bed_path)
    logger.info(f"Found {len(records)} primer records")

    # Stage 2: Validate
    logger.info("Validating primer names and pairs...")
    validated = validate_bed_records(records, fwd_suffix, rev_suffix)
    logger.success("Validation passed")

    # Stage 3: Parse reference
    logger.info(f"Parsing reference FASTA: {reference_path}")
    fasta_dict = parse_fasta(reference_path)
    logger.info(f"Loaded {len(fasta_dict)} sequence(s)")

    # Stage 4: Resplice
    logger.info("Resplicing primer combinations...")
    respliced = resplice_all_amplicons(validated, fwd_suffix, rev_suffix)
    logger.info(f"Generated {len(respliced)} respliced primers")

    # Stage 5: Extract sequences
    logger.info("Extracting primer sequences...")
    pairs = generate_primer_pairs(respliced, fasta_dict, fwd_suffix, rev_suffix)
    logger.info(f"Created {len(pairs)} primer pair(s)")

    # Stage 6: Generate output LazyFrames
    pairs_lf = _primer_pairs_to_dataframe(pairs)
    bed_lf = _respliced_to_bed_dataframe(respliced)

    return pairs_lf, bed_lf


def process_tsv_input(tsv_path: Path) -> pl.LazyFrame:
    """
    Process TSV input (validation and passthrough).

    Args:
        tsv_path: Path to primer TSV file

    Returns:
        Validated primer pairs LazyFrame
    """
    logger.info(f"Loading primer TSV: {tsv_path}")
    lf = validate_primer_tsv(tsv_path)

    # Collect to get row count for logging, then return lazy
    row_count = lf.collect().height
    logger.info(f"Loaded {row_count} primer pair(s)")

    return validate_primer_tsv(tsv_path)


def _configure_logging(verbosity: int) -> None:
    """Configure loguru logging based on verbosity level."""
    logger.remove()

    level = {
        0: "ERROR",
        1: "WARNING",
        2: "SUCCESS",
        3: "INFO",
        4: "DEBUG",
    }.get(min(verbosity, 4), "INFO")

    logger.add(
        sys.stderr,
        colorize=True,
        level=level,
        format="<green>{time:HH:mm:ss}</green> | <level>{level: <8}</level> | <level>{message}</level>",
    )


def parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # Input options (mutually exclusive)
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        "--input-bed",
        type=Path,
        help="Input BED file with primer coordinates",
    )
    input_group.add_argument(
        "--input-tsv",
        type=Path,
        help="Input TSV file with primer sequences (bypasses BED processing)",
    )

    # Reference (required for BED input)
    parser.add_argument(
        "--reference",
        type=Path,
        help="Reference FASTA file (required when using --input-bed)",
    )

    # Output options
    parser.add_argument(
        "--output-tsv",
        type=Path,
        default=Path("primer_pairs.tsv"),
        help="Output TSV file with primer pairs (default: primer_pairs.tsv)",
    )
    parser.add_argument(
        "--output-bed",
        type=Path,
        help="Output respliced BED file (optional)",
    )

    # Primer naming options
    parser.add_argument(
        "--fwd-suffix",
        type=str,
        default="_LEFT",
        help="Suffix for forward primers (default: _LEFT)",
    )
    parser.add_argument(
        "--rev-suffix",
        type=str,
        default="_RIGHT",
        help="Suffix for reverse primers (default: _RIGHT)",
    )

    # Verbosity
    parser.add_argument(
        "-v",
        "--verbose",
        action="count",
        default=3,
        help="Increase verbosity (default: INFO, -v for DEBUG)",
    )

    return parser.parse_args()


def validate_args(args: argparse.Namespace) -> None:
    """Validate command line arguments."""
    if args.input_bed:
        assert args.reference, "--reference is required when using --input-bed"
        assert args.input_bed.is_file(), f"BED file not found: {args.input_bed}"
        assert args.reference.is_file(), f"Reference not found: {args.reference}"

    if args.input_tsv:
        assert args.input_tsv.is_file(), f"TSV file not found: {args.input_tsv}"


def main() -> None:
    """Main entry point."""
    args = parse_args()
    _configure_logging(args.verbose)

    try:
        validate_args(args)
    except AssertionError as e:
        logger.error(f"{e}")
        sys.exit(1)

    try:
        if args.input_bed:
            pairs_lf, bed_lf = process_bed_input(
                args.input_bed,
                args.reference,
                args.fwd_suffix,
                args.rev_suffix,
            )

            if args.output_bed:
                logger.info(f"Writing respliced BED: {args.output_bed}")
                bed_lf.collect().write_csv(
                    args.output_bed,
                    separator="\t",
                    include_header=False,
                )
        else:
            pairs_lf = process_tsv_input(args.input_tsv)

        logger.info(f"Writing primer pairs TSV: {args.output_tsv}")
        collected = pairs_lf.collect()
        collected.write_csv(args.output_tsv, separator="\t")

        logger.success(f"Done! Prepared {collected.height} amplicon(s)")

    except (ValueError, AssertionError) as e:
        logger.error(f"{e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
