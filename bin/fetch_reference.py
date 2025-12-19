#!/usr/bin/env python3
"""Fetch and validate reference sequences from local paths or NCBI accessions.

This script resolves reference inputs that can be either:
- Local file paths (validated and normalized)
- NCBI accessions (fetched via Entrez API)

Sequences are normalized to contain only valid IUPAC nucleotide characters,
with invalid characters replaced by 'N' and a warning issued.

Usage:
    fetch_reference.py fasta NC_045512.2 --output reference.fasta
    fetch_reference.py genbank NC_045512.2 --output reference.gbk
    fetch_reference.py fasta /path/to/local.fasta --output reference.fasta

    fetch_reference.py --help
"""

from __future__ import annotations

import re
import time
from pathlib import Path
from typing import Annotated

import typer
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from rich.console import Console

# NCBI Entrez configuration
Entrez.email = "nrminor@wisc.edu"  # type: ignore[assignment]
Entrez.tool = "oneroof"  # type: ignore[assignment]

app = typer.Typer(
    name="fetch_reference",
    help="Fetch and validate reference sequences from local paths or NCBI accessions.",
    add_completion=False,
    rich_markup_mode="rich",
    no_args_is_help=True,
    context_settings={"help_option_names": ["-h", "--help"]},
)
console = Console()

# Valid IUPAC nucleotide characters (uppercase)
IUPAC_NUCLEOTIDES = set("ACGTUMRWSYKVHDBN")

# Regex patterns for NCBI accessions
# RefSeq: NC_045512.2, NM_001301717.2, NP_000509.1, etc.
# GenBank: MN908947.3, AB123456.1, etc.
REFSEQ_PATTERN = re.compile(r"^[A-Z]{2}_\d+(\.\d+)?$")
GENBANK_PATTERN = re.compile(r"^[A-Z]{1,2}\d{5,}(\.\d+)?$")


def is_ncbi_accession(value: str) -> bool:
    """Check if a string looks like an NCBI accession."""
    value = value.strip()
    return bool(REFSEQ_PATTERN.match(value) or GENBANK_PATTERN.match(value))


def is_local_path(value: str) -> bool:
    """Check if a string looks like a local file path."""
    # Contains path separators or starts with path-like characters
    if "/" in value or "\\" in value:
        return True
    if value.startswith((".", "~")):
        return True
    # Check if it exists as a file
    return Path(value).exists()


def normalize_sequence(seq: str) -> tuple[str, int]:
    """
    Normalize a sequence to contain only valid IUPAC nucleotide characters.

    Invalid characters are replaced with 'N'. Returns the normalized sequence
    and the count of characters that were replaced.

    Args:
        seq: The input sequence string

    Returns:
        Tuple of (normalized_sequence, replacement_count)
    """
    seq_upper = seq.upper()
    normalized = []
    replacements = 0

    for char in seq_upper:
        if char in IUPAC_NUCLEOTIDES:
            normalized.append(char)
        elif char.isspace():
            # Skip whitespace silently
            continue
        else:
            normalized.append("N")
            replacements += 1

    return "".join(normalized), replacements


def validate_fasta(filepath: Path) -> list[SeqRecord]:
    """
    Validate and parse a FASTA file.

    Args:
        filepath: Path to the FASTA file

    Returns:
        List of SeqRecord objects

    Raises:
        typer.Exit: If the file is invalid or empty
    """
    try:
        records = list(SeqIO.parse(filepath, "fasta"))
    except Exception as e:
        console.print(f"[red]Error:[/red] Failed to parse FASTA file: {e}")
        raise typer.Exit(1) from e

    if not records:
        console.print(f"[red]Error:[/red] FASTA file contains no sequences: {filepath}")
        raise typer.Exit(1)

    # Check for HTML error pages (common NCBI failure mode)
    first_seq = str(records[0].seq)[:100].lower()
    if "<html" in first_seq or "<!doctype" in first_seq:
        console.print(
            "[red]Error:[/red] File appears to be HTML, not FASTA. "
            "This may indicate a failed download.",
        )
        raise typer.Exit(1)

    return records


def validate_genbank(filepath: Path) -> list[SeqRecord]:
    """
    Validate and parse a GenBank file.

    Args:
        filepath: Path to the GenBank file

    Returns:
        List of SeqRecord objects

    Raises:
        typer.Exit: If the file is invalid or empty
    """
    try:
        records = list(SeqIO.parse(filepath, "genbank"))
    except Exception as e:
        console.print(f"[red]Error:[/red] Failed to parse GenBank file: {e}")
        raise typer.Exit(1) from e

    if not records:
        console.print(f"[red]Error:[/red] GenBank file contains no records: {filepath}")
        raise typer.Exit(1)

    return records


def normalize_and_write_fasta(
    records: list[SeqRecord],
    output: Path,
    source: str,
) -> None:
    """
    Normalize sequences and write to FASTA file.

    Args:
        records: List of SeqRecord objects
        output: Output file path
        source: Description of the source (for logging)
    """
    normalized_records = []
    total_replacements = 0

    for record in records:
        record_id = record.id or "unknown"
        seq_str = str(record.seq)
        normalized_seq, replacements = normalize_sequence(seq_str)
        total_replacements += replacements

        if len(normalized_seq) == 0:
            console.print(
                f"[red]Error:[/red] Sequence '{record_id}' is empty after normalization",
            )
            raise typer.Exit(1)

        normalized_record = SeqRecord(
            Seq(normalized_seq),
            id=record_id,
            name=record.name or "unknown",
            description=record.description or "",
        )
        normalized_records.append(normalized_record)

    if total_replacements > 0:
        console.print(
            f"[yellow]Warning:[/yellow] Replaced {total_replacements} invalid "
            f"character(s) with 'N' in sequences from {source}",
        )

    SeqIO.write(normalized_records, output, "fasta")
    console.print(
        f"[green]Success:[/green] Wrote {len(normalized_records)} sequence(s) to {output}",
    )


def normalize_and_write_genbank(
    records: list[SeqRecord],
    output: Path,
    source: str,
) -> None:
    """
    Normalize sequences and write to GenBank file.

    Args:
        records: List of SeqRecord objects
        output: Output file path
        source: Description of the source (for logging)
    """
    normalized_records = []
    total_replacements = 0

    for record in records:
        record_id = record.id or "unknown"
        seq_str = str(record.seq)
        normalized_seq, replacements = normalize_sequence(seq_str)
        total_replacements += replacements

        if len(normalized_seq) == 0:
            console.print(
                f"[red]Error:[/red] Sequence '{record_id}' is empty after normalization",
            )
            raise typer.Exit(1)

        # Create new record preserving annotations
        normalized_record = SeqRecord(
            Seq(normalized_seq),
            id=record_id,
            name=record.name or "unknown",
            description=record.description or "",
            annotations=record.annotations.copy(),
            features=record.features,
            dbxrefs=record.dbxrefs,
        )
        normalized_records.append(normalized_record)

    if total_replacements > 0:
        console.print(
            f"[yellow]Warning:[/yellow] Replaced {total_replacements} invalid "
            f"character(s) with 'N' in sequences from {source}",
        )

    SeqIO.write(normalized_records, output, "genbank")
    console.print(
        f"[green]Success:[/green] Wrote {len(normalized_records)} record(s) to {output}",
    )


class NCBIFetchError(Exception):
    """Error fetching data from NCBI."""


def _attempt_ncbi_fetch(accession: str, rettype: str, retmode: str) -> str:
    """
    Single attempt to fetch from NCBI.

    Returns content on success, raises NCBIFetchError on failure.
    """
    handle = Entrez.efetch(
        db="nucleotide",
        id=accession,
        rettype=rettype,
        retmode=retmode,
    )
    content = handle.read()
    handle.close()

    if isinstance(content, bytes):
        content = content.decode("utf-8")

    if not content or len(content.strip()) == 0:
        msg = "Empty response from NCBI"
        raise NCBIFetchError(msg)

    if "<error>" in content.lower() or "nothing has been found" in content.lower():
        msg = f"NCBI returned an error for accession: {accession}"
        raise NCBIFetchError(msg)

    return content


def fetch_from_ncbi(
    accession: str,
    rettype: str,
    retmode: str,
    max_retries: int = 3,
) -> str:
    """
    Fetch a sequence from NCBI Entrez with retry logic.

    Args:
        accession: NCBI accession number
        rettype: Return type ('fasta' or 'gb')
        retmode: Return mode ('text')
        max_retries: Maximum number of retry attempts

    Returns:
        The fetched content as a string

    Raises:
        typer.Exit: If fetch fails after all retries
    """
    console.print(f"[cyan]Fetching {accession} from NCBI ({rettype} format)...[/cyan]")

    last_error: Exception | None = None

    for attempt in range(1, max_retries + 1):
        result = _try_fetch_attempt(accession, rettype, retmode, attempt, max_retries)
        if isinstance(result, str):
            return result
        last_error = result

    console.print(
        f"[red]Error:[/red] Failed to fetch {accession} from NCBI "
        f"after {max_retries} attempts",
    )
    raise typer.Exit(1) from last_error


def _try_fetch_attempt(
    accession: str,
    rettype: str,
    retmode: str,
    attempt: int,
    max_retries: int,
) -> str | Exception:
    """
    Try a single fetch attempt, returning content or the exception.

    Returns the content string on success, or the Exception on failure.
    """
    try:
        return _attempt_ncbi_fetch(accession, rettype, retmode)
    except Exception as e:  # noqa: BLE001
        console.print(f"[yellow]Attempt {attempt}/{max_retries} failed:[/yellow] {e}")
        if attempt < max_retries:
            wait_time = 2**attempt  # Exponential backoff
            console.print(f"[dim]Retrying in {wait_time} seconds...[/dim]")
            time.sleep(wait_time)
        return e


@app.command("fasta")
def resolve_fasta(
    reference: Annotated[
        str,
        typer.Argument(
            help="Local file path or NCBI accession (e.g., NC_045512.2, MN908947.3)",
        ),
    ],
    output: Annotated[
        Path,
        typer.Option(
            "--output",
            "-o",
            help="Output FASTA file path",
        ),
    ],
) -> None:
    """
    Resolve a FASTA reference from a local path or NCBI accession.

    If REFERENCE is a local file path, validates and normalizes the sequences.
    If REFERENCE is an NCBI accession, fetches from Entrez and normalizes.
    """
    reference = reference.strip()

    if is_local_path(reference):
        # Local file path
        filepath = Path(reference)
        if not filepath.exists():
            console.print(
                f"[red]Error:[/red] File not found: {filepath}\n"
                f"If you meant to provide an NCBI accession, check the format "
                f"(e.g., NC_045512.2 or MN908947.3)",
            )
            raise typer.Exit(1)

        console.print(f"[cyan]Validating local FASTA file: {filepath}[/cyan]")
        records = validate_fasta(filepath)
        normalize_and_write_fasta(records, output, str(filepath))

    elif is_ncbi_accession(reference):
        # NCBI accession - fetch from Entrez
        content = fetch_from_ncbi(reference, rettype="fasta", retmode="text")

        # Write to temporary file for parsing
        temp_file = output.with_suffix(".tmp")
        temp_file.write_text(content)

        try:
            records = validate_fasta(temp_file)
            normalize_and_write_fasta(records, output, f"NCBI:{reference}")
        finally:
            temp_file.unlink(missing_ok=True)

    else:
        console.print(
            f"[red]Error:[/red] '{reference}' is neither a valid file path "
            f"nor a recognized NCBI accession.\n\n"
            f"[dim]Expected formats:[/dim]\n"
            f"  - Local path: /path/to/reference.fasta or ./reference.fasta\n"
            f"  - RefSeq accession: NC_045512.2, NM_001301717.2\n"
            f"  - GenBank accession: MN908947.3, AB123456.1",
        )
        raise typer.Exit(1)


@app.command("genbank")
def resolve_genbank(
    reference: Annotated[
        str,
        typer.Argument(
            help="Local file path or NCBI accession (e.g., NC_045512.2, MN908947.3)",
        ),
    ],
    output: Annotated[
        Path,
        typer.Option(
            "--output",
            "-o",
            help="Output GenBank file path",
        ),
    ],
) -> None:
    """
    Resolve a GenBank reference from a local path or NCBI accession.

    If REFERENCE is a local file path, validates and normalizes the sequences.
    If REFERENCE is an NCBI accession, fetches from Entrez and normalizes.
    """
    reference = reference.strip()

    if is_local_path(reference):
        # Local file path
        filepath = Path(reference)
        if not filepath.exists():
            console.print(
                f"[red]Error:[/red] File not found: {filepath}\n"
                f"If you meant to provide an NCBI accession, check the format "
                f"(e.g., NC_045512.2 or MN908947.3)",
            )
            raise typer.Exit(1)

        console.print(f"[cyan]Validating local GenBank file: {filepath}[/cyan]")
        records = validate_genbank(filepath)
        normalize_and_write_genbank(records, output, str(filepath))

    elif is_ncbi_accession(reference):
        # NCBI accession - fetch from Entrez
        content = fetch_from_ncbi(reference, rettype="gb", retmode="text")

        # Write to temporary file for parsing
        temp_file = output.with_suffix(".tmp")
        temp_file.write_text(content)

        try:
            records = validate_genbank(temp_file)
            normalize_and_write_genbank(records, output, f"NCBI:{reference}")
        finally:
            temp_file.unlink(missing_ok=True)

    else:
        console.print(
            f"[red]Error:[/red] '{reference}' is neither a valid file path "
            f"nor a recognized NCBI accession.\n\n"
            f"[dim]Expected formats:[/dim]\n"
            f"  - Local path: /path/to/reference.gbk or ./reference.gbk\n"
            f"  - RefSeq accession: NC_045512.2, NM_001301717.2\n"
            f"  - GenBank accession: MN908947.3, AB123456.1",
        )
        raise typer.Exit(1)


if __name__ == "__main__":
    app()
