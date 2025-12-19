"""
Haplotyping metrics extractor for OneRoof reporting.

Parses Devider output directory to extract haplotype phasing results.
Devider (https://github.com/bluenote-1577/devider) performs haplotype
reconstruction from long-read sequencing data.

This extractor is Nanopore-only, as Devider requires long reads for phasing.

Output directory structure from Devider:
    ${sample_id}_devider/
    ├── snp_haplotypes.fasta          # MSA of SNP alleles per haplotype
    ├── majority_vote_haplotypes.fasta # Base-level consensus sequences
    ├── ids.txt                        # Read-to-haplotype assignments
    ├── hap_info.txt                   # SNP positions and allele frequencies
    └── intermediate/                  # Debug files (not needed for metrics)

Example snp_haplotypes.fasta header:
    >Contig:OR483991.1,Range:ALL-ALL,Haplotype:0,Abundance:5.52,Depth:8.43

Example ids.txt line:
    Contig:X	Range:Y	Haplotype:0	read_id_1	read_id_2	read_id_3

Example hap_info.txt line:
    286	1:0.83	1:1.00	1:0.80
"""

import re
from pathlib import Path

from pydantic import BaseModel, Field


class HaplotypingExtractedMetrics(BaseModel):
    """Metrics extracted from Devider haplotype phasing output."""

    sample_id: str
    haplotypes_called: int = Field(
        ge=0,
        description="Number of distinct haplotypes identified",
    )
    reads_assigned: int = Field(ge=0, description="Total reads assigned to haplotypes")
    num_snps: int = Field(ge=0, description="Number of SNP positions used for phasing")
    haplotype_abundances: list[float] = Field(
        default_factory=list,
        description="Abundance percentage for each haplotype (may not sum to 100%)",
    )
    haplotype_depths: list[float] = Field(
        default_factory=list,
        description="Approximate depth for each haplotype",
    )


def parse_snp_haplotypes_fasta(fasta_path: Path) -> list[dict]:
    """
    Parse snp_haplotypes.fasta to extract haplotype metadata from headers.

    Each header contains comma-separated key:value pairs:
        >Contig:OR483991.1,Range:ALL-ALL,Haplotype:0,Abundance:5.52,Depth:8.43

    Args:
        fasta_path: Path to snp_haplotypes.fasta

    Returns:
        List of dicts with haplotype_id, abundance, and depth for each haplotype
    """
    haplotypes = []

    if not fasta_path.exists():
        return haplotypes

    with open(fasta_path) as f:
        for line in f:
            if not line.startswith(">"):
                continue

            # Parse header: >Key1:Val1,Key2:Val2,...
            header = line[1:].strip()
            header_dict: dict[str, str] = {}

            for part in header.split(","):
                if ":" in part:
                    key, value = part.split(":", 1)
                    header_dict[key] = value

            # Extract metrics we care about
            try:
                haplotype_id = int(header_dict.get("Haplotype", -1))
                abundance = float(header_dict.get("Abundance", 0.0))
                depth = float(header_dict.get("Depth", 0.0))

                if haplotype_id >= 0:
                    haplotypes.append(
                        {
                            "haplotype_id": haplotype_id,
                            "abundance": abundance,
                            "depth": depth,
                        },
                    )
            except (ValueError, TypeError):
                # Skip malformed headers
                continue

    # Sort by haplotype_id for consistent ordering
    haplotypes.sort(key=lambda h: h["haplotype_id"])

    return haplotypes


def count_assigned_reads(ids_file: Path) -> int:
    """
    Count total reads assigned to haplotypes from ids.txt.

    Each line has format:
        Contig:X	Range:Y	Haplotype:N	read_id_1	read_id_2	...

    The first 3 columns are metadata; columns 4+ are read IDs.

    Args:
        ids_file: Path to ids.txt

    Returns:
        Total count of assigned reads across all haplotypes
    """
    if not ids_file.exists():
        return 0

    total_reads = 0

    with open(ids_file) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            parts = line.split("\t")
            # Columns 0-2 are Contig, Range, Haplotype; 3+ are read IDs
            if len(parts) > 3:
                total_reads += len(parts) - 3

    return total_reads


def count_snps(hap_info_file: Path) -> int:
    """
    Count SNP positions from hap_info.txt.

    Each data line represents one SNP position:
        286	1:0.83	1:1.00	1:0.80

    The file may or may not have a header line. We count lines that
    start with a number (position).

    Args:
        hap_info_file: Path to hap_info.txt

    Returns:
        Number of SNP positions
    """
    if not hap_info_file.exists():
        return 0

    snp_count = 0

    with open(hap_info_file) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            # Check if line starts with a number (SNP position)
            # This handles both header and headerless formats
            first_field = line.split("\t")[0] if "\t" in line else line.split()[0]
            if re.match(r"^\d+$", first_field):
                snp_count += 1

    return snp_count


def extract(sample_id: str, devider_dir: Path | str) -> HaplotypingExtractedMetrics:
    """
    Extract haplotyping metrics from a Devider output directory.

    Args:
        sample_id: Sample identifier
        devider_dir: Path to the Devider output directory (e.g., sample1_devider/)

    Returns:
        Validated HaplotypingExtractedMetrics model
    """
    devider_dir = Path(devider_dir)

    # Define expected file paths
    snp_fasta = devider_dir / "snp_haplotypes.fasta"
    ids_file = devider_dir / "ids.txt"
    hap_info = devider_dir / "hap_info.txt"

    # Parse haplotypes from FASTA headers
    haplotypes = parse_snp_haplotypes_fasta(snp_fasta)

    # Count assigned reads
    reads_assigned = count_assigned_reads(ids_file)

    # Count SNPs
    num_snps = count_snps(hap_info)

    # Extract abundances and depths in haplotype order
    haplotype_abundances = [h["abundance"] for h in haplotypes]
    haplotype_depths = [h["depth"] for h in haplotypes]

    return HaplotypingExtractedMetrics(
        sample_id=sample_id,
        haplotypes_called=len(haplotypes),
        reads_assigned=reads_assigned,
        num_snps=num_snps,
        haplotype_abundances=haplotype_abundances,
        haplotype_depths=haplotype_depths,
    )
