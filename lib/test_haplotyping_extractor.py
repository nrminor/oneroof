"""Tests for the haplotyping metrics extractor."""

from pathlib import Path

import pytest
from reporting.extractors.haplotyping import (
    HaplotypingExtractedMetrics,
    count_assigned_reads,
    count_snps,
    extract,
    parse_snp_haplotypes_fasta,
)


@pytest.fixture
def sample_devider_dir(tmp_path: Path) -> Path:
    """
    Create a temporary Devider output directory with known values.

    Creates:
    - 3 haplotypes with abundances 50.0, 30.0, 20.0 and depths 10.0, 6.0, 4.0
    - 15 assigned reads (5 + 6 + 4)
    - 5 SNP positions
    """
    devider_dir = tmp_path / "sample1_devider"
    devider_dir.mkdir()

    # Create snp_haplotypes.fasta
    snp_fasta = devider_dir / "snp_haplotypes.fasta"
    snp_fasta.write_text(
        """>Contig:ref1,Range:ALL-ALL,Haplotype:0,Abundance:50.0,Depth:10.0
ACGTACGT
>Contig:ref1,Range:ALL-ALL,Haplotype:1,Abundance:30.0,Depth:6.0
ACGTTCGT
>Contig:ref1,Range:ALL-ALL,Haplotype:2,Abundance:20.0,Depth:4.0
ACGTGCGT
""",
    )

    # Create ids.txt with 15 reads total
    ids_file = devider_dir / "ids.txt"
    ids_file.write_text(
        """Contig:ref1\tRange:ALL-ALL\tHaplotype:0\tread1\tread2\tread3\tread4\tread5
Contig:ref1\tRange:ALL-ALL\tHaplotype:1\tread6\tread7\tread8\tread9\tread10\tread11
Contig:ref1\tRange:ALL-ALL\tHaplotype:2\tread12\tread13\tread14\tread15
""",
    )

    # Create hap_info.txt with 5 SNP positions
    hap_info = devider_dir / "hap_info.txt"
    hap_info.write_text(
        """100\t0:1.00\t1:1.00\t0:1.00
200\t1:0.95\t0:1.00\t1:0.90
300\t0:1.00\t0:1.00\t1:1.00
400\t1:1.00\t1:0.85\t0:1.00
500\t0:0.90\t1:1.00\t1:1.00
""",
    )

    return devider_dir


@pytest.fixture
def empty_devider_dir(tmp_path: Path) -> Path:
    """Create an empty Devider output directory (no files)."""
    devider_dir = tmp_path / "empty_devider"
    devider_dir.mkdir()
    return devider_dir


@pytest.fixture
def single_haplotype_dir(tmp_path: Path) -> Path:
    """Create a Devider directory with only one haplotype."""
    devider_dir = tmp_path / "single_devider"
    devider_dir.mkdir()

    snp_fasta = devider_dir / "snp_haplotypes.fasta"
    snp_fasta.write_text(
        """>Contig:ref1,Range:ALL-ALL,Haplotype:0,Abundance:100.0,Depth:50.0
ACGTACGT
""",
    )

    ids_file = devider_dir / "ids.txt"
    ids_file.write_text(
        """Contig:ref1\tRange:ALL-ALL\tHaplotype:0\tread1\tread2\tread3
""",
    )

    hap_info = devider_dir / "hap_info.txt"
    hap_info.write_text("")  # No SNPs needed for single haplotype

    return devider_dir


@pytest.fixture
def no_snps_dir(tmp_path: Path) -> Path:
    """Create a Devider directory with haplotypes but no SNPs."""
    devider_dir = tmp_path / "no_snps_devider"
    devider_dir.mkdir()

    snp_fasta = devider_dir / "snp_haplotypes.fasta"
    snp_fasta.write_text(
        """>Contig:ref1,Range:ALL-ALL,Haplotype:0,Abundance:60.0,Depth:12.0
ACGT
>Contig:ref1,Range:ALL-ALL,Haplotype:1,Abundance:40.0,Depth:8.0
ACGT
""",
    )

    ids_file = devider_dir / "ids.txt"
    ids_file.write_text(
        """Contig:ref1\tRange:ALL-ALL\tHaplotype:0\tread1\tread2
Contig:ref1\tRange:ALL-ALL\tHaplotype:1\tread3
""",
    )

    # Empty hap_info.txt
    hap_info = devider_dir / "hap_info.txt"
    hap_info.write_text("")

    return devider_dir


class TestParseSnpHaplotypesFasta:
    """Test the FASTA header parsing function."""

    def test_parses_multiple_haplotypes(self, sample_devider_dir: Path) -> None:
        """Test parsing multiple haplotypes from FASTA."""
        fasta_path = sample_devider_dir / "snp_haplotypes.fasta"
        haplotypes = parse_snp_haplotypes_fasta(fasta_path)

        assert len(haplotypes) == 3
        assert haplotypes[0]["haplotype_id"] == 0
        assert haplotypes[1]["haplotype_id"] == 1
        assert haplotypes[2]["haplotype_id"] == 2

    def test_extracts_abundances(self, sample_devider_dir: Path) -> None:
        """Test that abundances are correctly extracted."""
        fasta_path = sample_devider_dir / "snp_haplotypes.fasta"
        haplotypes = parse_snp_haplotypes_fasta(fasta_path)

        assert haplotypes[0]["abundance"] == 50.0
        assert haplotypes[1]["abundance"] == 30.0
        assert haplotypes[2]["abundance"] == 20.0

    def test_extracts_depths(self, sample_devider_dir: Path) -> None:
        """Test that depths are correctly extracted."""
        fasta_path = sample_devider_dir / "snp_haplotypes.fasta"
        haplotypes = parse_snp_haplotypes_fasta(fasta_path)

        assert haplotypes[0]["depth"] == 10.0
        assert haplotypes[1]["depth"] == 6.0
        assert haplotypes[2]["depth"] == 4.0

    def test_nonexistent_file_returns_empty(self, tmp_path: Path) -> None:
        """Test that nonexistent file returns empty list."""
        haplotypes = parse_snp_haplotypes_fasta(tmp_path / "nonexistent.fasta")
        assert haplotypes == []

    def test_empty_file_returns_empty(self, tmp_path: Path) -> None:
        """Test that empty file returns empty list."""
        empty_fasta = tmp_path / "empty.fasta"
        empty_fasta.write_text("")
        haplotypes = parse_snp_haplotypes_fasta(empty_fasta)
        assert haplotypes == []

    def test_sorts_by_haplotype_id(self, tmp_path: Path) -> None:
        """Test that haplotypes are sorted by ID."""
        fasta = tmp_path / "unsorted.fasta"
        fasta.write_text(
            """>Contig:ref1,Range:ALL-ALL,Haplotype:2,Abundance:20.0,Depth:4.0
ACGT
>Contig:ref1,Range:ALL-ALL,Haplotype:0,Abundance:50.0,Depth:10.0
ACGT
>Contig:ref1,Range:ALL-ALL,Haplotype:1,Abundance:30.0,Depth:6.0
ACGT
""",
        )
        haplotypes = parse_snp_haplotypes_fasta(fasta)

        assert [h["haplotype_id"] for h in haplotypes] == [0, 1, 2]


class TestCountAssignedReads:
    """Test the read counting function."""

    def test_counts_reads_correctly(self, sample_devider_dir: Path) -> None:
        """Test that reads are counted correctly."""
        ids_file = sample_devider_dir / "ids.txt"
        count = count_assigned_reads(ids_file)
        # 5 + 6 + 4 = 15 reads
        assert count == 15

    def test_nonexistent_file_returns_zero(self, tmp_path: Path) -> None:
        """Test that nonexistent file returns 0."""
        count = count_assigned_reads(tmp_path / "nonexistent.txt")
        assert count == 0

    def test_empty_file_returns_zero(self, tmp_path: Path) -> None:
        """Test that empty file returns 0."""
        empty_file = tmp_path / "empty.txt"
        empty_file.write_text("")
        count = count_assigned_reads(empty_file)
        assert count == 0

    def test_single_read_per_haplotype(self, tmp_path: Path) -> None:
        """Test counting with single read per haplotype."""
        ids_file = tmp_path / "single.txt"
        ids_file.write_text(
            """Contig:ref1\tRange:ALL\tHaplotype:0\tread1
Contig:ref1\tRange:ALL\tHaplotype:1\tread2
""",
        )
        count = count_assigned_reads(ids_file)
        assert count == 2


class TestCountSnps:
    """Test the SNP counting function."""

    def test_counts_snps_correctly(self, sample_devider_dir: Path) -> None:
        """Test that SNPs are counted correctly."""
        hap_info = sample_devider_dir / "hap_info.txt"
        count = count_snps(hap_info)
        assert count == 5

    def test_nonexistent_file_returns_zero(self, tmp_path: Path) -> None:
        """Test that nonexistent file returns 0."""
        count = count_snps(tmp_path / "nonexistent.txt")
        assert count == 0

    def test_empty_file_returns_zero(self, tmp_path: Path) -> None:
        """Test that empty file returns 0."""
        empty_file = tmp_path / "empty.txt"
        empty_file.write_text("")
        count = count_snps(empty_file)
        assert count == 0

    def test_handles_header_line(self, tmp_path: Path) -> None:
        """Test that non-numeric header lines are skipped."""
        hap_info = tmp_path / "with_header.txt"
        hap_info.write_text(
            """Position\tHap0\tHap1\tHap2
100\t0:1.00\t1:1.00\t0:1.00
200\t1:0.95\t0:1.00\t1:0.90
""",
        )
        count = count_snps(hap_info)
        assert count == 2  # Header line should be skipped


class TestExtract:
    """Test the main extract function."""

    def test_extraction_returns_valid_model(self, sample_devider_dir: Path) -> None:
        """Test that extraction returns a valid HaplotypingExtractedMetrics model."""
        metrics = extract("sample1", sample_devider_dir)
        assert isinstance(metrics, HaplotypingExtractedMetrics)

    def test_sample_id_set_correctly(self, sample_devider_dir: Path) -> None:
        """Test that sample_id is correctly set."""
        metrics = extract("my_sample_123", sample_devider_dir)
        assert metrics.sample_id == "my_sample_123"

    def test_haplotypes_called(self, sample_devider_dir: Path) -> None:
        """Test haplotype count."""
        metrics = extract("sample1", sample_devider_dir)
        assert metrics.haplotypes_called == 3

    def test_reads_assigned(self, sample_devider_dir: Path) -> None:
        """Test read assignment count."""
        metrics = extract("sample1", sample_devider_dir)
        assert metrics.reads_assigned == 15

    def test_num_snps(self, sample_devider_dir: Path) -> None:
        """Test SNP count."""
        metrics = extract("sample1", sample_devider_dir)
        assert metrics.num_snps == 5

    def test_haplotype_abundances(self, sample_devider_dir: Path) -> None:
        """Test haplotype abundances list."""
        metrics = extract("sample1", sample_devider_dir)
        assert metrics.haplotype_abundances == [50.0, 30.0, 20.0]

    def test_haplotype_depths(self, sample_devider_dir: Path) -> None:
        """Test haplotype depths list."""
        metrics = extract("sample1", sample_devider_dir)
        assert metrics.haplotype_depths == [10.0, 6.0, 4.0]

    def test_empty_directory(self, empty_devider_dir: Path) -> None:
        """Test with empty directory (no files)."""
        metrics = extract("empty", empty_devider_dir)
        assert metrics.haplotypes_called == 0
        assert metrics.reads_assigned == 0
        assert metrics.num_snps == 0
        assert metrics.haplotype_abundances == []
        assert metrics.haplotype_depths == []

    def test_single_haplotype(self, single_haplotype_dir: Path) -> None:
        """Test with single haplotype."""
        metrics = extract("single", single_haplotype_dir)
        assert metrics.haplotypes_called == 1
        assert metrics.reads_assigned == 3
        assert metrics.haplotype_abundances == [100.0]
        assert metrics.haplotype_depths == [50.0]

    def test_no_snps(self, no_snps_dir: Path) -> None:
        """Test with haplotypes but no SNPs."""
        metrics = extract("no_snps", no_snps_dir)
        assert metrics.haplotypes_called == 2
        assert metrics.num_snps == 0

    def test_model_serialization(self, sample_devider_dir: Path) -> None:
        """Test that metrics can be serialized to JSON."""
        metrics = extract("sample1", sample_devider_dir)
        json_str = metrics.model_dump_json()
        assert "sample1" in json_str
        assert "haplotypes_called" in json_str
        assert "50.0" in json_str


class TestEdgeCases:
    """Test edge cases and boundary conditions."""

    def test_path_as_string(self, sample_devider_dir: Path) -> None:
        """Test that string paths work correctly."""
        metrics = extract("sample1", str(sample_devider_dir))
        assert metrics.haplotypes_called == 3

    def test_malformed_fasta_header(self, tmp_path: Path) -> None:
        """Test handling of malformed FASTA headers."""
        devider_dir = tmp_path / "malformed_devider"
        devider_dir.mkdir()

        snp_fasta = devider_dir / "snp_haplotypes.fasta"
        snp_fasta.write_text(
            """>Contig:ref1,Range:ALL-ALL,Haplotype:0,Abundance:50.0,Depth:10.0
ACGT
>MalformedHeader
ACGT
>Contig:ref1,Range:ALL-ALL,Haplotype:1,Abundance:invalid,Depth:6.0
ACGT
""",
        )

        (devider_dir / "ids.txt").write_text("")
        (devider_dir / "hap_info.txt").write_text("")

        metrics = extract("malformed", devider_dir)
        # Should only parse the valid haplotype
        assert metrics.haplotypes_called == 1

    def test_whitespace_in_ids_file(self, tmp_path: Path) -> None:
        """Test handling of blank lines in ids.txt."""
        devider_dir = tmp_path / "whitespace_devider"
        devider_dir.mkdir()

        (devider_dir / "snp_haplotypes.fasta").write_text("")
        ids_file = devider_dir / "ids.txt"
        ids_file.write_text(
            """Contig:ref1\tRange:ALL\tHaplotype:0\tread1\tread2

Contig:ref1\tRange:ALL\tHaplotype:1\tread3

""",
        )
        (devider_dir / "hap_info.txt").write_text("")

        metrics = extract("whitespace", devider_dir)
        assert metrics.reads_assigned == 3
