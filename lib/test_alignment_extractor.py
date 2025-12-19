"""Tests for the alignment metrics extractor."""

from pathlib import Path

import pysam
import pytest
from reporting.extractors.alignment import (
    AlignmentExtractedMetrics,
    _count_reads_directly,
    extract,
)


@pytest.fixture
def sample_bam_file(tmp_path: Path) -> Path:
    """
    Create a temporary BAM file with known read counts.

    Creates a BAM with:
    - 80 mapped reads
    - 20 unmapped reads
    - Total: 100 reads
    - Mapping rate: 0.8
    """
    bam_path = tmp_path / "test.bam"

    # Create header with one reference sequence
    header = pysam.AlignmentHeader.from_dict(
        {
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [{"SN": "chr1", "LN": 10000}],
        },
    )

    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as outf:
        # Add 80 mapped reads
        for i in range(80):
            read = pysam.AlignedSegment()
            read.query_name = f"read_{i}"
            read.query_sequence = "ACGT" * 25  # 100bp read
            read.flag = 0  # Mapped, forward strand
            read.reference_id = 0  # chr1
            read.reference_start = i * 100
            read.mapping_quality = 60
            read.cigartuples = [(0, 100)]  # 100M
            read.query_qualities = pysam.qualitystring_to_array("I" * 100)
            outf.write(read)

        # Add 20 unmapped reads
        for i in range(20):
            read = pysam.AlignedSegment()
            read.query_name = f"unmapped_{i}"
            read.query_sequence = "ACGT" * 25
            read.flag = 4  # Unmapped
            read.reference_id = -1
            read.reference_start = -1
            read.mapping_quality = 0
            read.query_qualities = pysam.qualitystring_to_array("I" * 100)
            outf.write(read)

    # Index the BAM
    pysam.index(str(bam_path))

    return bam_path


@pytest.fixture
def all_mapped_bam(tmp_path: Path) -> Path:
    """Create a BAM file where all reads are mapped."""
    bam_path = tmp_path / "all_mapped.bam"

    header = pysam.AlignmentHeader.from_dict(
        {
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [{"SN": "chr1", "LN": 10000}],
        },
    )

    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as outf:
        for i in range(50):
            read = pysam.AlignedSegment()
            read.query_name = f"read_{i}"
            read.query_sequence = "ACGT" * 25
            read.flag = 0
            read.reference_id = 0
            read.reference_start = i * 100
            read.mapping_quality = 60
            read.cigartuples = [(0, 100)]
            read.query_qualities = pysam.qualitystring_to_array("I" * 100)
            outf.write(read)

    pysam.index(str(bam_path))
    return bam_path


@pytest.fixture
def all_unmapped_bam(tmp_path: Path) -> Path:
    """Create a BAM file where all reads are unmapped."""
    bam_path = tmp_path / "all_unmapped.bam"

    header = pysam.AlignmentHeader.from_dict(
        {
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [{"SN": "chr1", "LN": 10000}],
        },
    )

    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as outf:
        for i in range(30):
            read = pysam.AlignedSegment()
            read.query_name = f"unmapped_{i}"
            read.query_sequence = "ACGT" * 25
            read.flag = 4  # Unmapped
            read.reference_id = -1
            read.reference_start = -1
            read.mapping_quality = 0
            read.query_qualities = pysam.qualitystring_to_array("I" * 100)
            outf.write(read)

    pysam.index(str(bam_path))
    return bam_path


@pytest.fixture
def empty_bam(tmp_path: Path) -> Path:
    """Create an empty BAM file with header but no reads."""
    bam_path = tmp_path / "empty.bam"

    header = pysam.AlignmentHeader.from_dict(
        {
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [{"SN": "chr1", "LN": 10000}],
        },
    )

    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as _:
        pass  # Write no reads

    pysam.index(str(bam_path))
    return bam_path


class TestExtract:
    """Test the main extract function."""

    def test_extraction_returns_valid_model(self, sample_bam_file: Path) -> None:
        """Test that extraction returns a valid AlignmentExtractedMetrics model."""
        metrics = extract("test_sample", sample_bam_file)
        assert isinstance(metrics, AlignmentExtractedMetrics)

    def test_sample_id_set_correctly(self, sample_bam_file: Path) -> None:
        """Test that sample_id is correctly set."""
        metrics = extract("my_sample_123", sample_bam_file)
        assert metrics.sample_id == "my_sample_123"

    def test_read_counts(self, sample_bam_file: Path) -> None:
        """Test that read counts are correct."""
        metrics = extract("test", sample_bam_file)
        assert metrics.total_reads == 100
        assert metrics.mapped_reads == 80
        assert metrics.unmapped_reads == 20

    def test_mapping_rate_calculation(self, sample_bam_file: Path) -> None:
        """Test that mapping rate is correctly calculated."""
        metrics = extract("test", sample_bam_file)
        assert metrics.mapping_rate == 0.8

    def test_all_mapped_reads(self, all_mapped_bam: Path) -> None:
        """Test BAM with all reads mapped."""
        metrics = extract("test", all_mapped_bam)
        assert metrics.total_reads == 50
        assert metrics.mapped_reads == 50
        assert metrics.unmapped_reads == 0
        assert metrics.mapping_rate == 1.0

    def test_all_unmapped_reads(self, all_unmapped_bam: Path) -> None:
        """Test BAM with all reads unmapped."""
        metrics = extract("test", all_unmapped_bam)
        assert metrics.total_reads == 30
        assert metrics.mapped_reads == 0
        assert metrics.unmapped_reads == 30
        assert metrics.mapping_rate == 0.0

    def test_empty_bam(self, empty_bam: Path) -> None:
        """Test BAM with no reads."""
        metrics = extract("test", empty_bam)
        assert metrics.total_reads == 0
        assert metrics.mapped_reads == 0
        assert metrics.unmapped_reads == 0
        assert metrics.mapping_rate == 0.0

    def test_model_serialization(self, sample_bam_file: Path) -> None:
        """Test that metrics can be serialized to JSON."""
        metrics = extract("test", sample_bam_file)
        json_str = metrics.model_dump_json()
        assert "test" in json_str
        assert "mapped_reads" in json_str
        assert "80" in json_str
        assert "0.8" in json_str


class TestCountReadsDirectly:
    """Test the fallback direct read counting function."""

    def test_direct_count_matches_extract(self, sample_bam_file: Path) -> None:
        """Test that direct counting produces same results as extract."""
        metrics = extract("test", sample_bam_file)

        with pysam.AlignmentFile(str(sample_bam_file), "rb") as bam:
            mapped, unmapped = _count_reads_directly(bam)

        assert mapped == metrics.mapped_reads
        assert unmapped == metrics.unmapped_reads
