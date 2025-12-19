"""Tests for the consensus metrics extractor."""

from pathlib import Path

import pytest
from reporting.extractors.consensus import (
    ConsensusExtractedMetrics,
    extract,
)


@pytest.fixture
def sample_fasta_content() -> str:
    """
    Sample consensus FASTA with known values.

    Sequence: 100 bases total
    - 50 N bases (50% N, 50% completeness)
    - 20 G bases
    - 10 C bases
    - 10 A bases
    - 10 T bases
    - GC content of non-N: (20+10)/50 = 60%
    """
    # 50 Ns + 20 Gs + 10 Cs + 10 As + 10 Ts = 100 bases
    seq = "N" * 50 + "G" * 20 + "C" * 10 + "A" * 10 + "T" * 10
    return f">Consensus_test_sample\n{seq}\n"


@pytest.fixture
def sample_fasta_file(sample_fasta_content: str, tmp_path: Path) -> Path:
    """Create a temporary FASTA file with sample content."""
    fasta_file = tmp_path / "test.consensus.fa"
    fasta_file.write_text(sample_fasta_content)
    return fasta_file


@pytest.fixture
def perfect_fasta_file(tmp_path: Path) -> Path:
    """Create a FASTA file with no Ns (perfect sequence)."""
    # 100 bases: 30 G, 20 C, 25 A, 25 T = 50% GC
    seq = "G" * 30 + "C" * 20 + "A" * 25 + "T" * 25
    content = f">Consensus_perfect\n{seq}\n"
    fasta_file = tmp_path / "perfect.fa"
    fasta_file.write_text(content)
    return fasta_file


@pytest.fixture
def all_n_fasta_file(tmp_path: Path) -> Path:
    """Create a FASTA file with all Ns."""
    seq = "N" * 100
    content = f">Consensus_all_n\n{seq}\n"
    fasta_file = tmp_path / "all_n.fa"
    fasta_file.write_text(content)
    return fasta_file


@pytest.fixture
def empty_fasta_file(tmp_path: Path) -> Path:
    """Create an empty FASTA file."""
    fasta_file = tmp_path / "empty.fa"
    fasta_file.write_text("")
    return fasta_file


@pytest.fixture
def ambiguous_fasta_file(tmp_path: Path) -> Path:
    """Create a FASTA file with IUPAC ambiguity codes."""
    # 100 bases: 40 standard (ACGT), 50 N, 10 ambiguous (RYSWKM etc)
    seq = (
        "A" * 10
        + "C" * 10
        + "G" * 10
        + "T" * 10
        + "N" * 50
        + "R" * 3
        + "Y" * 3
        + "S" * 2
        + "W" * 2
    )
    content = f">Consensus_ambiguous\n{seq}\n"
    fasta_file = tmp_path / "ambiguous.fa"
    fasta_file.write_text(content)
    return fasta_file


@pytest.fixture
def multi_sequence_fasta_file(tmp_path: Path) -> Path:
    """Create a FASTA file with multiple sequences (multi-segment)."""
    content = """>Segment1
ACGTACGTNN
>Segment2
GGCCAATTNN
"""
    fasta_file = tmp_path / "multi.fa"
    fasta_file.write_text(content)
    return fasta_file


class TestExtract:
    """Test the main extract function."""

    def test_extraction_returns_valid_model(self, sample_fasta_file: Path) -> None:
        """Test that extraction returns a valid ConsensusExtractedMetrics model."""
        metrics = extract("test_sample", sample_fasta_file)
        assert isinstance(metrics, ConsensusExtractedMetrics)

    def test_sample_id_set_correctly(self, sample_fasta_file: Path) -> None:
        """Test that sample_id is correctly set."""
        metrics = extract("my_sample_123", sample_fasta_file)
        assert metrics.sample_id == "my_sample_123"

    def test_length(self, sample_fasta_file: Path) -> None:
        """Test total length calculation."""
        metrics = extract("test", sample_fasta_file)
        assert metrics.length == 100

    def test_n_count(self, sample_fasta_file: Path) -> None:
        """Test N base counting."""
        metrics = extract("test", sample_fasta_file)
        assert metrics.n_count == 50

    def test_n_percentage(self, sample_fasta_file: Path) -> None:
        """Test N percentage calculation."""
        metrics = extract("test", sample_fasta_file)
        assert metrics.n_percentage == 50.0

    def test_completeness(self, sample_fasta_file: Path) -> None:
        """Test completeness calculation (1 - n_count/length)."""
        metrics = extract("test", sample_fasta_file)
        assert metrics.completeness == 0.5

    def test_gc_content(self, sample_fasta_file: Path) -> None:
        """Test GC content calculation (of non-N bases)."""
        metrics = extract("test", sample_fasta_file)
        # (20 G + 10 C) / 50 non-N bases = 60%
        assert metrics.gc_content == 60.0

    def test_perfect_sequence(self, perfect_fasta_file: Path) -> None:
        """Test with perfect sequence (no Ns)."""
        metrics = extract("test", perfect_fasta_file)
        assert metrics.length == 100
        assert metrics.n_count == 0
        assert metrics.n_percentage == 0.0
        assert metrics.completeness == 1.0
        # (30 G + 20 C) / 100 = 50%
        assert metrics.gc_content == 50.0

    def test_all_n_sequence(self, all_n_fasta_file: Path) -> None:
        """Test with all-N sequence."""
        metrics = extract("test", all_n_fasta_file)
        assert metrics.length == 100
        assert metrics.n_count == 100
        assert metrics.n_percentage == 100.0
        assert metrics.completeness == 0.0
        # GC content undefined when all Ns, should be 0
        assert metrics.gc_content == 0.0

    def test_empty_file(self, empty_fasta_file: Path) -> None:
        """Test with empty FASTA file."""
        metrics = extract("test", empty_fasta_file)
        assert metrics.length == 0
        assert metrics.n_count == 0
        assert metrics.n_percentage == 0.0
        assert metrics.completeness == 0.0
        assert metrics.gc_content == 0.0

    def test_ambiguous_bases(self, ambiguous_fasta_file: Path) -> None:
        """Test counting of IUPAC ambiguity codes."""
        metrics = extract("test", ambiguous_fasta_file)
        assert metrics.length == 100
        assert metrics.ambiguous_bases == 10  # R*3 + Y*3 + S*2 + W*2

    def test_multi_sequence_fasta(self, multi_sequence_fasta_file: Path) -> None:
        """Test with multi-sequence FASTA (metrics combined)."""
        metrics = extract("test", multi_sequence_fasta_file)
        # Segment1: ACGTACGTNN (10 bases, 2 N)
        # Segment2: GGCCAATTNN (10 bases, 2 N)
        # Total: 20 bases, 4 N
        assert metrics.length == 20
        assert metrics.n_count == 4
        assert metrics.n_percentage == 20.0
        assert metrics.completeness == 0.8

    def test_model_serialization(self, sample_fasta_file: Path) -> None:
        """Test that metrics can be serialized to JSON."""
        metrics = extract("test", sample_fasta_file)
        json_str = metrics.model_dump_json()
        assert "test" in json_str
        assert "completeness" in json_str
        assert "gc_content" in json_str


class TestEdgeCases:
    """Test edge cases and boundary conditions."""

    def test_lowercase_sequence(self, tmp_path: Path) -> None:
        """Test that lowercase sequences are handled correctly."""
        content = ">test\nacgtnnnn\n"
        fasta_file = tmp_path / "lowercase.fa"
        fasta_file.write_text(content)

        metrics = extract("test", fasta_file)
        assert metrics.length == 8
        assert metrics.n_count == 4
        assert metrics.completeness == 0.5

    def test_mixed_case_sequence(self, tmp_path: Path) -> None:
        """Test that mixed case sequences are handled correctly."""
        content = ">test\nAcGtNnNn\n"
        fasta_file = tmp_path / "mixed.fa"
        fasta_file.write_text(content)

        metrics = extract("test", fasta_file)
        assert metrics.length == 8
        assert metrics.n_count == 4

    def test_multiline_sequence(self, tmp_path: Path) -> None:
        """Test that multiline sequences are handled correctly."""
        content = ">test\nACGT\nACGT\nNNNN\n"
        fasta_file = tmp_path / "multiline.fa"
        fasta_file.write_text(content)

        metrics = extract("test", fasta_file)
        assert metrics.length == 12
        assert metrics.n_count == 4
