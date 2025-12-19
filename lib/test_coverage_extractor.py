"""Tests for the coverage metrics extractor."""

from pathlib import Path

import polars as pl
import pytest

from reporting.extractors.coverage import (
    CoverageMetrics,
    compute_coverage_stats,
    extract,
    load_coverage_bed,
)


@pytest.fixture
def sample_bed_content() -> str:
    """Sample BED content with known coverage values."""
    # 5 regions, 100bp each = 500bp total
    # Depths: 0, 50, 100, 150, 200
    # Weighted mean: (0*100 + 50*100 + 100*100 + 150*100 + 200*100) / 500 = 100
    # At 1x: 400/500 = 0.8
    # At 10x: 400/500 = 0.8
    # At 100x: 300/500 = 0.6
    return "\n".join(
        [
            "chr1\t0\t100\t0",
            "chr1\t100\t200\t50",
            "chr1\t200\t300\t100",
            "chr1\t300\t400\t150",
            "chr1\t400\t500\t200",
        ]
    )


@pytest.fixture
def sample_bed_file(sample_bed_content: str, tmp_path: Path) -> Path:
    """Create a temporary BED file with sample content."""
    bed_file = tmp_path / "test.per-base.bed"
    bed_file.write_text(sample_bed_content)
    return bed_file


@pytest.fixture
def uniform_bed_file(tmp_path: Path) -> Path:
    """Create a BED file with uniform coverage."""
    content = "\n".join(
        [
            "chr1\t0\t100\t50",
            "chr1\t100\t200\t50",
            "chr1\t200\t300\t50",
        ]
    )
    bed_file = tmp_path / "uniform.bed"
    bed_file.write_text(content)
    return bed_file


@pytest.fixture
def empty_bed_file(tmp_path: Path) -> Path:
    """Create an empty BED file."""
    bed_file = tmp_path / "empty.bed"
    bed_file.write_text("")
    return bed_file


class TestLoadCoverageBed:
    """Test BED file loading."""

    def test_load_valid_bed(self, sample_bed_file: Path) -> None:
        """Test loading a valid BED file."""
        df = load_coverage_bed(sample_bed_file)
        assert len(df) == 5
        assert df.shape == (5, 4)

    def test_correct_column_names(self, sample_bed_file: Path) -> None:
        """Test that loaded DataFrame has correct column names."""
        df = load_coverage_bed(sample_bed_file)
        assert list(df.columns) == ["chrom", "start", "end", "depth"]

    def test_correct_dtypes(self, sample_bed_file: Path) -> None:
        """Test that columns have correct data types."""
        df = load_coverage_bed(sample_bed_file)
        assert df["chrom"].dtype == pl.Utf8
        assert df["start"].dtype == pl.Int64
        assert df["end"].dtype == pl.Int64
        assert df["depth"].dtype == pl.Int64

    def test_empty_file_returns_empty_df(self, empty_bed_file: Path) -> None:
        """Test that empty file returns empty DataFrame."""
        df = load_coverage_bed(empty_bed_file)
        assert len(df) == 0


class TestComputeCoverageStats:
    """Test coverage statistics computation."""

    def test_mean_coverage_calculation(self, sample_bed_file: Path) -> None:
        """Test mean coverage with known values."""
        df = load_coverage_bed(sample_bed_file)
        stats = compute_coverage_stats(df)
        # (0*100 + 50*100 + 100*100 + 150*100 + 200*100) / 500 = 100
        assert stats["mean_coverage"] == 100.0

    def test_total_bases(self, sample_bed_file: Path) -> None:
        """Test total bases calculation."""
        df = load_coverage_bed(sample_bed_file)
        stats = compute_coverage_stats(df)
        assert stats["total_bases"] == 500

    def test_genome_coverage_at_1x(self, sample_bed_file: Path) -> None:
        """Test genome coverage at 1x threshold."""
        df = load_coverage_bed(sample_bed_file)
        stats = compute_coverage_stats(df)
        # 400 bases have depth >= 1 out of 500
        assert stats["genome_coverage_at_1x"] == 0.8

    def test_genome_coverage_at_10x(self, sample_bed_file: Path) -> None:
        """Test genome coverage at 10x threshold."""
        df = load_coverage_bed(sample_bed_file)
        stats = compute_coverage_stats(df)
        # 400 bases have depth >= 10 out of 500
        assert stats["genome_coverage_at_10x"] == 0.8

    def test_genome_coverage_at_100x(self, sample_bed_file: Path) -> None:
        """Test genome coverage at 100x threshold."""
        df = load_coverage_bed(sample_bed_file)
        stats = compute_coverage_stats(df)
        # 300 bases have depth >= 100 out of 500
        assert stats["genome_coverage_at_100x"] == 0.6

    def test_min_max_coverage(self, sample_bed_file: Path) -> None:
        """Test min and max coverage values."""
        df = load_coverage_bed(sample_bed_file)
        stats = compute_coverage_stats(df)
        assert stats["min_coverage"] == 0
        assert stats["max_coverage"] == 200

    def test_uniform_coverage(self, uniform_bed_file: Path) -> None:
        """Test with uniform coverage across all positions."""
        df = load_coverage_bed(uniform_bed_file)
        stats = compute_coverage_stats(df)
        assert stats["mean_coverage"] == 50.0
        assert stats["median_coverage"] == 50.0
        assert stats["min_coverage"] == 50
        assert stats["max_coverage"] == 50
        assert stats["genome_coverage_at_1x"] == 1.0
        assert stats["genome_coverage_at_10x"] == 1.0
        assert stats["genome_coverage_at_100x"] == 0.0

    def test_empty_dataframe(self) -> None:
        """Test with empty DataFrame returns zeros."""
        df = pl.DataFrame(
            {"chrom": [], "start": [], "end": [], "depth": []},
            schema={
                "chrom": pl.Utf8,
                "start": pl.Int64,
                "end": pl.Int64,
                "depth": pl.Int64,
            },
        )
        stats = compute_coverage_stats(df)
        assert stats["total_bases"] == 0
        assert stats["mean_coverage"] == 0.0
        assert stats["genome_coverage_at_1x"] == 0.0


class TestExtract:
    """Test the main extract function."""

    def test_full_extraction(self, sample_bed_file: Path) -> None:
        """Test full extraction pipeline."""
        metrics = extract("test_sample", sample_bed_file)
        assert isinstance(metrics, CoverageMetrics)

    def test_sample_id_set_correctly(self, sample_bed_file: Path) -> None:
        """Test that sample_id is correctly set."""
        metrics = extract("my_sample_123", sample_bed_file)
        assert metrics.sample_id == "my_sample_123"

    def test_metrics_values(self, sample_bed_file: Path) -> None:
        """Test that extracted metrics have expected values."""
        metrics = extract("test", sample_bed_file)
        assert metrics.total_bases == 500
        assert metrics.mean_coverage == 100.0
        assert metrics.genome_coverage_at_1x == 0.8
        assert metrics.genome_coverage_at_10x == 0.8
        assert metrics.genome_coverage_at_100x == 0.6
        assert metrics.min_coverage == 0
        assert metrics.max_coverage == 200

    def test_model_serialization(self, sample_bed_file: Path) -> None:
        """Test that metrics can be serialized to JSON."""
        metrics = extract("test", sample_bed_file)
        json_str = metrics.model_dump_json()
        assert "test" in json_str
        assert "mean_coverage" in json_str
        assert "100.0" in json_str
