"""Tests for the assemble_report module."""

import json
from pathlib import Path

import pytest
from assemble_report import (
    DEFAULT_QC_THRESHOLDS,
    compute_summary,
    determine_qc_status,
    load_sample_metrics,
)
from reporting.schema import QCStatus

# --- Fixtures ---


@pytest.fixture
def sample_coverage_metrics() -> dict:
    """Sample coverage metrics with good values."""
    return {
        "sample_id": "sample1",
        "total_bases": 29903,
        "mean_coverage": 1250.5,
        "median_coverage": 1180.0,
        "genome_coverage_at_1x": 0.998,
        "genome_coverage_at_10x": 0.965,
        "genome_coverage_at_100x": 0.92,
        "min_coverage": 0,
        "max_coverage": 3500,
    }


@pytest.fixture
def sample_alignment_metrics() -> dict:
    """Sample alignment metrics with good values."""
    return {
        "sample_id": "sample1",
        "total_reads": 50000,
        "mapped_reads": 48000,
        "unmapped_reads": 2000,
        "mapping_rate": 0.96,
    }


@pytest.fixture
def sample_alignment_metrics_warn() -> dict:
    """Sample alignment metrics with marginal values (warn)."""
    return {
        "sample_id": "sample2",
        "total_reads": 800,
        "mapped_reads": 750,
        "unmapped_reads": 50,
        "mapping_rate": 0.94,
    }


@pytest.fixture
def sample_alignment_metrics_fail() -> dict:
    """Sample alignment metrics with poor values (fail)."""
    return {
        "sample_id": "sample3",
        "total_reads": 80,
        "mapped_reads": 50,
        "unmapped_reads": 30,
        "mapping_rate": 0.625,
    }


@pytest.fixture
def sample_coverage_metrics_warn() -> dict:
    """Sample coverage metrics with marginal values (warn)."""
    return {
        "sample_id": "sample2",
        "total_bases": 29903,
        "mean_coverage": 450.2,
        "median_coverage": 380.0,
        "genome_coverage_at_1x": 0.95,
        "genome_coverage_at_10x": 0.85,
        "genome_coverage_at_100x": 0.65,
        "min_coverage": 0,
        "max_coverage": 1200,
    }


@pytest.fixture
def sample_coverage_metrics_fail() -> dict:
    """Sample coverage metrics with poor values (fail)."""
    return {
        "sample_id": "sample3",
        "total_bases": 29903,
        "mean_coverage": 25.0,
        "median_coverage": 15.0,
        "genome_coverage_at_1x": 0.75,
        "genome_coverage_at_10x": 0.45,
        "genome_coverage_at_100x": 0.05,
        "min_coverage": 0,
        "max_coverage": 150,
    }


@pytest.fixture
def metrics_dir_single(tmp_path: Path, sample_coverage_metrics: dict) -> Path:
    """Create a temp directory with a single sample's metrics."""
    metrics_file = tmp_path / "sample1_coverage_metrics.json"
    metrics_file.write_text(json.dumps(sample_coverage_metrics))
    return tmp_path


@pytest.fixture
def metrics_dir_multiple(
    tmp_path: Path,
    sample_coverage_metrics: dict,
    sample_coverage_metrics_warn: dict,
    sample_coverage_metrics_fail: dict,
) -> Path:
    """Create a temp directory with multiple samples' metrics."""
    (tmp_path / "sample1_coverage_metrics.json").write_text(
        json.dumps(sample_coverage_metrics),
    )
    (tmp_path / "sample2_coverage_metrics.json").write_text(
        json.dumps(sample_coverage_metrics_warn),
    )
    (tmp_path / "sample3_coverage_metrics.json").write_text(
        json.dumps(sample_coverage_metrics_fail),
    )
    return tmp_path


@pytest.fixture
def metrics_dir_empty(tmp_path: Path) -> Path:
    """Create an empty temp directory."""
    return tmp_path


# --- TestLoadSampleMetrics ---


class TestLoadSampleMetrics:
    """Tests for load_sample_metrics function."""

    def test_load_single_sample(
        self,
        metrics_dir_single: Path,
        sample_coverage_metrics: dict,
    ) -> None:
        """Test loading a single sample's metrics."""
        samples = load_sample_metrics(metrics_dir_single)

        assert len(samples) == 1
        assert "sample1" in samples
        assert samples["sample1"]["coverage"] == sample_coverage_metrics

    def test_load_multiple_samples(self, metrics_dir_multiple: Path) -> None:
        """Test loading multiple samples' metrics."""
        samples = load_sample_metrics(metrics_dir_multiple)

        assert len(samples) == 3
        assert "sample1" in samples
        assert "sample2" in samples
        assert "sample3" in samples

    def test_load_empty_directory(self, metrics_dir_empty: Path) -> None:
        """Test loading from an empty directory returns empty dict."""
        samples = load_sample_metrics(metrics_dir_empty)

        assert samples == {}

    def test_merging_metrics_from_same_sample(self, tmp_path: Path) -> None:
        """Test that multiple metric types for same sample are merged."""
        # Create coverage metrics
        coverage = {
            "sample_id": "sample1",
            "mean_coverage": 100.0,
            "genome_coverage_at_10x": 0.95,
        }
        (tmp_path / "sample1_coverage_metrics.json").write_text(json.dumps(coverage))

        # Create consensus metrics (simulated - would come from Phase 2)
        consensus = {
            "sample_id": "sample1",
            "length": 29903,
            "n_count": 50,
            "completeness": 0.998,
        }
        (tmp_path / "sample1_consensus_metrics.json").write_text(json.dumps(consensus))

        samples = load_sample_metrics(tmp_path)

        assert len(samples) == 1
        assert "sample1" in samples
        assert "coverage" in samples["sample1"]
        assert "consensus" in samples["sample1"]
        assert samples["sample1"]["coverage"]["mean_coverage"] == 100.0
        assert samples["sample1"]["consensus"]["completeness"] == 0.998

    def test_sample_id_extracted_correctly(self, metrics_dir_single: Path) -> None:
        """Test that sample_id is correctly extracted from metrics."""
        samples = load_sample_metrics(metrics_dir_single)

        assert samples["sample1"]["sample_id"] == "sample1"


# --- TestDetermineQCStatus ---


class TestDetermineQCStatus:
    """Tests for determine_qc_status function."""

    def test_pass_case_above_all_thresholds(
        self,
        sample_coverage_metrics: dict,
        sample_alignment_metrics: dict,
    ) -> None:
        """Test that sample above all thresholds gets PASS status."""
        sample = {
            "coverage": sample_coverage_metrics,
            "alignment": sample_alignment_metrics,
        }
        status, notes = determine_qc_status(sample, DEFAULT_QC_THRESHOLDS)

        assert status == QCStatus.PASS
        assert notes == []

    def test_warn_case_coverage_between_thresholds(
        self,
        sample_coverage_metrics_warn: dict,
        sample_alignment_metrics: dict,
    ) -> None:
        """Test that sample with coverage between thresholds gets WARN status."""
        sample = {
            "coverage": sample_coverage_metrics_warn,
            "alignment": sample_alignment_metrics,
        }
        status, notes = determine_qc_status(sample, DEFAULT_QC_THRESHOLDS)

        assert status == QCStatus.WARN
        assert len(notes) == 1
        assert "Marginal coverage" in notes[0]

    def test_warn_case_low_read_count(
        self,
        sample_coverage_metrics: dict,
        sample_alignment_metrics_warn: dict,
    ) -> None:
        """Test that sample with low read count gets WARN status."""
        sample = {
            "coverage": sample_coverage_metrics,
            "alignment": sample_alignment_metrics_warn,
        }
        status, notes = determine_qc_status(sample, DEFAULT_QC_THRESHOLDS)

        assert status == QCStatus.WARN
        assert len(notes) == 1
        assert "Low read count" in notes[0]

    def test_fail_case_below_coverage_threshold(
        self,
        sample_coverage_metrics_fail: dict,
        sample_alignment_metrics: dict,
    ) -> None:
        """Test that sample below coverage thresholds gets FAIL status."""
        sample = {
            "coverage": sample_coverage_metrics_fail,
            "alignment": sample_alignment_metrics,
        }
        status, notes = determine_qc_status(sample, DEFAULT_QC_THRESHOLDS)

        assert status == QCStatus.FAIL
        assert any("Low coverage" in note for note in notes)

    def test_fail_case_very_low_read_count(
        self,
        sample_coverage_metrics: dict,
        sample_alignment_metrics_fail: dict,
    ) -> None:
        """Test that sample with very low read count gets FAIL status."""
        sample = {
            "coverage": sample_coverage_metrics,
            "alignment": sample_alignment_metrics_fail,
        }
        status, notes = determine_qc_status(sample, DEFAULT_QC_THRESHOLDS)

        assert status == QCStatus.FAIL
        assert any("Very low read count" in note for note in notes)

    def test_custom_thresholds(
        self,
        sample_coverage_metrics_warn: dict,
        sample_alignment_metrics: dict,
    ) -> None:
        """Test that custom thresholds are respected."""
        sample = {
            "coverage": sample_coverage_metrics_warn,
            "alignment": sample_alignment_metrics,
        }

        # With lower thresholds, this sample should pass
        lenient_thresholds = {
            "coverage_pass": 0.80,
            "coverage_warn": 0.60,
            "min_reads_warn": 100,
            "min_reads_fail": 10,
        }
        status, notes = determine_qc_status(sample, lenient_thresholds)

        assert status == QCStatus.PASS
        assert notes == []

    def test_completeness_check_when_consensus_present(
        self,
        sample_alignment_metrics: dict,
    ) -> None:
        """Test that completeness is checked when consensus metrics present."""
        sample = {
            "coverage": {
                "genome_coverage_at_10x": 0.99,  # Good coverage
            },
            "alignment": sample_alignment_metrics,
            "consensus": {
                "completeness": 0.92,  # Between warn (0.90) and pass (0.98) thresholds
                "n_percentage": 2.0,  # Below warn threshold
            },
        }
        status, notes = determine_qc_status(sample, DEFAULT_QC_THRESHOLDS)

        assert status == QCStatus.WARN
        assert any("completeness" in note.lower() for note in notes)

    def test_n_percentage_warn(self, sample_alignment_metrics: dict) -> None:
        """Test that high N percentage triggers WARN."""
        sample = {
            "coverage": {"genome_coverage_at_10x": 0.99},
            "alignment": sample_alignment_metrics,
            "consensus": {
                "completeness": 0.99,
                "n_percentage": 7.0,  # Between warn (5%) and fail (10%)
            },
        }
        status, notes = determine_qc_status(sample, DEFAULT_QC_THRESHOLDS)

        assert status == QCStatus.WARN
        assert any("High N bases" in note for note in notes)

    def test_n_percentage_fail(self, sample_alignment_metrics: dict) -> None:
        """Test that excessive N percentage triggers FAIL."""
        sample = {
            "coverage": {"genome_coverage_at_10x": 0.99},
            "alignment": sample_alignment_metrics,
            "consensus": {
                "completeness": 0.99,
                "n_percentage": 15.0,  # Above fail (10%)
            },
        }
        status, notes = determine_qc_status(sample, DEFAULT_QC_THRESHOLDS)

        assert status == QCStatus.FAIL
        assert any("Excessive N bases" in note for note in notes)

    def test_fail_overrides_warn(self, sample_alignment_metrics: dict) -> None:
        """Test that FAIL status takes precedence over WARN."""
        sample = {
            "coverage": {
                "genome_coverage_at_10x": 0.50,  # Below fail threshold
            },
            "alignment": sample_alignment_metrics,
            "consensus": {
                "completeness": 0.92,  # Between warn and pass (would be WARN)
                "n_percentage": 2.0,
            },
        }
        status, notes = determine_qc_status(sample, DEFAULT_QC_THRESHOLDS)

        assert status == QCStatus.FAIL
        assert len(notes) == 2  # Both issues noted

    def test_empty_sample_fails(self) -> None:
        """Test that sample with no metrics fails."""
        sample: dict = {}
        status, notes = determine_qc_status(sample, DEFAULT_QC_THRESHOLDS)

        assert status == QCStatus.FAIL
        # Should fail on read count (checked first) and coverage
        assert any("read count" in note.lower() for note in notes)
        assert any("coverage" in note.lower() for note in notes)


# --- TestComputeSummary ---


class TestComputeSummary:
    """Tests for compute_summary function."""

    def test_all_passing_samples(self, sample_coverage_metrics: dict) -> None:
        """Test summary with all passing samples."""
        samples = {
            "sample1": {
                "qc_status": QCStatus.PASS,
                "coverage": sample_coverage_metrics,
            },
            "sample2": {
                "qc_status": QCStatus.PASS,
                "coverage": sample_coverage_metrics,
            },
        }
        summary = compute_summary(samples)

        assert summary.sample_count == 2
        assert summary.samples_pass == 2
        assert summary.samples_warn == 0
        assert summary.samples_fail == 0

    def test_mixed_pass_warn_fail(
        self,
        sample_coverage_metrics: dict,
        sample_coverage_metrics_warn: dict,
        sample_coverage_metrics_fail: dict,
    ) -> None:
        """Test summary with mixed QC statuses."""
        samples = {
            "sample1": {
                "qc_status": QCStatus.PASS,
                "coverage": sample_coverage_metrics,
            },
            "sample2": {
                "qc_status": QCStatus.WARN,
                "coverage": sample_coverage_metrics_warn,
            },
            "sample3": {
                "qc_status": QCStatus.FAIL,
                "coverage": sample_coverage_metrics_fail,
            },
        }
        summary = compute_summary(samples)

        assert summary.sample_count == 3
        assert summary.samples_pass == 1
        assert summary.samples_warn == 1
        assert summary.samples_fail == 1

    def test_mean_coverage_calculation(
        self,
        sample_coverage_metrics: dict,
        sample_coverage_metrics_warn: dict,
    ) -> None:
        """Test that mean coverage is calculated correctly."""
        samples = {
            "sample1": {
                "qc_status": QCStatus.PASS,
                "coverage": sample_coverage_metrics,  # mean_coverage: 1250.5
            },
            "sample2": {
                "qc_status": QCStatus.WARN,
                "coverage": sample_coverage_metrics_warn,  # mean_coverage: 450.2
            },
        }
        summary = compute_summary(samples)

        expected_mean = (1250.5 + 450.2) / 2
        assert summary.mean_coverage_depth == pytest.approx(expected_mean)

    def test_mean_genome_coverage_calculation(
        self,
        sample_coverage_metrics: dict,
        sample_coverage_metrics_warn: dict,
    ) -> None:
        """Test that mean genome coverage is calculated correctly."""
        samples = {
            "sample1": {
                "qc_status": QCStatus.PASS,
                "coverage": sample_coverage_metrics,  # genome_coverage_at_10x: 0.965
            },
            "sample2": {
                "qc_status": QCStatus.WARN,
                "coverage": sample_coverage_metrics_warn,  # genome_coverage_at_10x: 0.85
            },
        }
        summary = compute_summary(samples)

        expected_mean = (0.965 + 0.85) / 2
        assert summary.mean_genome_coverage == pytest.approx(expected_mean)

    def test_empty_samples(self) -> None:
        """Test summary with no samples."""
        samples: dict = {}
        summary = compute_summary(samples)

        assert summary.sample_count == 0
        assert summary.samples_pass == 0
        assert summary.samples_warn == 0
        assert summary.samples_fail == 0
        assert summary.mean_coverage_depth == 0.0
        assert summary.mean_genome_coverage == 0.0

    def test_variant_totals_accumulated(self) -> None:
        """Test that variant totals are accumulated across samples."""
        samples = {
            "sample1": {
                "qc_status": QCStatus.PASS,
                "coverage": {"mean_coverage": 100, "genome_coverage_at_10x": 0.95},
                "variants": {"total_called": 25},
            },
            "sample2": {
                "qc_status": QCStatus.PASS,
                "coverage": {"mean_coverage": 100, "genome_coverage_at_10x": 0.95},
                "variants": {"total_called": 30},
            },
        }
        summary = compute_summary(samples)

        assert summary.total_variants_called == 55

    def test_string_qc_status_handled(self) -> None:
        """Test that string QC status values are handled (for JSON roundtrip)."""
        samples = {
            "sample1": {
                "qc_status": "pass",  # String instead of enum
                "coverage": {"mean_coverage": 100, "genome_coverage_at_10x": 0.95},
            },
            "sample2": {
                "qc_status": "warn",
                "coverage": {"mean_coverage": 100, "genome_coverage_at_10x": 0.95},
            },
            "sample3": {
                "qc_status": "fail",
                "coverage": {"mean_coverage": 100, "genome_coverage_at_10x": 0.95},
            },
        }
        summary = compute_summary(samples)

        assert summary.samples_pass == 1
        assert summary.samples_warn == 1
        assert summary.samples_fail == 1
