"""Tests for the QC dashboard visualizations."""

from pathlib import Path

import polars as pl
import pytest
from reporting.visualizations.qc_dashboard import (
    completeness_distribution,
    coverage_distribution,
    prepare_distribution_data,
    prepare_qc_status_data,
    qc_scatter,
    qc_status_summary,
)


@pytest.fixture
def multi_sample_metrics() -> dict[str, dict]:
    """
    Sample metrics dict with multiple samples having different QC statuses.
    """
    return {
        "sample1": {
            "qc_status": "pass",
            "coverage": {"mean_coverage": 150.0},
            "consensus": {"completeness": 0.98},
        },
        "sample2": {
            "qc_status": "pass",
            "coverage": {"mean_coverage": 120.0},
            "consensus": {"completeness": 0.95},
        },
        "sample3": {
            "qc_status": "warn",
            "coverage": {"mean_coverage": 80.0},
            "consensus": {"completeness": 0.85},
        },
        "sample4": {
            "qc_status": "fail",
            "coverage": {"mean_coverage": 10.0},
            "consensus": {"completeness": 0.50},
        },
        "sample5": {
            "qc_status": "fail",
            "coverage": {"mean_coverage": 5.0},
            "consensus": {"completeness": 0.20},
        },
    }


@pytest.fixture
def single_sample_metrics() -> dict[str, dict]:
    """Single sample metrics."""
    return {
        "sample1": {
            "qc_status": "pass",
            "coverage": {"mean_coverage": 150.0},
            "consensus": {"completeness": 0.98},
        },
    }


@pytest.fixture
def all_pass_metrics() -> dict[str, dict]:
    """All samples passing."""
    return {
        "sample1": {
            "qc_status": "pass",
            "coverage": {"mean_coverage": 150.0},
            "consensus": {"completeness": 0.98},
        },
        "sample2": {
            "qc_status": "pass",
            "coverage": {"mean_coverage": 120.0},
            "consensus": {"completeness": 0.95},
        },
    }


@pytest.fixture
def all_fail_metrics() -> dict[str, dict]:
    """All samples failing."""
    return {
        "sample1": {
            "qc_status": "fail",
            "coverage": {"mean_coverage": 5.0},
            "consensus": {"completeness": 0.10},
        },
        "sample2": {
            "qc_status": "fail",
            "coverage": {"mean_coverage": 2.0},
            "consensus": {"completeness": 0.05},
        },
    }


@pytest.fixture
def empty_metrics() -> dict[str, dict]:
    """Empty metrics dict."""
    return {}


@pytest.fixture
def missing_data_metrics() -> dict[str, dict]:
    """Metrics with missing coverage/consensus data."""
    return {
        "sample1": {
            "qc_status": "fail",
        },
        "sample2": {
            "qc_status": "warn",
            "coverage": {"mean_coverage": 50.0},
            # No consensus
        },
    }


class TestPrepareQcStatusData:
    """Test the QC status data preparation function."""

    def test_counts_all_statuses(self, multi_sample_metrics: dict[str, dict]) -> None:
        """Test that all status types are counted."""
        df = prepare_qc_status_data(multi_sample_metrics)

        assert len(df) == 3  # pass, warn, fail

        pass_count = df.filter(pl.col("status") == "pass")["count"].item()
        warn_count = df.filter(pl.col("status") == "warn")["count"].item()
        fail_count = df.filter(pl.col("status") == "fail")["count"].item()

        assert pass_count == 2
        assert warn_count == 1
        assert fail_count == 2

    def test_all_pass(self, all_pass_metrics: dict[str, dict]) -> None:
        """Test with all passing samples."""
        df = prepare_qc_status_data(all_pass_metrics)

        assert len(df) == 1  # Only pass status has counts
        assert df["status"].item() == "pass"
        assert df["count"].item() == 2

    def test_all_fail(self, all_fail_metrics: dict[str, dict]) -> None:
        """Test with all failing samples."""
        df = prepare_qc_status_data(all_fail_metrics)

        assert len(df) == 1
        assert df["status"].item() == "fail"
        assert df["count"].item() == 2

    def test_empty_metrics(self, empty_metrics: dict[str, dict]) -> None:
        """Test with empty metrics."""
        df = prepare_qc_status_data(empty_metrics)
        assert len(df) == 0


class TestPrepareDistributionData:
    """Test the distribution data preparation function."""

    def test_extracts_all_fields(self, multi_sample_metrics: dict[str, dict]) -> None:
        """Test that all fields are extracted."""
        df = prepare_distribution_data(multi_sample_metrics)

        assert len(df) == 5
        assert set(df.columns) == {
            "sample_id",
            "mean_coverage",
            "completeness",
            "qc_status",
        }

    def test_completeness_as_percentage(
        self,
        single_sample_metrics: dict[str, dict],
    ) -> None:
        """Test that completeness is converted to percentage."""
        df = prepare_distribution_data(single_sample_metrics)

        # 0.98 should become 98.0
        assert df["completeness"].item() == pytest.approx(98.0)

    def test_handles_missing_data(self, missing_data_metrics: dict[str, dict]) -> None:
        """Test handling of missing coverage/consensus data."""
        df = prepare_distribution_data(missing_data_metrics)

        assert len(df) == 2

        # Sample with no coverage/consensus should have 0 values
        sample1 = df.filter(pl.col("sample_id") == "sample1")
        assert sample1["mean_coverage"].item() == 0
        assert sample1["completeness"].item() == 0

    def test_empty_metrics(self, empty_metrics: dict[str, dict]) -> None:
        """Test with empty metrics."""
        df = prepare_distribution_data(empty_metrics)
        assert len(df) == 0


class TestQcStatusSummary:
    """Test the QC status summary chart generation."""

    def test_generates_html_file(
        self,
        multi_sample_metrics: dict[str, dict],
        tmp_path: Path,
    ) -> None:
        """Test that HTML file is generated."""
        output_path = tmp_path / "qc_status"
        paths = qc_status_summary(multi_sample_metrics, output_path)

        assert len(paths) == 1
        assert paths[0].suffix == ".html"
        assert paths[0].exists()

    def test_html_contains_expected_content(
        self,
        multi_sample_metrics: dict[str, dict],
        tmp_path: Path,
    ) -> None:
        """Test that generated HTML contains expected content."""
        output_path = tmp_path / "qc_status"
        paths = qc_status_summary(multi_sample_metrics, output_path)

        content = paths[0].read_text()
        assert "pass" in content.lower()
        assert "fail" in content.lower()

    def test_empty_metrics_returns_empty_list(
        self,
        empty_metrics: dict[str, dict],
        tmp_path: Path,
    ) -> None:
        """Test that empty metrics returns empty list."""
        output_path = tmp_path / "qc_status"
        paths = qc_status_summary(empty_metrics, output_path)
        assert paths == []

    def test_single_sample(
        self,
        single_sample_metrics: dict[str, dict],
        tmp_path: Path,
    ) -> None:
        """Test that single sample works."""
        output_path = tmp_path / "qc_status"
        paths = qc_status_summary(single_sample_metrics, output_path)

        assert len(paths) == 1
        assert paths[0].exists()


class TestCoverageDistribution:
    """Test the coverage distribution chart generation."""

    def test_generates_html_file(
        self,
        multi_sample_metrics: dict[str, dict],
        tmp_path: Path,
    ) -> None:
        """Test that HTML file is generated."""
        output_path = tmp_path / "coverage_dist"
        paths = coverage_distribution(multi_sample_metrics, output_path)

        assert len(paths) == 1
        assert paths[0].suffix == ".html"
        assert paths[0].exists()

    def test_empty_metrics_returns_empty_list(
        self,
        empty_metrics: dict[str, dict],
        tmp_path: Path,
    ) -> None:
        """Test that empty metrics returns empty list."""
        output_path = tmp_path / "coverage_dist"
        paths = coverage_distribution(empty_metrics, output_path)
        assert paths == []

    def test_single_sample_generates_bar(
        self,
        single_sample_metrics: dict[str, dict],
        tmp_path: Path,
    ) -> None:
        """Test that single sample generates a bar chart instead of histogram."""
        output_path = tmp_path / "coverage_dist"
        paths = coverage_distribution(single_sample_metrics, output_path)

        assert len(paths) == 1
        assert paths[0].exists()
        # Should contain sample_id since it's a bar chart
        content = paths[0].read_text()
        assert "sample1" in content


class TestCompletenessDistribution:
    """Test the completeness distribution chart generation."""

    def test_generates_html_file(
        self,
        multi_sample_metrics: dict[str, dict],
        tmp_path: Path,
    ) -> None:
        """Test that HTML file is generated."""
        output_path = tmp_path / "completeness_dist"
        paths = completeness_distribution(multi_sample_metrics, output_path)

        assert len(paths) == 1
        assert paths[0].suffix == ".html"
        assert paths[0].exists()

    def test_empty_metrics_returns_empty_list(
        self,
        empty_metrics: dict[str, dict],
        tmp_path: Path,
    ) -> None:
        """Test that empty metrics returns empty list."""
        output_path = tmp_path / "completeness_dist"
        paths = completeness_distribution(empty_metrics, output_path)
        assert paths == []

    def test_single_sample_generates_bar(
        self,
        single_sample_metrics: dict[str, dict],
        tmp_path: Path,
    ) -> None:
        """Test that single sample generates a bar chart instead of histogram."""
        output_path = tmp_path / "completeness_dist"
        paths = completeness_distribution(single_sample_metrics, output_path)

        assert len(paths) == 1
        assert paths[0].exists()


class TestQcScatter:
    """Test the QC scatter plot generation."""

    def test_generates_html_file(
        self,
        multi_sample_metrics: dict[str, dict],
        tmp_path: Path,
    ) -> None:
        """Test that HTML file is generated."""
        output_path = tmp_path / "qc_scatter"
        paths = qc_scatter(multi_sample_metrics, output_path)

        assert len(paths) == 1
        assert paths[0].suffix == ".html"
        assert paths[0].exists()

    def test_html_contains_expected_content(
        self,
        multi_sample_metrics: dict[str, dict],
        tmp_path: Path,
    ) -> None:
        """Test that generated HTML contains expected content."""
        output_path = tmp_path / "qc_scatter"
        paths = qc_scatter(multi_sample_metrics, output_path)

        content = paths[0].read_text()
        assert "sample1" in content
        assert "Coverage" in content

    def test_empty_metrics_returns_empty_list(
        self,
        empty_metrics: dict[str, dict],
        tmp_path: Path,
    ) -> None:
        """Test that empty metrics returns empty list."""
        output_path = tmp_path / "qc_scatter"
        paths = qc_scatter(empty_metrics, output_path)
        assert paths == []

    def test_single_sample(
        self,
        single_sample_metrics: dict[str, dict],
        tmp_path: Path,
    ) -> None:
        """Test that single sample works."""
        output_path = tmp_path / "qc_scatter"
        paths = qc_scatter(single_sample_metrics, output_path)

        assert len(paths) == 1
        assert paths[0].exists()

    def test_custom_title(
        self,
        multi_sample_metrics: dict[str, dict],
        tmp_path: Path,
    ) -> None:
        """Test that custom title is used."""
        output_path = tmp_path / "qc_scatter"
        paths = qc_scatter(multi_sample_metrics, output_path, title="Custom QC Scatter")

        content = paths[0].read_text()
        assert "Custom QC Scatter" in content
