"""Tests for the variant summary visualizations."""

from pathlib import Path

import polars as pl
import pytest
from reporting.visualizations.variant_summary import (
    prepare_variant_effect_data,
    prepare_variant_type_data,
    variant_effect_bar,
    variant_type_bar,
)


@pytest.fixture
def sample_metrics() -> dict[str, dict]:
    """
    Sample metrics dict with variant data for multiple samples.

    Creates 3 samples with different variant profiles.
    """
    return {
        "sample1": {
            "variants": {
                "snps": 10,
                "insertions": 2,
                "deletions": 3,
                "mnps": 1,
                "by_effect": {
                    "missense_variant": 5,
                    "synonymous_variant": 4,
                    "stop_gained": 1,
                },
            },
        },
        "sample2": {
            "variants": {
                "snps": 20,
                "insertions": 5,
                "deletions": 2,
                "mnps": 0,
                "by_effect": {
                    "missense_variant": 12,
                    "synonymous_variant": 8,
                    "frameshift_variant": 2,
                },
            },
        },
        "sample3": {
            "variants": {
                "snps": 5,
                "insertions": 1,
                "deletions": 1,
                "mnps": 0,
                "by_effect": {
                    "synonymous_variant": 3,
                    "intergenic_region": 2,
                },
            },
        },
    }


@pytest.fixture
def single_sample_metrics() -> dict[str, dict]:
    """Single sample with variant data."""
    return {
        "sample1": {
            "variants": {
                "snps": 10,
                "insertions": 2,
                "deletions": 3,
                "mnps": 1,
                "by_effect": {
                    "missense_variant": 5,
                    "synonymous_variant": 4,
                },
            },
        },
    }


@pytest.fixture
def empty_metrics() -> dict[str, dict]:
    """Empty metrics dict."""
    return {}


@pytest.fixture
def no_variants_metrics() -> dict[str, dict]:
    """Metrics with no variant data."""
    return {
        "sample1": {
            "coverage": {"mean_coverage": 100},
        },
        "sample2": {
            "variants": {},
        },
    }


class TestPrepareVariantTypeData:
    """Test the variant type data preparation function."""

    def test_prepares_data_for_multiple_samples(
        self,
        sample_metrics: dict[str, dict],
    ) -> None:
        """Test that data is prepared correctly for multiple samples."""
        df = prepare_variant_type_data(sample_metrics)

        assert len(df) == 12  # 3 samples * 4 mutation types
        assert set(df["sample_id"].unique().to_list()) == {
            "sample1",
            "sample2",
            "sample3",
        }
        assert set(df["mutation_type"].unique().to_list()) == {
            "SNP",
            "insertion",
            "deletion",
            "MNP",
        }

    def test_correct_counts(self, sample_metrics: dict[str, dict]) -> None:
        """Test that counts are extracted correctly."""
        df = prepare_variant_type_data(sample_metrics)

        sample1_snps = df.filter(
            (pl.col("sample_id") == "sample1") & (pl.col("mutation_type") == "SNP"),
        )["count"].item()
        assert sample1_snps == 10

        sample2_insertions = df.filter(
            (pl.col("sample_id") == "sample2")
            & (pl.col("mutation_type") == "insertion"),
        )["count"].item()
        assert sample2_insertions == 5

    def test_empty_metrics_returns_empty_df(
        self,
        empty_metrics: dict[str, dict],
    ) -> None:
        """Test that empty metrics returns empty DataFrame."""
        df = prepare_variant_type_data(empty_metrics)
        assert len(df) == 0

    def test_no_variants_returns_empty_df(
        self,
        no_variants_metrics: dict[str, dict],
    ) -> None:
        """Test that samples without variant data are skipped."""
        df = prepare_variant_type_data(no_variants_metrics)
        assert len(df) == 0


class TestPrepareVariantEffectData:
    """Test the variant effect data preparation function."""

    def test_prepares_data_for_multiple_samples(
        self,
        sample_metrics: dict[str, dict],
    ) -> None:
        """Test that data is prepared correctly for multiple samples."""
        df = prepare_variant_effect_data(sample_metrics)

        # sample1: 3 effects, sample2: 3 effects, sample3: 2 effects = 8 rows
        assert len(df) == 8
        assert set(df["sample_id"].unique().to_list()) == {
            "sample1",
            "sample2",
            "sample3",
        }

    def test_correct_effect_counts(self, sample_metrics: dict[str, dict]) -> None:
        """Test that effect counts are extracted correctly."""
        df = prepare_variant_effect_data(sample_metrics)

        sample1_missense = df.filter(
            (pl.col("sample_id") == "sample1")
            & (pl.col("effect_type") == "missense_variant"),
        )["count"].item()
        assert sample1_missense == 5

    def test_empty_metrics_returns_empty_df(
        self,
        empty_metrics: dict[str, dict],
    ) -> None:
        """Test that empty metrics returns empty DataFrame."""
        df = prepare_variant_effect_data(empty_metrics)
        assert len(df) == 0

    def test_no_effects_returns_empty_df(
        self,
        no_variants_metrics: dict[str, dict],
    ) -> None:
        """Test that samples without effect data are skipped."""
        df = prepare_variant_effect_data(no_variants_metrics)
        assert len(df) == 0


class TestVariantTypeBar:
    """Test the variant type bar chart generation."""

    def test_generates_html_file(
        self,
        sample_metrics: dict[str, dict],
        tmp_path: Path,
    ) -> None:
        """Test that HTML file is generated."""
        output_path = tmp_path / "variant_type"
        paths = variant_type_bar(sample_metrics, output_path)

        assert len(paths) == 1
        assert paths[0].suffix == ".html"
        assert paths[0].exists()

    def test_html_contains_expected_content(
        self,
        sample_metrics: dict[str, dict],
        tmp_path: Path,
    ) -> None:
        """Test that generated HTML contains expected content."""
        output_path = tmp_path / "variant_type"
        paths = variant_type_bar(sample_metrics, output_path)

        content = paths[0].read_text()
        assert "sample1" in content
        assert "SNP" in content

    def test_empty_metrics_returns_empty_list(
        self,
        empty_metrics: dict[str, dict],
        tmp_path: Path,
    ) -> None:
        """Test that empty metrics returns empty list."""
        output_path = tmp_path / "variant_type"
        paths = variant_type_bar(empty_metrics, output_path)
        assert paths == []

    def test_single_sample(
        self,
        single_sample_metrics: dict[str, dict],
        tmp_path: Path,
    ) -> None:
        """Test that single sample works correctly."""
        output_path = tmp_path / "variant_type"
        paths = variant_type_bar(single_sample_metrics, output_path)

        assert len(paths) == 1
        assert paths[0].exists()

    def test_custom_title(
        self,
        sample_metrics: dict[str, dict],
        tmp_path: Path,
    ) -> None:
        """Test that custom title is used."""
        output_path = tmp_path / "variant_type"
        paths = variant_type_bar(
            sample_metrics,
            output_path,
            title="Custom Variant Title",
        )

        content = paths[0].read_text()
        assert "Custom Variant Title" in content


class TestVariantEffectBar:
    """Test the variant effect bar chart generation."""

    def test_generates_html_file(
        self,
        sample_metrics: dict[str, dict],
        tmp_path: Path,
    ) -> None:
        """Test that HTML file is generated."""
        output_path = tmp_path / "variant_effect"
        paths = variant_effect_bar(sample_metrics, output_path)

        assert len(paths) == 1
        assert paths[0].suffix == ".html"
        assert paths[0].exists()

    def test_html_contains_expected_content(
        self,
        sample_metrics: dict[str, dict],
        tmp_path: Path,
    ) -> None:
        """Test that generated HTML contains expected content."""
        output_path = tmp_path / "variant_effect"
        paths = variant_effect_bar(sample_metrics, output_path)

        content = paths[0].read_text()
        assert "sample1" in content
        assert "missense_variant" in content

    def test_empty_metrics_returns_empty_list(
        self,
        empty_metrics: dict[str, dict],
        tmp_path: Path,
    ) -> None:
        """Test that empty metrics returns empty list."""
        output_path = tmp_path / "variant_effect"
        paths = variant_effect_bar(empty_metrics, output_path)
        assert paths == []

    def test_single_sample(
        self,
        single_sample_metrics: dict[str, dict],
        tmp_path: Path,
    ) -> None:
        """Test that single sample works correctly."""
        output_path = tmp_path / "variant_effect"
        paths = variant_effect_bar(single_sample_metrics, output_path)

        assert len(paths) == 1
        assert paths[0].exists()

    def test_custom_title(
        self,
        sample_metrics: dict[str, dict],
        tmp_path: Path,
    ) -> None:
        """Test that custom title is used."""
        output_path = tmp_path / "variant_effect"
        paths = variant_effect_bar(
            sample_metrics,
            output_path,
            title="Custom Effect Title",
        )

        content = paths[0].read_text()
        assert "Custom Effect Title" in content
