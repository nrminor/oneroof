"""Tests for the amplicon efficiency visualizations."""

from pathlib import Path

import polars as pl
import pytest
from reporting.visualizations.amplicon_efficiency import (
    amplicon_dropout_scatter,
    amplicon_heatmap,
    amplicon_ranking_bar,
    load_amplicon_summary,
    prepare_amplicon_stats,
    prepare_heatmap_data,
)


@pytest.fixture
def amplicon_summary_tsv(tmp_path: Path) -> Path:
    """Create a sample amplicon_summary.tsv with varied performance."""
    content = "\n".join(
        [
            "sample_name\tamplicon_name\tstart_pos\tend_pos\treads",
            # amp1: good performer (high reads, no dropouts)
            "sampleA\tamp1\t100\t200\t500",
            "sampleB\tamp1\t100\t200\t450",
            "sampleC\tamp1\t100\t200\t480",
            # amp2: moderate performer
            "sampleA\tamp2\t300\t400\t100",
            "sampleB\tamp2\t300\t400\t80",
            "sampleC\tamp2\t300\t400\t90",
            # amp3: poor performer (low reads, some dropouts)
            "sampleA\tamp3\t500\t600\t10",
            "sampleB\tamp3\t500\t600\t0",  # dropout
            "sampleC\tamp3\t500\t600\t5",
            # amp4: very poor (all dropouts)
            "sampleA\tamp4\t700\t800\t0",
            "sampleB\tamp4\t700\t800\t0",
            "sampleC\tamp4\t700\t800\t0",
        ],
    )
    tsv_file = tmp_path / "amplicon_summary.tsv"
    tsv_file.write_text(content)
    return tsv_file


@pytest.fixture
def single_sample_summary(tmp_path: Path) -> Path:
    """Create an amplicon summary with only one sample."""
    content = "\n".join(
        [
            "sample_name\tamplicon_name\tstart_pos\tend_pos\treads",
            "sampleA\tamp1\t100\t200\t500",
            "sampleA\tamp2\t300\t400\t100",
            "sampleA\tamp3\t500\t600\t0",
        ],
    )
    tsv_file = tmp_path / "amplicon_summary.tsv"
    tsv_file.write_text(content)
    return tsv_file


@pytest.fixture
def empty_summary(tmp_path: Path) -> Path:
    """Create an empty amplicon summary (header only)."""
    content = "sample_name\tamplicon_name\tstart_pos\tend_pos\treads"
    tsv_file = tmp_path / "amplicon_summary.tsv"
    tsv_file.write_text(content)
    return tsv_file


class TestLoadAmpliconSummary:
    """Test loading amplicon summary TSV."""

    def test_load_valid_summary(self, amplicon_summary_tsv: Path) -> None:
        """Test loading a valid amplicon summary."""
        df = load_amplicon_summary(amplicon_summary_tsv)
        assert len(df) == 12  # 4 amplicons × 3 samples
        assert list(df.columns) == [
            "sample_name",
            "amplicon_name",
            "start_pos",
            "end_pos",
            "reads",
        ]

    def test_load_empty_summary(self, empty_summary: Path) -> None:
        """Test loading an empty summary returns empty DataFrame."""
        df = load_amplicon_summary(empty_summary)
        assert len(df) == 0


class TestPrepareAmpliconStats:
    """Test computing per-amplicon statistics."""

    def test_computes_median_reads(self, amplicon_summary_tsv: Path) -> None:
        """Test that median reads are computed correctly."""
        data = load_amplicon_summary(amplicon_summary_tsv)
        stats = prepare_amplicon_stats(data)

        # amp1 has reads [500, 450, 480], median = 480
        amp1 = stats.filter(pl.col("amplicon_name") == "amp1")
        assert amp1["median_reads"][0] == 480.0

    def test_computes_dropout_rate(self, amplicon_summary_tsv: Path) -> None:
        """Test that dropout rate is computed correctly."""
        data = load_amplicon_summary(amplicon_summary_tsv)
        stats = prepare_amplicon_stats(data)

        # amp3 has 1 dropout out of 3 samples = 33.3%
        amp3 = stats.filter(pl.col("amplicon_name") == "amp3")
        assert abs(amp3["dropout_rate"][0] - 1 / 3) < 0.01

        # amp4 has all dropouts = 100%
        amp4 = stats.filter(pl.col("amplicon_name") == "amp4")
        assert amp4["dropout_rate"][0] == 1.0

    def test_assigns_performance_tiers(self, amplicon_summary_tsv: Path) -> None:
        """Test that performance tiers are assigned based on median."""
        data = load_amplicon_summary(amplicon_summary_tsv)
        stats = prepare_amplicon_stats(data)

        # amp1 should be "good" (highest median)
        amp1 = stats.filter(pl.col("amplicon_name") == "amp1")
        assert amp1["performance_tier"][0] == "good"

        # amp4 should be "poor" (zero reads)
        amp4 = stats.filter(pl.col("amplicon_name") == "amp4")
        assert amp4["performance_tier"][0] == "poor"

    def test_empty_data_returns_empty_stats(self) -> None:
        """Test that empty input returns empty stats DataFrame."""
        empty_df = pl.DataFrame(
            schema={
                "sample_name": pl.Utf8,
                "amplicon_name": pl.Utf8,
                "reads": pl.Int64,
            },
        )
        stats = prepare_amplicon_stats(empty_df)
        assert len(stats) == 0

    def test_sample_count(self, amplicon_summary_tsv: Path) -> None:
        """Test that sample count is correct."""
        data = load_amplicon_summary(amplicon_summary_tsv)
        stats = prepare_amplicon_stats(data)

        # Each amplicon should have 3 samples
        assert stats["sample_count"].to_list() == [3, 3, 3, 3]


class TestAmpliconRankingBar:
    """Test amplicon ranking bar chart generation."""

    def test_generates_html_file(
        self,
        amplicon_summary_tsv: Path,
        tmp_path: Path,
    ) -> None:
        """Test that HTML file is generated."""
        output_path = tmp_path / "ranking"
        paths = amplicon_ranking_bar(
            amplicon_summary_tsv,
            output_path,
            formats=["html"],
        )

        assert len(paths) == 1
        assert paths[0].suffix == ".html"
        assert paths[0].exists()

    def test_html_contains_amplicon_names(
        self,
        amplicon_summary_tsv: Path,
        tmp_path: Path,
    ) -> None:
        """Test that HTML contains amplicon names."""
        output_path = tmp_path / "ranking"
        paths = amplicon_ranking_bar(
            amplicon_summary_tsv,
            output_path,
            formats=["html"],
        )

        content = paths[0].read_text()
        assert "amp1" in content
        assert "amp2" in content
        assert "amp3" in content
        assert "amp4" in content

    def test_empty_summary_returns_empty_list(
        self,
        empty_summary: Path,
        tmp_path: Path,
    ) -> None:
        """Test that empty summary returns no files."""
        output_path = tmp_path / "ranking"
        paths = amplicon_ranking_bar(empty_summary, output_path, formats=["html"])
        assert paths == []

    def test_single_sample(self, single_sample_summary: Path, tmp_path: Path) -> None:
        """Test chart generation with single sample."""
        output_path = tmp_path / "ranking"
        paths = amplicon_ranking_bar(
            single_sample_summary,
            output_path,
            formats=["html"],
        )

        assert len(paths) == 1
        assert paths[0].exists()


class TestAmpliconDropoutScatter:
    """Test amplicon dropout scatter plot generation."""

    def test_generates_html_file(
        self,
        amplicon_summary_tsv: Path,
        tmp_path: Path,
    ) -> None:
        """Test that HTML file is generated."""
        output_path = tmp_path / "dropout"
        paths = amplicon_dropout_scatter(
            amplicon_summary_tsv,
            output_path,
            formats=["html"],
        )

        assert len(paths) == 1
        assert paths[0].suffix == ".html"
        assert paths[0].exists()

    def test_html_contains_amplicon_names(
        self,
        amplicon_summary_tsv: Path,
        tmp_path: Path,
    ) -> None:
        """Test that HTML contains amplicon names."""
        output_path = tmp_path / "dropout"
        paths = amplicon_dropout_scatter(
            amplicon_summary_tsv,
            output_path,
            formats=["html"],
        )

        content = paths[0].read_text()
        assert "amp1" in content
        assert "amp4" in content

    def test_empty_summary_returns_empty_list(
        self,
        empty_summary: Path,
        tmp_path: Path,
    ) -> None:
        """Test that empty summary returns no files."""
        output_path = tmp_path / "dropout"
        paths = amplicon_dropout_scatter(empty_summary, output_path, formats=["html"])
        assert paths == []

    def test_interactive_chart(
        self,
        amplicon_summary_tsv: Path,
        tmp_path: Path,
    ) -> None:
        """Test that chart is interactive (contains selection params)."""
        output_path = tmp_path / "dropout"
        paths = amplicon_dropout_scatter(
            amplicon_summary_tsv,
            output_path,
            formats=["html"],
        )

        content = paths[0].read_text()
        # Interactive charts have selection/zoom params
        assert "param" in content.lower() or "selection" in content.lower()


class TestPrepareHeatmapData:
    """Test heatmap data preparation."""

    def test_adds_log_reads_column(self, amplicon_summary_tsv: Path) -> None:
        """Test that log_reads column is added."""
        data_lf = pl.scan_csv(amplicon_summary_tsv, separator="\t")
        result = prepare_heatmap_data(data_lf)

        assert "log_reads" in result.columns

    def test_log_reads_handles_zeros(self, amplicon_summary_tsv: Path) -> None:
        """Test that log1p handles zero reads correctly."""
        data_lf = pl.scan_csv(amplicon_summary_tsv, separator="\t")
        result = prepare_heatmap_data(data_lf)

        # Rows with reads=0 should have log_reads=0 (log1p(0) = 0)
        zero_rows = result.filter(pl.col("reads") == 0)
        assert zero_rows["log_reads"].to_list() == [0.0] * len(zero_rows)

    def test_preserves_original_reads(self, amplicon_summary_tsv: Path) -> None:
        """Test that original reads column is preserved."""
        data_lf = pl.scan_csv(amplicon_summary_tsv, separator="\t")
        result = prepare_heatmap_data(data_lf)

        assert "reads" in result.columns
        # Check a known value
        amp1_a = result.filter(
            (pl.col("amplicon_name") == "amp1") & (pl.col("sample_name") == "sampleA"),
        )
        assert amp1_a["reads"][0] == 500


class TestAmpliconHeatmap:
    """Test amplicon heatmap generation."""

    def test_generates_html_file(
        self,
        amplicon_summary_tsv: Path,
        tmp_path: Path,
    ) -> None:
        """Test that HTML file is generated."""
        output_path = tmp_path / "heatmap"
        paths = amplicon_heatmap(amplicon_summary_tsv, output_path, formats=["html"])

        assert len(paths) == 1
        assert paths[0].suffix == ".html"
        assert paths[0].exists()

    def test_html_contains_sample_names(
        self,
        amplicon_summary_tsv: Path,
        tmp_path: Path,
    ) -> None:
        """Test that HTML contains sample names."""
        output_path = tmp_path / "heatmap"
        paths = amplicon_heatmap(amplicon_summary_tsv, output_path, formats=["html"])

        content = paths[0].read_text()
        assert "sampleA" in content
        assert "sampleB" in content
        assert "sampleC" in content

    def test_html_contains_amplicon_names(
        self,
        amplicon_summary_tsv: Path,
        tmp_path: Path,
    ) -> None:
        """Test that HTML contains amplicon names."""
        output_path = tmp_path / "heatmap"
        paths = amplicon_heatmap(amplicon_summary_tsv, output_path, formats=["html"])

        content = paths[0].read_text()
        assert "amp1" in content
        assert "amp2" in content
        assert "amp3" in content
        assert "amp4" in content

    def test_empty_summary_returns_empty_list(
        self,
        empty_summary: Path,
        tmp_path: Path,
    ) -> None:
        """Test that empty summary returns no files."""
        output_path = tmp_path / "heatmap"
        paths = amplicon_heatmap(empty_summary, output_path, formats=["html"])
        assert paths == []

    def test_single_sample(self, single_sample_summary: Path, tmp_path: Path) -> None:
        """Test heatmap generation with single sample."""
        output_path = tmp_path / "heatmap"
        paths = amplicon_heatmap(single_sample_summary, output_path, formats=["html"])

        assert len(paths) == 1
        assert paths[0].exists()
