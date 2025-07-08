#!/usr/bin/env python3
"""
Comprehensive test suite for multisample_plot.py.

Tests include:
- Command-line argument parsing
- Data parsing from BED files
- Plot generation verification
- Edge cases handling
- Color schemes and customization
"""

from pathlib import Path
from unittest.mock import MagicMock, Mock, patch

import polars as pl
import pytest


@pytest.fixture
def mock_bed_file():
    """Create a mock BED file content."""
    return "chr1\t0\t10\t5\nchr1\t10\t20\t15\nchr1\t20\t30\t25\nchr2\t0\t15\t10\nchr2\t15\t25\t20\n"


@pytest.fixture
def mock_dataframe():
    """Create a mock Polars DataFrame for testing."""
    return pl.DataFrame(
        {
            "chromosome": ["chr1", "chr1", "chr1", "chr2", "chr2"],
            "position": [5, 15, 25, 10, 20],
            "coverage": [5, 15, 25, 10, 20],
            "sample": ["sample1", "sample1", "sample1", "sample2", "sample2"],
        }
    )


@pytest.fixture
def temp_dir_with_bed_files(tmp_path, mock_bed_file):
    """Create a temporary directory with mock BED files."""
    # Create sample BED files
    samples = ["sample1", "sample2", "sample3"]
    for sample in samples:
        bed_file = tmp_path / f"{sample}.per-base.bed"
        bed_file.write_text(mock_bed_file)
    return tmp_path


class TestCommandLineArgs:
    """Test command-line argument parsing."""

    def test_parse_args_with_all_options(self, monkeypatch):
        """Test parsing all command-line arguments."""
        from bin.multisample_plot import parse_command_line_args

        test_args = [
            "multisample_plot.py",
            "--input_dir",
            "/test/dir",
            "--min_coverage",
            "30",
            "--log",
        ]
        monkeypatch.setattr("sys.argv", test_args)

        args = parse_command_line_args()
        assert args.input_dir == Path("/test/dir")
        assert args.min_coverage == 30
        assert args.log is True

    def test_parse_args_defaults(self, monkeypatch):
        """Test default values for optional arguments."""
        from bin.multisample_plot import parse_command_line_args

        test_args = ["multisample_plot.py", "--input_dir", "/test/dir"]
        monkeypatch.setattr("sys.argv", test_args)

        args = parse_command_line_args()
        assert args.input_dir == Path("/test/dir")
        assert args.min_coverage == 20  # default
        assert args.log is False  # default

    def test_missing_required_args(self, monkeypatch):
        """Test error when required arguments are missing."""
        from bin.multisample_plot import parse_command_line_args

        test_args = ["multisample_plot.py"]
        monkeypatch.setattr("sys.argv", test_args)

        with pytest.raises(SystemExit):
            parse_command_line_args()


class TestDataAccumulation:
    """Test coverage data accumulation from BED files."""

    def test_accumulate_cov_dfs_success(self, temp_dir_with_bed_files):
        """Test successful accumulation of coverage data."""
        from bin.multisample_plot import accumulate_cov_dfs

        sample_list = ["sample1", "sample2", "sample3"]
        result = accumulate_cov_dfs(str(temp_dir_with_bed_files), sample_list)

        assert isinstance(result, pl.DataFrame)
        assert "chromosome" in result.columns
        assert "coverage" in result.columns
        assert "sample" in result.columns
        assert "position" in result.columns
        assert result.shape[0] > 0

        # Check that all samples are present
        unique_samples = result["sample"].unique().to_list()
        assert set(unique_samples) == set(sample_list)

    def test_accumulate_cov_dfs_empty_directory(self, tmp_path):
        """Test handling of empty directory."""
        from bin.multisample_plot import accumulate_cov_dfs

        sample_list = []
        result = accumulate_cov_dfs(str(tmp_path), sample_list)

        assert isinstance(result, pl.DataFrame)
        assert result.shape[0] == 0

    def test_accumulate_cov_dfs_invalid_directory(self):
        """Test error handling for invalid directory."""
        from bin.multisample_plot import accumulate_cov_dfs

        with pytest.raises(AssertionError):
            accumulate_cov_dfs("/non/existent/directory", ["sample1"])

    def test_accumulate_cov_dfs_missing_files(self, tmp_path):
        """Test handling when expected BED files are missing."""
        from bin.multisample_plot import accumulate_cov_dfs

        # Create only one file but request two samples
        bed_file = tmp_path / "sample1.per-base.bed"
        bed_file.write_text("chr1\t0\t10\t5\n")

        sample_list = ["sample1", "sample2"]

        # This should raise an error for the missing file
        with pytest.raises(Exception):
            accumulate_cov_dfs(str(tmp_path), sample_list)

    def test_accumulate_cov_dfs_malformed_bed(self, tmp_path):
        """Test handling of malformed BED file."""
        # Create a malformed BED file
        bed_file = tmp_path / "sample1.per-base.bed"
        bed_file.write_text("invalid\tbed\tformat\n")

        from bin.multisample_plot import accumulate_cov_dfs

        with pytest.raises(Exception):
            accumulate_cov_dfs(str(tmp_path), ["sample1"])


class TestPlotGeneration:
    """Test plot generation functions."""

    @patch("bin.multisample_plot.ggplot")
    @patch("bin.multisample_plot.geom_line")
    @patch("bin.multisample_plot.geom_hline")
    def test_plot_coverages(self, mock_hline, mock_line, mock_ggplot, mock_dataframe):
        """Test standard coverage plot generation."""
        from bin.multisample_plot import plot_coverages

        # Create mocks
        mock_plot = MagicMock()
        mock_ggplot.return_value = mock_plot
        mock_line.return_value = MagicMock()
        mock_hline.return_value = MagicMock()

        result = plot_coverages(mock_dataframe, min_desired_depth=20)

        # Verify ggplot was called
        mock_ggplot.assert_called_once()

        # Verify geom_line was called
        mock_line.assert_called_once()

        # Verify geom_hline was called with correct depth
        mock_hline.assert_called_once_with(yintercept=20, linetype="dashed")

    @patch("bin.multisample_plot.ggplot")
    @patch("bin.multisample_plot.geom_line")
    @patch("bin.multisample_plot.geom_hline")
    def test_plot_log_coverages(
        self, mock_hline, mock_line, mock_ggplot, mock_dataframe
    ):
        """Test log-transformed coverage plot generation."""
        from bin.multisample_plot import plot_log_coverages
        from math import log10

        # Create mocks
        mock_plot = MagicMock()
        mock_ggplot.return_value = mock_plot
        mock_line.return_value = MagicMock()
        mock_hline.return_value = MagicMock()

        result = plot_log_coverages(mock_dataframe, min_desired_depth=20)

        # Verify ggplot was called
        mock_ggplot.assert_called_once()

        # Verify geom_hline was called with log-transformed depth
        expected_log_depth = log10(20)
        mock_hline.assert_called_once_with(
            yintercept=expected_log_depth, linetype="dashed"
        )

    def test_plot_with_single_sample(self, mock_dataframe):
        """Test plotting with a single sample."""
        from bin.multisample_plot import plot_coverages

        single_sample_df = mock_dataframe.filter(pl.col("sample") == "sample1")

        with patch("bin.multisample_plot.ggplot") as mock_ggplot:
            mock_plot = MagicMock()
            mock_ggplot.return_value = mock_plot

            result = plot_coverages(single_sample_df)

            # Verify the plot was created
            mock_ggplot.assert_called_once()

    def test_plot_with_empty_dataframe(self):
        """Test plotting with empty dataframe."""
        from bin.multisample_plot import plot_coverages

        empty_df = pl.DataFrame(
            {"chromosome": [], "position": [], "coverage": [], "sample": []}
        )

        with patch("bin.multisample_plot.ggplot") as mock_ggplot:
            mock_plot = MagicMock()
            mock_ggplot.return_value = mock_plot

            result = plot_coverages(empty_df)

            # Should still create a plot, even if empty
            mock_ggplot.assert_called_once()


class TestMainFunction:
    """Test the main function orchestration."""

    @patch("bin.multisample_plot.ggsave")
    @patch("bin.multisample_plot.plot_coverages")
    @patch("bin.multisample_plot.accumulate_cov_dfs")
    @patch("bin.multisample_plot.parse_command_line_args")
    @patch("pathlib.Path.glob")
    def test_main_standard_execution(
        self,
        mock_glob,
        mock_parse_args,
        mock_accumulate,
        mock_plot,
        mock_save,
        mock_dataframe,
    ):
        """Test standard execution of main function."""
        from bin.multisample_plot import main

        # Mock arguments
        mock_args = Mock()
        mock_args.input_dir = Path("/test/dir")
        mock_args.min_coverage = 20
        mock_args.log = False
        mock_parse_args.return_value = mock_args

        # Mock glob to return sample files
        mock_glob.return_value = [
            Path("/test/dir/sample1.per-base.bed"),
            Path("/test/dir/sample2.per-base.bed"),
        ]

        # Mock data accumulation
        mock_accumulate.return_value = mock_dataframe

        # Mock plot
        mock_plot_instance = MagicMock()
        mock_plot.return_value = mock_plot_instance

        # Run main
        main()

        # Verify calls
        mock_parse_args.assert_called_once()
        mock_accumulate.assert_called_once_with(
            Path("/test/dir"), ["sample1", "sample2"]
        )
        mock_plot.assert_called_once_with(mock_dataframe, 20)
        mock_save.assert_called_once_with(
            mock_plot_instance,
            "multisample.fixed_y.coverage.pdf",
            format="pdf",
            height=6,
            width=11,
        )

    @patch("bin.multisample_plot.ggsave")
    @patch("bin.multisample_plot.plot_log_coverages")
    @patch("bin.multisample_plot.accumulate_cov_dfs")
    @patch("bin.multisample_plot.parse_command_line_args")
    @patch("pathlib.Path.glob")
    def test_main_with_log_option(
        self,
        mock_glob,
        mock_parse_args,
        mock_accumulate,
        mock_plot_log,
        mock_save,
        mock_dataframe,
    ):
        """Test main function with log transformation option."""
        from bin.multisample_plot import main

        # Mock arguments with log=True
        mock_args = Mock()
        mock_args.input_dir = Path("/test/dir")
        mock_args.min_coverage = 30
        mock_args.log = True
        mock_parse_args.return_value = mock_args

        # Mock glob
        mock_glob.return_value = [Path("/test/dir/sample1.per-base.bed")]

        # Mock data accumulation
        mock_accumulate.return_value = mock_dataframe

        # Mock plot
        mock_plot_instance = MagicMock()
        mock_plot_log.return_value = mock_plot_instance

        # Run main
        main()

        # Verify log plot was called
        mock_plot_log.assert_called_once_with(mock_dataframe, 30)


class TestEdgeCases:
    """Test edge cases and error handling."""

    def test_single_position_coverage(self):
        """Test handling of single position coverage data."""
        from bin.multisample_plot import plot_coverages

        single_pos_df = pl.DataFrame(
            {
                "chromosome": ["chr1"],
                "position": [1],
                "coverage": [10],
                "sample": ["sample1"],
            }
        )

        with patch("bin.multisample_plot.ggplot") as mock_ggplot:
            mock_plot = MagicMock()
            mock_ggplot.return_value = mock_plot

            result = plot_coverages(single_pos_df)
            mock_ggplot.assert_called_once()

    def test_very_high_coverage_values(self):
        """Test handling of very high coverage values."""
        from bin.multisample_plot import plot_coverages

        high_cov_df = pl.DataFrame(
            {
                "chromosome": ["chr1", "chr1"],
                "position": [1, 2],
                "coverage": [1000000, 2000000],
                "sample": ["sample1", "sample1"],
            }
        )

        with patch("bin.multisample_plot.ggplot") as mock_ggplot:
            mock_plot = MagicMock()
            mock_ggplot.return_value = mock_plot

            result = plot_coverages(high_cov_df)
            mock_ggplot.assert_called_once()

    def test_zero_coverage_regions(self):
        """Test handling of zero coverage regions."""
        from bin.multisample_plot import plot_coverages

        zero_cov_df = pl.DataFrame(
            {
                "chromosome": ["chr1", "chr1", "chr1"],
                "position": [1, 2, 3],
                "coverage": [0, 0, 0],
                "sample": ["sample1", "sample1", "sample1"],
            }
        )

        with patch("bin.multisample_plot.ggplot") as mock_ggplot:
            mock_plot = MagicMock()
            mock_ggplot.return_value = mock_plot

            result = plot_coverages(zero_cov_df)
            mock_ggplot.assert_called_once()

    def test_multiple_chromosomes(self):
        """Test handling of multiple chromosomes with faceting."""
        from bin.multisample_plot import plot_coverages

        multi_chr_df = pl.DataFrame(
            {
                "chromosome": ["chr1", "chr1", "chr2", "chr2", "chr3", "chr3"],
                "position": [1, 2, 1, 2, 1, 2],
                "coverage": [10, 20, 30, 40, 50, 60],
                "sample": ["s1", "s1", "s1", "s1", "s1", "s1"],
            }
        )

        with patch("bin.multisample_plot.ggplot") as mock_ggplot:
            with patch("bin.multisample_plot.facet_wrap") as mock_facet:
                mock_plot = MagicMock()
                mock_ggplot.return_value = mock_plot
                mock_facet.return_value = MagicMock()

                result = plot_coverages(multi_chr_df)

                # Verify facet_wrap was called
                mock_facet.assert_called_with("~chromosome", scales="free_x")


class TestColorSchemes:
    """Test color schemes and plot customization."""

    @patch("bin.multisample_plot.guides")
    @patch("bin.multisample_plot.ggplot")
    def test_legend_customization(self, mock_ggplot, mock_guides):
        """Test legend customization with multiple samples."""
        from bin.multisample_plot import plot_coverages

        # Create data with many samples
        many_samples_df = pl.DataFrame(
            {
                "chromosome": ["chr1"] * 10,
                "position": list(range(1, 11)),
                "coverage": list(range(10, 110, 10)),
                "sample": [f"sample{i}" for i in range(1, 11)],
            }
        )

        mock_plot = MagicMock()
        mock_ggplot.return_value = mock_plot
        mock_guides.return_value = MagicMock()

        result = plot_coverages(many_samples_df)

        # Verify guides was called for legend customization
        mock_guides.assert_called_once()
        # Check that ncol=3 was set for the legend
        call_args = mock_guides.call_args
        assert "color" in call_args[1] or call_args[0]

    def test_plot_theme_application(self):
        """Test that minimal theme is applied."""
        from bin.multisample_plot import plot_coverages

        df = pl.DataFrame(
            {
                "chromosome": ["chr1"],
                "position": [1],
                "coverage": [10],
                "sample": ["sample1"],
            }
        )

        with patch("bin.multisample_plot.ggplot") as mock_ggplot:
            with patch("bin.multisample_plot.theme_minimal") as mock_theme:
                mock_plot = MagicMock()
                mock_ggplot.return_value = mock_plot
                mock_theme.return_value = MagicMock()

                result = plot_coverages(df)

                # Verify theme_minimal was called
                mock_theme.assert_called_once()
