#!/usr/bin/env python3
"""
Comprehensive test suite for plot_coverage.py.

Tests include:
- Command-line argument parsing
- BED file parsing and data transformation
- Plot generation with different configurations
- Edge cases (empty data, single contig, multiple contigs)
- Coverage percentage calculations
- File output verification
"""

from pathlib import Path
from unittest.mock import MagicMock, Mock, patch

import polars as pl
import pytest


@pytest.fixture
def mock_bed_content():
    """Create mock BED file content for testing."""
    return """chr1\t0\t100\t50
chr1\t100\t200\t25
chr1\t200\t300\t10
chr1\t300\t400\t5
chr2\t0\t50\t30
chr2\t50\t100\t15"""


@pytest.fixture
def mock_coverage_lazyframe():
    """Create a mock LazyFrame for coverage data."""
    df = pl.DataFrame(
        {
            "chromosome": ["chr1", "chr1", "chr1", "chr2", "chr2"],
            "position": [50, 150, 250, 25, 75],
            "coverage": [50, 25, 10, 30, 15],
            "reference length": [400, 400, 400, 100, 100],
            "passes_depth_cutoff": [True, True, True, True, True],
        }
    )
    return df.lazy()


@pytest.fixture
def temp_bed_file(tmp_path, mock_bed_content):
    """Create a temporary BED file."""
    bed_file = tmp_path / "test.per-base.bed"
    bed_file.write_text(mock_bed_content)
    return bed_file


class TestCommandLineArgs:
    """Test command-line argument parsing."""

    def test_parse_args_all_options(self, monkeypatch):
        """Test parsing all command-line arguments."""
        from bin.plot_coverage import parse_command_line_args

        test_args = [
            "plot_coverage.py",
            "--input",
            "/path/to/coverage.bed",
            "--label",
            "test_sample",
            "--depth",
            "30",
        ]
        monkeypatch.setattr("sys.argv", test_args)

        args = parse_command_line_args()
        assert args.input == Path("/path/to/coverage.bed")
        assert args.label == "test_sample"
        assert args.depth == 30

    def test_parse_args_defaults(self, monkeypatch):
        """Test default values for optional arguments."""
        from bin.plot_coverage import parse_command_line_args

        test_args = ["plot_coverage.py", "--input", "/path/to/coverage.bed"]
        monkeypatch.setattr("sys.argv", test_args)

        args = parse_command_line_args()
        assert args.input == Path("/path/to/coverage.bed")
        assert args.label == "sample"  # default
        assert args.depth == 20  # default

    def test_missing_required_input(self, monkeypatch):
        """Test error when required input is missing."""
        from bin.plot_coverage import parse_command_line_args

        test_args = ["plot_coverage.py"]
        monkeypatch.setattr("sys.argv", test_args)

        with pytest.raises(SystemExit):
            parse_command_line_args()


class TestBedFileProcessing:
    """Test BED file parsing and data transformation."""

    def test_adapt_bed_to_df_success(self, temp_bed_file):
        """Test successful BED file adaptation to LazyFrame."""
        from bin.plot_coverage import adapt_bed_to_df

        result_lf = adapt_bed_to_df(temp_bed_file, depth=20)
        result_df = result_lf.collect()

        assert isinstance(result_lf, pl.LazyFrame)
        assert "chromosome" in result_df.columns
        assert "position" in result_df.columns
        assert "coverage" in result_df.columns
        assert "reference length" in result_df.columns
        assert "passes_depth_cutoff" in result_df.columns

        # Check that positions were exploded correctly
        assert result_df.shape[0] > 0

        # Verify reference length calculation
        chr1_data = result_df.filter(pl.col("chromosome") == "chr1")
        assert chr1_data["reference length"][0] == 400  # 0-400 for chr1

    def test_adapt_bed_to_df_depth_cutoff(self, temp_bed_file):
        """Test depth cutoff determination."""
        from bin.plot_coverage import adapt_bed_to_df

        # Test with high depth threshold
        result_lf = adapt_bed_to_df(temp_bed_file, depth=100)
        result_df = result_lf.collect()

        # Check that chromosomes with max coverage < 100 are marked as not passing
        chr1_data = (
            result_df.filter(pl.col("chromosome") == "chr1")
            .select("passes_depth_cutoff")
            .unique()
        )
        chr2_data = (
            result_df.filter(pl.col("chromosome") == "chr2")
            .select("passes_depth_cutoff")
            .unique()
        )

        # chr1 max coverage is 50, chr2 max is 30, both < 100
        assert not chr1_data["passes_depth_cutoff"][0]
        assert not chr2_data["passes_depth_cutoff"][0]

    def test_adapt_bed_empty_file(self, tmp_path):
        """Test handling of empty BED file."""
        from bin.plot_coverage import adapt_bed_to_df

        empty_bed = tmp_path / "empty.bed"
        empty_bed.write_text("")

        result_lf = adapt_bed_to_df(empty_bed, depth=20)
        result_df = result_lf.collect()

        assert result_df.shape[0] == 0

    def test_position_explosion(self, tmp_path):
        """Test that positions are correctly exploded from ranges."""
        from bin.plot_coverage import adapt_bed_to_df

        # Create a simple BED file with known range
        bed_file = tmp_path / "range_test.bed"
        bed_file.write_text("chr1\t0\t3\t10\n")  # 3 positions: 0, 1, 2

        result_lf = adapt_bed_to_df(bed_file, depth=5)
        result_df = result_lf.collect()

        assert result_df.shape[0] == 3
        positions = sorted(result_df["position"].to_list())
        assert positions == [0, 1, 2]


class TestPlotConstruction:
    """Test plot construction and customization."""

    @patch("bin.plot_coverage.ggplot")
    @patch("bin.plot_coverage.geom_rect")
    @patch("bin.plot_coverage.geom_line")
    @patch("bin.plot_coverage.geom_hline")
    def test_construct_plot_basic(
        self, mock_hline, mock_line, mock_rect, mock_ggplot, mock_coverage_lazyframe
    ):
        """Test basic plot construction."""
        from bin.plot_coverage import construct_plot

        mock_plot = MagicMock()
        mock_ggplot.return_value = mock_plot

        result = construct_plot(mock_coverage_lazyframe, "test_sample", 20)

        # Verify plot components were called
        mock_ggplot.assert_called_once()
        mock_rect.assert_called_once()  # For low coverage regions
        mock_line.assert_called_once()
        mock_hline.assert_called_once_with(yintercept=20, linetype="dashed")

    def test_construct_plot_low_coverage_detection(self):
        """Test detection of low coverage regions."""
        from bin.plot_coverage import construct_plot

        # Create data with low coverage regions
        df = pl.DataFrame(
            {
                "chromosome": ["chr1"] * 10,
                "position": list(range(10)),
                "coverage": [5, 5, 5, 25, 25, 5, 5, 25, 25, 25],
                "reference length": [10] * 10,
                "passes_depth_cutoff": [True] * 10,
            }
        )

        with patch("bin.plot_coverage.ggplot") as mock_ggplot:
            with patch("bin.plot_coverage.geom_rect") as mock_rect:
                mock_plot = MagicMock()
                mock_ggplot.return_value = mock_plot

                result = construct_plot(df.lazy(), "test", 20)

                # Verify geom_rect was called for low coverage visualization
                mock_rect.assert_called_once()

                # Check the data passed to geom_rect
                rect_call_args = mock_rect.call_args
                low_cov_data = rect_call_args[1]["data"]

                # Should have detected the low coverage regions
                assert low_cov_data.shape[0] > 0

    def test_construct_plot_all_high_coverage(self):
        """Test plot when all positions have high coverage."""
        from bin.plot_coverage import construct_plot

        # Create data with all high coverage
        df = pl.DataFrame(
            {
                "chromosome": ["chr1"] * 5,
                "position": list(range(5)),
                "coverage": [100, 150, 200, 250, 300],
                "reference length": [5] * 5,
                "passes_depth_cutoff": [True] * 5,
            }
        )

        with patch("bin.plot_coverage.ggplot") as mock_ggplot:
            mock_plot = MagicMock()
            mock_ggplot.return_value = mock_plot

            # This should raise assertion error as no low coverage regions exist
            with pytest.raises(AssertionError):
                construct_plot(df.lazy(), "test", 20)


class TestPlotFinishing:
    """Test plot finishing with different contig configurations."""

    def test_finish_plot_single_contig(self):
        """Test plot finishing for single contig."""
        from bin.plot_coverage import finish_plot

        mock_plot = MagicMock()

        with patch("bin.plot_coverage.facet_wrap") as mock_facet:
            mock_facet.return_value = MagicMock()

            result = finish_plot(mock_plot, 1)

            # For single contig, should return single plot
            assert not isinstance(result, tuple)
            mock_facet.assert_called_once_with("~chromosome", scales="free_x")

    def test_finish_plot_multiple_contigs(self):
        """Test plot finishing for multiple contigs."""
        from bin.plot_coverage import finish_plot

        mock_plot = MagicMock()

        with patch("bin.plot_coverage.facet_wrap") as mock_facet:
            mock_facet.return_value = MagicMock()

            result = finish_plot(mock_plot, 3)

            # For multiple contigs, should return tuple of two plots
            assert isinstance(result, tuple)
            assert len(result) == 2

            # Verify facet_wrap was called twice with different scales
            assert mock_facet.call_count == 2
            calls = mock_facet.call_args_list
            assert calls[0][0] == ("~chromosome",)
            assert calls[0][1] == {"scales": "free_x"}
            assert calls[1][0] == ("~chromosome",)
            assert calls[1][1] == {"scales": "free"}

    def test_finish_plot_invalid_contig_count(self):
        """Test error handling for invalid contig count."""
        from bin.plot_coverage import finish_plot

        mock_plot = MagicMock()

        with pytest.raises(AssertionError):
            finish_plot(mock_plot, 0)

        with pytest.raises(AssertionError):
            finish_plot(mock_plot, -1)


class TestCoveragePercentageCalculation:
    """Test coverage percentage calculations."""

    def test_compute_perc_cov_basic(self, mock_coverage_lazyframe):
        """Test basic percentage coverage calculation."""
        from bin.plot_coverage import compute_perc_cov

        result_lf = compute_perc_cov(mock_coverage_lazyframe, "test_sample", 2, 20)
        result_df = result_lf.collect()

        assert "sample id" in result_df.columns
        assert "chromosome" in result_df.columns
        assert "reference length" in result_df.columns
        assert "proportion ≥ 20X coverage" in result_df.columns

        # Check sample ID was set correctly
        assert all(result_df["sample id"] == "test_sample")

        # Should have one row per chromosome
        assert result_df.shape[0] == 2

    def test_compute_perc_cov_with_low_coverage(self):
        """Test percentage calculation with low coverage regions."""
        from bin.plot_coverage import compute_perc_cov

        # Create data with mixed coverage
        df = pl.DataFrame(
            {
                "chromosome": ["chr1"] * 10,
                "position": list(range(10)),
                "coverage": [5, 5, 25, 25, 25, 25, 25, 5, 5, 5],  # 5/10 positions >= 20
                "reference length": [10] * 10,
                "passes_depth_cutoff": [True] * 10,
            }
        )

        result_lf = compute_perc_cov(df.lazy(), "test", 1, 20)
        result_df = result_lf.collect()

        # Should calculate 50% coverage (5 out of 10 positions)
        proportion = result_df["proportion ≥ 20X coverage"][0]
        assert proportion == pytest.approx(0.5)

    def test_compute_perc_cov_all_below_cutoff(self):
        """Test when entire chromosome is below cutoff."""
        from bin.plot_coverage import compute_perc_cov

        df = pl.DataFrame(
            {
                "chromosome": ["chr1"] * 5,
                "position": list(range(5)),
                "coverage": [10, 10, 10, 10, 10],  # All below 20
                "reference length": [5] * 5,
                "passes_depth_cutoff": [False] * 5,  # Doesn't pass cutoff
            }
        )

        result_lf = compute_perc_cov(df.lazy(), "test", 1, 20)
        result_df = result_lf.collect()

        # Should be 0% coverage
        proportion = result_df["proportion ≥ 20X coverage"][0]
        assert proportion == 0.0

    def test_compute_perc_cov_assertion_error(self):
        """Test assertion error when row count doesn't match contig count."""
        from bin.plot_coverage import compute_perc_cov

        df = pl.DataFrame(
            {
                "chromosome": ["chr1", "chr2"],
                "position": [1, 1],
                "coverage": [20, 20],
                "reference length": [1, 1],
                "passes_depth_cutoff": [True, True],
            }
        )

        # Expect 3 contigs but only have 2
        with pytest.raises(AssertionError):
            compute_perc_cov(df.lazy(), "test", 3, 20).collect()


class TestMainFunction:
    """Test main function orchestration."""

    @patch("bin.plot_coverage.compute_perc_cov")
    @patch("bin.plot_coverage.ggsave")
    @patch("bin.plot_coverage.finish_plot")
    @patch("bin.plot_coverage.construct_plot")
    @patch("bin.plot_coverage.adapt_bed_to_df")
    @patch("bin.plot_coverage.parse_command_line_args")
    @patch("pathlib.Path.is_file")
    def test_main_single_contig(
        self,
        mock_is_file,
        mock_parse_args,
        mock_adapt,
        mock_construct,
        mock_finish,
        mock_save,
        mock_compute,
        mock_coverage_lazyframe,
    ):
        """Test main function with single contig."""
        from bin.plot_coverage import main

        # Mock setup
        mock_args = Mock()
        mock_args.input = Path("/test/coverage.bed")
        mock_args.label = "test_sample"
        mock_args.depth = 20
        mock_parse_args.return_value = mock_args
        mock_is_file.return_value = True

        mock_adapt.return_value = mock_coverage_lazyframe
        mock_plot = MagicMock()
        mock_construct.return_value = mock_plot
        mock_finish.return_value = mock_plot  # Single plot for single contig

        mock_perc_lf = MagicMock()
        mock_perc_df = MagicMock()
        mock_perc_df.write_csv = MagicMock()
        mock_perc_lf.collect.return_value = mock_perc_df
        mock_compute.return_value = mock_perc_lf

        # Run main
        main()

        # Verify calls
        mock_parse_args.assert_called_once()
        mock_adapt.assert_called_once_with(Path("/test/coverage.bed"), 20)
        mock_construct.assert_called_once()
        mock_finish.assert_called_once()

        # Should save single plot
        mock_save.assert_called_once_with(
            mock_plot, "test_sample.coverage.pdf", format="pdf", height=6, width=11
        )

        # Verify coverage percentage was computed and saved
        mock_compute.assert_called_once()
        mock_perc_df.write_csv.assert_called_once_with(
            "test_sample.passing_cov.tsv", separator="\t"
        )

    @patch("bin.plot_coverage.compute_perc_cov")
    @patch("bin.plot_coverage.ggsave")
    @patch("bin.plot_coverage.finish_plot")
    @patch("bin.plot_coverage.construct_plot")
    @patch("bin.plot_coverage.adapt_bed_to_df")
    @patch("bin.plot_coverage.parse_command_line_args")
    @patch("pathlib.Path.is_file")
    def test_main_multiple_contigs(
        self,
        mock_is_file,
        mock_parse_args,
        mock_adapt,
        mock_construct,
        mock_finish,
        mock_save,
        mock_compute,
        mock_coverage_lazyframe,
    ):
        """Test main function with multiple contigs."""
        from bin.plot_coverage import main

        # Mock setup
        mock_args = Mock()
        mock_args.input = Path("/test/coverage.bed")
        mock_args.label = "multi_sample"
        mock_args.depth = 30
        mock_parse_args.return_value = mock_args
        mock_is_file.return_value = True

        # Mock multiple contigs
        multi_contig_lf = pl.DataFrame(
            {
                "chromosome": ["chr1", "chr2", "chr3"],
                "position": [1, 1, 1],
                "coverage": [20, 20, 20],
                "reference length": [1, 1, 1],
                "passes_depth_cutoff": [True, True, True],
            }
        ).lazy()

        mock_adapt.return_value = multi_contig_lf
        mock_plot = MagicMock()
        mock_construct.return_value = mock_plot

        # Return tuple of plots for multiple contigs
        mock_plot1 = MagicMock()
        mock_plot2 = MagicMock()
        mock_finish.return_value = (mock_plot1, mock_plot2)

        mock_perc_lf = MagicMock()
        mock_perc_df = MagicMock()
        mock_perc_df.write_csv = MagicMock()
        mock_perc_lf.collect.return_value = mock_perc_df
        mock_compute.return_value = mock_perc_lf

        # Run main
        main()

        # Should save two plots
        assert mock_save.call_count == 2
        calls = mock_save.call_args_list

        # Check both plots were saved with correct names
        assert calls[0][0][0] == mock_plot1
        assert calls[0][0][1] == "multi_sample.fixed_y.coverage.pdf"
        assert calls[1][0][0] == mock_plot2
        assert calls[1][0][1] == "multi_sample.free_scales.coverage.pdf"

    def test_main_file_not_exists(self, monkeypatch):
        """Test error when input file doesn't exist."""
        from bin.plot_coverage import main

        test_args = ["plot_coverage.py", "--input", "/non/existent/file.bed"]
        monkeypatch.setattr("sys.argv", test_args)

        with pytest.raises(AssertionError):
            main()


class TestEdgeCases:
    """Test edge cases and error conditions."""

    def test_empty_chromosome(self, tmp_path):
        """Test handling of chromosome with no data."""
        from bin.plot_coverage import adapt_bed_to_df

        # BED file with empty region
        bed_file = tmp_path / "empty_chr.bed"
        bed_file.write_text("chr1\t0\t0\t0\n")  # Empty range

        result_lf = adapt_bed_to_df(bed_file, 20)
        result_df = result_lf.collect()

        # Should handle empty range gracefully
        assert result_df.shape[0] == 0

    def test_very_large_positions(self, tmp_path):
        """Test handling of very large genomic positions."""
        from bin.plot_coverage import adapt_bed_to_df

        bed_file = tmp_path / "large_pos.bed"
        bed_file.write_text("chr1\t1000000000\t1000000010\t50\n")

        result_lf = adapt_bed_to_df(bed_file, 20)
        result_df = result_lf.collect()

        assert result_df.shape[0] == 10  # 10 positions
        assert result_df["position"].min() == 1000000000
        assert result_df["position"].max() == 1000000009

    def test_malformed_bed_file(self, tmp_path):
        """Test handling of malformed BED file."""
        from bin.plot_coverage import adapt_bed_to_df

        bed_file = tmp_path / "malformed.bed"
        bed_file.write_text("chr1\tnot_a_number\t100\t50\n")

        with pytest.raises(Exception):
            adapt_bed_to_df(bed_file, 20).collect()

    def test_unicode_chromosome_names(self, tmp_path):
        """Test handling of unicode in chromosome names."""
        from bin.plot_coverage import adapt_bed_to_df

        bed_file = tmp_path / "unicode.bed"
        bed_file.write_text("chrÜ\t0\t10\t20\n")

        result_lf = adapt_bed_to_df(bed_file, 20)
        result_df = result_lf.collect()

        assert result_df.shape[0] == 10
        assert result_df["chromosome"][0] == "chrÜ"
