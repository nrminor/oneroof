"""Tests for the amplicon summary generator."""

# Import the module's internal functions for testing
import sys
from pathlib import Path

import polars as pl
import pytest

sys.path.insert(0, str(Path(__file__).parent))

from generate_amplicon_summary import (
    _create_summary,
    _load_positions_from_tsv,
    _parse_stats_files,
)


@pytest.fixture
def primer_pairs_tsv(tmp_path: Path) -> Path:
    """Create a primer_pairs.tsv with 3 amplicons."""
    content = "\n".join(
        [
            "amplicon_name\tfwd_sequence\trev_sequence\tchrom\tamplicon_start\tamplicon_end",
            "amp1\tACGT\tTGCA\tchr1\t100\t200",
            "amp2\tGGGG\tCCCC\tchr1\t300\t400",
            "amp3\tAAAA\tTTTT\tchr1\t500\t600",
        ],
    )
    tsv_file = tmp_path / "primer_pairs.tsv"
    tsv_file.write_text(content)
    return tsv_file


@pytest.fixture
def stats_all_amplicons(tmp_path: Path) -> Path:
    """Create stats file where all amplicons have reads."""
    # seqkit stats output format
    content = "\n".join(
        [
            "file\tformat\ttype\tnum_seqs\tsum_len\tmin_len\tavg_len\tmax_len",
            "sampleA.amp1.fasta.gz\tFASTA\tDNA\t100\t10000\t80\t100\t120",
            "sampleA.amp2.fasta.gz\tFASTA\tDNA\t150\t15000\t80\t100\t120",
            "sampleA.amp3.fasta.gz\tFASTA\tDNA\t200\t20000\t80\t100\t120",
        ],
    )
    stats_file = tmp_path / "stats_sampleA.tsv"
    stats_file.write_text(content)
    return stats_file


@pytest.fixture
def stats_with_dropout(tmp_path: Path) -> Path:
    """Create stats file where amp3 has dropped out (no file = no row)."""
    content = "\n".join(
        [
            "file\tformat\ttype\tnum_seqs\tsum_len\tmin_len\tavg_len\tmax_len",
            "sampleB.amp1.fasta.gz\tFASTA\tDNA\t50\t5000\t80\t100\t120",
            "sampleB.amp2.fasta.gz\tFASTA\tDNA\t75\t7500\t80\t100\t120",
            # amp3 is missing - this is a dropout
        ],
    )
    stats_file = tmp_path / "stats_sampleB.tsv"
    stats_file.write_text(content)
    return stats_file


@pytest.fixture
def stats_multiple_samples(tmp_path: Path) -> tuple[Path, Path]:
    """Create stats files for two samples with different dropout patterns."""
    # Sample A: has all amplicons
    content_a = "\n".join(
        [
            "file\tformat\ttype\tnum_seqs\tsum_len\tmin_len\tavg_len\tmax_len",
            "sampleA.amp1.fasta.gz\tFASTA\tDNA\t100\t10000\t80\t100\t120",
            "sampleA.amp2.fasta.gz\tFASTA\tDNA\t150\t15000\t80\t100\t120",
            "sampleA.amp3.fasta.gz\tFASTA\tDNA\t200\t20000\t80\t100\t120",
        ],
    )
    stats_a = tmp_path / "stats_sampleA.tsv"
    stats_a.write_text(content_a)

    # Sample B: amp3 dropped out
    content_b = "\n".join(
        [
            "file\tformat\ttype\tnum_seqs\tsum_len\tmin_len\tavg_len\tmax_len",
            "sampleB.amp1.fasta.gz\tFASTA\tDNA\t50\t5000\t80\t100\t120",
            "sampleB.amp2.fasta.gz\tFASTA\tDNA\t75\t7500\t80\t100\t120",
        ],
    )
    stats_b = tmp_path / "stats_sampleB.tsv"
    stats_b.write_text(content_b)

    return stats_a, stats_b


class TestLoadPositionsFromTsv:
    """Test primer pairs TSV loading."""

    def test_load_valid_tsv(self, primer_pairs_tsv: Path) -> None:
        """Test loading a valid primer pairs TSV."""
        lf = _load_positions_from_tsv(str(primer_pairs_tsv))
        df = lf.collect()
        assert len(df) == 3
        assert list(df.columns) == ["amplicon_name", "start_pos", "end_pos"]

    def test_amplicon_names(self, primer_pairs_tsv: Path) -> None:
        """Test that amplicon names are correctly extracted."""
        lf = _load_positions_from_tsv(str(primer_pairs_tsv))
        df = lf.collect()
        names = df["amplicon_name"].to_list()
        assert names == ["amp1", "amp2", "amp3"]

    def test_positions(self, primer_pairs_tsv: Path) -> None:
        """Test that positions are correctly extracted."""
        lf = _load_positions_from_tsv(str(primer_pairs_tsv))
        df = lf.collect()
        assert df["start_pos"].to_list() == [100, 300, 500]
        assert df["end_pos"].to_list() == [200, 400, 600]


class TestParseStatsFiles:
    """Test stats file parsing."""

    def test_parse_single_file(self, stats_all_amplicons: Path) -> None:
        """Test parsing a single stats file."""
        lf = _parse_stats_files(str(stats_all_amplicons))
        df = lf.collect()
        assert len(df) == 3

    def test_extract_sample_name(self, stats_all_amplicons: Path) -> None:
        """Test that sample name is correctly extracted from filename."""
        lf = _parse_stats_files(str(stats_all_amplicons))
        df = lf.collect()
        assert df["sample_name"].unique().to_list() == ["sampleA"]

    def test_extract_amplicon_name(self, stats_all_amplicons: Path) -> None:
        """Test that amplicon name is correctly extracted from filename."""
        lf = _parse_stats_files(str(stats_all_amplicons))
        df = lf.collect()
        amplicons = sorted(df["amplicon_name"].to_list())
        assert amplicons == ["amp1", "amp2", "amp3"]

    def test_extract_reads(self, stats_all_amplicons: Path) -> None:
        """Test that read counts are correctly extracted."""
        lf = _parse_stats_files(str(stats_all_amplicons))
        df = lf.collect().sort("amplicon_name")
        assert df["reads"].to_list() == [100, 150, 200]


class TestCreateSummary:
    """Test the summary creation with dropout handling."""

    def test_all_amplicons_present(
        self,
        primer_pairs_tsv: Path,
        stats_all_amplicons: Path,
        tmp_path: Path,
    ) -> None:
        """Test when all amplicons have reads."""
        positions_lf = _load_positions_from_tsv(str(primer_pairs_tsv))
        stats_lf = _parse_stats_files(str(stats_all_amplicons))
        output_path = tmp_path / "summary.tsv"

        _create_summary(stats_lf, positions_lf, str(output_path))

        result = pl.read_csv(output_path, separator="\t")
        assert len(result) == 3  # 1 sample × 3 amplicons
        assert result["reads"].to_list() == [100, 150, 200]

    def test_dropout_filled_with_zero(
        self,
        primer_pairs_tsv: Path,
        stats_with_dropout: Path,
        tmp_path: Path,
    ) -> None:
        """Test that dropout amplicons are filled with reads=0."""
        positions_lf = _load_positions_from_tsv(str(primer_pairs_tsv))
        stats_lf = _parse_stats_files(str(stats_with_dropout))
        output_path = tmp_path / "summary.tsv"

        _create_summary(stats_lf, positions_lf, str(output_path))

        result = pl.read_csv(output_path, separator="\t")
        # Should have 3 rows (all amplicons), even though only 2 had reads
        assert len(result) == 3

        # Check amp3 has reads=0
        amp3_row = result.filter(pl.col("amplicon_name") == "amp3")
        assert len(amp3_row) == 1
        assert amp3_row["reads"][0] == 0

    def test_multiple_samples_with_different_dropouts(
        self,
        primer_pairs_tsv: Path,
        stats_multiple_samples: tuple[Path, Path],
        tmp_path: Path,
    ) -> None:
        """Test with multiple samples having different dropout patterns."""
        stats_a, stats_b = stats_multiple_samples
        positions_lf = _load_positions_from_tsv(str(primer_pairs_tsv))

        # Parse both stats files using glob pattern
        stats_lf = _parse_stats_files(str(tmp_path / "stats_*.tsv"))
        output_path = tmp_path / "summary.tsv"

        _create_summary(stats_lf, positions_lf, str(output_path))

        result = pl.read_csv(output_path, separator="\t")

        # Should have 6 rows: 2 samples × 3 amplicons
        assert len(result) == 6

        # Sample A should have all reads
        sample_a = result.filter(pl.col("sample_name") == "sampleA").sort(
            "amplicon_name",
        )
        assert sample_a["reads"].to_list() == [100, 150, 200]

        # Sample B should have amp3 = 0
        sample_b = result.filter(pl.col("sample_name") == "sampleB").sort(
            "amplicon_name",
        )
        assert sample_b["reads"].to_list() == [50, 75, 0]

    def test_positions_preserved(
        self,
        primer_pairs_tsv: Path,
        stats_with_dropout: Path,
        tmp_path: Path,
    ) -> None:
        """Test that positions are preserved for all amplicons including dropouts."""
        positions_lf = _load_positions_from_tsv(str(primer_pairs_tsv))
        stats_lf = _parse_stats_files(str(stats_with_dropout))
        output_path = tmp_path / "summary.tsv"

        _create_summary(stats_lf, positions_lf, str(output_path))

        result = pl.read_csv(output_path, separator="\t")

        # Check that amp3 (the dropout) still has its positions
        amp3_row = result.filter(pl.col("amplicon_name") == "amp3")
        assert str(amp3_row["start_pos"][0]) == "500"
        assert str(amp3_row["end_pos"][0]) == "600"

    def test_output_sorted(
        self,
        primer_pairs_tsv: Path,
        stats_multiple_samples: tuple[Path, Path],
        tmp_path: Path,
    ) -> None:
        """Test that output is sorted by sample_name, then amplicon_name."""
        positions_lf = _load_positions_from_tsv(str(primer_pairs_tsv))
        stats_lf = _parse_stats_files(str(tmp_path / "stats_*.tsv"))
        output_path = tmp_path / "summary.tsv"

        _create_summary(stats_lf, positions_lf, str(output_path))

        result = pl.read_csv(output_path, separator="\t")

        # Check sorting
        expected_order = [
            ("sampleA", "amp1"),
            ("sampleA", "amp2"),
            ("sampleA", "amp3"),
            ("sampleB", "amp1"),
            ("sampleB", "amp2"),
            ("sampleB", "amp3"),
        ]
        actual_order = list(
            zip(result["sample_name"].to_list(), result["amplicon_name"].to_list()),
        )
        assert actual_order == expected_order


class TestEdgeCases:
    """Test edge cases."""

    def test_single_sample_all_dropouts_except_one(
        self,
        primer_pairs_tsv: Path,
        tmp_path: Path,
    ) -> None:
        """Test when a sample has only one amplicon with reads."""
        content = "\n".join(
            [
                "file\tformat\ttype\tnum_seqs\tsum_len\tmin_len\tavg_len\tmax_len",
                "lonely.amp2.fasta.gz\tFASTA\tDNA\t10\t1000\t80\t100\t120",
            ],
        )
        stats_file = tmp_path / "stats_lonely.tsv"
        stats_file.write_text(content)

        positions_lf = _load_positions_from_tsv(str(primer_pairs_tsv))
        stats_lf = _parse_stats_files(str(stats_file))
        output_path = tmp_path / "summary.tsv"

        _create_summary(stats_lf, positions_lf, str(output_path))

        result = pl.read_csv(output_path, separator="\t")

        # Should still have 3 rows
        assert len(result) == 3

        # amp1 and amp3 should be 0, amp2 should be 10
        result_sorted = result.sort("amplicon_name")
        assert result_sorted["reads"].to_list() == [0, 10, 0]
