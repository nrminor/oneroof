#!/usr/bin/env python3
"""
Comprehensive test suite for validate_primer_bed.py

Tests the validation and normalization of BED files containing primer information,
focusing on biological correctness and edge cases.
"""

import pytest
from pathlib import Path
import sys
import os

# Add the bin directory to the path so we can import the module
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "bin"))

from validate_primer_bed import (
    parse_command_line_args,
    check_for_suffixes,
    check_for_pairs,
    orient_primer_coords,
    normalize_bed_lines,
)


class TestValidatePrimerBed:
    """Test suite for primer BED validation functionality"""

    @pytest.fixture
    def valid_bed_file(self, tmp_path):
        """Create a valid BED file with properly formatted primers"""
        bed_content = """chr1\t100\t120\tamplicon1_LEFT\t60\t+
chr1\t500\t520\tamplicon1_RIGHT\t60\t-
chr1\t1000\t1020\tamplicon2_LEFT\t60\t+
chr1\t1400\t1420\tamplicon2_RIGHT\t60\t-"""

        bed_file = tmp_path / "valid.bed"
        bed_file.write_text(bed_content)
        return bed_file

    @pytest.fixture
    def reversed_coords_bed_file(self, tmp_path):
        """Create a BED file with reversed coordinates"""
        bed_content = """chr1\t120\t100\tamplicon1_LEFT\t60\t+
chr1\t520\t500\tamplicon1_RIGHT\t60\t-
chr1\t1020\t1000\tamplicon2_LEFT\t60\t-
chr1\t1420\t1400\tamplicon2_RIGHT\t60\t+"""

        bed_file = tmp_path / "reversed.bed"
        bed_file.write_text(bed_content)
        return bed_file

    @pytest.fixture
    def missing_suffix_bed_file(self, tmp_path):
        """Create a BED file with missing primer suffixes"""
        bed_content = """chr1\t100\t120\tamplicon1\t60\t+
chr1\t500\t520\tamplicon1_RIGHT\t60\t-"""

        bed_file = tmp_path / "missing_suffix.bed"
        bed_file.write_text(bed_content)
        return bed_file

    @pytest.fixture
    def unpaired_bed_file(self, tmp_path):
        """Create a BED file with unpaired primers"""
        bed_content = """chr1\t100\t120\tamplicon1_LEFT\t60\t+
chr1\t1000\t1020\tamplicon2_LEFT\t60\t+
chr1\t1400\t1420\tamplicon2_RIGHT\t60\t-"""

        bed_file = tmp_path / "unpaired.bed"
        bed_file.write_text(bed_content)
        return bed_file

    @pytest.fixture
    def incomplete_bed_file(self, tmp_path):
        """Create a BED file with fewer than 6 columns"""
        bed_content = """chr1\t100\t120\tamplicon1_LEFT\t60
chr1\t500\t520\tamplicon1_RIGHT"""

        bed_file = tmp_path / "incomplete.bed"
        bed_file.write_text(bed_content)
        return bed_file

    def test_parse_command_line_args_defaults(self):
        """Test default command line arguments"""
        test_args = ["validate_primer_bed.py", "-i", "test.bed"]

        with pytest.MonkeyPatch.context() as m:
            m.setattr(sys, "argv", test_args)
            args = parse_command_line_args()

            assert args.input_bed == Path("test.bed")
            assert args.output_prefix == "validated"
            assert args.fwd_suffix == "_LEFT"
            assert args.rev_suffix == "_RIGHT"

    def test_parse_command_line_args_custom(self):
        """Test custom command line arguments"""
        test_args = [
            "validate_primer_bed.py",
            "-i",
            "input.bed",
            "-o",
            "output",
            "-f",
            "_FWD",
            "-r",
            "_REV",
        ]

        with pytest.MonkeyPatch.context() as m:
            m.setattr(sys, "argv", test_args)
            args = parse_command_line_args()

            assert args.output_prefix == "output"
            assert args.fwd_suffix == "_FWD"
            assert args.rev_suffix == "_REV"

    def test_check_for_suffixes_valid(self):
        """Test suffix checking with valid primer names"""
        row = ["chr1", "100", "120", "amplicon1_LEFT", "60", "+"]
        result = check_for_suffixes(row, "_LEFT", "_RIGHT")
        assert result is None

        row = ["chr1", "500", "520", "amplicon1_RIGHT", "60", "-"]
        result = check_for_suffixes(row, "_LEFT", "_RIGHT")
        assert result is None

    def test_check_for_suffixes_invalid(self):
        """Test suffix checking with invalid primer names"""
        row = ["chr1", "100", "120", "amplicon1", "60", "+"]
        result = check_for_suffixes(row, "_LEFT", "_RIGHT")
        assert result == "amplicon1"

    def test_check_for_suffixes_insufficient_columns(self):
        """Test suffix checking with insufficient columns"""
        row = ["chr1", "100", "120", "amplicon1_LEFT", "60"]  # Missing strand

        with pytest.raises(AssertionError) as excinfo:
            check_for_suffixes(row)
        assert "fewer than the required 6 columns" in str(excinfo.value)

    def test_check_for_pairs_valid(self):
        """Test pair checking with valid primer pairs"""
        rows = [
            ["chr1", "100", "120", "amplicon1_LEFT", "60", "+"],
            ["chr1", "500", "520", "amplicon1_RIGHT", "60", "-"],
            ["chr1", "1000", "1020", "amplicon2_LEFT", "60", "+"],
            ["chr1", "1400", "1420", "amplicon2_RIGHT", "60", "-"],
        ]

        singletons = check_for_pairs(rows)
        assert len(singletons) == 0

    def test_check_for_pairs_singletons(self):
        """Test pair checking with unpaired primers"""
        rows = [
            ["chr1", "100", "120", "amplicon1_LEFT", "60", "+"],
            ["chr1", "1000", "1020", "amplicon2_LEFT", "60", "+"],
            ["chr1", "1400", "1420", "amplicon2_RIGHT", "60", "-"],
        ]

        singletons = check_for_pairs(rows)
        assert len(singletons) == 1
        assert "amplicon1" in singletons

    def test_check_for_pairs_with_spike_ins(self):
        """Test pair checking with spike-in primers"""
        rows = [
            ["chr1", "100", "120", "amplicon1_LEFT-1", "60", "+"],
            ["chr1", "110", "130", "amplicon1_LEFT-2", "60", "+"],
            ["chr1", "500", "520", "amplicon1_RIGHT-1", "60", "-"],
            ["chr1", "510", "530", "amplicon1_RIGHT-2", "60", "-"],
        ]

        singletons = check_for_pairs(rows)
        assert len(singletons) == 0  # All should be paired

    def test_orient_primer_coords_no_change(self):
        """Test coordinate orientation when already correct"""
        row = ["chr1", "100", "120", "amplicon1_LEFT", "60", "+"]
        result = orient_primer_coords(row, 0)
        assert result == row

    def test_orient_primer_coords_reversed_positive(self):
        """Test coordinate orientation with reversed positive strand"""
        row = ["chr1", "120", "100", "amplicon1_LEFT", "60", "+"]
        result = orient_primer_coords(row, 0)

        assert result[1] == "100"
        assert result[2] == "120"
        assert result[5] == "-"

    def test_orient_primer_coords_reversed_negative(self):
        """Test coordinate orientation with reversed negative strand"""
        row = ["chr1", "520", "500", "amplicon1_RIGHT", "60", "-"]
        result = orient_primer_coords(row, 0)

        assert result[1] == "500"
        assert result[2] == "520"
        assert result[5] == "+"

    def test_orient_primer_coords_invalid_strand(self):
        """Test coordinate orientation with invalid strand character"""
        row = ["chr1", "120", "100", "amplicon1_LEFT", "60", "?"]

        with pytest.raises(ValueError) as excinfo:
            orient_primer_coords(row, 0)
        assert "unsupported value encountered in the strand column" in str(
            excinfo.value
        )

    def test_normalize_bed_lines_valid(self, valid_bed_file, tmp_path):
        """Test normalization of valid BED file"""
        output_prefix = "test_output"
        normalize_bed_lines(valid_bed_file, output_prefix)

        output_file = Path(f"{output_prefix}.bed")
        assert output_file.exists()

        # Check content
        lines = output_file.read_text().strip().split("\n")
        assert len(lines) == 4

        # Verify format preserved
        for line in lines:
            columns = line.split("\t")
            assert len(columns) == 6

    def test_normalize_bed_lines_reversed_coords(
        self, reversed_coords_bed_file, tmp_path
    ):
        """Test normalization fixes reversed coordinates"""
        os.chdir(tmp_path)
        output_prefix = "fixed"
        normalize_bed_lines(reversed_coords_bed_file, output_prefix)

        output_file = Path(f"{output_prefix}.bed")
        lines = output_file.read_text().strip().split("\n")

        # Check all coordinates are properly oriented
        for line in lines:
            columns = line.split("\t")
            start = int(columns[1])
            stop = int(columns[2])
            assert start < stop

    def test_normalize_bed_lines_missing_suffix(
        self, missing_suffix_bed_file, tmp_path
    ):
        """Test normalization fails with missing suffixes"""
        with pytest.raises(AssertionError) as excinfo:
            normalize_bed_lines(missing_suffix_bed_file)
        assert "Invalid primer label(s)" in str(excinfo.value)

    def test_normalize_bed_lines_unpaired(self, unpaired_bed_file, tmp_path):
        """Test normalization fails with unpaired primers"""
        with pytest.raises(AssertionError) as excinfo:
            normalize_bed_lines(unpaired_bed_file)
        assert "do not appear in more than one primer" in str(excinfo.value)

    def test_normalize_bed_lines_incomplete(self, incomplete_bed_file, tmp_path):
        """Test normalization fails with incomplete BED lines"""
        with pytest.raises(AssertionError) as excinfo:
            normalize_bed_lines(incomplete_bed_file)
        assert "fewer than the required 6 columns" in str(excinfo.value)

    def test_normalize_bed_lines_empty_file(self, tmp_path):
        """Test normalization with empty BED file"""
        empty_file = tmp_path / "empty.bed"
        empty_file.write_text("")

        os.chdir(tmp_path)
        normalize_bed_lines(empty_file, "empty_output")

        output_file = Path("empty_output.bed")
        assert output_file.exists()
        assert output_file.read_text() == ""

    def test_normalize_bed_lines_custom_suffixes(self, tmp_path):
        """Test normalization with custom primer suffixes"""
        bed_content = """chr1\t100\t120\tamplicon1_FWD\t60\t+
chr1\t500\t520\tamplicon1_REV\t60\t-"""

        bed_file = tmp_path / "custom.bed"
        bed_file.write_text(bed_content)

        os.chdir(tmp_path)
        normalize_bed_lines(bed_file, "custom_output", "_FWD", "_REV")

        assert Path("custom_output.bed").exists()

    def test_normalize_bed_lines_extra_columns(self, tmp_path):
        """Test BED file with extra columns beyond the required 6"""
        bed_content = """chr1\t100\t120\tamplicon1_LEFT\t60\t+\textra1\textra2
chr1\t500\t520\tamplicon1_RIGHT\t60\t-\textra1\textra2"""

        bed_file = tmp_path / "extra_cols.bed"
        bed_file.write_text(bed_content)

        os.chdir(tmp_path)
        normalize_bed_lines(bed_file, "extra_output")

        output_file = Path("extra_output.bed")
        lines = output_file.read_text().strip().split("\n")

        # Extra columns should be preserved
        for line in lines:
            columns = line.split("\t")
            assert len(columns) == 8

    def test_normalize_bed_lines_whitespace_handling(self, tmp_path):
        """Test handling of various whitespace in BED file"""
        bed_content = """chr1\t100\t120\tamplicon1_LEFT\t60\t+
chr1\t500\t520\tamplicon1_RIGHT\t60\t-\t
  chr1\t1000\t1020\tamplicon2_LEFT\t60\t+
chr1\t1400\t1420\tamplicon2_RIGHT\t60\t-"""

        bed_file = tmp_path / "whitespace.bed"
        bed_file.write_text(bed_content)

        os.chdir(tmp_path)
        normalize_bed_lines(bed_file, "whitespace_output")

        output_file = Path("whitespace_output.bed")
        lines = output_file.read_text().strip().split("\n")
        assert len(lines) == 4

    def test_complex_amplicon_names(self, tmp_path):
        """Test handling of complex amplicon names with special characters"""
        bed_content = """chr1\t100\t120\tamplicon1.v2_LEFT\t60\t+
chr1\t500\t520\tamplicon1.v2_RIGHT\t60\t-
chr1\t1000\t1020\tamplicon-2_LEFT\t60\t+
chr1\t1400\t1420\tamplicon-2_RIGHT\t60\t-"""

        bed_file = tmp_path / "complex_names.bed"
        bed_file.write_text(bed_content)

        os.chdir(tmp_path)
        normalize_bed_lines(bed_file, "complex_output")

        assert Path("complex_output.bed").exists()

    def test_mixed_case_chromosomes(self, tmp_path):
        """Test handling of mixed case chromosome names"""
        bed_content = """Chr1\t100\t120\tamplicon1_LEFT\t60\t+
CHR1\t500\t520\tamplicon1_RIGHT\t60\t-
chrX\t1000\t1020\tamplicon2_LEFT\t60\t+
chrx\t1400\t1420\tamplicon2_RIGHT\t60\t-"""

        bed_file = tmp_path / "mixed_case.bed"
        bed_file.write_text(bed_content)

        os.chdir(tmp_path)
        normalize_bed_lines(bed_file, "mixed_case_output")

        output_file = Path("mixed_case_output.bed")
        lines = output_file.read_text().strip().split("\n")

        # Original case should be preserved
        assert "Chr1" in lines[0]
        assert "CHR1" in lines[1]
        assert "chrX" in lines[2]
        assert "chrx" in lines[3]
