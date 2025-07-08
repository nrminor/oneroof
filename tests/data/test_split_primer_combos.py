#!/usr/bin/env python3
"""
Comprehensive test suite for split_primer_combos.py

Tests the splitting of BED files containing multiple primer combinations
into individual BED files per combination.
"""

import pytest
import polars as pl
from pathlib import Path
import sys
import os

# Add the bin directory to the path so we can import the module
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "bin"))

from split_primer_combos import parse_command_line_args, main


class TestSplitPrimerCombos:
    """Test suite for primer combination splitting functionality"""

    @pytest.fixture
    def basic_bed_file(self, tmp_path):
        """Create a basic BED file with multiple primer combinations"""
        bed_content = """chr1\t100\t120\tamplicon1_LEFT_splice1\t60\t+
chr1\t500\t520\tamplicon1_RIGHT_splice1\t60\t-
chr1\t100\t120\tamplicon1_LEFT_splice2\t60\t+
chr1\t510\t530\tamplicon1_RIGHT_splice2\t60\t-
chr1\t1000\t1020\tamplicon2_LEFT_splice1\t60\t+
chr1\t1400\t1420\tamplicon2_RIGHT_splice1\t60\t-"""

        bed_file = tmp_path / "combined.bed"
        bed_file.write_text(bed_content)
        return bed_file

    @pytest.fixture
    def complex_bed_file(self, tmp_path):
        """Create a complex BED file with multiple amplicons and splice variants"""
        bed_content = """chr1\t100\t120\tamplicon1_LEFT_splice1\t60\t+
chr1\t500\t520\tamplicon1_RIGHT_splice1\t60\t-
chr1\t100\t120\tamplicon1_LEFT_splice2\t60\t+
chr1\t510\t530\tamplicon1_RIGHT_splice2\t60\t-
chr1\t110\t130\tamplicon1_LEFT_splice3\t60\t+
chr1\t520\t540\tamplicon1_RIGHT_splice3\t60\t-
chr2\t2000\t2020\tamplicon2_LEFT_splice1\t60\t+
chr2\t2400\t2420\tamplicon2_RIGHT_splice1\t60\t-
chr3\t3000\t3020\tamplicon3_LEFT_splice1\t60\t+
chr3\t3400\t3420\tamplicon3_RIGHT_splice1\t60\t-"""

        bed_file = tmp_path / "complex.bed"
        bed_file.write_text(bed_content)
        return bed_file

    @pytest.fixture
    def custom_suffix_bed_file(self, tmp_path):
        """Create a BED file with custom primer suffixes"""
        bed_content = """chr1\t100\t120\tamplicon1_FWD_splice1\t60\t+
chr1\t500\t520\tamplicon1_REV_splice1\t60\t-
chr1\t100\t120\tamplicon1_FWD_splice2\t60\t+
chr1\t510\t530\tamplicon1_REV_splice2\t60\t-"""

        bed_file = tmp_path / "custom_suffix.bed"
        bed_file.write_text(bed_content)
        return bed_file

    @pytest.fixture
    def invalid_bed_file(self, tmp_path):
        """Create an invalid BED file with unpaired primers"""
        bed_content = """chr1\t100\t120\tamplicon1_LEFT_splice1\t60\t+
chr1\t500\t520\tamplicon1_RIGHT_splice1\t60\t-
chr1\t1000\t1020\tamplicon2_LEFT_splice1\t60\t+"""  # Missing RIGHT primer

        bed_file = tmp_path / "invalid.bed"
        bed_file.write_text(bed_content)
        return bed_file

    def test_parse_command_line_args_basic(self):
        """Test basic command line argument parsing"""
        test_args = ["split_primer_combos.py", "-i", "test.bed"]

        with pytest.MonkeyPatch.context() as m:
            m.setattr(sys, "argv", test_args)
            bed_path, fwd_suffix, rev_suffix = parse_command_line_args()

            assert bed_path == Path("test.bed")
            assert fwd_suffix == "_LEFT"  # Default
            assert rev_suffix == "_RIGHT"  # Default

    def test_parse_command_line_args_custom_suffixes(self):
        """Test command line argument parsing with custom suffixes"""
        test_args = [
            "split_primer_combos.py",
            "-i",
            "test.bed",
            "-f",
            "_FWD",
            "-r",
            "_REV",
        ]

        with pytest.MonkeyPatch.context() as m:
            m.setattr(sys, "argv", test_args)
            bed_path, fwd_suffix, rev_suffix = parse_command_line_args()

            assert fwd_suffix == "_FWD"
            assert rev_suffix == "_REV"

    def test_split_basic_bed_file(self, basic_bed_file, tmp_path):
        """Test splitting a basic BED file into individual combinations"""
        # Change to temp directory to capture output files
        original_cwd = Path.cwd()
        os.chdir(tmp_path)

        try:
            test_args = ["split_primer_combos.py", "-i", str(basic_bed_file)]

            with pytest.MonkeyPatch.context() as m:
                m.setattr(sys, "argv", test_args)
                main()

            # Check that the correct files were created
            assert Path("amplicon1_splice1.bed").exists()
            assert Path("amplicon1_splice2.bed").exists()
            assert Path("amplicon2_splice1.bed").exists()

            # Verify content of one output file
            df = pl.read_csv("amplicon1_splice1.bed", separator="\t", has_header=False)
            assert df.shape[0] == 2  # Should have 2 primers (LEFT and RIGHT)

        finally:
            os.chdir(original_cwd)

    def test_split_complex_bed_file(self, complex_bed_file, tmp_path):
        """Test splitting a complex BED file with multiple amplicons"""
        original_cwd = Path.cwd()
        os.chdir(tmp_path)

        try:
            test_args = ["split_primer_combos.py", "-i", str(complex_bed_file)]

            with pytest.MonkeyPatch.context() as m:
                m.setattr(sys, "argv", test_args)
                main()

            # Check all expected files
            assert Path("amplicon1_splice1.bed").exists()
            assert Path("amplicon1_splice2.bed").exists()
            assert Path("amplicon1_splice3.bed").exists()
            assert Path("amplicon2_splice1.bed").exists()
            assert Path("amplicon3_splice1.bed").exists()

            # Verify each file has exactly 2 primers
            for bed_file in tmp_path.glob("amplicon*.bed"):
                df = pl.read_csv(bed_file, separator="\t", has_header=False)
                assert df.shape[0] == 2

        finally:
            os.chdir(original_cwd)

    def test_split_custom_suffix_bed(self, custom_suffix_bed_file, tmp_path):
        """Test splitting with custom primer suffixes"""
        original_cwd = Path.cwd()
        os.chdir(tmp_path)

        try:
            test_args = [
                "split_primer_combos.py",
                "-i",
                str(custom_suffix_bed_file),
                "-f",
                "_FWD",
                "-r",
                "_REV",
            ]

            with pytest.MonkeyPatch.context() as m:
                m.setattr(sys, "argv", test_args)
                main()

            assert Path("amplicon1_splice1.bed").exists()
            assert Path("amplicon1_splice2.bed").exists()

        finally:
            os.chdir(original_cwd)

    def test_invalid_bed_file_assertion(self, invalid_bed_file, tmp_path):
        """Test that invalid BED file with unpaired primers raises assertion"""
        original_cwd = Path.cwd()
        os.chdir(tmp_path)

        try:
            test_args = ["split_primer_combos.py", "-i", str(invalid_bed_file)]

            with pytest.MonkeyPatch.context() as m:
                m.setattr(sys, "argv", test_args)

                with pytest.raises(AssertionError) as excinfo:
                    main()

                assert "Problematic splicing occurred" in str(excinfo.value)

        finally:
            os.chdir(original_cwd)

    def test_empty_bed_file(self, tmp_path):
        """Test behavior with empty BED file"""
        bed_file = tmp_path / "empty.bed"
        bed_file.write_text("")

        original_cwd = Path.cwd()
        os.chdir(tmp_path)

        try:
            test_args = ["split_primer_combos.py", "-i", str(bed_file)]

            with pytest.MonkeyPatch.context() as m:
                m.setattr(sys, "argv", test_args)
                main()

            # No files should be created
            bed_files = list(tmp_path.glob("*.bed"))
            assert len(bed_files) == 1  # Only the input file

        finally:
            os.chdir(original_cwd)

    def test_bed_file_with_duplicate_combinations(self, tmp_path):
        """Test BED file with duplicate primer combinations"""
        bed_content = """chr1\t100\t120\tamplicon1_LEFT_splice1\t60\t+
chr1\t500\t520\tamplicon1_RIGHT_splice1\t60\t-
chr1\t100\t120\tamplicon1_LEFT_splice1\t60\t+
chr1\t500\t520\tamplicon1_RIGHT_splice1\t60\t-"""

        bed_file = tmp_path / "duplicates.bed"
        bed_file.write_text(bed_content)

        original_cwd = Path.cwd()
        os.chdir(tmp_path)

        try:
            test_args = ["split_primer_combos.py", "-i", str(bed_file)]

            with pytest.MonkeyPatch.context() as m:
                m.setattr(sys, "argv", test_args)

                with pytest.raises(AssertionError):
                    main()  # Should fail due to 4 primers with same name

        finally:
            os.chdir(original_cwd)

    def test_bed_file_with_mixed_chromosomes(self, tmp_path):
        """Test BED file where same amplicon spans different chromosomes"""
        bed_content = """chr1\t100\t120\tamplicon1_LEFT_splice1\t60\t+
chr2\t500\t520\tamplicon1_RIGHT_splice1\t60\t-"""

        bed_file = tmp_path / "mixed_chr.bed"
        bed_file.write_text(bed_content)

        original_cwd = Path.cwd()
        os.chdir(tmp_path)

        try:
            test_args = ["split_primer_combos.py", "-i", str(bed_file)]

            with pytest.MonkeyPatch.context() as m:
                m.setattr(sys, "argv", test_args)
                main()

            # Should still create the file despite different chromosomes
            assert Path("amplicon1_splice1.bed").exists()

            # Verify content preserved original chromosomes
            df = pl.read_csv("amplicon1_splice1.bed", separator="\t", has_header=False)
            chrs = df.column(0).to_list()
            assert "chr1" in chrs
            assert "chr2" in chrs

        finally:
            os.chdir(original_cwd)

    def test_preserve_all_bed_columns(self, tmp_path):
        """Test that all BED columns are preserved in output"""
        bed_content = """chr1\t100\t120\tamplicon1_LEFT_splice1\t60\t+
chr1\t500\t520\tamplicon1_RIGHT_splice1\t60\t-"""

        bed_file = tmp_path / "preserve.bed"
        bed_file.write_text(bed_content)

        original_cwd = Path.cwd()
        os.chdir(tmp_path)

        try:
            test_args = ["split_primer_combos.py", "-i", str(bed_file)]

            with pytest.MonkeyPatch.context() as m:
                m.setattr(sys, "argv", test_args)
                main()

            df = pl.read_csv("amplicon1_splice1.bed", separator="\t", has_header=False)

            # Check all 6 columns are present
            assert df.shape[1] == 6

            # Verify specific values
            assert df.row(0) == ("chr1", 100, 120, "amplicon1_splice1", 60, "+")
            assert df.row(1) == ("chr1", 500, 520, "amplicon1_splice1", 60, "-")

        finally:
            os.chdir(original_cwd)

    def test_special_characters_in_amplicon_names(self, tmp_path):
        """Test handling of special characters in amplicon names"""
        bed_content = """chr1\t100\t120\tamplicon1.v2_LEFT_splice1\t60\t+
chr1\t500\t520\tamplicon1.v2_RIGHT_splice1\t60\t-"""

        bed_file = tmp_path / "special.bed"
        bed_file.write_text(bed_content)

        original_cwd = Path.cwd()
        os.chdir(tmp_path)

        try:
            test_args = ["split_primer_combos.py", "-i", str(bed_file)]

            with pytest.MonkeyPatch.context() as m:
                m.setattr(sys, "argv", test_args)
                main()

            # File should be created with special characters preserved
            assert Path("amplicon1.v2_splice1.bed").exists()

        finally:
            os.chdir(original_cwd)
