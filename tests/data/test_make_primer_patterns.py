#!/usr/bin/env python3
"""
Comprehensive test suite for make_primer_patterns.py

Tests the generation of regular expression patterns from primer FASTA files,
focusing on biological correctness and edge cases.
"""

import pytest
from pathlib import Path
import sys
import os

# Add the bin directory to the path so we can import the module
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "bin"))

from make_primer_patterns import generate_regex_patterns, parse_command_line_args


class TestMakePrimerPatterns:
    """Test suite for primer pattern generation"""

    @pytest.fixture
    def valid_primer_fasta(self, tmp_path):
        """Create a valid primer FASTA file with bedtools getfasta format"""
        fasta_file = tmp_path / "primers.fasta"
        content = """>MN908947.3:0-16
ATCGATCGATCGATCG
>MN908947.3:400-416
GCTAGCTAGCTAGCTA
"""
        fasta_file.write_text(content)
        return str(fasta_file)

    @pytest.fixture
    def invalid_primer_fasta_single(self, tmp_path):
        """Create an invalid FASTA with only one sequence"""
        fasta_file = tmp_path / "single_primer.fasta"
        content = """>MN908947.3:0-16
ATCGATCGATCGATCG
"""
        fasta_file.write_text(content)
        return str(fasta_file)

    @pytest.fixture
    def invalid_primer_fasta_three(self, tmp_path):
        """Create an invalid FASTA with three sequences"""
        fasta_file = tmp_path / "three_primers.fasta"
        content = """>MN908947.3:0-16
ATCGATCGATCGATCG
>MN908947.3:400-416
GCTAGCTAGCTAGCTA
>MN908947.3:800-816
TTTTTTTTTTTTTTTT
"""
        fasta_file.write_text(content)
        return str(fasta_file)

    @pytest.fixture
    def reversed_coords_fasta(self, tmp_path):
        """Create a FASTA where reverse primer appears before forward primer"""
        fasta_file = tmp_path / "reversed.fasta"
        content = """>MN908947.3:400-416
GCTAGCTAGCTAGCTA
>MN908947.3:0-16
ATCGATCGATCGATCG
"""
        fasta_file.write_text(content)
        return str(fasta_file)

    def test_generate_regex_patterns_valid(self, valid_primer_fasta, tmp_path):
        """Test generation of regex patterns with valid input"""
        os.chdir(tmp_path)
        label = "test_patterns"

        generate_regex_patterns(valid_primer_fasta, label, r"^(.*?)", r"^(.*?)")

        output_file = Path(f"{label}.txt")
        assert output_file.exists()

        lines = output_file.read_text().strip().split("\n")
        assert len(lines) == 2
        assert lines[0] == r"^(.*?)ATCGATCGATCGATCG"
        assert lines[1] == r"GCTAGCTAGCTAGCTA^(.*?)"

    def test_custom_patterns(self, valid_primer_fasta, tmp_path):
        """Test with custom forward and reverse patterns"""
        os.chdir(tmp_path)
        label = "custom_patterns"

        generate_regex_patterns(
            valid_primer_fasta,
            label,
            r".*?",  # Different forward pattern
            r"(.*)",  # Different reverse pattern
        )

        output_file = Path(f"{label}.txt")
        lines = output_file.read_text().strip().split("\n")
        assert lines[0] == r".*?ATCGATCGATCGATCG"
        assert lines[1] == r"GCTAGCTAGCTAGCTA(.*)"

    def test_invalid_single_sequence(self, invalid_primer_fasta_single):
        """Test that single sequence FASTA raises assertion error"""
        with pytest.raises(AssertionError) as excinfo:
            generate_regex_patterns(
                invalid_primer_fasta_single, "test", r"^(.*?)", r"^(.*?)"
            )
        assert "does not contain exactly two sequences" in str(excinfo.value)

    def test_invalid_three_sequences(self, invalid_primer_fasta_three):
        """Test that three sequences FASTA raises assertion error"""
        with pytest.raises(AssertionError) as excinfo:
            generate_regex_patterns(
                invalid_primer_fasta_three, "test", r"^(.*?)", r"^(.*?)"
            )
        assert "does not contain exactly two sequences" in str(excinfo.value)

    def test_reversed_coordinates_warning(self, reversed_coords_fasta, tmp_path):
        """Test that reversed coordinates trigger a warning"""
        os.chdir(tmp_path)
        with pytest.warns(UserWarning, match="bedtools getfasta"):
            generate_regex_patterns(reversed_coords_fasta, "test", r"^(.*?)", r"^(.*?)")

    def test_empty_fasta(self, tmp_path):
        """Test behavior with empty FASTA file"""
        fasta_file = tmp_path / "empty.fasta"
        fasta_file.write_text("")

        with pytest.raises(AssertionError):
            generate_regex_patterns(str(fasta_file), "test", r"^(.*?)", r"^(.*?)")

    def test_fasta_with_empty_sequences(self, tmp_path):
        """Test FASTA with headers but no sequences"""
        fasta_file = tmp_path / "empty_seqs.fasta"
        content = """>MN908947.3:0-16
>MN908947.3:400-416
"""
        fasta_file.write_text(content)

        with pytest.raises(AssertionError):
            generate_regex_patterns(str(fasta_file), "test", r"^(.*?)", r"^(.*?)")

    def test_special_characters_in_sequences(self, tmp_path):
        """Test handling of sequences with special regex characters"""
        fasta_file = tmp_path / "special.fasta"
        content = """>chr1:0-16
ATCG[AT]CG*TCG+ATC
>chr1:400-416
GCT(AGC)TAG{CTA}GC
"""
        fasta_file.write_text(content)
        os.chdir(tmp_path)

        generate_regex_patterns(str(fasta_file), "special", r"^(.*?)", r"^(.*?)")

        output_file = Path("special.txt")
        lines = output_file.read_text().strip().split("\n")
        # Special characters should be preserved as-is
        assert lines[0] == r"^(.*?)ATCG[AT]CG*TCG+ATC"
        assert lines[1] == r"GCT(AGC)TAG{CTA}GC^(.*?)"

    def test_multiline_sequences(self, tmp_path):
        """Test FASTA with sequences spanning multiple lines"""
        fasta_file = tmp_path / "multiline.fasta"
        content = """>chr1:0-48
ATCGATCGATCGATCG
ATCGATCGATCGATCG
ATCGATCGATCGATCG
>chr1:400-448
GCTAGCTAGCTAGCTA
GCTAGCTAGCTAGCTA
GCTAGCTAGCTAGCTA
"""
        fasta_file.write_text(content)
        os.chdir(tmp_path)

        generate_regex_patterns(str(fasta_file), "multiline", r"^(.*?)", r"^(.*?)")

        output_file = Path("multiline.txt")
        lines = output_file.read_text().strip().split("\n")
        # Sequences should be concatenated
        assert lines[0] == r"^(.*?)ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
        assert lines[1] == r"GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA^(.*?)"

    def test_parse_command_line_args(self):
        """Test command line argument parsing"""
        test_args = [
            "make_primer_patterns.py",
            "-i",
            "test.fasta",
            "-o",
            "output_prefix",
            "-f",
            r".*?",
            "-r",
            r"(.*)",
        ]

        with pytest.MonkeyPatch.context() as m:
            m.setattr(sys, "argv", test_args)
            args = parse_command_line_args()

            assert args.input_fasta == Path("test.fasta")
            assert args.output_prefix == "output_prefix"
            assert args.forward_pattern == r".*?"
            assert args.reverse_pattern == r"(.*)"

    def test_default_arguments(self):
        """Test default command line arguments"""
        test_args = ["make_primer_patterns.py", "-i", "test.fasta"]

        with pytest.MonkeyPatch.context() as m:
            m.setattr(sys, "argv", test_args)
            args = parse_command_line_args()

            assert args.output_prefix == "primer_patterns"
            assert args.forward_pattern == r"^(.*?)"
            assert args.reverse_pattern == r"^(.*?)"

    def test_nonstandard_header_format(self, tmp_path):
        """Test FASTA with non-bedtools header format"""
        fasta_file = tmp_path / "nonstandard.fasta"
        content = """>primer1_forward
ATCGATCGATCGATCG
>primer2_reverse
GCTAGCTAGCTAGCTA
"""
        fasta_file.write_text(content)
        os.chdir(tmp_path)

        # Should work but might trigger coordinate check
        generate_regex_patterns(str(fasta_file), "nonstandard", r"^(.*?)", r"^(.*?)")

        output_file = Path("nonstandard.txt")
        assert output_file.exists()
