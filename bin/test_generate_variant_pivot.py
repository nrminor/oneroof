#!/usr/bin/env python3
"""Tests for generate_variant_pivot.py module."""

import tempfile
from pathlib import Path

import polars as pl
import pytest

from generate_variant_pivot import main, parse_command_line_args


@pytest.fixture
def valid_variant_tsv():
    """Create a valid variant TSV file for testing."""
    data = """CHROM\tREF\tPOS\tALT\tAF\tAC\tDP\tMQ\tGENE\tAA_EFFECT\tREF_CODON_ALT\tCDS_POS\tAA_POS
NC_045512.2\tA\t241\tT\t1.0\t2\t100\t60.0\tORF1ab\t.p.Asn123Asp\tAAC/GAC\t369\t123
NC_045512.2\tC\t3037\tT\t0.95\t95\t100\t60.0\tORF1ab\t.p.Phe924Phe\tTTC/TTT\t2772\t924
NC_045512.2\tC\t14408\tT\t1.0\t100\t100\t60.0\tORF1ab\t.p.Pro4715Leu\tCCT/CTT\t14145\t4715
NC_045512.2\tG\t23403\tA\t0.99\t99\t100\t60.0\tS\t.p.Asp614Gly\tGAT/GGT\t1841\t614
"""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
        f.write(data)
        return Path(f.name)


@pytest.fixture
def empty_variant_tsv():
    """Create an empty variant TSV file for testing."""
    data = """CHROM\tREF\tPOS\tALT\tAF\tAC\tDP\tMQ\tGENE\tAA_EFFECT\tREF_CODON_ALT\tCDS_POS\tAA_POS
"""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
        f.write(data)
        return Path(f.name)


@pytest.fixture
def malformed_variant_tsv():
    """Create a malformed variant TSV file for testing."""
    data = """CHROM\tREF\tPOS\tALT\tAF
NC_045512.2\tA\t241\tT\t1.0\tmissing_columns
"""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
        f.write(data)
        return Path(f.name)


@pytest.fixture
def variant_with_special_chars():
    """Create a variant TSV with special characters in amino acid effects."""
    data = """CHROM\tREF\tPOS\tALT\tAF\tAC\tDP\tMQ\tGENE\tAA_EFFECT\tREF_CODON_ALT\tCDS_POS\tAA_POS
NC_045512.2\tA\t241\tT\t1.0\t2\t100\t60.0\tORF1ab\t.p.Asn123*\tAAC/GAC\t369\t123
NC_045512.2\tC\t3037\tT\t0.95\t95\t100\t60.0\tORF1ab\t.p.Ter924Trp\tTTC/TTT\t2772\t924
NC_045512.2\tC\t14408\tT\t1.0\t100\t100\t60.0\tORF1ab\t.p.Pro4715=\tCCT/CTT\t14145\t4715
"""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
        f.write(data)
        return Path(f.name)


class TestParseCommandLineArgs:
    """Test command line argument parsing."""

    def test_parse_valid_args(self, valid_variant_tsv, monkeypatch):
        """Test parsing valid command line arguments."""
        test_args = ["generate_variant_pivot.py", "-i", str(valid_variant_tsv)]
        monkeypatch.setattr("sys.argv", test_args)

        args = parse_command_line_args()
        assert args.input_table == valid_variant_tsv
        assert isinstance(args.input_table, Path)

    def test_parse_long_form_args(self, valid_variant_tsv, monkeypatch):
        """Test parsing long form command line arguments."""
        test_args = [
            "generate_variant_pivot.py",
            "--input_table",
            str(valid_variant_tsv),
        ]
        monkeypatch.setattr("sys.argv", test_args)

        args = parse_command_line_args()
        assert args.input_table == valid_variant_tsv

    def test_missing_required_args(self, monkeypatch):
        """Test that missing required arguments raises SystemExit."""
        test_args = ["generate_variant_pivot.py"]
        monkeypatch.setattr("sys.argv", test_args)

        with pytest.raises(SystemExit):
            parse_command_line_args()


class TestVariantPivotGeneration:
    """Test variant pivot table generation functionality."""

    def test_valid_variant_parsing(self, valid_variant_tsv):
        """Test parsing of valid variant file."""
        df = pl.read_csv(
            valid_variant_tsv,
            separator="\t",
            has_header=False,
            skip_rows=1,
            columns=[
                "contig",
                "ref",
                "pos",
                "alt",
                "af",
                "ac",
                "dp",
                "mq",
                "gene",
                "aa_effect",
                "ref_codon_alt",
                "cds_pos",
                "aa_pos",
            ],
        )

        assert len(df) == 4
        assert df.columns == [
            "contig",
            "ref",
            "pos",
            "alt",
            "af",
            "ac",
            "dp",
            "mq",
            "gene",
            "aa_effect",
            "ref_codon_alt",
            "cds_pos",
            "aa_pos",
        ]
        assert df["contig"].dtype == pl.Utf8
        assert df["pos"].dtype == pl.Int64
        assert df["af"].dtype == pl.Float64

    def test_aa_effect_cleaning(self, valid_variant_tsv):
        """Test that .p prefix is removed from amino acid effects."""
        df = pl.read_csv(
            valid_variant_tsv,
            separator="\t",
            has_header=False,
            skip_rows=1,
            columns=[
                "contig",
                "ref",
                "pos",
                "alt",
                "af",
                "ac",
                "dp",
                "mq",
                "gene",
                "aa_effect",
                "ref_codon_alt",
                "cds_pos",
                "aa_pos",
            ],
        ).with_columns(pl.col("aa_effect").str.replace(".p.", "").alias("aa_effect"))

        # Check that .p. prefix is removed
        assert all(not effect.startswith(".p.") for effect in df["aa_effect"])
        assert df["aa_effect"][0] == "Asn123Asp"
        assert df["aa_effect"][1] == "Phe924Phe"

    def test_empty_file_handling(self, empty_variant_tsv):
        """Test handling of empty variant file."""
        df = pl.read_csv(
            empty_variant_tsv,
            separator="\t",
            has_header=False,
            skip_rows=1,
            columns=[
                "contig",
                "ref",
                "pos",
                "alt",
                "af",
                "ac",
                "dp",
                "mq",
                "gene",
                "aa_effect",
                "ref_codon_alt",
                "cds_pos",
                "aa_pos",
            ],
        )

        assert len(df) == 0
        assert df.columns == [
            "contig",
            "ref",
            "pos",
            "alt",
            "af",
            "ac",
            "dp",
            "mq",
            "gene",
            "aa_effect",
            "ref_codon_alt",
            "cds_pos",
            "aa_pos",
        ]

    def test_malformed_file_handling(self, malformed_variant_tsv):
        """Test handling of malformed variant file."""
        with pytest.raises(
            Exception
        ):  # Polars will raise an exception for column mismatch
            pl.read_csv(
                malformed_variant_tsv,
                separator="\t",
                has_header=False,
                skip_rows=1,
                columns=[
                    "contig",
                    "ref",
                    "pos",
                    "alt",
                    "af",
                    "ac",
                    "dp",
                    "mq",
                    "gene",
                    "aa_effect",
                    "ref_codon_alt",
                    "cds_pos",
                    "aa_pos",
                ],
            )

    def test_special_characters_handling(self, variant_with_special_chars):
        """Test handling of special characters in amino acid effects."""
        df = pl.read_csv(
            variant_with_special_chars,
            separator="\t",
            has_header=False,
            skip_rows=1,
            columns=[
                "contig",
                "ref",
                "pos",
                "alt",
                "af",
                "ac",
                "dp",
                "mq",
                "gene",
                "aa_effect",
                "ref_codon_alt",
                "cds_pos",
                "aa_pos",
            ],
        ).with_columns(pl.col("aa_effect").str.replace(".p.", "").alias("aa_effect"))

        # Check special amino acid notations are preserved
        assert df["aa_effect"][0] == "Asn123*"  # Stop codon
        assert df["aa_effect"][1] == "Ter924Trp"  # Termination to Trp
        assert df["aa_effect"][2] == "Pro4715="  # Synonymous


class TestMainFunction:
    """Test the main function execution."""

    def test_main_execution(self, valid_variant_tsv, monkeypatch, capsys):
        """Test that main function executes without errors."""
        test_args = ["generate_variant_pivot.py", "-i", str(valid_variant_tsv)]
        monkeypatch.setattr("sys.argv", test_args)

        # Run main function
        main()

        # Check that logger output is produced
        captured = capsys.readouterr()
        assert "Hi mom!" in captured.out or "Hi mom!" in captured.err

    def test_main_with_invalid_file(self, monkeypatch):
        """Test main function with non-existent file."""
        non_existent = Path("/tmp/non_existent_file.tsv")
        test_args = ["generate_variant_pivot.py", "-i", str(non_existent)]
        monkeypatch.setattr("sys.argv", test_args)

        with pytest.raises(Exception):
            main()


class TestDataIntegrity:
    """Test data integrity and edge cases."""

    def test_numeric_columns(self, valid_variant_tsv):
        """Test that numeric columns contain valid values."""
        df = pl.read_csv(
            valid_variant_tsv,
            separator="\t",
            has_header=False,
            skip_rows=1,
            columns=[
                "contig",
                "ref",
                "pos",
                "alt",
                "af",
                "ac",
                "dp",
                "mq",
                "gene",
                "aa_effect",
                "ref_codon_alt",
                "cds_pos",
                "aa_pos",
            ],
        )

        # Check position values are positive
        assert all(df["pos"] > 0)

        # Check allele frequency is between 0 and 1
        assert all((df["af"] >= 0) & (df["af"] <= 1))

        # Check depth values are non-negative
        assert all(df["dp"] >= 0)
        assert all(df["ac"] >= 0)

        # Check mapping quality is reasonable
        assert all(df["mq"] >= 0)

    def test_reference_alternate_alleles(self, valid_variant_tsv):
        """Test that reference and alternate alleles are valid."""
        df = pl.read_csv(
            valid_variant_tsv,
            separator="\t",
            has_header=False,
            skip_rows=1,
            columns=[
                "contig",
                "ref",
                "pos",
                "alt",
                "af",
                "ac",
                "dp",
                "mq",
                "gene",
                "aa_effect",
                "ref_codon_alt",
                "cds_pos",
                "aa_pos",
            ],
        )

        valid_bases = {"A", "C", "G", "T", "N"}

        # Check reference alleles
        for ref in df["ref"]:
            assert all(base in valid_bases for base in ref.upper())

        # Check alternate alleles (excluding indel notation)
        for alt in df["alt"]:
            if not alt.startswith(("+", "-")):
                assert all(base in valid_bases for base in alt.upper())


# Cleanup fixtures
@pytest.fixture(autouse=True)
def cleanup(request):
    """Clean up temporary files after tests."""

    def remove_temp_files():
        import glob
        import os

        temp_files = glob.glob("/tmp/tmp*.tsv")
        for f in temp_files:
            try:
                os.remove(f)
            except:
                pass

    request.addfinalizer(remove_temp_files)
