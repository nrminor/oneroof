#!/usr/bin/env python3
"""
Test module for concat_consensus.py

Tests the functionality of concatenating multiple consensus FASTA files into a single output file.
"""

import sys
from pathlib import Path
from textwrap import dedent
from unittest.mock import patch

import pytest
from Bio import SeqIO

# Add the bin directory to the Python path to import the module
sys.path.insert(0, str(Path(__file__).parent))
from concat_consensus import main


@pytest.fixture
def temp_dir(tmp_path):
    """Create a temporary directory for test files."""
    return tmp_path


@pytest.fixture
def sample_fasta_files(temp_dir):
    """Create sample consensus FASTA files for testing."""
    # Sample 1: Single sequence
    sample1_content = dedent("""\
        >contig1
        ATCGATCGATCG
        ATCGATCGATCG
        """)
    sample1_path = temp_dir / "sample1.consensus.fasta"
    sample1_path.write_text(sample1_content)

    # Sample 2: Multiple sequences
    sample2_content = dedent("""\
        >contig1
        GCTAGCTAGCTA
        >contig2
        TTTAAATTTAAA
        """)
    sample2_path = temp_dir / "sample2.consensus.fasta"
    sample2_path.write_text(sample2_content)

    # Sample 3: Empty sequences
    sample3_content = dedent("""\
        >contig1

        >contig2
        AAAA
        """)
    sample3_path = temp_dir / "sample3.consensus.fasta"
    sample3_path.write_text(sample3_content)

    return temp_dir


@pytest.fixture
def empty_fasta_file(temp_dir):
    """Create an empty FASTA file."""
    empty_path = temp_dir / "empty.consensus.fasta"
    empty_path.write_text("")
    return temp_dir


@pytest.fixture
def malformed_fasta_file(temp_dir):
    """Create a malformed FASTA file."""
    malformed_content = "This is not a valid FASTA file\nNo headers here"
    malformed_path = temp_dir / "malformed.consensus.fasta"
    malformed_path.write_text(malformed_content)
    return temp_dir


@pytest.fixture
def single_fasta_file(temp_dir):
    """Create a single consensus FASTA file."""
    content = dedent("""\
        >contig1
        ATCGATCGATCGATCG
        >contig2
        GCTAGCTAGCTAGCTA
        """)
    path = temp_dir / "single.consensus.fasta"
    path.write_text(content)
    return temp_dir


class TestConcatConsensus:
    """Test class for concat_consensus.py functionality."""

    def test_main_with_valid_inputs(self, sample_fasta_files, monkeypatch):
        """Test main function with valid input files."""
        # Change to the temp directory
        monkeypatch.chdir(sample_fasta_files)

        # Run the main function
        main()

        # Check that output file was created
        output_file = Path("all_sample_consensus.fasta")
        assert output_file.exists()

        # Parse and verify the output
        records = list(SeqIO.parse(output_file, "fasta"))
        assert len(records) == 3

        # Check sample1: should have concatenated sequence
        sample1_record = next(r for r in records if r.id == "sample1")
        assert str(sample1_record.seq) == "ATCGATCGATCGATCGATCGATCG"

        # Check sample2: should have both sequences concatenated
        sample2_record = next(r for r in records if r.id == "sample2")
        assert str(sample2_record.seq) == "GCTAGCTAGCTATTTAAATTTAAA"

        # Check sample3: should have empty string + "AAAA"
        sample3_record = next(r for r in records if r.id == "sample3")
        assert str(sample3_record.seq) == "AAAA"

    def test_main_with_single_file(self, single_fasta_file, monkeypatch):
        """Test main function with a single input file."""
        monkeypatch.chdir(single_fasta_file)

        main()

        output_file = Path("all_sample_consensus.fasta")
        assert output_file.exists()

        records = list(SeqIO.parse(output_file, "fasta"))
        assert len(records) == 1
        assert records[0].id == "single"
        assert str(records[0].seq) == "ATCGATCGATCGATCGGCTAGCTAGCTAGCTA"

    def test_main_with_empty_fasta(self, empty_fasta_file, monkeypatch):
        """Test main function with an empty FASTA file."""
        monkeypatch.chdir(empty_fasta_file)

        main()

        output_file = Path("all_sample_consensus.fasta")
        assert output_file.exists()

        records = list(SeqIO.parse(output_file, "fasta"))
        assert len(records) == 1
        assert records[0].id == "empty"
        assert str(records[0].seq) == ""

    def test_main_with_no_consensus_files(self, temp_dir, monkeypatch):
        """Test main function when no consensus files are found."""
        monkeypatch.chdir(temp_dir)

        # Should raise AssertionError
        with pytest.raises(AssertionError) as excinfo:
            main()

        assert "Please double check that the working directory" in str(excinfo.value)

    def test_main_with_malformed_fasta(self, malformed_fasta_file, monkeypatch):
        """Test main function with a malformed FASTA file."""
        monkeypatch.chdir(malformed_fasta_file)

        # This should still work - BioPython will just not parse any sequences
        main()

        output_file = Path("all_sample_consensus.fasta")
        assert output_file.exists()

        records = list(SeqIO.parse(output_file, "fasta"))
        assert len(records) == 1
        assert records[0].id == "malformed"
        assert str(records[0].seq) == ""  # No valid sequences parsed

    def test_file_naming_with_different_extensions(self, temp_dir, monkeypatch):
        """Test that .consensus.fasta and .consensus.fa files are processed."""
        monkeypatch.chdir(temp_dir)

        # Create files with different extensions
        (temp_dir / "sample.consensus.fasta").write_text(">seq\nATCG")
        (temp_dir / "other.fasta").write_text(">seq\nGCTA")  # Should NOT match
        (temp_dir / "another.consensus.fa").write_text(
            ">seq\nTTTT",
        )  # Should match - .fa is valid

        main()

        output_file = Path("all_sample_consensus.fasta")
        records = list(SeqIO.parse(output_file, "fasta"))

        # Both .consensus.fasta and .consensus.fa files should be processed
        assert len(records) == 2
        record_ids = {r.id for r in records}
        assert record_ids == {"sample", "another.consensus.fa"}

    def test_sample_name_extraction(self, temp_dir, monkeypatch):
        """Test correct extraction of sample names from file paths."""
        monkeypatch.chdir(temp_dir)

        # Create files with various naming patterns
        files = {
            "simple.consensus.fasta": "simple",
            "sample_with_underscores.consensus.fasta": "sample_with_underscores",
            "sample.with.dots.consensus.fasta": "sample.with.dots",
            "sample-with-dashes.consensus.fasta": "sample-with-dashes",
        }

        for filename, expected_id in files.items():
            (temp_dir / filename).write_text(">seq\nATCG")

        main()

        output_file = Path("all_sample_consensus.fasta")
        records = list(SeqIO.parse(output_file, "fasta"))

        # Check that all sample names were extracted correctly
        record_ids = {r.id for r in records}
        expected_ids = set(files.values())
        assert record_ids == expected_ids

    def test_large_sequences(self, temp_dir, monkeypatch):
        """Test handling of large sequences."""
        monkeypatch.chdir(temp_dir)

        # Create a file with a large sequence
        large_seq = "A" * 10000 + "T" * 10000 + "C" * 10000 + "G" * 10000
        content = f">large_contig\n{large_seq}"
        (temp_dir / "large.consensus.fasta").write_text(content)

        main()

        output_file = Path("all_sample_consensus.fasta")
        records = list(SeqIO.parse(output_file, "fasta"))

        assert len(records) == 1
        assert len(str(records[0].seq)) == 40000
        assert str(records[0].seq) == large_seq

    def test_multiline_fasta_format(self, temp_dir, monkeypatch):
        """Test handling of FASTA files with sequences split across multiple lines."""
        monkeypatch.chdir(temp_dir)

        # Create a multiline FASTA file
        content = dedent("""\
            >contig1
            ATCGATCGATCG
            ATCGATCGATCG
            ATCGATCGATCG
            >contig2
            GCTAGCTAGCTA
            GCTAGCTAGCTA
            """)
        (temp_dir / "multiline.consensus.fasta").write_text(content)

        main()

        output_file = Path("all_sample_consensus.fasta")
        records = list(SeqIO.parse(output_file, "fasta"))

        assert len(records) == 1
        expected_seq = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGGCTAGCTAGCTAGCTAGCTAGCTA"
        assert str(records[0].seq) == expected_seq

    def test_special_characters_in_sequences(self, temp_dir, monkeypatch):
        """Test handling of sequences with special characters (lowercase, gaps, etc.)."""
        monkeypatch.chdir(temp_dir)

        # Create a file with mixed case and gap characters
        content = dedent("""\
            >contig1
            ATCGatcg
            >contig2
            NNNN----NNNN
            >contig3
            RYSWKMrykswm
            """)
        (temp_dir / "special.consensus.fasta").write_text(content)

        main()

        output_file = Path("all_sample_consensus.fasta")
        records = list(SeqIO.parse(output_file, "fasta"))

        assert len(records) == 1
        # BioPython preserves the original characters
        assert str(records[0].seq) == "ATCGatcgNNNN----NNNNRYSWKMrykswm"

    def test_output_file_overwrite(self, temp_dir, monkeypatch):
        """Test that the output file is overwritten if it already exists."""
        monkeypatch.chdir(temp_dir)

        # Create an existing output file
        output_file = Path("all_sample_consensus.fasta")
        output_file.write_text(">old_data\nOLDOLDOLD")

        # Create a new consensus file
        (temp_dir / "new.consensus.fasta").write_text(">new\nNEWNEWNEW")

        main()

        # Check that the old content was overwritten
        records = list(SeqIO.parse(output_file, "fasta"))
        assert len(records) == 1
        assert records[0].id == "new"
        assert str(records[0].seq) == "NEWNEWNEW"

    def test_permission_error_handling(self, temp_dir, monkeypatch):
        """Test handling of permission errors when writing output file."""
        monkeypatch.chdir(temp_dir)

        (temp_dir / "sample.consensus.fasta").write_text(">seq\nATCG")

        # Mock the open function to raise a PermissionError when writing
        original_open = open

        def mock_open_func(file, mode="r", *args, **kwargs):
            if "w" in mode and "all_sample_consensus.fasta" in str(file):
                raise PermissionError("Permission denied")
            return original_open(file, mode, *args, **kwargs)

        with patch("builtins.open", side_effect=mock_open_func):
            with pytest.raises(PermissionError):
                main()

    def test_glob_pattern_matching(self, temp_dir, monkeypatch):
        """Test that the glob pattern correctly matches consensus files."""
        monkeypatch.chdir(temp_dir)

        # Create files that should and shouldn't match
        # Pattern is *.consensus.fa* so matches .fasta, .fa, .fastq etc.
        matching_files = [
            "sample1.consensus.fasta",
            "sample2.consensus.fasta",
            "SAMPLE3.consensus.fasta",  # Test case sensitivity
            "sample4.consensus.fa",  # .fa extension also matches
        ]

        non_matching_files = [
            "sample_consensus.fasta",  # Missing dot before consensus
            "sample.consensus",  # Missing extension after .consensus
            "consensus.fasta",  # Missing sample name prefix
        ]

        for filename in matching_files:
            (temp_dir / filename).write_text(f">{filename}\nATCG")

        for filename in non_matching_files:
            (temp_dir / filename).write_text(f">{filename}\nGCTA")

        main()

        output_file = Path("all_sample_consensus.fasta")
        records = list(SeqIO.parse(output_file, "fasta"))

        # Only matching files should be processed
        assert len(records) == len(matching_files)
        record_ids = {r.id for r in records}
        # Note: .consensus.fasta is stripped, but .consensus.fa is not
        expected_ids = {"sample1", "sample2", "SAMPLE3", "sample4.consensus.fa"}
        assert record_ids == expected_ids


@pytest.mark.parametrize(
    "num_files,num_seqs_per_file",
    [
        (10, 5),
        (50, 2),
        (100, 1),
    ],
)
def test_performance_with_many_files(
    temp_dir,
    monkeypatch,
    num_files,
    num_seqs_per_file,
):
    """Test performance with many input files."""
    monkeypatch.chdir(temp_dir)

    # Create many consensus files
    for i in range(num_files):
        content = ""
        for j in range(num_seqs_per_file):
            content += f">contig{j}\nATCGATCG\n"
        (temp_dir / f"sample{i}.consensus.fasta").write_text(content)

    main()

    output_file = Path("all_sample_consensus.fasta")
    records = list(SeqIO.parse(output_file, "fasta"))

    assert len(records) == num_files
    # Each record should have all sequences concatenated
    for record in records:
        assert len(str(record.seq)) == num_seqs_per_file * 8  # "ATCGATCG" = 8 chars


def test_unicode_handling(temp_dir, monkeypatch):
    """Test handling of unicode characters in file names and content."""
    monkeypatch.chdir(temp_dir)

    # Create a file with unicode in the name (if filesystem supports it)
    try:
        unicode_filename = "sample_é.consensus.fasta"
        (temp_dir / unicode_filename).write_text(">seq\nATCG", encoding="utf-8")

        main()

        output_file = Path("all_sample_consensus.fasta")
        records = list(SeqIO.parse(output_file, "fasta"))

        assert len(records) == 1
        assert records[0].id == "sample_é"
    except (UnicodeEncodeError, OSError):
        # Skip test if filesystem doesn't support unicode filenames
        pytest.skip("Filesystem doesn't support unicode filenames")


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
