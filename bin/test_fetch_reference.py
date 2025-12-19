"""Tests for fetch_reference.py module."""

import tempfile
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
from fetch_reference import (
    GENBANK_PATTERN,
    REFSEQ_PATTERN,
    app,
    is_local_path,
    is_ncbi_accession,
    normalize_sequence,
    validate_fasta,
    validate_genbank,
)
from typer.testing import CliRunner

runner = CliRunner()


# =============================================================================
# Test Fixtures
# =============================================================================


@pytest.fixture
def valid_fasta_file():
    """Create a valid FASTA file for testing."""
    content = """>NC_045512.2 Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1
ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCT
GTTCTCTAAACGAACTTTAAAATCTGTGTGGCTGTCACTCGGCTGCATGCTTAGTGCACT
"""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
        f.write(content)
        return Path(f.name)


@pytest.fixture
def fasta_with_invalid_chars():
    """Create a FASTA file with invalid characters for testing."""
    content = """>test_seq Invalid characters
ATCGXYZATCG123ATCG
"""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
        f.write(content)
        return Path(f.name)


@pytest.fixture
def valid_genbank_file():
    """Create a valid GenBank file for testing."""
    # GenBank format requires proper spacing and molecule_type annotation
    content = """LOCUS       NC_045512                100 bp    RNA     linear   VRL 18-JUL-2020
DEFINITION  Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1.
ACCESSION   NC_045512
VERSION     NC_045512.2
KEYWORDS    RefSeq.
SOURCE      Severe acute respiratory syndrome coronavirus 2 (SARS-CoV-2)
  ORGANISM  Severe acute respiratory syndrome coronavirus 2
            Viruses; Coronaviridae; Betacoronavirus.
FEATURES             Location/Qualifiers
     source          1..100
                     /organism="SARS-CoV-2"
                     /mol_type="genomic RNA"
ORIGIN
        1 attaaaggtt tataccttcc caggtaacaa accaaccaac tttcgatctc ttgtagatct
       61 gttctctaaa cgaactttaa aatctgtgtg gctgtcactc
//
"""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".gbk", delete=False) as f:
        f.write(content)
        return Path(f.name)


@pytest.fixture
def empty_fasta_file():
    """Create an empty FASTA file for testing."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
        f.write("")
        return Path(f.name)


@pytest.fixture
def fasta_header_only():
    """Create a FASTA file with header but no sequence."""
    content = """>test_seq No sequence here
"""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
        f.write(content)
        return Path(f.name)


@pytest.fixture
def html_error_file():
    """Create a file with HTML content starting with a FASTA-like header."""
    # This simulates a failed download that might look like FASTA initially
    content = """>error
<!DOCTYPE html><html><body>Error</body></html>
"""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
        f.write(content)
        return Path(f.name)


# =============================================================================
# Test Accession Pattern Detection
# =============================================================================


class TestAccessionPatterns:
    """Tests for NCBI accession pattern detection."""

    def test_refseq_pattern_valid(self):
        """Test valid RefSeq accession patterns."""
        valid_refseq = [
            "NC_045512.2",
            "NC_045512",
            "NM_001301717.2",
            "NP_000509.1",
            "XM_123456.1",
            "XP_789012.3",
        ]
        for acc in valid_refseq:
            assert REFSEQ_PATTERN.match(acc), f"Should match RefSeq: {acc}"

    def test_refseq_pattern_invalid(self):
        """Test invalid RefSeq patterns."""
        invalid_refseq = [
            "NC045512.2",  # Missing underscore
            "nc_045512.2",  # Lowercase
            "N_045512.2",  # Single letter prefix
            "NCX_045512.2",  # Three letter prefix
        ]
        for acc in invalid_refseq:
            assert not REFSEQ_PATTERN.match(acc), f"Should not match RefSeq: {acc}"

    def test_genbank_pattern_valid(self):
        """Test valid GenBank accession patterns."""
        valid_genbank = [
            "MN908947.3",
            "MN908947",
            "AB123456.1",
            "U12345.1",
            "AF123456789.1",
        ]
        for acc in valid_genbank:
            assert GENBANK_PATTERN.match(acc), f"Should match GenBank: {acc}"

    def test_genbank_pattern_invalid(self):
        """Test invalid GenBank patterns."""
        invalid_genbank = [
            "mn908947.3",  # Lowercase
            "123456.1",  # No letter prefix
            "ABCD.1",  # No numbers
        ]
        for acc in invalid_genbank:
            assert not GENBANK_PATTERN.match(acc), f"Should not match GenBank: {acc}"


class TestIsNcbiAccession:
    """Tests for is_ncbi_accession function."""

    def test_refseq_accessions(self):
        """Test RefSeq accessions are recognized."""
        assert is_ncbi_accession("NC_045512.2")
        assert is_ncbi_accession("NM_001301717.2")
        assert is_ncbi_accession("  NC_045512.2  ")  # With whitespace

    def test_genbank_accessions(self):
        """Test GenBank accessions are recognized."""
        assert is_ncbi_accession("MN908947.3")
        assert is_ncbi_accession("AB123456.1")

    def test_file_paths_not_accessions(self):
        """Test file paths are not recognized as accessions."""
        assert not is_ncbi_accession("/path/to/file.fasta")
        assert not is_ncbi_accession("./reference.fasta")
        assert not is_ncbi_accession("reference.fasta")


class TestIsLocalPath:
    """Tests for is_local_path function."""

    def test_absolute_paths(self):
        """Test absolute paths are recognized."""
        assert is_local_path("/path/to/file.fasta")
        assert is_local_path("/home/user/reference.gbk")

    def test_relative_paths(self):
        """Test relative paths are recognized."""
        assert is_local_path("./reference.fasta")
        assert is_local_path("../data/reference.gbk")

    def test_home_paths(self):
        """Test home directory paths are recognized."""
        assert is_local_path("~/reference.fasta")

    def test_existing_file(self, valid_fasta_file):
        """Test existing files are recognized as paths."""
        assert is_local_path(str(valid_fasta_file))

    def test_accessions_not_paths(self):
        """Test accessions are not recognized as paths."""
        # Note: This depends on the file not existing
        assert not is_local_path("NC_045512.2")
        assert not is_local_path("MN908947.3")


# =============================================================================
# Test Sequence Normalization
# =============================================================================


class TestNormalizeSequence:
    """Tests for sequence normalization."""

    def test_valid_sequence_unchanged(self):
        """Test valid IUPAC sequences are unchanged."""
        seq = "ACGTUMRWSYKVHDBN"
        normalized, replacements = normalize_sequence(seq)
        assert normalized == seq
        assert replacements == 0

    def test_lowercase_uppercased(self):
        """Test lowercase bases are uppercased."""
        seq = "acgt"
        normalized, replacements = normalize_sequence(seq)
        assert normalized == "ACGT"
        assert replacements == 0

    def test_invalid_chars_replaced(self):
        """Test invalid characters are replaced with N."""
        # Use characters that are definitely not IUPAC: X, Z, and numbers
        seq = "ATCGXZATCG"
        normalized, replacements = normalize_sequence(seq)
        assert normalized == "ATCGNNATCG"
        assert replacements == 2

    def test_numbers_replaced(self):
        """Test numbers are replaced with N."""
        seq = "ATCG123ATCG"
        normalized, replacements = normalize_sequence(seq)
        assert normalized == "ATCGNNNATCG"
        assert replacements == 3

    def test_whitespace_stripped(self):
        """Test whitespace is stripped silently."""
        seq = "ATCG ATCG\nATCG\tATCG"
        normalized, replacements = normalize_sequence(seq)
        assert normalized == "ATCGATCGATCGATCG"
        assert replacements == 0


# =============================================================================
# Test File Validation
# =============================================================================


class TestValidateFasta:
    """Tests for FASTA file validation."""

    def test_valid_fasta(self, valid_fasta_file):
        """Test valid FASTA file is parsed correctly."""
        records = validate_fasta(valid_fasta_file)
        assert len(records) == 1
        assert records[0].id == "NC_045512.2"

    def test_empty_fasta_raises(self, empty_fasta_file):
        """Test empty FASTA file raises error."""
        # Empty file should result in no records, which triggers exit
        # typer.Exit raises click.exceptions.Exit
        from click.exceptions import Exit

        with pytest.raises(Exit) as exc_info:
            validate_fasta(empty_fasta_file)
        assert exc_info.value.exit_code == 1

    def test_html_error_raises(self, html_error_file):
        """Test HTML content (failed download) raises error."""
        # HTML content in sequence triggers detection
        from click.exceptions import Exit

        with pytest.raises(Exit) as exc_info:
            validate_fasta(html_error_file)
        assert exc_info.value.exit_code == 1

    def test_header_only_fasta(self, fasta_header_only):
        """Test FASTA with header but empty sequence."""
        # Biopython parses this but sequence is empty
        records = validate_fasta(fasta_header_only)
        # Empty sequence should still be valid (Biopython accepts it)
        assert len(records) == 1


class TestValidateGenbank:
    """Tests for GenBank file validation."""

    def test_valid_genbank(self, valid_genbank_file):
        """Test valid GenBank file is parsed correctly."""
        records = validate_genbank(valid_genbank_file)
        assert len(records) == 1
        # GenBank VERSION line includes version suffix
        assert records[0].id == "NC_045512.2"


# =============================================================================
# Test CLI Commands
# =============================================================================


class TestFastaCommand:
    """Tests for the fasta CLI command."""

    def test_local_file_success(self, valid_fasta_file):
        """Test resolving a local FASTA file."""
        with tempfile.NamedTemporaryFile(suffix=".fasta", delete=False) as out:
            output_path = Path(out.name)

        result = runner.invoke(
            app,
            ["fasta", str(valid_fasta_file), "--output", str(output_path)],
        )
        assert result.exit_code == 0
        assert output_path.exists()
        assert "Success" in result.stdout

    def test_local_file_not_found(self):
        """Test error when local file doesn't exist."""
        result = runner.invoke(
            app,
            ["fasta", "/nonexistent/file.fasta", "--output", "out.fasta"],
        )
        assert result.exit_code == 1
        assert "not found" in result.stdout.lower() or "Error" in result.stdout

    def test_invalid_input_format(self):
        """Test error for invalid input (not path or accession)."""
        result = runner.invoke(
            app,
            ["fasta", "not_valid_input", "--output", "out.fasta"],
        )
        assert result.exit_code == 1
        assert "neither" in result.stdout.lower() or "Error" in result.stdout

    def test_normalization_warning(self, fasta_with_invalid_chars):
        """Test warning is shown when characters are normalized."""
        with tempfile.NamedTemporaryFile(suffix=".fasta", delete=False) as out:
            output_path = Path(out.name)

        result = runner.invoke(
            app,
            ["fasta", str(fasta_with_invalid_chars), "--output", str(output_path)],
        )
        assert result.exit_code == 0
        assert "Warning" in result.stdout
        assert "replaced" in result.stdout.lower()


class TestGenbankCommand:
    """Tests for the genbank CLI command."""

    def test_local_file_success(self, valid_genbank_file):
        """Test resolving a local GenBank file."""
        with tempfile.NamedTemporaryFile(suffix=".gbk", delete=False) as out:
            output_path = Path(out.name)

        result = runner.invoke(
            app,
            ["genbank", str(valid_genbank_file), "--output", str(output_path)],
        )
        assert result.exit_code == 0
        assert output_path.exists()
        assert "Success" in result.stdout


# =============================================================================
# Test NCBI Fetching (Mocked)
# =============================================================================


class TestNcbiFetching:
    """Tests for NCBI Entrez fetching (mocked to avoid network calls)."""

    @patch("fetch_reference.Entrez.efetch")
    def test_fetch_fasta_from_ncbi(self, mock_efetch):
        """Test fetching FASTA from NCBI accession."""
        mock_handle = MagicMock()
        mock_handle.read.return_value = """>NC_045512.2 Test sequence
ATCGATCGATCG
"""
        mock_handle.close = MagicMock()
        mock_efetch.return_value = mock_handle

        with tempfile.NamedTemporaryFile(suffix=".fasta", delete=False) as out:
            output_path = Path(out.name)

        result = runner.invoke(
            app,
            ["fasta", "NC_045512.2", "--output", str(output_path)],
        )

        assert result.exit_code == 0
        assert output_path.exists()
        mock_efetch.assert_called_once()

    @patch("fetch_reference.Entrez.efetch")
    def test_fetch_genbank_from_ncbi(self, mock_efetch):
        """Test fetching GenBank from NCBI accession."""
        mock_handle = MagicMock()
        # Proper GenBank format with molecule_type
        mock_handle.read.return_value = """LOCUS       NC_045512                 12 bp    RNA     linear   VRL 01-JAN-2020
DEFINITION  Test sequence.
ACCESSION   NC_045512
VERSION     NC_045512.2
FEATURES             Location/Qualifiers
     source          1..12
                     /mol_type="genomic RNA"
ORIGIN
        1 atcgatcgat cg
//
"""
        mock_handle.close = MagicMock()
        mock_efetch.return_value = mock_handle

        with tempfile.NamedTemporaryFile(suffix=".gbk", delete=False) as out:
            output_path = Path(out.name)

        result = runner.invoke(
            app,
            ["genbank", "NC_045512.2", "--output", str(output_path)],
        )

        assert result.exit_code == 0
        assert output_path.exists()
        mock_efetch.assert_called_once()

    @patch("fetch_reference.Entrez.efetch")
    def test_fetch_retry_on_failure(self, mock_efetch):
        """Test retry logic on NCBI fetch failure."""
        # First two calls fail, third succeeds
        mock_efetch.side_effect = [
            Exception("Network error"),
            Exception("Timeout"),
            MagicMock(
                read=MagicMock(return_value=">test\nATCG\n"),
                close=MagicMock(),
            ),
        ]

        with tempfile.NamedTemporaryFile(suffix=".fasta", delete=False) as out:
            output_path = Path(out.name)

        result = runner.invoke(
            app,
            ["fasta", "NC_045512.2", "--output", str(output_path)],
        )

        assert result.exit_code == 0
        assert mock_efetch.call_count == 3
