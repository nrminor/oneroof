"""Tests for the metagenomics metrics extractor."""

from pathlib import Path

import pytest
from reporting.extractors.metagenomics import (
    MetagenomicsExtractedMetrics,
    MetagenomicsHit,
    extract,
    load_sylph_profile,
)


@pytest.fixture
def sample_sylph_content() -> str:
    """
    Sample Sylph profile TSV with known values.

    Contains 5 hits with taxonomic abundances summing to 80%,
    meaning 20% unknown fraction.
    """
    header = "Sample_file\tGenome_file\tTaxonomic_abundance\tSequence_abundance\tAdjusted_ANI\tEff_cov\tMedian_cov\tMean_cov_geq1\tContainment_ind\tNaive_ANI\tContig_name"
    rows = [
        "sample.fq\tgenome_A.fa\t35.5\t38.2\t98.5\t150.0\t145\t152.3\t5000/5500\t98.2\tSpecies_A chromosome",
        "sample.fq\tgenome_B.fa\t25.0\t27.1\t97.2\t100.0\t98\t101.5\t4000/4500\t97.0\tSpecies_B chromosome",
        "sample.fq\tgenome_C.fa\t12.5\t13.8\t96.8\t50.0\t48\t51.2\t2000/2500\t96.5\tSpecies_C chromosome",
        "sample.fq\tgenome_D.fa\t5.0\t5.5\t95.5\t20.0\t19\t20.8\t800/1000\t95.2\tSpecies_D chromosome",
        "sample.fq\tgenome_E.fa\t2.0\t2.2\t94.0\t8.0\t7\t8.3\t300/500\t93.8\tSpecies_E chromosome",
    ]
    return header + "\n" + "\n".join(rows)


@pytest.fixture
def sample_sylph_file(sample_sylph_content: str, tmp_path: Path) -> Path:
    """Create a temporary Sylph profile TSV file."""
    tsv_file = tmp_path / "sample_sylph.tsv"
    tsv_file.write_text(sample_sylph_content)
    return tsv_file


@pytest.fixture
def empty_sylph_file(tmp_path: Path) -> Path:
    """Create an empty Sylph profile file."""
    tsv_file = tmp_path / "empty_sylph.tsv"
    tsv_file.write_text("")
    return tsv_file


@pytest.fixture
def header_only_sylph_file(tmp_path: Path) -> Path:
    """Create a Sylph profile with header but no data."""
    header = "Sample_file\tGenome_file\tTaxonomic_abundance\tSequence_abundance\tAdjusted_ANI\tEff_cov"
    tsv_file = tmp_path / "header_only_sylph.tsv"
    tsv_file.write_text(header + "\n")
    return tsv_file


@pytest.fixture
def full_abundance_sylph_file(tmp_path: Path) -> Path:
    """Create a Sylph profile where abundances sum to 100%."""
    header = "Sample_file\tGenome_file\tTaxonomic_abundance\tSequence_abundance\tAdjusted_ANI\tEff_cov"
    rows = [
        "sample.fq\tgenome_A.fa\t60.0\t62.0\t98.5\t150.0",
        "sample.fq\tgenome_B.fa\t40.0\t38.0\t97.2\t100.0",
    ]
    tsv_file = tmp_path / "full_abundance_sylph.tsv"
    tsv_file.write_text(header + "\n" + "\n".join(rows))
    return tsv_file


@pytest.fixture
def single_hit_sylph_file(tmp_path: Path) -> Path:
    """Create a Sylph profile with a single hit."""
    header = "Sample_file\tGenome_file\tTaxonomic_abundance\tSequence_abundance\tAdjusted_ANI\tEff_cov"
    row = "sample.fq\tgenome_A.fa\t50.0\t52.0\t99.0\t200.0"
    tsv_file = tmp_path / "single_hit_sylph.tsv"
    tsv_file.write_text(header + "\n" + row)
    return tsv_file


class TestLoadSylphProfile:
    """Test Sylph profile loading."""

    def test_load_valid_profile(self, sample_sylph_file: Path) -> None:
        """Test loading a valid Sylph profile."""
        df = load_sylph_profile(sample_sylph_file)
        assert len(df) == 5

    def test_column_names_standardized(self, sample_sylph_file: Path) -> None:
        """Test that column names are standardized to lowercase."""
        df = load_sylph_profile(sample_sylph_file)
        assert "taxonomic_abundance" in df.columns
        assert "adjusted_ani" in df.columns
        assert "genome_file" in df.columns

    def test_empty_file_returns_empty_df(self, empty_sylph_file: Path) -> None:
        """Test that empty file returns empty DataFrame."""
        df = load_sylph_profile(empty_sylph_file)
        assert len(df) == 0

    def test_header_only_returns_empty_df(self, header_only_sylph_file: Path) -> None:
        """Test that header-only file returns empty DataFrame."""
        df = load_sylph_profile(header_only_sylph_file)
        assert len(df) == 0


class TestExtract:
    """Test the main extract function."""

    def test_extraction_returns_valid_model(self, sample_sylph_file: Path) -> None:
        """Test that extraction returns a valid MetagenomicsExtractedMetrics model."""
        metrics = extract("test_sample", sample_sylph_file)
        assert isinstance(metrics, MetagenomicsExtractedMetrics)

    def test_sample_id_set_correctly(self, sample_sylph_file: Path) -> None:
        """Test that sample_id is correctly set."""
        metrics = extract("my_sample_123", sample_sylph_file)
        assert metrics.sample_id == "my_sample_123"

    def test_total_hits(self, sample_sylph_file: Path) -> None:
        """Test total hits count."""
        metrics = extract("test", sample_sylph_file)
        assert metrics.total_hits == 5

    def test_total_abundance(self, sample_sylph_file: Path) -> None:
        """Test total abundance calculation."""
        metrics = extract("test", sample_sylph_file)
        # 35.5 + 25.0 + 12.5 + 5.0 + 2.0 = 80.0
        assert metrics.total_abundance == 80.0

    def test_unknown_fraction(self, sample_sylph_file: Path) -> None:
        """Test unknown fraction calculation."""
        metrics = extract("test", sample_sylph_file)
        # Total abundance is 80%, so unknown is 20% = 0.2
        assert metrics.unknown_fraction == pytest.approx(0.2)

    def test_top_hits_default_count(self, sample_sylph_file: Path) -> None:
        """Test that default top_n=5 returns 5 hits."""
        metrics = extract("test", sample_sylph_file)
        assert len(metrics.top_hits) == 5

    def test_top_hits_custom_count(self, sample_sylph_file: Path) -> None:
        """Test custom top_n parameter."""
        metrics = extract("test", sample_sylph_file, top_n=3)
        assert len(metrics.top_hits) == 3

    def test_top_hits_sorted_by_abundance(self, sample_sylph_file: Path) -> None:
        """Test that top hits are sorted by abundance descending."""
        metrics = extract("test", sample_sylph_file)
        abundances = [hit.relative_abundance for hit in metrics.top_hits]
        assert abundances == sorted(abundances, reverse=True)

    def test_top_hit_values(self, sample_sylph_file: Path) -> None:
        """Test that top hit has correct values."""
        metrics = extract("test", sample_sylph_file)
        top_hit = metrics.top_hits[0]
        assert top_hit.taxon == "genome_A.fa"
        assert top_hit.ani == 98.5
        assert top_hit.relative_abundance == 35.5

    def test_full_abundance_zero_unknown(self, full_abundance_sylph_file: Path) -> None:
        """Test that 100% abundance gives 0% unknown."""
        metrics = extract("test", full_abundance_sylph_file)
        assert metrics.total_abundance == 100.0
        assert metrics.unknown_fraction == 0.0

    def test_single_hit(self, single_hit_sylph_file: Path) -> None:
        """Test profile with single hit."""
        metrics = extract("test", single_hit_sylph_file)
        assert metrics.total_hits == 1
        assert len(metrics.top_hits) == 1
        assert metrics.top_hits[0].taxon == "genome_A.fa"
        assert metrics.unknown_fraction == pytest.approx(0.5)

    def test_empty_profile(self, empty_sylph_file: Path) -> None:
        """Test empty profile."""
        metrics = extract("test", empty_sylph_file)
        assert metrics.total_hits == 0
        assert len(metrics.top_hits) == 0
        assert metrics.total_abundance == 0.0
        assert metrics.unknown_fraction == 1.0

    def test_header_only_profile(self, header_only_sylph_file: Path) -> None:
        """Test header-only profile."""
        metrics = extract("test", header_only_sylph_file)
        assert metrics.total_hits == 0
        assert len(metrics.top_hits) == 0
        assert metrics.unknown_fraction == 1.0

    def test_model_serialization(self, sample_sylph_file: Path) -> None:
        """Test that metrics can be serialized to JSON."""
        metrics = extract("test", sample_sylph_file)
        json_str = metrics.model_dump_json()
        assert "test" in json_str
        assert "top_hits" in json_str
        assert "unknown_fraction" in json_str
        assert "genome_A.fa" in json_str


class TestMetagenomicsHit:
    """Test the MetagenomicsHit model."""

    def test_hit_creation(self) -> None:
        """Test creating a MetagenomicsHit."""
        hit = MetagenomicsHit(
            taxon="E. coli",
            ani=98.5,
            relative_abundance=25.0,
        )
        assert hit.taxon == "E. coli"
        assert hit.ani == 98.5
        assert hit.relative_abundance == 25.0

    def test_hit_serialization(self) -> None:
        """Test MetagenomicsHit JSON serialization."""
        hit = MetagenomicsHit(
            taxon="E. coli",
            ani=98.5,
            relative_abundance=25.0,
        )
        json_str = hit.model_dump_json()
        assert "E. coli" in json_str
        assert "98.5" in json_str
