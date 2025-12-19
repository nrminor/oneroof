"""Tests for the variant metrics extractor."""

from pathlib import Path

import pytest
from reporting.extractors.variants import (
    VariantExtractedMetrics,
    classify_mutation_type,
    extract,
    load_variant_effects_tsv,
)


@pytest.fixture
def sample_tsv_content() -> str:
    """
    Sample SnpSift extractFields TSV content with known values.

    Contains 5 variants:
    - 3 SNPs (2 consensus AF>=0.5, 1 subclonal AF<0.5)
    - 1 insertion (consensus)
    - 1 deletion (subclonal)

    Effects:
    - 2 missense_variant
    - 2 synonymous_variant
    - 1 frameshift_variant
    """
    header = "CHROM\tREF\tPOS\tALT\tAF\tAC\tDP\tGEN[0].REF_DP\tGEN[0].ALT_DP\tGEN[0].ALT_FREQ\tMQ\tANN[0].GENE\tANN[0].EFFECT\tANN[0].HGVS_P\tANN[0].CDS_POS\tANN[0].AA_POS"
    rows = [
        "chr1\tA\t100\tG\t0.8\t80\t100\t20\t80\t0.8\t60\tgeneA\tmissense_variant\tp.Asp100Gly\t298\t100",
        "chr1\tC\t200\tT\t0.6\t60\t100\t40\t60\t0.6\t60\tgeneA\tsynonymous_variant\tp.Ala200=\t600\t200",
        "chr1\tG\t300\tA\t0.3\t30\t100\t70\t30\t0.3\t60\tgeneB\tmissense_variant\tp.Glu300Lys\t898\t300",
        "chr1\tA\t400\tATG\t0.7\t70\t100\t30\t70\t0.7\t60\tgeneB\tframeshift_variant\t.\t.\t.",
        "chr1\tGAT\t500\tG\t0.4\t40\t100\t60\t40\t0.4\t60\tgeneC\tsynonymous_variant\tp.Val500=\t1498\t500",
    ]
    return header + "\n" + "\n".join(rows)


@pytest.fixture
def sample_tsv_file(sample_tsv_content: str, tmp_path: Path) -> Path:
    """Create a temporary TSV file with sample content."""
    tsv_file = tmp_path / "test_variant_effects.tsv"
    tsv_file.write_text(sample_tsv_content)
    return tsv_file


@pytest.fixture
def empty_tsv_file(tmp_path: Path) -> Path:
    """Create an empty TSV file (no header, no data)."""
    tsv_file = tmp_path / "empty.tsv"
    tsv_file.write_text("")
    return tsv_file


@pytest.fixture
def header_only_tsv_file(tmp_path: Path) -> Path:
    """Create a TSV file with header but no data rows."""
    header = "CHROM\tREF\tPOS\tALT\tAF\tAC\tDP\tGEN[0].REF_DP\tGEN[0].ALT_DP\tGEN[0].ALT_FREQ\tMQ\tANN[0].GENE\tANN[0].EFFECT\tANN[0].HGVS_P\tANN[0].CDS_POS\tANN[0].AA_POS"
    tsv_file = tmp_path / "header_only.tsv"
    tsv_file.write_text(header + "\n")
    return tsv_file


@pytest.fixture
def all_consensus_tsv_file(tmp_path: Path) -> Path:
    """Create a TSV file where all variants are at consensus level (AF >= 0.5)."""
    header = "CHROM\tREF\tPOS\tALT\tAF\tAC\tDP\tGEN[0].REF_DP\tGEN[0].ALT_DP\tGEN[0].ALT_FREQ\tMQ\tANN[0].GENE\tANN[0].EFFECT\tANN[0].HGVS_P\tANN[0].CDS_POS\tANN[0].AA_POS"
    rows = [
        "chr1\tA\t100\tG\t0.9\t90\t100\t10\t90\t0.9\t60\tgeneA\tmissense_variant\tp.Asp100Gly\t298\t100",
        "chr1\tC\t200\tT\t0.8\t80\t100\t20\t80\t0.8\t60\tgeneA\tsynonymous_variant\tp.Ala200=\t600\t200",
        "chr1\tG\t300\tA\t0.5\t50\t100\t50\t50\t0.5\t60\tgeneB\tmissense_variant\tp.Glu300Lys\t898\t300",
    ]
    tsv_file = tmp_path / "all_consensus.tsv"
    tsv_file.write_text(header + "\n" + "\n".join(rows))
    return tsv_file


class TestClassifyMutationType:
    """Test mutation type classification."""

    def test_snp(self) -> None:
        """Test SNP classification (single base substitution)."""
        assert classify_mutation_type("A", "G") == "SNP"
        assert classify_mutation_type("C", "T") == "SNP"

    def test_insertion(self) -> None:
        """Test insertion classification (alt longer than ref)."""
        assert classify_mutation_type("A", "AT") == "insertion"
        assert classify_mutation_type("G", "GAT") == "insertion"

    def test_deletion(self) -> None:
        """Test deletion classification (ref longer than alt)."""
        assert classify_mutation_type("AT", "A") == "deletion"
        assert classify_mutation_type("GAT", "G") == "deletion"

    def test_mnp(self) -> None:
        """Test MNP classification (multi-nucleotide polymorphism)."""
        assert classify_mutation_type("AT", "GC") == "MNP"
        assert classify_mutation_type("AAA", "GGG") == "MNP"


class TestLoadVariantEffectsTsv:
    """Test TSV file loading."""

    def test_load_valid_tsv(self, sample_tsv_file: Path) -> None:
        """Test loading a valid TSV file."""
        df = load_variant_effects_tsv(sample_tsv_file)
        assert len(df) == 5

    def test_column_names_standardized(self, sample_tsv_file: Path) -> None:
        """Test that column names are standardized (lowercase, no brackets)."""
        df = load_variant_effects_tsv(sample_tsv_file)
        # Original: GEN[0].REF_DP -> gen_0_ref_dp
        assert "gen_0_ref_dp" in df.columns
        # Original: ANN[0].EFFECT -> ann_0_effect
        assert "ann_0_effect" in df.columns

    def test_empty_file_returns_empty_df(self, empty_tsv_file: Path) -> None:
        """Test that empty file returns empty DataFrame."""
        df = load_variant_effects_tsv(empty_tsv_file)
        assert len(df) == 0

    def test_header_only_returns_empty_df(self, header_only_tsv_file: Path) -> None:
        """Test that header-only file returns empty DataFrame."""
        df = load_variant_effects_tsv(header_only_tsv_file)
        assert len(df) == 0


class TestExtract:
    """Test the main extract function."""

    def test_extraction_returns_valid_model(self, sample_tsv_file: Path) -> None:
        """Test that extraction returns a valid VariantExtractedMetrics model."""
        metrics = extract("test_sample", sample_tsv_file)
        assert isinstance(metrics, VariantExtractedMetrics)

    def test_sample_id_set_correctly(self, sample_tsv_file: Path) -> None:
        """Test that sample_id is correctly set."""
        metrics = extract("my_sample_123", sample_tsv_file)
        assert metrics.sample_id == "my_sample_123"

    def test_total_called(self, sample_tsv_file: Path) -> None:
        """Test total variant count."""
        metrics = extract("test", sample_tsv_file)
        assert metrics.total_called == 5

    def test_consensus_vs_subclonal(self, sample_tsv_file: Path) -> None:
        """Test consensus vs subclonal classification with default threshold 0.5."""
        metrics = extract("test", sample_tsv_file, consensus_threshold=0.5)
        # AF values: 0.8, 0.6, 0.3, 0.7, 0.4
        # Consensus (>=0.5): 0.8, 0.6, 0.7 = 3
        # Subclonal (<0.5): 0.3, 0.4 = 2
        assert metrics.consensus_variants == 3
        assert metrics.subclonal_variants == 2

    def test_custom_consensus_threshold(self, sample_tsv_file: Path) -> None:
        """Test with custom consensus threshold."""
        metrics = extract("test", sample_tsv_file, consensus_threshold=0.7)
        # AF values: 0.8, 0.6, 0.3, 0.7, 0.4
        # Consensus (>=0.7): 0.8, 0.7 = 2
        # Subclonal (<0.7): 0.6, 0.3, 0.4 = 3
        assert metrics.consensus_variants == 2
        assert metrics.subclonal_variants == 3

    def test_mutation_type_counts(self, sample_tsv_file: Path) -> None:
        """Test mutation type counting."""
        metrics = extract("test", sample_tsv_file)
        # 3 SNPs: A>G, C>T, G>A
        # 1 insertion: A>ATG
        # 1 deletion: GAT>G
        assert metrics.snps == 3
        assert metrics.insertions == 1
        assert metrics.deletions == 1
        assert metrics.mnps == 0

    def test_effect_counting(self, sample_tsv_file: Path) -> None:
        """Test counting by effect type."""
        metrics = extract("test", sample_tsv_file)
        assert metrics.by_effect["missense_variant"] == 2
        assert metrics.by_effect["synonymous_variant"] == 2
        assert metrics.by_effect["frameshift_variant"] == 1

    def test_empty_file(self, empty_tsv_file: Path) -> None:
        """Test with empty file."""
        metrics = extract("test", empty_tsv_file)
        assert metrics.total_called == 0
        assert metrics.consensus_variants == 0
        assert metrics.subclonal_variants == 0
        assert metrics.snps == 0
        assert metrics.by_effect == {}

    def test_header_only_file(self, header_only_tsv_file: Path) -> None:
        """Test with header-only file."""
        metrics = extract("test", header_only_tsv_file)
        assert metrics.total_called == 0
        assert metrics.consensus_variants == 0
        assert metrics.by_effect == {}

    def test_all_consensus_variants(self, all_consensus_tsv_file: Path) -> None:
        """Test file where all variants are at consensus level."""
        metrics = extract("test", all_consensus_tsv_file, consensus_threshold=0.5)
        assert metrics.total_called == 3
        assert metrics.consensus_variants == 3
        assert metrics.subclonal_variants == 0

    def test_model_serialization(self, sample_tsv_file: Path) -> None:
        """Test that metrics can be serialized to JSON."""
        metrics = extract("test", sample_tsv_file)
        json_str = metrics.model_dump_json()
        assert "test" in json_str
        assert "total_called" in json_str
        assert "missense_variant" in json_str
