"""Tests for the MultiQC custom content generators."""

from pathlib import Path
from typing import Any

import pytest
import yaml
from reporting.multiqc import (
    generate_amplicon_efficiency_tsv,
    generate_amplicon_heatmap_tsv,
    generate_variant_bargraph_tsv,
)


def parse_multiqc_tsv_header(content: str) -> dict[str, Any]:
    """
    Parse and validate the YAML header from a MultiQC custom content TSV file.

    This function extracts the comment lines (starting with #), strips the #
    prefix, and parses them as YAML. If the YAML is invalid, this will raise
    an exception - ensuring our tests catch malformed headers.

    Args:
        content: The full content of a MultiQC TSV file

    Returns:
        Parsed YAML configuration dict

    Raises:
        yaml.YAMLError: If the header contains invalid YAML
    """
    lines = content.strip().split("\n")
    header_lines = []
    for line in lines:
        if line.startswith("#"):
            # Strip the "# " prefix
            header_lines.append(line[2:] if line.startswith("# ") else line[1:])
        else:
            break

    yaml_str = "\n".join(header_lines)
    return yaml.safe_load(yaml_str)


# --- Fixtures for variant data ---


@pytest.fixture
def samples_with_variants() -> dict[str, dict]:
    """Sample data with variant metrics."""
    return {
        "sample1": {
            "variants": {
                "snps": 10,
                "insertions": 2,
                "deletions": 3,
                "mnps": 1,
            },
        },
        "sample2": {
            "variants": {
                "snps": 25,
                "insertions": 5,
                "deletions": 0,
                "mnps": 0,
            },
        },
        "sample3": {
            "variants": {
                "snps": 5,
                "insertions": 1,
                "deletions": 1,
                "mnps": 2,
            },
        },
    }


@pytest.fixture
def samples_no_variants() -> dict[str, dict]:
    """Sample data with no variant metrics."""
    return {
        "sample1": {"coverage": {"mean_coverage": 100}},
        "sample2": {"coverage": {"mean_coverage": 200}},
    }


@pytest.fixture
def samples_zero_variants() -> dict[str, dict]:
    """Sample data with variant metrics but all zeros."""
    return {
        "sample1": {
            "variants": {
                "snps": 0,
                "insertions": 0,
                "deletions": 0,
                "mnps": 0,
            },
        },
    }


@pytest.fixture
def samples_mixed_variants() -> dict[str, dict]:
    """Sample data with some samples having variants, some not."""
    return {
        "sample1": {
            "variants": {
                "snps": 10,
                "insertions": 2,
                "deletions": 3,
                "mnps": 1,
            },
        },
        "sample2": {"coverage": {"mean_coverage": 200}},  # No variants key
        "sample3": {
            "variants": {
                "snps": 0,
                "insertions": 0,
                "deletions": 0,
                "mnps": 0,
            },
        },
    }


# --- Fixtures for amplicon data ---


@pytest.fixture
def amplicon_summary_tsv(tmp_path: Path) -> Path:
    """Create a sample amplicon_summary.tsv with varied performance."""
    content = "\n".join(
        [
            "sample_name\tamplicon_name\tstart_pos\tend_pos\treads",
            # amp1: good performer (high reads, no dropouts)
            "sampleA\tamp1\t100\t200\t500",
            "sampleB\tamp1\t100\t200\t450",
            "sampleC\tamp1\t100\t200\t480",
            # amp2: moderate performer
            "sampleA\tamp2\t300\t400\t100",
            "sampleB\tamp2\t300\t400\t80",
            "sampleC\tamp2\t300\t400\t90",
            # amp3: poor performer (low reads, some dropouts)
            "sampleA\tamp3\t500\t600\t10",
            "sampleB\tamp3\t500\t600\t0",  # dropout
            "sampleC\tamp3\t500\t600\t5",
        ],
    )
    tsv_file = tmp_path / "amplicon_summary.tsv"
    tsv_file.write_text(content)
    return tsv_file


@pytest.fixture
def empty_amplicon_summary(tmp_path: Path) -> Path:
    """Create an empty amplicon summary (header only)."""
    content = "sample_name\tamplicon_name\tstart_pos\tend_pos\treads"
    tsv_file = tmp_path / "amplicon_summary.tsv"
    tsv_file.write_text(content)
    return tsv_file


@pytest.fixture
def single_sample_amplicon_summary(tmp_path: Path) -> Path:
    """Create an amplicon summary with only one sample."""
    content = "\n".join(
        [
            "sample_name\tamplicon_name\tstart_pos\tend_pos\treads",
            "sampleA\tamp1\t100\t200\t500",
            "sampleA\tamp2\t300\t400\t100",
            "sampleA\tamp3\t500\t600\t10",
        ],
    )
    tsv_file = tmp_path / "amplicon_summary.tsv"
    tsv_file.write_text(content)
    return tsv_file


# --- Tests for generate_variant_bargraph_tsv ---


class TestGenerateVariantBargraphTsv:
    """Tests for the variant bar graph TSV generator."""

    def test_generates_valid_tsv(
        self, samples_with_variants: dict, tmp_path: Path
    ) -> None:
        """Test that a valid TSV file is generated with correct structure."""
        output_path = tmp_path / "variants_mqc.tsv"
        result = generate_variant_bargraph_tsv(samples_with_variants, output_path)

        assert result == output_path
        assert output_path.exists()

        content = output_path.read_text()
        lines = content.strip().split("\n")

        # Parse and validate YAML header - this will raise if YAML is invalid
        header_config = parse_multiqc_tsv_header(content)
        assert header_config["id"] == "oneroof_variants"
        assert header_config["plot_type"] == "bargraph"
        assert "section_name" in header_config
        assert "pconfig" in header_config

        # Find the data lines (after header comments)
        data_lines = [line for line in lines if not line.startswith("#")]
        assert len(data_lines) == 4  # header + 3 samples

        # Check header row
        assert data_lines[0] == "Sample\tSNPs\tInsertions\tDeletions\tMNPs"

        # Check data rows contain expected samples
        data_content = "\n".join(data_lines[1:])
        assert "sample1" in data_content
        assert "sample2" in data_content
        assert "sample3" in data_content

    def test_returns_none_when_no_variants(
        self, samples_no_variants: dict, tmp_path: Path
    ) -> None:
        """Test that None is returned when no samples have variant data."""
        output_path = tmp_path / "variants_mqc.tsv"
        result = generate_variant_bargraph_tsv(samples_no_variants, output_path)

        assert result is None
        assert not output_path.exists()

    def test_returns_none_when_all_zero_variants(
        self, samples_zero_variants: dict, tmp_path: Path
    ) -> None:
        """Test that None is returned when all variant counts are zero."""
        output_path = tmp_path / "variants_mqc.tsv"
        result = generate_variant_bargraph_tsv(samples_zero_variants, output_path)

        assert result is None
        assert not output_path.exists()

    def test_filters_samples_without_variants(
        self, samples_mixed_variants: dict, tmp_path: Path
    ) -> None:
        """Test that samples without variants or with zero variants are filtered."""
        output_path = tmp_path / "variants_mqc.tsv"
        result = generate_variant_bargraph_tsv(samples_mixed_variants, output_path)

        assert result == output_path
        content = output_path.read_text()

        # Only sample1 should be included (has non-zero variants)
        assert "sample1" in content
        assert "sample2" not in content  # No variants key
        assert "sample3" not in content  # All zeros

    def test_correct_variant_counts(
        self, samples_with_variants: dict, tmp_path: Path
    ) -> None:
        """Test that variant counts are correctly written."""
        output_path = tmp_path / "variants_mqc.tsv"
        generate_variant_bargraph_tsv(samples_with_variants, output_path)

        content = output_path.read_text()
        lines = [line for line in content.split("\n") if not line.startswith("#")]

        # Find sample1 row and verify counts
        sample1_line = next(line for line in lines if line.startswith("sample1"))
        parts = sample1_line.split("\t")
        assert parts == ["sample1", "10", "2", "3", "1"]  # SNPs, Ins, Del, MNPs


# --- Tests for generate_amplicon_efficiency_tsv ---


class TestGenerateAmpliconEfficiencyTsv:
    """Tests for the amplicon efficiency table TSV generator."""

    def test_generates_valid_tsv(
        self, amplicon_summary_tsv: Path, tmp_path: Path
    ) -> None:
        """Test that a valid TSV file is generated with correct structure."""
        output_path = tmp_path / "amplicon_efficiency_mqc.tsv"
        result = generate_amplicon_efficiency_tsv(amplicon_summary_tsv, output_path)

        assert result == output_path
        assert output_path.exists()

        content = output_path.read_text()
        lines = content.strip().split("\n")

        # Parse and validate YAML header - this will raise if YAML is invalid
        header_config = parse_multiqc_tsv_header(content)
        assert header_config["id"] == "oneroof_amplicon_efficiency"
        assert header_config["plot_type"] == "table"
        assert "section_name" in header_config
        assert "pconfig" in header_config
        assert "headers" in header_config
        # Verify nested config is properly structured (not a string)
        assert isinstance(header_config["headers"], dict)
        assert "performance_tier" in header_config["headers"]

        # Find the data lines
        data_lines = [line for line in lines if not line.startswith("#")]
        assert len(data_lines) == 4  # header + 3 amplicons

        # Check header row
        header = data_lines[0]
        assert "Amplicon" in header
        assert "median_reads" in header
        assert "dropout_pct" in header
        assert "performance_tier" in header

    def test_returns_none_for_empty_data(
        self, empty_amplicon_summary: Path, tmp_path: Path
    ) -> None:
        """Test that None is returned for empty amplicon summary."""
        output_path = tmp_path / "amplicon_efficiency_mqc.tsv"
        result = generate_amplicon_efficiency_tsv(empty_amplicon_summary, output_path)

        assert result is None
        assert not output_path.exists()

    def test_returns_none_for_missing_file(self, tmp_path: Path) -> None:
        """Test that None is returned when input file doesn't exist."""
        nonexistent = tmp_path / "nonexistent.tsv"
        output_path = tmp_path / "amplicon_efficiency_mqc.tsv"
        result = generate_amplicon_efficiency_tsv(nonexistent, output_path)

        assert result is None
        assert not output_path.exists()

    def test_performance_tier_assignment(
        self, amplicon_summary_tsv: Path, tmp_path: Path
    ) -> None:
        """Test that performance tiers are correctly assigned."""
        output_path = tmp_path / "amplicon_efficiency_mqc.tsv"
        generate_amplicon_efficiency_tsv(amplicon_summary_tsv, output_path)

        content = output_path.read_text()
        data_lines = [
            line for line in content.split("\n") if not line.startswith("#") and line
        ]

        # amp1 has highest reads (480 median) -> should be "good"
        amp1_line = next(line for line in data_lines if line.startswith("amp1"))
        assert "good" in amp1_line

        # amp3 has lowest reads (5 median) -> should be "poor"
        amp3_line = next(line for line in data_lines if line.startswith("amp3"))
        assert "poor" in amp3_line

    def test_dropout_rate_calculation(
        self, amplicon_summary_tsv: Path, tmp_path: Path
    ) -> None:
        """Test that dropout rates are correctly calculated."""
        output_path = tmp_path / "amplicon_efficiency_mqc.tsv"
        generate_amplicon_efficiency_tsv(amplicon_summary_tsv, output_path)

        content = output_path.read_text()
        data_lines = [
            line for line in content.split("\n") if not line.startswith("#") and line
        ]

        # amp3 has 1 dropout out of 3 samples = 33.3%
        amp3_line = next(line for line in data_lines if line.startswith("amp3"))
        parts = amp3_line.split("\t")
        dropout_pct = float(parts[2])  # dropout_pct column
        assert 33.0 <= dropout_pct <= 34.0  # ~33.3%

        # amp1 has no dropouts = 0%
        amp1_line = next(line for line in data_lines if line.startswith("amp1"))
        parts = amp1_line.split("\t")
        dropout_pct = float(parts[2])
        assert dropout_pct == 0.0

    def test_sorted_by_median_reads_descending(
        self, amplicon_summary_tsv: Path, tmp_path: Path
    ) -> None:
        """Test that amplicons are sorted by median reads (descending)."""
        output_path = tmp_path / "amplicon_efficiency_mqc.tsv"
        generate_amplicon_efficiency_tsv(amplicon_summary_tsv, output_path)

        content = output_path.read_text()
        data_lines = [
            line
            for line in content.split("\n")
            if not line.startswith("#") and line and not line.startswith("Amplicon")
        ]

        # Extract amplicon names in order
        amplicon_order = [line.split("\t")[0] for line in data_lines]

        # amp1 (highest) should be first, amp3 (lowest) should be last
        assert amplicon_order[0] == "amp1"
        assert amplicon_order[-1] == "amp3"


# --- Tests for generate_amplicon_heatmap_tsv ---


class TestGenerateAmpliconHeatmapTsv:
    """Tests for the amplicon heatmap TSV generator."""

    def test_generates_valid_tsv(
        self, amplicon_summary_tsv: Path, tmp_path: Path
    ) -> None:
        """Test that a valid TSV file is generated with correct structure."""
        output_path = tmp_path / "amplicon_heatmap_mqc.tsv"
        result = generate_amplicon_heatmap_tsv(amplicon_summary_tsv, output_path)

        assert result == output_path
        assert output_path.exists()

        content = output_path.read_text()
        lines = content.strip().split("\n")

        # Parse and validate YAML header - this will raise if YAML is invalid
        header_config = parse_multiqc_tsv_header(content)
        assert header_config["id"] == "oneroof_amplicon_heatmap"
        assert header_config["plot_type"] == "heatmap"
        assert "section_name" in header_config
        assert "pconfig" in header_config
        # Verify nested config is properly structured (not a string)
        assert isinstance(header_config["pconfig"], dict)

        # Find the data lines
        data_lines = [line for line in lines if not line.startswith("#")]
        assert len(data_lines) == 4  # header + 3 samples

    def test_returns_none_for_empty_data(
        self, empty_amplicon_summary: Path, tmp_path: Path
    ) -> None:
        """Test that None is returned for empty amplicon summary."""
        output_path = tmp_path / "amplicon_heatmap_mqc.tsv"
        result = generate_amplicon_heatmap_tsv(empty_amplicon_summary, output_path)

        assert result is None
        assert not output_path.exists()

    def test_returns_none_for_missing_file(self, tmp_path: Path) -> None:
        """Test that None is returned when input file doesn't exist."""
        nonexistent = tmp_path / "nonexistent.tsv"
        output_path = tmp_path / "amplicon_heatmap_mqc.tsv"
        result = generate_amplicon_heatmap_tsv(nonexistent, output_path)

        assert result is None
        assert not output_path.exists()

    def test_amplicons_sorted_by_genomic_position(
        self, amplicon_summary_tsv: Path, tmp_path: Path
    ) -> None:
        """Test that amplicon columns are sorted by genomic position."""
        output_path = tmp_path / "amplicon_heatmap_mqc.tsv"
        generate_amplicon_heatmap_tsv(amplicon_summary_tsv, output_path)

        content = output_path.read_text()
        data_lines = [line for line in content.split("\n") if not line.startswith("#")]

        # Get header row (amplicon column order)
        header = data_lines[0].split("\t")
        amplicon_cols = header[1:]  # Skip "Sample" column

        # amp1 (start=100) should be first, amp3 (start=500) should be last
        assert amplicon_cols == ["amp1", "amp2", "amp3"]

    def test_samples_sorted_alphabetically(
        self, amplicon_summary_tsv: Path, tmp_path: Path
    ) -> None:
        """Test that sample rows are sorted alphabetically."""
        output_path = tmp_path / "amplicon_heatmap_mqc.tsv"
        generate_amplicon_heatmap_tsv(amplicon_summary_tsv, output_path)

        content = output_path.read_text()
        data_lines = [
            line
            for line in content.split("\n")
            if not line.startswith("#") and line and not line.startswith("Sample")
        ]

        # Extract sample names in order
        sample_order = [line.split("\t")[0] for line in data_lines]

        assert sample_order == ["sampleA", "sampleB", "sampleC"]

    def test_correct_read_counts_in_matrix(
        self, amplicon_summary_tsv: Path, tmp_path: Path
    ) -> None:
        """Test that read counts are correctly placed in the matrix."""
        output_path = tmp_path / "amplicon_heatmap_mqc.tsv"
        generate_amplicon_heatmap_tsv(amplicon_summary_tsv, output_path)

        content = output_path.read_text()
        data_lines = [
            line
            for line in content.split("\n")
            if not line.startswith("#") and line and not line.startswith("Sample")
        ]

        # sampleA row: amp1=500, amp2=100, amp3=10
        sample_a_line = next(line for line in data_lines if line.startswith("sampleA"))
        parts = sample_a_line.split("\t")
        assert parts == ["sampleA", "500", "100", "10"]

        # sampleB row: amp1=450, amp2=80, amp3=0 (dropout)
        sample_b_line = next(line for line in data_lines if line.startswith("sampleB"))
        parts = sample_b_line.split("\t")
        assert parts == ["sampleB", "450", "80", "0"]

    def test_single_sample_works(
        self, single_sample_amplicon_summary: Path, tmp_path: Path
    ) -> None:
        """Test that heatmap works with a single sample."""
        output_path = tmp_path / "amplicon_heatmap_mqc.tsv"
        result = generate_amplicon_heatmap_tsv(
            single_sample_amplicon_summary, output_path
        )

        assert result == output_path
        content = output_path.read_text()
        data_lines = [
            line
            for line in content.split("\n")
            if not line.startswith("#") and line  # Filter empty lines
        ]

        # Should have header + 1 sample row
        assert len(data_lines) == 2
