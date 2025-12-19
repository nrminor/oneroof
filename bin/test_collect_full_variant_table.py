"""Tests for collect_full_variant_table.py module."""

import polars as pl
import pytest

from collect_full_variant_table import (
    DEFAULT_CONSENSUS_THRESHOLD,
    FINAL_COLUMNS,
    VARIANT_EFFECTS_SUFFIX,
    _parse_hgvs_to_short,
    _three_letter_to_one,
    add_cross_sample_metrics,
    add_derived_columns,
    discover_input_files,
    is_consensus_expr,
    load_and_concat,
    main,
    mutation_type_expr,
    parse_command_line_args,
    select_and_sort,
    variant_id_expr,
)


SNPSIFT_HEADER = "CHROM\tREF\tPOS\tALT\tAF\tAC\tDP\tGEN[0].REF_DP\tGEN[0].ALT_DP\tGEN[0].ALT_FREQ\tMQ\tANN[0].GENE\tANN[0].EFFECT\tANN[0].HGVS_P\tANN[0].CDS_POS\tANN[0].AA_POS"


@pytest.fixture
def sample_variant_data():
    """Sample variant data rows matching SnpSift output format."""
    return [
        "NC_045512.2\tA\t241\tT\t0.98\t98\t100\t2\t98\t0.98\t60\t5UTR\tnon_coding_transcript_variant\t\t\t",
        "NC_045512.2\tC\t3037\tT\t0.95\t95\t100\t5\t95\t0.95\t60\tORF1ab\tsynonymous_variant\tp.Phe924Phe\t2772\t924",
        "NC_045512.2\tA\t23403\tG\t0.99\t99\t100\t1\t99\t0.99\t60\tS\tmissense_variant\tp.Asp614Gly\t1841\t614",
        "NC_045512.2\tG\t28881\tA\t0.45\t45\t100\t55\t45\t0.45\t55\tN\tmissense_variant\tp.Arg203Lys\t608\t203",
    ]


@pytest.fixture
def input_dir_with_samples(sample_variant_data, tmp_path):
    """Create a temporary directory with multiple sample variant files."""
    for sample_id in ["sample_01", "sample_02", "sample_03"]:
        file_path = tmp_path / f"{sample_id}{VARIANT_EFFECTS_SUFFIX}"
        with open(file_path, "w") as f:
            f.write(SNPSIFT_HEADER + "\n")
            for row in sample_variant_data:
                f.write(row + "\n")
    return tmp_path


@pytest.fixture
def input_dir_single_sample(sample_variant_data, tmp_path):
    """Create a temporary directory with a single sample variant file."""
    file_path = tmp_path / f"sample_01{VARIANT_EFFECTS_SUFFIX}"
    with open(file_path, "w") as f:
        f.write(SNPSIFT_HEADER + "\n")
        for row in sample_variant_data:
            f.write(row + "\n")
    return tmp_path


@pytest.fixture
def empty_input_dir(tmp_path):
    """Create an empty temporary directory."""
    return tmp_path


class TestAminoAcidConversion:
    """Test amino acid conversion utilities."""

    def test_three_letter_to_one_standard(self):
        """Test standard amino acid conversions."""
        assert _three_letter_to_one("Asp") == "D"
        assert _three_letter_to_one("Gly") == "G"
        assert _three_letter_to_one("Phe") == "F"
        assert _three_letter_to_one("Arg") == "R"
        assert _three_letter_to_one("Lys") == "K"

    def test_three_letter_to_one_case_insensitive(self):
        """Test case insensitivity."""
        assert _three_letter_to_one("ASP") == "D"
        assert _three_letter_to_one("asp") == "D"
        assert _three_letter_to_one("Asp") == "D"

    def test_three_letter_to_one_stop_codon(self):
        """Test stop codon conversion."""
        assert _three_letter_to_one("Ter") == "*"

    def test_three_letter_to_one_unknown(self):
        """Test unknown amino acid returns X."""
        assert _three_letter_to_one("Xyz") == "X"
        assert _three_letter_to_one("???") == "X"


class TestHgvsParsing:
    """Test HGVS notation parsing."""

    def test_parse_hgvs_missense(self):
        """Test parsing missense mutations."""
        assert _parse_hgvs_to_short("p.Asp614Gly") == "D614G"
        assert _parse_hgvs_to_short("p.Arg203Lys") == "R203K"
        assert _parse_hgvs_to_short("p.Asn501Tyr") == "N501Y"

    def test_parse_hgvs_synonymous(self):
        """Test parsing synonymous mutations."""
        assert _parse_hgvs_to_short("p.Phe924=") == "F924F"
        assert _parse_hgvs_to_short("p.Asp614=") == "D614D"

    def test_parse_hgvs_stop_codon(self):
        """Test parsing stop codon mutations."""
        assert _parse_hgvs_to_short("p.Gln493Ter") == "Q493*"

    def test_parse_hgvs_none(self):
        """Test parsing None input."""
        assert _parse_hgvs_to_short(None) is None

    def test_parse_hgvs_empty(self):
        """Test parsing empty string."""
        assert _parse_hgvs_to_short("") is None

    def test_parse_hgvs_invalid(self):
        """Test parsing invalid HGVS notation."""
        assert _parse_hgvs_to_short("invalid") is None
        assert _parse_hgvs_to_short("p.123") is None


class TestDiscoverInputFiles:
    """Test input file discovery."""

    def test_discover_multiple_files(self, input_dir_with_samples):
        """Test discovering multiple sample files."""
        files = discover_input_files(input_dir_with_samples)
        assert len(files) == 3
        sample_ids = [f[0] for f in files]
        assert "sample_01" in sample_ids
        assert "sample_02" in sample_ids
        assert "sample_03" in sample_ids

    def test_discover_single_file(self, input_dir_single_sample):
        """Test discovering a single sample file."""
        files = discover_input_files(input_dir_single_sample)
        assert len(files) == 1
        assert files[0][0] == "sample_01"

    def test_discover_empty_dir_raises(self, empty_input_dir):
        """Test that empty directory raises assertion error."""
        with pytest.raises(AssertionError):
            discover_input_files(empty_input_dir)


class TestPolarsExpressions:
    """Test Polars expression builders."""

    @pytest.fixture
    def sample_df(self):
        """Create a sample DataFrame for expression testing."""
        return pl.DataFrame(
            {
                "chrom": ["NC_045512.2", "NC_045512.2", "NC_045512.2"],
                "pos": [241, 3037, 23403],
                "ref": ["A", "C", "A"],
                "alt": ["T", "T", "G"],
                "af": [0.98, 0.45, 0.75],
                "gene": ["5UTR", "ORF1ab", "S"],
                "hgvs_p": ["", "p.Phe924Phe", "p.Asp614Gly"],
            }
        )

    def test_variant_id_expr(self, sample_df):
        """Test variant ID expression."""
        result = sample_df.with_columns(variant_id_expr())
        assert result["variant_id"][0] == "NC_045512.2:241:A>T"
        assert result["variant_id"][1] == "NC_045512.2:3037:C>T"
        assert result["variant_id"][2] == "NC_045512.2:23403:A>G"

    def test_mutation_type_expr_snp(self, sample_df):
        """Test mutation type expression for SNPs."""
        result = sample_df.with_columns(mutation_type_expr())
        assert all(result["mutation_type"] == "SNP")

    def test_mutation_type_expr_indels(self):
        """Test mutation type expression for indels."""
        df = pl.DataFrame(
            {
                "ref": ["A", "AT", "A"],
                "alt": ["AG", "A", "ATG"],
            }
        )
        result = df.with_columns(mutation_type_expr())
        assert result["mutation_type"][0] == "insertion"
        assert result["mutation_type"][1] == "deletion"
        assert result["mutation_type"][2] == "insertion"

    def test_mutation_type_expr_mnp(self):
        """Test mutation type expression for MNPs."""
        df = pl.DataFrame(
            {
                "ref": ["AT", "GGG"],
                "alt": ["GC", "AAA"],
            }
        )
        result = df.with_columns(mutation_type_expr())
        assert all(result["mutation_type"] == "MNP")

    def test_is_consensus_expr(self, sample_df):
        """Test consensus threshold expression."""
        result = sample_df.with_columns(is_consensus_expr(0.8))
        assert result["is_consensus"][0] is True  # 0.98 >= 0.8
        assert result["is_consensus"][1] is False  # 0.45 < 0.8
        assert result["is_consensus"][2] is False  # 0.75 < 0.8

    def test_is_consensus_expr_custom_threshold(self, sample_df):
        """Test consensus expression with custom threshold."""
        result = sample_df.with_columns(is_consensus_expr(0.5))
        assert result["is_consensus"][0] is True  # 0.98 >= 0.5
        assert result["is_consensus"][1] is False  # 0.45 < 0.5
        assert result["is_consensus"][2] is True  # 0.75 >= 0.5


class TestLoadAndConcat:
    """Test file loading and concatenation."""

    def test_load_and_concat_multiple_samples(self, input_dir_with_samples):
        """Test loading and concatenating multiple sample files."""
        files = discover_input_files(input_dir_with_samples)
        lf = load_and_concat(files)
        df = lf.collect()

        # 3 samples × 4 variants each = 12 rows
        assert len(df) == 12

        # Check sample_id column was added
        assert "sample_id" in df.columns
        assert df["sample_id"].n_unique() == 3

    def test_load_and_concat_single_sample(self, input_dir_single_sample):
        """Test loading a single sample file."""
        files = discover_input_files(input_dir_single_sample)
        lf = load_and_concat(files)
        df = lf.collect()

        assert len(df) == 4
        assert df["sample_id"].n_unique() == 1
        assert df["sample_id"][0] == "sample_01"


class TestDerivedColumns:
    """Test derived column generation."""

    def test_add_derived_columns(self, input_dir_single_sample):
        """Test adding derived columns."""
        files = discover_input_files(input_dir_single_sample)
        lf = load_and_concat(files)
        result = add_derived_columns(lf, DEFAULT_CONSENSUS_THRESHOLD).collect()

        assert "variant_id" in result.columns
        assert "mutation_type" in result.columns
        assert "is_consensus" in result.columns
        assert "aa_change" in result.columns


class TestCrossSampleMetrics:
    """Test cross-sample metric computation."""

    def test_sample_count_shared_variants(self, input_dir_with_samples):
        """Test sample count for shared variants."""
        files = discover_input_files(input_dir_with_samples)
        lf = load_and_concat(files)
        lf = add_derived_columns(lf, DEFAULT_CONSENSUS_THRESHOLD)
        result = add_cross_sample_metrics(lf).collect()

        # All variants appear in all 3 samples
        assert all(result["sample_count"] == 3)
        assert all(result["is_shared"])

    def test_sample_count_unique_variants(self, input_dir_single_sample):
        """Test sample count for unique variants."""
        files = discover_input_files(input_dir_single_sample)
        lf = load_and_concat(files)
        lf = add_derived_columns(lf, DEFAULT_CONSENSUS_THRESHOLD)
        result = add_cross_sample_metrics(lf).collect()

        # All variants appear in only 1 sample
        assert all(result["sample_count"] == 1)
        assert not any(result["is_shared"])


class TestFinalOutput:
    """Test final output generation."""

    def test_select_and_sort(self, input_dir_with_samples):
        """Test final column selection and sorting."""
        files = discover_input_files(input_dir_with_samples)
        result = (
            load_and_concat(files)
            .pipe(add_derived_columns, DEFAULT_CONSENSUS_THRESHOLD)
            .pipe(add_cross_sample_metrics)
            .pipe(select_and_sort)
            .collect()
        )

        assert list(result.columns) == FINAL_COLUMNS

        # Check sorting: should be by chrom, pos, sample_id
        positions = result["pos"].to_list()
        assert positions == sorted(positions)

    def test_output_files_created(self, input_dir_with_samples, tmp_path):
        """Test that output files are created."""
        output_basename = tmp_path / "variants"

        files = discover_input_files(input_dir_with_samples)
        result = (
            load_and_concat(files)
            .pipe(add_derived_columns, DEFAULT_CONSENSUS_THRESHOLD)
            .pipe(add_cross_sample_metrics)
            .pipe(select_and_sort)
            .collect()
        )

        # Write outputs
        result.write_csv(output_basename.with_suffix(".tsv"), separator="\t")
        result.write_parquet(output_basename.with_suffix(".parquet"))

        assert output_basename.with_suffix(".tsv").exists()
        assert output_basename.with_suffix(".parquet").exists()


class TestCommandLineArgs:
    """Test command line argument parsing."""

    def test_parse_required_args(self, input_dir_with_samples, tmp_path, monkeypatch):
        """Test parsing required arguments."""
        output_path = tmp_path / "output"
        test_args = [
            "collect_full_variant_table.py",
            "--input-dir",
            str(input_dir_with_samples),
            "--output",
            str(output_path),
        ]
        monkeypatch.setattr("sys.argv", test_args)

        args = parse_command_line_args()
        assert args.input_dir == input_dir_with_samples
        assert args.output == output_path
        assert args.consensus_threshold == DEFAULT_CONSENSUS_THRESHOLD

    def test_parse_custom_threshold(
        self, input_dir_with_samples, tmp_path, monkeypatch
    ):
        """Test parsing custom consensus threshold."""
        output_path = tmp_path / "output"
        test_args = [
            "collect_full_variant_table.py",
            "-i",
            str(input_dir_with_samples),
            "-o",
            str(output_path),
            "-t",
            "0.9",
        ]
        monkeypatch.setattr("sys.argv", test_args)

        args = parse_command_line_args()
        assert args.consensus_threshold == 0.9

    def test_missing_required_args(self, monkeypatch):
        """Test that missing required arguments raises SystemExit."""
        test_args = ["collect_full_variant_table.py"]
        monkeypatch.setattr("sys.argv", test_args)

        with pytest.raises(SystemExit):
            parse_command_line_args()


class TestMainFunction:
    """Test main function execution."""

    def test_main_execution(self, input_dir_with_samples, tmp_path, monkeypatch):
        """Test that main function executes successfully."""
        output_path = tmp_path / "variants"
        test_args = [
            "collect_full_variant_table.py",
            "--input-dir",
            str(input_dir_with_samples),
            "--output",
            str(output_path),
        ]
        monkeypatch.setattr("sys.argv", test_args)

        main()

        assert output_path.with_suffix(".tsv").exists()
        assert output_path.with_suffix(".parquet").exists()

        # Verify output content
        df = pl.read_csv(output_path.with_suffix(".tsv"), separator="\t")
        assert len(df) == 12  # 3 samples × 4 variants
        assert list(df.columns) == FINAL_COLUMNS

    def test_main_with_empty_dir(self, empty_input_dir, tmp_path, monkeypatch):
        """Test main function with empty input directory."""
        output_path = tmp_path / "variants"
        test_args = [
            "collect_full_variant_table.py",
            "--input-dir",
            str(empty_input_dir),
            "--output",
            str(output_path),
        ]
        monkeypatch.setattr("sys.argv", test_args)

        with pytest.raises(AssertionError):
            main()
