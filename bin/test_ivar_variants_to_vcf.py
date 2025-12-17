#!/usr/bin/env python3
"""Tests for ivar_variants_to_vcf.py module."""

import gzip
import tempfile
from pathlib import Path

import polars as pl
import pytest
from typer.testing import CliRunner

from ivar_variants_to_vcf import (
    ConversionConfig,
    FilterType,
    IvarVariant,
    VariantType,
    VcfVariant,
    app,
    calculate_strand_bias_pvalue,
    create_filter_expr,
    create_sample_info_expr,
    create_strand_bias_expr,
    determine_variant_type_expr,
    generate_vcf_header,
    transform_ivar_to_vcf,
    transform_ref_alt_expr,
)

runner = CliRunner()


@pytest.fixture
def valid_ivar_tsv():
    """Create a valid iVar TSV file for testing."""
    data = """REGION\tPOS\tREF\tALT\tREF_DP\tREF_RV\tREF_QUAL\tALT_DP\tALT_RV\tALT_QUAL\tALT_FREQ\tTOTAL_DP\tPVAL\tPASS\tGFF_FEATURE\tREF_CODON\tREF_AA\tALT_CODON\tALT_AA
NC_045512.2\t241\tC\tT\t3\t1\t37.0\t997\t499\t37.0\t0.997\t1000\t0.0\tTRUE\tgene-ORF1ab\tAAC\tN\tAAT\tN
NC_045512.2\t3037\tC\tT\t50\t25\t35.0\t950\t475\t36.0\t0.95\t1000\t0.0\tTRUE\tgene-ORF1ab\tTTC\tF\tTTT\tF
NC_045512.2\t10000\tA\t+ACGT\t100\t50\t38.0\t900\t450\t38.0\t0.9\t1000\t0.0\tTRUE\tgene-ORF1ab\tAAA\tK\tAAAA\tK
NC_045512.2\t15000\tACGT\t-CGT\t200\t100\t36.0\t800\t400\t37.0\t0.8\t1000\t0.0\tTRUE\tgene-ORF1ab\tACGT\tT\tA\t-
NC_045512.2\t20000\tG\tA\t500\t250\t30.0\t500\t250\t15.0\t0.5\t1000\t0.05\tFALSE\tgene-S\tGAT\tD\tAAT\tN
"""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
        f.write(data)
        return Path(f.name)


@pytest.fixture
def empty_ivar_tsv():
    """Create an empty iVar TSV file for testing."""
    data = """REGION\tPOS\tREF\tALT\tREF_DP\tREF_RV\tREF_QUAL\tALT_DP\tALT_RV\tALT_QUAL\tALT_FREQ\tTOTAL_DP\tPVAL\tPASS\tGFF_FEATURE\tREF_CODON\tREF_AA\tALT_CODON\tALT_AA
"""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
        f.write(data)
        return Path(f.name)


@pytest.fixture
def malformed_ivar_tsv():
    """Create a malformed iVar TSV file for testing."""
    data = """REGION\tPOS\tREF
NC_045512.2\t241\tC\tmissing_columns
"""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
        f.write(data)
        return Path(f.name)


@pytest.fixture
def strand_bias_ivar_tsv():
    """Create an iVar TSV file with strand bias for testing."""
    data = """REGION\tPOS\tREF\tALT\tREF_DP\tREF_RV\tREF_QUAL\tALT_DP\tALT_RV\tALT_QUAL\tALT_FREQ\tTOTAL_DP\tPVAL\tPASS\tGFF_FEATURE\tREF_CODON\tREF_AA\tALT_CODON\tALT_AA
NC_045512.2\t1000\tA\tG\t100\t90\t37.0\t900\t50\t37.0\t0.9\t1000\t0.0\tTRUE\tgene-ORF1ab\tAAA\tK\tGAA\tE
NC_045512.2\t2000\tC\tT\t200\t10\t36.0\t800\t790\t36.0\t0.8\t1000\t0.0\tTRUE\tgene-ORF1ab\tCCC\tP\tTCC\tS
"""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
        f.write(data)
        return Path(f.name)


@pytest.fixture
def reference_fasta():
    """Create a reference FASTA file for testing."""
    data = """>NC_045512.2 Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1
ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCTGTTCTCTAAA
CGAACTTTAAAATCTGTGTGGCTGTCACTCGGCTGCATGCTTAGTGCACTCACGCAGTATAATTAATAAC
"""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
        f.write(data)
        return Path(f.name)


@pytest.fixture
def output_vcf_path():
    """Create a temporary output VCF path."""
    return Path(tempfile.mktemp(suffix=".vcf"))


@pytest.fixture
def output_vcf_gz_path():
    """Create a temporary output gzipped VCF path."""
    return Path(tempfile.mktemp(suffix=".vcf.gz"))


class TestPydanticModels:
    """Test Pydantic model validation."""

    def test_ivar_variant_valid(self):
        """Test valid IvarVariant creation."""
        variant = IvarVariant(
            region="NC_045512.2",
            pos=241,
            ref="C",
            alt="T",
            ref_dp=10,
            ref_rv=5,
            ref_qual=37.0,
            alt_dp=990,
            alt_rv=495,
            alt_qual=37.0,
            alt_freq=0.99,
            total_dp=1000,
            pval=0.0,
            pass_test=True,
        )
        assert variant.pos == 241
        assert variant.alt_freq == 0.99

    def test_ivar_variant_invalid_position(self):
        """Test IvarVariant with invalid position."""
        with pytest.raises(ValueError):
            IvarVariant(
                region="NC_045512.2",
                pos=0,  # Invalid: must be > 0
                ref="C",
                alt="T",
                ref_dp=10,
                ref_rv=5,
                ref_qual=37.0,
                alt_dp=990,
                alt_rv=495,
                alt_qual=37.0,
                alt_freq=0.99,
                total_dp=1000,
                pval=0.0,
                pass_test=True,
            )

    def test_ivar_variant_invalid_freq(self):
        """Test IvarVariant with invalid frequency."""
        with pytest.raises(ValueError):
            IvarVariant(
                region="NC_045512.2",
                pos=241,
                ref="C",
                alt="T",
                ref_dp=10,
                ref_rv=5,
                ref_qual=37.0,
                alt_dp=990,
                alt_rv=495,
                alt_qual=37.0,
                alt_freq=1.5,  # Invalid: must be <= 1
                total_dp=1000,
                pval=0.0,
                pass_test=True,
            )

    def test_vcf_variant_valid(self):
        """Test valid VcfVariant creation."""
        variant = VcfVariant(
            chrom="NC_045512.2",
            pos=241,
            ref="C",
            alt="T",
            filter="PASS",
            info="TYPE=SNP",
            format="GT:DP",
            sample="1:1000",
        )
        assert variant.pos == 241
        assert variant.filter == "PASS"

    def test_vcf_variant_invalid_allele(self):
        """Test VcfVariant with invalid allele."""
        with pytest.raises(ValueError):
            VcfVariant(
                chrom="NC_045512.2",
                pos=241,
                ref="X",  # Invalid nucleotide
                alt="T",
                filter="PASS",
                info="TYPE=SNP",
                format="GT:DP",
                sample="1:1000",
            )

    def test_conversion_config_valid(self, valid_ivar_tsv, output_vcf_path):
        """Test valid ConversionConfig creation."""
        config = ConversionConfig(
            file_in=valid_ivar_tsv,
            file_out=output_vcf_path,
            freq_threshold=0.05,
            bad_qual_threshold=20.0,
        )
        assert config.file_in == valid_ivar_tsv
        assert config.freq_threshold == 0.05

    def test_conversion_config_invalid_threshold(self, valid_ivar_tsv, output_vcf_path):
        """Test ConversionConfig with invalid threshold."""
        with pytest.raises(ValueError):
            ConversionConfig(
                file_in=valid_ivar_tsv,
                file_out=output_vcf_path,
                freq_threshold=1.5,  # Invalid: must be <= 1
            )


class TestPureFunctions:
    """Test pure transformation functions."""

    def test_calculate_strand_bias_pvalue(self):
        """Test strand bias p-value calculation."""
        # No strand bias
        pval1 = calculate_strand_bias_pvalue(100, 50, 900, 450)
        assert pval1 > 0.05

        # Significant strand bias
        pval2 = calculate_strand_bias_pvalue(100, 90, 900, 50)
        assert pval2 < 0.05

    def test_determine_variant_type_expr(self):
        """Test variant type determination."""
        df = pl.DataFrame({"ALT": ["T", "+ACGT", "-CGT", "A"]})

        result = df.select(determine_variant_type_expr().alias("TYPE"))
        assert result["TYPE"][0] == VariantType.SNP.value
        assert result["TYPE"][1] == VariantType.INS.value
        assert result["TYPE"][2] == VariantType.DEL.value
        assert result["TYPE"][3] == VariantType.SNP.value

    def test_transform_ref_alt_expr(self):
        """Test REF/ALT transformation for indels."""
        df = pl.DataFrame(
            {"REF": ["A", "A", "ACGT", "G"], "ALT": ["T", "+CGTC", "-CGT", "C"]}
        )

        ref_expr, alt_expr = transform_ref_alt_expr()
        result = df.select([ref_expr.alias("NEW_REF"), alt_expr.alias("NEW_ALT")])

        # SNP: unchanged
        assert result["NEW_REF"][0] == "A"
        assert result["NEW_ALT"][0] == "T"

        # Insertion: REF unchanged, ALT = REF + inserted sequence
        assert result["NEW_REF"][1] == "A"
        assert result["NEW_ALT"][1] == "ACGTC"

        # Deletion: REF = REF + deleted sequence, ALT = REF
        assert result["NEW_REF"][2] == "ACGTCGT"
        assert result["NEW_ALT"][2] == "ACGT"

    def test_create_sample_info_expr(self):
        """Test sample information expression."""
        df = pl.DataFrame(
            {
                "TOTAL_DP": [1000],
                "REF_DP": [10],
                "REF_RV": [5],
                "REF_QUAL": [37.0],
                "ALT_DP": [990],
                "ALT_RV": [495],
                "ALT_QUAL": [37.0],
                "ALT_FREQ": [0.99],
            }
        )

        result = df.select(create_sample_info_expr().alias("SAMPLE"))
        expected = "1:1000:10:5:37.0:990:495:37.0:0.99"
        assert result["SAMPLE"][0] == expected


class TestTransformations:
    """Test data transformation functions."""

    def test_transform_ivar_to_vcf_basic(self, valid_ivar_tsv, output_vcf_path):
        """Test basic iVar to VCF transformation."""
        config = ConversionConfig(
            file_in=valid_ivar_tsv,
            file_out=output_vcf_path,
        )

        ivar_lf = pl.scan_csv(str(valid_ivar_tsv), separator="\t")
        vcf_lf = transform_ivar_to_vcf(ivar_lf, config)
        vcf_df = vcf_lf.collect()

        assert len(vcf_df) == 5
        assert "CHROM" in vcf_df.columns
        assert "POS" in vcf_df.columns
        assert "REF" in vcf_df.columns
        assert "ALT" in vcf_df.columns
        assert "FILTER" in vcf_df.columns
        assert "INFO" in vcf_df.columns
        assert "FORMAT" in vcf_df.columns
        assert "SAMPLE" in vcf_df.columns

    def test_transform_with_filters(self, valid_ivar_tsv, output_vcf_path):
        """Test transformation with filters applied."""
        config = ConversionConfig(
            file_in=valid_ivar_tsv,
            file_out=output_vcf_path,
            pass_only=True,
            freq_threshold=0.9,
        )

        ivar_lf = pl.scan_csv(str(valid_ivar_tsv), separator="\t")
        vcf_lf = transform_ivar_to_vcf(ivar_lf, config)
        vcf_df = vcf_lf.collect()

        # Should filter out the variant with PASS=FALSE and those below 0.9 freq
        assert len(vcf_df) == 3

    def test_indel_transformation(self, valid_ivar_tsv, output_vcf_path):
        """Test correct transformation of insertions and deletions."""
        config = ConversionConfig(
            file_in=valid_ivar_tsv,
            file_out=output_vcf_path,
        )

        ivar_lf = pl.scan_csv(str(valid_ivar_tsv), separator="\t")
        vcf_lf = transform_ivar_to_vcf(ivar_lf, config)
        vcf_df = vcf_lf.collect()

        # Check insertion (position 10000)
        ins_row = vcf_df.filter(pl.col("POS") == 10000).row(0, named=True)
        assert ins_row["REF"] == "A"
        assert ins_row["ALT"] == "AACGT"
        assert "TYPE=INS" in ins_row["INFO"]

        # Check deletion (position 15000)
        del_row = vcf_df.filter(pl.col("POS") == 15000).row(0, named=True)
        assert del_row["REF"] == "ACGTCGT"
        assert del_row["ALT"] == "ACGT"
        assert "TYPE=DEL" in del_row["INFO"]


class TestFilterExpressions:
    """Test filter expression creation."""

    def test_create_filter_expr_all_pass(self, output_vcf_path):
        """Test filter expression when all variants pass."""
        config = ConversionConfig(
            file_in=Path("dummy.tsv"),
            file_out=output_vcf_path,
            bad_qual_threshold=10.0,
            ignore_strand_bias=True,
        )

        df = pl.DataFrame(
            {
                "PASS": [True, True],
                "ALT_QUAL": [30.0, 40.0],
            }
        )

        filter_expr = create_filter_expr(config)
        result = df.select(filter_expr.alias("FILTER"))

        assert all(result["FILTER"] == FilterType.PASS.value)

    def test_create_filter_expr_quality_fail(self, output_vcf_path):
        """Test filter expression with quality failures."""
        config = ConversionConfig(
            file_in=Path("dummy.tsv"),
            file_out=output_vcf_path,
            bad_qual_threshold=30.0,
            ignore_strand_bias=True,
        )

        df = pl.DataFrame(
            {
                "PASS": [True, True],
                "ALT_QUAL": [20.0, 40.0],
            }
        )

        filter_expr = create_filter_expr(config)
        result = df.select(filter_expr.alias("FILTER"))

        assert result["FILTER"][0] == FilterType.BAD_QUALITY.value
        assert result["FILTER"][1] == FilterType.PASS.value


class TestVCFGeneration:
    """Test VCF file generation."""

    def test_generate_vcf_header_basic(self, output_vcf_path):
        """Test basic VCF header generation."""
        config = ConversionConfig(
            file_in=Path("dummy.tsv"),
            file_out=output_vcf_path,
        )

        headers = generate_vcf_header(config)

        assert "##fileformat=VCFv4.2" in headers
        assert "##source=iVar" in headers
        assert any("##INFO=<ID=TYPE" in h for h in headers)
        assert any("##FILTER=<ID=PASS" in h for h in headers)
        assert any("##FORMAT=<ID=GT" in h for h in headers)

    def test_generate_vcf_header_with_reference(self, output_vcf_path, reference_fasta):
        """Test VCF header generation with reference."""
        config = ConversionConfig(
            file_in=Path("dummy.tsv"),
            file_out=output_vcf_path,
            ref_fasta=reference_fasta,
        )

        headers = generate_vcf_header(config)

        # Should include contig information
        assert any("##contig=<ID=NC_045512.2" in h for h in headers)


class TestCLI:
    """Test command-line interface."""

    def test_convert_basic(self, valid_ivar_tsv, output_vcf_path):
        """Test basic convert command."""
        result = runner.invoke(
            app, ["convert", str(valid_ivar_tsv), str(output_vcf_path)]
        )

        assert result.exit_code == 0
        assert output_vcf_path.exists()

        # Check that all_hap file is also created
        all_hap_path = output_vcf_path.parent / f"{output_vcf_path.stem}_all_hap.vcf"
        assert all_hap_path.exists()

    def test_convert_with_options(self, valid_ivar_tsv, output_vcf_path):
        """Test convert command with options."""
        result = runner.invoke(
            app,
            [
                "convert",
                str(valid_ivar_tsv),
                str(output_vcf_path),
                "--pass-only",
                "--allele-freq",
                "0.9",
                "--min-quality",
                "30",
                "--ignore-strand-bias",
            ],
        )

        assert result.exit_code == 0
        assert output_vcf_path.exists()

    def test_convert_gzipped_output(self, valid_ivar_tsv, output_vcf_gz_path):
        """Test convert command with gzipped output."""
        result = runner.invoke(
            app, ["convert", str(valid_ivar_tsv), str(output_vcf_gz_path)]
        )

        assert result.exit_code == 0
        assert output_vcf_gz_path.exists()

        # Verify it's actually gzipped
        with gzip.open(output_vcf_gz_path, "rt") as f:
            first_line = f.readline()
            assert first_line.startswith("##fileformat=VCFv4")

    def test_convert_missing_input(self, output_vcf_path):
        """Test convert command with missing input file."""
        result = runner.invoke(
            app, ["convert", "non_existent.tsv", str(output_vcf_path)]
        )

        assert result.exit_code == 1

    def test_stats_command(self, valid_ivar_tsv):
        """Test stats command."""
        result = runner.invoke(app, ["stats", str(valid_ivar_tsv)])

        assert result.exit_code == 0
        assert "Total variants:" in result.output
        assert "Variant Types:" in result.output
        assert "SNPs:" in result.output
        assert "Insertions:" in result.output
        assert "Deletions:" in result.output

    def test_validate_command(self, valid_ivar_tsv, output_vcf_path):
        """Test validate command."""
        # First create a VCF file
        runner.invoke(app, ["convert", str(valid_ivar_tsv), str(output_vcf_path)])

        # Then validate it
        result = runner.invoke(app, ["validate", str(output_vcf_path)])

        assert result.exit_code == 0
        assert "File Statistics:" in result.output
        assert "All required VCF columns present" in result.output


class TestEdgeCases:
    """Test edge cases and error conditions."""

    def test_empty_file_handling(self, empty_ivar_tsv, output_vcf_path):
        """Test handling of empty input file."""
        result = runner.invoke(
            app, ["convert", str(empty_ivar_tsv), str(output_vcf_path)]
        )

        assert result.exit_code == 0
        assert output_vcf_path.exists()

        # Check VCF has headers but no data
        with open(output_vcf_path) as f:
            lines = f.readlines()
            header_lines = [l for l in lines if l.startswith("##")]
            assert len(header_lines) > 0

            # Should have column header line
            assert any(l.startswith("#CHROM") for l in lines)

    def test_malformed_file_handling(self, malformed_ivar_tsv, output_vcf_path):
        """Test handling of malformed input file."""
        result = runner.invoke(
            app, ["convert", str(malformed_ivar_tsv), str(output_vcf_path)]
        )

        assert result.exit_code == 1

    def test_strand_bias_detection(self, strand_bias_ivar_tsv, output_vcf_path):
        """Test strand bias detection."""
        ConversionConfig(
            file_in=strand_bias_ivar_tsv,
            file_out=output_vcf_path,
            ignore_strand_bias=False,
        )

        ivar_lf = pl.scan_csv(str(strand_bias_ivar_tsv), separator="\t")

        # Add strand bias detection columns
        df_with_bias = ivar_lf.with_columns(
            create_strand_bias_expr().alias("has_strand_bias")
        ).collect()

        # At least one variant should have strand bias
        assert any(df_with_bias["has_strand_bias"])

    def test_consecutive_variants_handling(self, output_vcf_path):
        """Test handling of consecutive variants."""
        # Create test data with consecutive positions
        data = """REGION\tPOS\tREF\tALT\tREF_DP\tREF_RV\tREF_QUAL\tALT_DP\tALT_RV\tALT_QUAL\tALT_FREQ\tTOTAL_DP\tPVAL\tPASS\tGFF_FEATURE\tREF_CODON\tREF_AA\tALT_CODON\tALT_AA
NC_045512.2\t1000\tA\tG\t100\t50\t37.0\t900\t450\t37.0\t0.9\t1000\t0.0\tTRUE\tgene-ORF1ab\tAAA\tK\tGAA\tE
NC_045512.2\t1001\tC\tT\t100\t50\t36.0\t900\t450\t36.0\t0.9\t1000\t0.0\tTRUE\tgene-ORF1ab\tCCC\tP\tTCC\tS
NC_045512.2\t1002\tG\tA\t100\t50\t38.0\t900\t450\t38.0\t0.9\t1000\t0.0\tTRUE\tgene-ORF1ab\tGGG\tG\tAGG\tR
"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
            f.write(data)
            consecutive_file = Path(f.name)

        try:
            result = runner.invoke(
                app,
                [
                    "convert",
                    str(consecutive_file),
                    str(output_vcf_path),
                    "--no-merge",  # Disable merging to test basic handling
                ],
            )

            assert result.exit_code == 0
            assert output_vcf_path.exists()
        finally:
            consecutive_file.unlink()


# Cleanup fixtures
@pytest.fixture(autouse=True)
def cleanup(request):
    """Clean up temporary files after tests."""

    def remove_temp_files():
        import glob
        import os

        # Clean up temporary TSV files
        for pattern in [
            "/tmp/tmp*.tsv",
            "/tmp/tmp*.vcf",
            "/tmp/tmp*.vcf.gz",
            "/tmp/tmp*.fasta",
        ]:
            for f in glob.glob(pattern):
                try:
                    os.remove(f)
                except:
                    pass

    request.addfinalizer(remove_temp_files)
