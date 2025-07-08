#!/usr/bin/env python3
"""
Comprehensive test suite for resplice_primers.py

Tests the resplicing of spike-in primers into all possible amplicon combinations,
focusing on biological correctness and edge cases.
"""

import pytest
import polars as pl
import sys
import os

# Add the bin directory to the path so we can import the module
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "bin"))

from resplice_primers import (
    check_bed_existence,
    check_idx_delims,
    partition_by_amplicon,
    assign_new_indices,
    normalize_indices,
    resolve_primer_names,
    resplice_primers,
    finalize_primer_pairings,
    _convertable_to_int,
)


class TestResplicePrimers:
    """Test suite for primer resplicing functionality"""

    @pytest.fixture
    def basic_bed_df(self):
        """Create a basic BED DataFrame with standard primers"""
        data = {
            "Ref": ["chr1", "chr1", "chr1", "chr1"],
            "Start Position": [100, 500, 1000, 1400],
            "Stop Position": [120, 520, 1020, 1420],
            "NAME": [
                "amplicon1_LEFT",
                "amplicon1_RIGHT",
                "amplicon2_LEFT",
                "amplicon2_RIGHT",
            ],
            "INDEX": [60, 60, 60, 60],
            "SENSE": ["+", "-", "+", "-"],
        }
        return pl.DataFrame(data)

    @pytest.fixture
    def spike_in_bed_df(self):
        """Create a BED DataFrame with spike-in primers"""
        data = {
            "Ref": ["chr1"] * 6,
            "Start Position": [100, 110, 500, 510, 1000, 1400],
            "Stop Position": [120, 130, 520, 530, 1020, 1420],
            "NAME": [
                "amplicon1_LEFT-1",
                "amplicon1_LEFT-2",
                "amplicon1_RIGHT-1",
                "amplicon1_RIGHT-2",
                "amplicon2_LEFT",
                "amplicon2_RIGHT",
            ],
            "INDEX": [60] * 6,
            "SENSE": ["+", "+", "-", "-", "+", "-"],
        }
        return pl.DataFrame(data)

    @pytest.fixture
    def complex_spike_bed_df(self):
        """Create a complex BED DataFrame with multiple spike-ins"""
        data = {
            "Ref": ["chr1"] * 9,
            "Start Position": [100, 110, 120, 500, 510, 520, 1000, 1010, 1400],
            "Stop Position": [120, 130, 140, 520, 530, 540, 1020, 1030, 1420],
            "NAME": [
                "amplicon1_LEFT-1",
                "amplicon1_LEFT-2",
                "amplicon1_LEFT-3",
                "amplicon1_RIGHT-1",
                "amplicon1_RIGHT-2",
                "amplicon1_RIGHT-3",
                "amplicon2_LEFT-1",
                "amplicon2_LEFT-2",
                "amplicon2_RIGHT",
            ],
            "INDEX": [60] * 9,
            "SENSE": ["+", "+", "+", "-", "-", "-", "+", "+", "-"],
        }
        return pl.DataFrame(data)

    def test_check_bed_existence_valid(self, tmp_path):
        """Test bed existence check with valid file"""
        bed_file = tmp_path / "test.bed"
        bed_file.write_text("chr1\t100\t120\tprimer1\t60\t+\n")

        # Should not raise an exception
        check_bed_existence(str(bed_file))

    def test_check_bed_existence_invalid(self, tmp_path):
        """Test bed existence check with non-existent file"""
        with pytest.raises(SystemExit):
            check_bed_existence(str(tmp_path / "nonexistent.bed"))

    def test_check_idx_delims_no_warning(self, basic_bed_df, caplog):
        """Test that no warning is issued for primers without multiple delimiters"""
        check_idx_delims(basic_bed_df, "-", -1)
        assert "primer names contained more than one" not in caplog.text

    def test_check_idx_delims_with_warning(self, caplog):
        """Test warning for primers with multiple delimiter symbols"""
        df = pl.DataFrame(
            {
                "NAME": ["amplicon1-spike-1_LEFT", "amplicon2_RIGHT"],
                "Ref": ["chr1", "chr1"],
                "Start Position": [100, 500],
                "Stop Position": [120, 520],
                "INDEX": [60, 60],
                "SENSE": ["+", "-"],
            }
        )

        check_idx_delims(df, "-", -1)
        assert "primer names contained more than one" in caplog.text

    def test_partition_by_amplicon(self, spike_in_bed_df):
        """Test partitioning primers by amplicon"""
        partitions = partition_by_amplicon(spike_in_bed_df, "_LEFT", "_RIGHT")

        assert len(partitions) == 2  # Two amplicons

        # Check that primers are grouped correctly
        amplicon_names = set()
        for partition in partitions:
            amplicon = partition.select("Amplicon").unique().item()
            amplicon_names.add(amplicon)

        assert "amplicon1" in amplicon_names
        assert "amplicon2" in amplicon_names

    def test_assign_new_indices(self):
        """Test assignment of new indices to primers"""
        df = pl.DataFrame(
            {
                "Ref": ["chr1", "chr1"],
                "Start Position": [100, 110],
                "Stop Position": [120, 130],
                "ORIG_NAME": ["amplicon1_LEFT-1", "amplicon1_LEFT-2"],
                "NAME": ["amplicon1_LEFT", "amplicon1_LEFT"],
                "INDEX": [60, 60],
                "SENSE": ["+", "+"],
                "Amplicon": ["amplicon1", "amplicon1"],
            }
        )

        result = assign_new_indices(df, "_LEFT", "_RIGHT", "-", -1)

        names = result.select("NAME").to_series().to_list()
        assert names[0] == "amplicon1_LEFT-1"
        assert names[1] == "amplicon1_LEFT-2"

    def test_normalize_indices(self, spike_in_bed_df):
        """Test normalization of indices across amplicons"""
        partitions = partition_by_amplicon(spike_in_bed_df, "_LEFT", "_RIGHT")
        normalized = normalize_indices(partitions, "_LEFT", "_RIGHT")

        # Check that all amplicons are present
        assert len(normalized) == 2
        assert "amplicon1" in normalized
        assert "amplicon2" in normalized

    def test_convertable_to_int(self):
        """Test integer conversion checking"""
        assert _convertable_to_int("1") is True
        assert _convertable_to_int("123") is True
        assert _convertable_to_int("-1") is True
        assert _convertable_to_int("abc") is False
        assert _convertable_to_int("1.5") is False
        assert _convertable_to_int(None) is False

    def test_resolve_primer_names_basic(self):
        """Test basic primer name resolution"""
        fwd_primers = ["amplicon1_LEFT-1", "amplicon1_LEFT-2"]
        rev_primers = ["amplicon1_RIGHT-1", "amplicon1_RIGHT-2"]

        old_names, new_names = resolve_primer_names(fwd_primers, rev_primers)

        # Should create 4 combinations
        assert len(old_names) == 8  # 4 forward + 4 reverse
        assert len(new_names) == 8

        # Check that splice indices are added
        assert all("_splice" in name for name in new_names)

    def test_resolve_primer_names_asymmetric(self):
        """Test primer resolution with different numbers of forward/reverse primers"""
        fwd_primers = ["amplicon1_LEFT-1", "amplicon1_LEFT-2", "amplicon1_LEFT-3"]
        rev_primers = ["amplicon1_RIGHT-1", "amplicon1_RIGHT-2"]

        old_names, new_names = resolve_primer_names(fwd_primers, rev_primers)

        # Should create 6 combinations (3x2)
        assert len(old_names) == 12  # 6 forward + 6 reverse
        assert len(new_names) == 12

    def test_resolve_primer_names_invalid_index(self):
        """Test error handling for invalid primer indices"""
        fwd_primers = ["amplicon1_LEFT-abc"]  # Invalid index
        rev_primers = ["amplicon1_RIGHT-1"]

        with pytest.raises(SystemExit):
            resolve_primer_names(fwd_primers, rev_primers)

    def test_resplice_primers_standard_pairs(self, basic_bed_df):
        """Test resplicing with standard primer pairs (no spike-ins)"""
        partitions = partition_by_amplicon(basic_bed_df, "_LEFT", "_RIGHT")
        normalized = normalize_indices(partitions, "_LEFT", "_RIGHT")

        result = resplice_primers(normalized, "_LEFT", "_RIGHT")

        # Should return same number of amplicons
        assert len(result) == 2

        # Each amplicon should still have 2 primers
        for df in result:
            assert df.shape[0] == 2

    def test_resplice_primers_with_spikes(self, spike_in_bed_df):
        """Test resplicing with spike-in primers"""
        partitions = partition_by_amplicon(spike_in_bed_df, "_LEFT", "_RIGHT")
        normalized = normalize_indices(partitions, "_LEFT", "_RIGHT")

        result = resplice_primers(normalized, "_LEFT", "_RIGHT")

        # Check that amplicon1 has been expanded to all combinations
        amplicon1_found = False
        for df in result:
            if "amplicon1" in df.select("Amplicon").unique().item():
                amplicon1_found = True
                # 2 forward x 2 reverse = 4 combinations
                assert df.shape[0] == 8  # 4 forward + 4 reverse primers

        assert amplicon1_found

    def test_resplice_primers_single_primer_error(self):
        """Test error handling for amplicon with single primer"""
        df = pl.DataFrame(
            {
                "Ref": ["chr1"],
                "Start Position": [100],
                "Stop Position": [120],
                "NAME": ["amplicon1_LEFT"],
                "INDEX": [60],
                "SENSE": ["+"],
                "Amplicon": ["amplicon1"],
            }
        )

        amplicon_dfs = {"amplicon1": df}

        with pytest.raises(ValueError):
            resplice_primers(amplicon_dfs, "_LEFT", "_RIGHT")

    def test_finalize_primer_pairings(self):
        """Test finalization of primer pairings"""
        # Create a valid respliced dataframe
        df1 = pl.DataFrame(
            {
                "Ref": ["chr1", "chr1"],
                "Start Position": [100, 500],
                "Stop Position": [120, 520],
                "NAME": ["amplicon1_LEFT_splice1", "amplicon1_RIGHT_splice1"],
                "ORIG_NAME": ["amplicon1_LEFT-1", "amplicon1_RIGHT-1"],
                "INDEX": [60, 60],
                "SENSE": ["+", "-"],
                "Amplicon": ["amplicon1", "amplicon1"],
            }
        )

        result = finalize_primer_pairings([df1], "_LEFT", "_RIGHT")

        assert result.shape[0] == 2
        assert "Amplicon" not in result.columns
        assert "ORIG_NAME" not in result.columns

    def test_finalize_primer_pairings_missing_reverse(self, caplog):
        """Test handling of amplicon missing reverse primer"""
        df = pl.DataFrame(
            {
                "Ref": ["chr1", "chr1"],
                "Start Position": [100, 110],
                "Stop Position": [120, 130],
                "NAME": ["amplicon1_LEFT_splice1", "amplicon1_LEFT_splice2"],
                "ORIG_NAME": ["amplicon1_LEFT-1", "amplicon1_LEFT-2"],
                "INDEX": [60, 60],
                "SENSE": ["+", "+"],
                "Amplicon": ["amplicon1", "amplicon1"],
            }
        )

        result = finalize_primer_pairings([df], "_LEFT", "_RIGHT")

        # Should return empty dataframe as no valid pairs
        assert result.shape[0] == 0
        assert "Incorrect splicing occurred" in caplog.text

    def test_complex_spike_in_scenario(self, complex_spike_bed_df):
        """Test complex scenario with multiple spike-ins per primer position"""
        partitions = partition_by_amplicon(complex_spike_bed_df, "_LEFT", "_RIGHT")
        normalized = normalize_indices(partitions, "_LEFT", "_RIGHT")
        result = resplice_primers(normalized, "_LEFT", "_RIGHT")

        # Find amplicon1 result
        for df in result:
            if "amplicon1" in df.select("Amplicon").unique().item():
                # 3 forward x 3 reverse = 9 combinations
                assert df.shape[0] == 18  # 9 forward + 9 reverse primers

                # Check that all have splice indices
                names = df.select("NAME").to_series().to_list()
                assert all("_splice" in name for name in names)

    def test_custom_delimiter(self):
        """Test with custom delimiter symbol"""
        df = pl.DataFrame(
            {
                "Ref": ["chr1", "chr1"],
                "Start Position": [100, 110],
                "Stop Position": [120, 130],
                "NAME": ["amplicon1_LEFT_1", "amplicon1_LEFT_2"],
                "INDEX": [60, 60],
                "SENSE": ["+", "+"],
            }
        )

        check_idx_delims(df, "_", -1)  # Using underscore as delimiter

    def test_edge_case_empty_dataframe(self):
        """Test handling of empty dataframe"""
        df = pl.DataFrame(
            {
                "Ref": [],
                "Start Position": [],
                "Stop Position": [],
                "NAME": [],
                "INDEX": [],
                "SENSE": [],
            }
        )

        partitions = partition_by_amplicon(df, "_LEFT", "_RIGHT")
        assert len(partitions) == 0

    def test_primer_name_with_special_chars(self):
        """Test handling of primer names with special characters"""
        df = pl.DataFrame(
            {
                "Ref": ["chr1", "chr1"],
                "Start Position": [100, 500],
                "Stop Position": [120, 520],
                "NAME": ["amplicon1.v2_LEFT-1", "amplicon1.v2_RIGHT-1"],
                "INDEX": [60, 60],
                "SENSE": ["+", "-"],
            }
        )

        partitions = partition_by_amplicon(df, "_LEFT", "_RIGHT")
        assert len(partitions) == 1
