# OneRoof Pipeline Utility Scripts

This directory contains utility scripts that support various operations in the OneRoof bioinformatics pipeline. These scripts handle tasks ranging from consensus sequence concatenation to primer validation, variant processing, and quality control visualization.

All Python scripts include PEP 723 inline dependencies and can be run directly with `uv run` without additional setup.

## Scripts Overview

### concat_consensus.py
**Purpose**: Concatenates consensus sequences from multiple FASTA files into a single multi-sample FASTA file.

**Key Features**:
- Automatically discovers all `.consensus.fasta` files in the current directory
- Concatenates sequences from each sample into a single sequence
- Outputs a unified `all_sample_consensus.fasta` file

**Usage Example**:
```bash
uv run concat_consensus.py
```

### file_watcher.py
**Purpose**: Monitors remote directories for new files and automatically transfers them via SSH.

**Key Features**:
- Watches remote paths for files matching specified patterns
- Supports credential-based authentication
- Validates file integrity using hash comparison
- Configurable watch duration

**Usage Example**:
```bash
uv run file_watcher.py --config credentials.yaml
```

### generate_variant_pivot.py
**Purpose**: Processes variant data from annotated VCF files into a pivot table format.

**Key Features**:
- Parses VCF fields extracted with SnpSift
- Creates pivot tables for variant analysis
- Uses Polars for efficient data processing

**Usage Example**:
```bash
uv run generate_variant_pivot.py --input_table variants.tsv
```

### ivar_variants_to_vcf.py
**Purpose**: Converts iVar TSV variant files to standard VCF format with advanced filtering and validation.

**Key Features**:
- Beautiful CLI interface with Typer and Rich
- Functional programming approach using Polars
- Pydantic validation for data integrity
- Fisher's exact test for strand bias filtering
- Support for both pass-only and all-variants output
- Configurable variant frequency thresholds

**Usage Example**:
```bash
uv run ivar_variants_to_vcf.py input.tsv output.vcf --ref reference.fasta --pass-only
```

### make_primer_patterns.py
**Purpose**: Generates primer pattern files from FASTA sequences of primers.

**Key Features**:
- Creates pattern files for forward and reverse primers
- Supports custom pattern prefixes/suffixes
- Handles primer sequences for amplicon-based sequencing

**Usage Example**:
```bash
uv run make_primer_patterns.py -i primers.fasta -o primer_patterns
```

### multisample_plot.py
**Purpose**: Creates multi-sample coverage plots from BED files.

**Key Features**:
- Visualizes coverage across multiple samples and chromosomes
- Log-scale transformation for better visualization
- Faceted plots by chromosome
- Customizable minimum coverage thresholds
- Sample lookup support for metadata integration

**Usage Example**:
```bash
uv run multisample_plot.py --input_dir bed_files/ --sample_lookup metadata.json --min_coverage 20
```

### plot_coverage.py
**Purpose**: Generates coverage plots from Mosdepth output files.

**Key Features**:
- Single-sample coverage visualization
- Chromosome-level faceting
- Customizable plot labels
- Highlights low-coverage regions

**Usage Example**:
```bash
uv run plot_coverage.py --input sample.mosdepth.global.dist.txt --label MySample
```

### resplice_primers.py
**Purpose**: Finds all possible combinations of spike-in primers in an amplicon scheme.

**Key Features**:
- Identifies primer combinations for complex amplicon schemes
- Handles spike-in primer variants
- Configurable primer naming conventions
- Outputs all valid amplicon combinations

**Usage Example**:
```bash
uv run resplice_primers.py -i primers_with_spikein.bed -o respliced_amplicons
```

### slack_alerts.py
**Purpose**: Sends Slack notifications about sample coverage statistics.

**Key Features**:
- Analyzes depth of coverage TSV files
- Reports samples meeting coverage thresholds
- Integrates with Slack API for notifications
- Customizable coverage depth requirements

**Usage Example**:
```bash
uv run slack_alerts.py --input_tsv_dir coverage_data/ --depth 20 --run_label "Run_001"
```

### split_primer_combos.py
**Purpose**: Splits a BED file containing multiple primer combinations into separate files for each combination.

**Key Features**:
- Separates complex primer schemes into individual combinations
- Maintains primer pair relationships
- Configurable primer naming suffixes

**Usage Example**:
```bash
uv run split_primer_combos.py -i combined_primers.bed -f _LEFT -r _RIGHT
```

### validate_primer_bed.py
**Purpose**: Validates and normalizes BED files containing primer coordinates.

**Key Features**:
- Checks for correct primer naming conventions
- Validates coordinate orientation
- Ensures proper primer pairing
- Outputs normalized BED file

**Usage Example**:
```bash
uv run validate_primer_bed.py -i primers.bed -o validated_primers
```

## Test Files

The directory includes test files for several scripts:
- `test_concat_consensus.py`
- `test_file_watcher.py`
- `test_generate_variant_pivot.py`
- `test_ivar_variants_to_vcf.py`
- `test_slack_alerts.py`

To run tests, use pytest:
```bash
pytest test_*.py
```

## Other Files

- `__init__.py` and `__main__.py`: Python package initialization files
- `resplice_primers.rs`: Rust implementation of primer resplicing (alternative to Python version)

## Dependencies

All scripts use PEP 723 inline script dependencies, which means you can run them directly with `uv run` without manual dependency installation. Common dependencies include:
- `polars` or `polars-lts-cpu`: High-performance dataframe library
- `biopython`: Biological sequence handling
- `loguru`: Advanced logging
- `plotnine`: Grammar of graphics plotting
- `typer` and `rich`: Beautiful CLI interfaces
- `pydantic`: Data validation

## Notes

- Most scripts include comprehensive help messages accessible via `--help`
- Scripts follow functional programming patterns where appropriate
- Error handling and logging are implemented throughout
- File paths should typically be absolute paths when provided as arguments
