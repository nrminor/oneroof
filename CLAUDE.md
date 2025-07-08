# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

OneRoof is a Nextflow-based bioinformatics pipeline for base-calling, variant-calling, and consensus-calling of amplicon sequencing data. It supports both Nanopore (pod5/BAM/FASTQ) and Illumina (paired-end FASTQ) data, with particular focus on SARS-CoV-2 and H5N1 influenza genomic surveillance.

## Key Commands

### Development Environment Setup
```bash
# For environments with conda dependencies (full pipeline)
pixi install --frozen
pixi shell --frozen

# For PyPI-only environments (Python development)
uv venv
source .venv/bin/activate  # or .venv\Scripts\activate on Windows
uv sync
```

### Running the Pipeline
```bash
# Nanopore data from raw POD5s
nextflow run . \
  --pod5_dir my_pod5_dir \
  --primer_bed my_primers.bed \
  --refseq my_ref.fasta \
  --ref_gbk my_ref.gbk \
  --kit "SQK-NBD114-24"

# Illumina data
nextflow run . \
  --illumina_fastq_dir my_illumina_reads/ \
  --primer_bed my_primers.bed \
  --refseq my_ref.fasta \
  --ref_gbk my_ref.gbk

# Run without containers (requires pixi environment)
nextflow run . -profile containerless [options]
```

### Code Quality & Testing
```bash
# Python linting and formatting
ruff check . --exit-zero --fix --unsafe-fixes
ruff format .

# Run Python tests (using uv for speed)
uv run pytest bin/test_*.py
# Or run tests with tox for multiple environments
tox

# Build documentation
just docs

# IMPORTANT: Modifying README.md
# The README.md in the project root is generated from docs/index.qmd
# NEVER edit README.md directly - it will be overwritten
# Always edit docs/index.qmd and re-render:
just make-readme  # or: just docs

# Docker operations
just docker-build
just docker-push
```

## Architecture

### Directory Structure
- `main.nf` - Main workflow entry point that orchestrates platform-specific workflows
- `workflows/` - Platform-specific workflows (nanopore.nf, illumina.nf)
- `subworkflows/` - Reusable workflow components (alignment, variant_calling, primer_handling, etc.)
- `modules/` - Individual process definitions for tools (dorado, minimap2, ivar, etc.)
- `bin/` - Python utility scripts with PEP 723 inline dependencies (fully portable with uv)
- `conf/` - Configuration files for different platforms and tools

### Key Workflow Components

1. **Data Ingestion** - Handles multiple input formats (pod5, BAM, FASTQ) with optional remote file watching
2. **Primer Handling** - Validates primers, trims reads, and ensures complete amplicons
3. **Alignment & Variant Calling** - Platform-specific alignment and variant calling using minimap2 and ivar/bcftools
4. **Quality Control** - FastQC, MultiQC, and custom coverage plotting
5. **Consensus Generation** - Creates consensus sequences with configurable frequency thresholds
6. **Optional Features** - Metagenomics (Sylph), phylogenetics (Nextclade), haplotyping (Devider)

### Technology Stack
- **Workflow Engine**: Nextflow DSL2
- **Container Support**: Docker, Singularity/Apptainer
- **Environment Management**: Pixi (combines conda and PyPI dependencies), UV (fast Python package management)
- **Languages**: Nextflow (Groovy), Python 3.10+
- **Key Tools**: Dorado (basecalling), minimap2 (alignment), ivar/bcftools (variants), FastQC/MultiQC (QC)

### Configuration Philosophy
- Parameters are primarily set via command line arguments
- Platform-specific configs (nanopore.config, illumina.config) are auto-loaded based on input data type
- Container profiles (docker, singularity, apptainer, containerless) control execution environment
- Advanced users can modify nextflow.config for fine-tuning

### Important Parameters
- `--pod5_batch_size`: Controls GPU memory usage during basecalling
- `--min_variant_frequency`: Platform-specific defaults (0.05 for Illumina, 0.10 for Nanopore)
- `--downsample_to`: Manages computational resources by limiting coverage depth
- `--model`: Nanopore basecalling model (defaults to sup@latest)

## Dependency Management

### Python Package Management
- **Always use `uv` instead of `pip`** for any Python package installation - it's significantly faster and more reliable
- **Use `uv` for PyPI-only environments**: When working with Python scripts that only need PyPI dependencies
- **Use `pixi` for mixed environments**: When conda dependencies are required (e.g., for the full pipeline)
- **Script execution**: Always use `uv run` instead of `python3` to execute Python scripts
  ```bash
  # Good - uses inline dependencies from PEP 723 headers
  uv run bin/some_script.py

  # Avoid - doesn't guarantee dependencies
  python3 bin/some_script.py
  ```
- **Portable scripts**: All scripts in `bin/` include PEP 723 inline dependencies, making them fully portable with uv
- **Benefits**: This approach eliminates dependency hell in Python by ensuring consistent, reproducible environments

### Testing Infrastructure
- **Comprehensive test coverage**: Python scripts in `bin/` have extensive test coverage using pytest
- **Test execution**: Tests can be run quickly with UV for PyPI-only environments
  ```bash
  # Run all tests
  uv run pytest bin/test_*.py

  # Run specific test
  uv run pytest bin/test_specific_module.py
  ```
- **CI/CD**: The continuous integration pipeline uses UV instead of pip for improved speed and reliability
- **Test organization**: Test files follow the pattern `test_*.py` and are colocated with the scripts they test

## Development Notes

1. **Testing**: Python scripts have comprehensive test coverage; Nextflow workflow tests are planned for future implementation
2. **GPU Requirements**: Nanopore basecalling requires CUDA-capable GPUs
3. **Memory Management**: Use `--low_memory` flag for resource-constrained environments
4. **Slack Integration**: Optional alerts can be configured for pipeline completion
5. **Dependency Management**: Always use `uv` for Python operations to ensure fast, reliable dependency resolution
