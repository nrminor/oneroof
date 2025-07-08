# OneRoof Test Suite

This directory contains the nf-test based testing framework for the OneRoof pipeline. The test suite is designed to validate pipeline functionality with minimal intrusion into the existing codebase.

## Overview

The test suite uses [nf-test](https://www.nf-test.com/), the official testing framework for Nextflow pipelines. Tests are organized hierarchically to match the pipeline structure:

```
tests/
├── data/                 # Test data files
├── modules/              # Tests for individual processes
├── subworkflows/         # Tests for sub-workflows (future)
├── workflows/            # Tests for platform-specific workflows
└── pipelines/            # End-to-end pipeline tests
```

## Installation

nf-test is included in the project's Pixi environment. To access it:

```bash
# Install all project dependencies including nf-test
pixi install --frozen

# Activate the Pixi environment
pixi shell --frozen

# Verify nf-test is available
nf-test version
```

No separate installation is required as nf-test is managed alongside all other project dependencies through Pixi.

## Running Tests

All test commands assume you're working within the Pixi environment:

```bash
# First, activate the Pixi environment if not already active
pixi shell --frozen
```

### Run all tests
```bash
just test
# Or directly: nf-test test
```

### Run specific test categories
```bash
# Run only module tests
just test-modules
# Or: nf-test test --tag modules

# Run only workflow tests
just test-workflows
# Or: nf-test test --tag workflows

# Run only pipeline tests
just test-pipeline
# Or: nf-test test --tag pipeline
```

### Run specific test files
```bash
# Test a specific module using just
just test-file tests/modules/minimap2_align.nf.test

# Or directly with nf-test
nf-test test tests/modules/minimap2_align.nf.test
```

### Update test snapshots
```bash
# Update all snapshots after intentional changes
just test-update
# Or: nf-test test --update-snapshot
```

### Clean test outputs
```bash
# Remove all test artifacts
just test-clean
```

## Test Structure

### Module Tests
Located in `tests/modules/`, these tests validate individual processes:
- `minimap2_align.nf.test` - Tests read alignment functionality
- `ivar_variants.nf.test` - Tests variant calling
- `ivar_consensus.nf.test` - Tests consensus sequence generation

### Workflow Tests
Located in `tests/workflows/`, these tests validate platform-specific workflows:
- `illumina.nf.test` - Tests the Illumina paired-end workflow
- `nanopore.nf.test` - Tests the Nanopore workflow (without basecalling)

### Pipeline Tests
Located in `tests/pipelines/`, these tests validate end-to-end functionality:
- `main.nf.test` - Tests the complete pipeline with various input combinations

## Test Data

The `tests/data/` directory contains minimal test datasets:
- Reference genome (1.5kb synthetic sequence)
- Primer BED files
- Small FASTQ files (10 reads each)
- Pre-aligned BAM files
- Sample metadata files

All test data is internally consistent and designed to exercise key pipeline features while maintaining small file sizes for fast test execution.

## CI/CD Integration

Tests are automatically run on:
- Push to main/dev branches
- Pull requests
- Manual workflow dispatch

The GitHub Actions workflow is defined in `.github/workflows/test.yml`.

## Writing New Tests

When adding new functionality to the pipeline:

1. Create test data if needed in `tests/data/`
2. Write a test file following the naming convention: `<module_name>.nf.test`
3. Include appropriate tags for test organization
4. Run the test locally before committing
5. Update snapshots if output changes are expected

Example test structure:
```groovy
nextflow_process {
    name "Test Process Name"
    script "../../../modules/process_name.nf"
    process "PROCESS_NAME"
    tag "modules"
    tag "process_name"

    test("basic functionality") {
        when {
            process {
                """
                input[0] = file("${projectDir}/tests/data/test_file.txt")
                """
            }
        }
        then {
            assert process.success
            assert path(process.out.output[0]).exists()
        }
    }
}
```

## Troubleshooting

### Common Issues

1. **Tests fail with "file not found"**
   - Ensure test data files exist in `tests/data/`
   - Check file paths use `${projectDir}` prefix

2. **Snapshot mismatches**
   - Review changes with `nf-test test --verbose`
   - Update snapshots if changes are intentional

3. **Resource errors**
   - The test profile limits resources; adjust in `nextflow.config` if needed

### Debug Mode

Run tests with increased verbosity:
```bash
nf-test test --verbose --debug
```

## Python Script Testing

### Overview

The `bin/` directory contains Python utility scripts that support the pipeline's data processing and analysis. These scripts have their own test suite using pytest and tox, independent of the Nextflow tests.

### Test Files

Python test files follow the naming convention `test_<script_name>.py` and are located in the `bin/` directory alongside the scripts they test:
- `test_concat_consensus.py` - Tests consensus sequence concatenation
- `test_file_watcher.py` - Tests remote file watching functionality
- `test_generate_variant_pivot.py` - Tests variant pivot table generation
- `test_ivar_variants_to_vcf.py` - Tests ivar to VCF conversion
- `test_slack_alerts.py` - Tests Slack notification functionality

### Setting Up the Test Environment

#### Quick Setup with UV (Recommended for Python-only testing)

For rapid Python testing without conda dependencies, use UV:

```bash
# Install UV if not already available
curl -LsSf https://astral.sh/uv/install.sh | sh

# Create a virtual environment and install dependencies
uv venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate
uv pip install -e .
uv pip install pytest pytest-cov pytest-mock pytest-xdist tox
```

#### Full Setup with Pixi

For the complete environment including all bioinformatics tools:

```bash
pixi install --frozen
pixi shell --frozen
```

### Running Python Tests

#### Direct pytest execution

```bash
# Run all Python tests
pytest bin/ -v

# Run with coverage report
pytest bin/ --cov=bin --cov-report=html --cov-report=term-missing

# Run a specific test file
pytest bin/test_ivar_variants_to_vcf.py -v

# Run tests in parallel
pytest bin/ -n auto -v

# Run only tests matching a pattern
pytest bin/ -k "variant" -v

# Run tests with specific markers
pytest bin/ -m "not slow" -v
```

#### Using tox for multi-version testing

```bash
# Run tests across all supported Python versions (3.10, 3.11, 3.12, 3.13)
tox

# Run tests for a specific Python version
tox -e py312

# Run linting checks
tox -e lint

# Run formatting
tox -e format

# Run specific tox environments in parallel
tox -p auto
```

#### Using just recipes

```bash
# Run all Python tests
just py-test

# Run with coverage
just py-test-cov

# Run tests for a specific module
just py-test-module ivar_variants_to_vcf

# Run with tox
just py-test-tox

# Run in parallel
just py-test-parallel

# Clean test artifacts
just py-test-clean
```

#### CI/CD Testing

Python tests are automatically run in GitHub Actions for all supported Python versions:

```yaml
# Triggered on push to main/dev/experimental branches and PRs
# Tests run with pytest and tox
# Coverage reports uploaded to Codecov
```

### Writing Python Tests

When adding new Python scripts or modifying existing ones:

1. Create a test file named `test_<script_name>.py` in the `bin/` directory
2. Use pytest fixtures for common test data
3. Mock external dependencies (file I/O, network calls, etc.)
4. Aim for high test coverage (>80%)
5. Use descriptive test names that explain what is being tested

Example test structure:

```python
import pytest
from unittest.mock import Mock, patch
from your_script import function_to_test

class TestYourFunction:
    """Test suite for your_function."""

    @pytest.fixture
    def sample_data(self):
        """Provide sample test data."""
        return {"key": "value"}

    def test_basic_functionality(self, sample_data):
        """Test that function works with valid input."""
        result = function_to_test(sample_data)
        assert result == expected_output

    @patch('your_script.external_dependency')
    def test_with_mocked_dependency(self, mock_dep):
        """Test function with mocked external calls."""
        mock_dep.return_value = "mocked_result"
        result = function_to_test()
        assert result == "expected_with_mock"
```

### Test Configuration

Test configuration is defined in `pyproject.toml`:

- **pytest settings**: Minimum version, test paths, markers, coverage options
- **coverage settings**: Source paths, omit patterns, exclusion rules
- **tox settings**: Python versions, test environments, dependencies

## Future Enhancements

- Sub-workflow tests for primer handling, alignment, and QC
- Performance benchmarking tests
- Integration tests with real sequencing data
- Automated test data generation scripts

## Environment Management

This test suite is designed to work seamlessly with the OneRoof Pixi environment. All test dependencies, including nf-test itself, are managed through the project's `pyproject.toml` and Pixi configuration. This ensures consistent testing environments across different machines and CI/CD systems.

When updating test dependencies:
1. Add them to `pyproject.toml` in the project root
2. Run `pixi install --frozen` to update the environment
3. Commit both `pyproject.toml` and `pixi.lock` to version control
