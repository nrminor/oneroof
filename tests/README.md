# OneRoof Test Suite

The test suite is designed to validate pipeline functionality with minimal intrusion into the existing codebase.

## Installation

```bash
# Install all project dependencies 
pixi install --frozen

# Activate the Pixi environment
pixi shell --frozen
```

## Running Tests

All test commands assume you're working within the Pixi environment:

```bash
# First, activate the Pixi environment if not already active
pixi shell --frozen
```

### Run all tests 
```bash 
just test-all
```

### Run all tests per sequencing type
```bash
just test-illumina 
just test-nanopore
```

### Run specific Illumina tests 
```bash
just test-illumina-bad-primers
just test-illumina-missing-fastq
just test-illumina-with-metagenomics
just test-illumina-with-phylo
just test-illumina-with-primers     
just test-illumina-without-primers
```

### Run specific Nanopore tests
```bash
just test-nanopore-bad-primers
just test-nanopore-missing-input
just test-nanopore-with-haplo
just test-nanopore-with-metagenomics
just test-nanopore-with-phylo
just test-nanopore-with-primers    
just test-nanopore-without-primers
```

### Clean test outputs
```bash
# Remove all test artifacts
just clean-logs
```

## Test Structure

### Illumina Tests
Located in `conf/illumina_tests/`

### Nanopore Tests
Located in `conf/nanopore/`

## Test Data

The `tests/data/` directory contains minimal test datasets

All test data is internally consistent and designed to exercise key pipeline features while maintaining small file sizes for fast test execution.


## Writing New Tests

When adding new functionality to the pipeline:

1. Create test data if needed in `tests/data/`
2. Write a test config file in the correct directory under `conf/`
3. Add the test as a new process in `nextflow.config`
4. Create a new just recipe for the test in the `justfile`
5. Run the test locally before committing

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
