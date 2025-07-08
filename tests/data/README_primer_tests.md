# Primer-Related Test Suite

This directory contains comprehensive pytest test modules for the primer-related scripts in the OneRoof bioinformatics pipeline.

## Test Modules

### `test_make_primer_patterns.py`
Tests for `make_primer_patterns.py` which generates regex patterns from primer FASTA files.

**Key test cases:**
- Valid primer FASTA processing (bedtools getfasta format)
- Custom regex pattern generation
- Invalid FASTA handling (single sequence, multiple sequences)
- Coordinate validation and warnings
- Special characters in sequences
- Multi-line sequence handling
- Command-line argument parsing

### `test_resplice_primers.py`
Tests for `resplice_primers.py` which finds all possible combinations of spike-in primers.

**Key test cases:**
- BED file existence validation
- Primer index delimiter checking
- Amplicon partitioning
- Index assignment and normalization
- Primer name resolution for all combinations
- Complex spike-in scenarios
- Error handling for unpaired primers
- Custom delimiter support

### `test_split_primer_combos.py`
Tests for `split_primer_combos.py` which splits combined BED files into individual primer combinations.

**Key test cases:**
- Basic BED file splitting
- Complex multi-amplicon splitting
- Custom primer suffix handling
- Invalid/unpaired primer detection
- Empty file handling
- Duplicate combination handling
- Mixed chromosome support
- Column preservation

### `test_validate_primer_bed.py`
Tests for `validate_primer_bed.py` which validates and normalizes primer BED files.

**Key test cases:**
- Primer suffix validation
- Primer pair checking
- Coordinate orientation correction
- BED line normalization
- Missing suffix detection
- Unpaired primer detection
- Incomplete BED line handling
- Custom suffix support
- Extra column preservation

## Running the Tests

### Run all primer tests:
```bash
pytest tests/data/test_*.py -v
```

### Run a specific test module:
```bash
pytest tests/data/test_make_primer_patterns.py -v
```

### Run with coverage:
```bash
pytest tests/data/test_*.py --cov=bin --cov-report=html
```

### Use the test runner script:
```bash
python tests/data/run_primer_tests.py
```

## Test Design Principles

1. **Biological Correctness**: Tests focus on ensuring primer handling follows biological conventions
2. **Edge Case Coverage**: Comprehensive testing of malformed inputs, missing data, and boundary conditions
3. **Fixture-Based**: Reusable test data created through pytest fixtures
4. **Isolation**: Each test is independent and uses temporary directories for file operations
5. **Realistic Scenarios**: Test cases based on real-world primer design and spike-in scenarios

## Dependencies

The test suite requires:
- pytest
- polars (for dataframe-based tests)
- loguru (for logging tests)
- Python 3.10+

## Adding New Tests

When adding new test cases:
1. Use descriptive test names that explain what is being tested
2. Create fixtures for reusable test data
3. Test both successful operations and error conditions
4. Include edge cases and boundary conditions
5. Document complex test scenarios with comments
6. Ensure tests clean up any created files
