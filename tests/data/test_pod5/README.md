# Test POD5 Directory

This directory would contain POD5 files for testing Nanopore basecalling.
Since POD5 files are binary format and require specific tools to create,
this directory serves as a placeholder for test purposes.

In actual tests, you would either:
1. Use small real POD5 files
2. Mock the basecalling step
3. Skip POD5 tests and use FASTQ directly

For pipeline testing, the test_nanopore.fastq file in the parent directory
can be used to test the downstream processing steps.
