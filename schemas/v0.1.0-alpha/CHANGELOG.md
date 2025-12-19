# Schema Changelog

## v0.1.0-alpha (Initial Alpha)

Initial schema version for OneRoof reporting. Alpha version indicates
the schema is experimental and breaking changes should be expected.

### Top-level Structure

- `schema_version`: Version string (e.g., "0.1.0-alpha")
- `generated_at`: ISO 8601 timestamp
- `run_metadata`: Pipeline run configuration and metadata
- `summary`: Aggregate statistics across all samples
- `samples`: Per-sample metrics keyed by sample_id

### Sample Metrics

Each sample includes:

- `sample_id`: Unique identifier
- `qc_status`: One of "pass", "warn", "fail"
- `qc_notes`: Human-readable QC issue descriptions
- `alignment`: Read mapping statistics
- `variants`: Variant calling results (optional)
- `consensus`: Consensus sequence metrics (optional)
- `metagenomics`: Sylph profiling results (optional)
- `haplotyping`: Devider haplotype phasing results (optional, Nanopore only)

### QC Thresholds (defaults)

- Coverage: pass ≥95%, warn ≥80%, fail <80% genome at ≥10x
- Completeness: pass ≥98%, warn ≥90%, fail <90%
- N percentage: pass <1%, warn <5%, fail ≥5%
- Mapped reads: pass ≥1000, warn ≥100, fail <100
