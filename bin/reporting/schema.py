"""
Pydantic models for OneRoof report JSON schema.

These models define the canonical data format consumed by:
- MultiQC custom content generation
- Static visualization scripts
- Future web dashboards

The schema is versioned (see OneRoofReport.schema_version) to support
backward compatibility as the format evolves.
"""

from datetime import datetime
from enum import Enum
from typing import Optional

from pydantic import BaseModel, Field


class QCStatus(str, Enum):
    """Quality control status for a sample."""

    PASS = "pass"
    WARN = "warn"
    FAIL = "fail"


class Platform(str, Enum):
    """Sequencing platform identifier."""

    ONT = "ont"
    ILLUMINA = "illumina"


class AlignmentMetrics(BaseModel):
    """Metrics from read alignment."""

    total_reads: int = Field(description="Total reads in input")
    mapped_reads: int = Field(description="Reads aligned to reference")
    mapping_rate: float = Field(ge=0, le=1, description="Fraction of reads mapped")
    mean_coverage: float = Field(ge=0, description="Mean depth across reference")
    median_coverage: Optional[float] = Field(default=None, ge=0)
    genome_coverage_at_1x: float = Field(
        ge=0, le=1, description="Fraction of genome with ≥1x"
    )
    genome_coverage_at_10x: float = Field(
        ge=0, le=1, description="Fraction of genome with ≥10x"
    )
    genome_coverage_at_100x: float = Field(
        ge=0, le=1, description="Fraction of genome with ≥100x"
    )


class VariantMetrics(BaseModel):
    """Metrics from variant calling."""

    total_called: int = Field(ge=0, description="Total variants called")
    consensus_variants: int = Field(
        ge=0, description="Variants above consensus threshold"
    )
    subclonal_variants: int = Field(
        ge=0, description="Variants below consensus threshold"
    )
    snps: int = Field(ge=0)
    insertions: int = Field(ge=0)
    deletions: int = Field(ge=0)
    by_effect: dict[str, int] = Field(
        default_factory=dict, description="Counts by SnpEff effect"
    )


class ConsensusMetrics(BaseModel):
    """Metrics from consensus sequence."""

    length: int = Field(gt=0, description="Total consensus length")
    n_count: int = Field(ge=0, description="Number of N bases")
    n_percentage: float = Field(ge=0, le=100, description="Percentage of N bases")
    completeness: float = Field(ge=0, le=1, description="1 - (n_count / length)")
    ambiguous_bases: int = Field(default=0, ge=0, description="Non-ACGTN bases")


class MetagenomicsHit(BaseModel):
    """Single hit from metagenomic profiling."""

    taxon: str
    taxid: Optional[int] = None
    ani: Optional[float] = Field(default=None, ge=0, le=100)
    relative_abundance: float = Field(ge=0, le=1)


class MetagenomicsMetrics(BaseModel):
    """Metrics from Sylph metagenomic profiling."""

    top_hits: list[MetagenomicsHit] = Field(default_factory=list)
    unknown_fraction: float = Field(ge=0, le=1, description="Fraction not classified")
    primary_on_target: bool = Field(
        description="Whether top hit matches expected target"
    )


class HaplotypingMetrics(BaseModel):
    """Metrics from Devider haplotype phasing (Nanopore only)."""

    haplotypes_called: int = Field(ge=0, description="Number of distinct haplotypes")
    reads_assigned: int = Field(ge=0, description="Reads successfully assigned")
    reads_unassigned: int = Field(ge=0, description="Reads not assigned")
    assignment_rate: float = Field(ge=0, le=1, description="Fraction of reads assigned")
    num_snps: int = Field(ge=0, description="SNP positions used for phasing")
    haplotype_abundances: list[float] = Field(
        default_factory=list, description="Abundance % per haplotype"
    )
    haplotype_depths: list[float] = Field(
        default_factory=list, description="Depth per haplotype"
    )


class SampleMetrics(BaseModel):
    """All metrics for a single sample."""

    sample_id: str
    qc_status: QCStatus
    qc_notes: list[str] = Field(
        default_factory=list, description="Human-readable QC notes"
    )

    alignment: AlignmentMetrics
    variants: Optional[VariantMetrics] = None
    consensus: Optional[ConsensusMetrics] = None
    metagenomics: Optional[MetagenomicsMetrics] = None
    haplotyping: Optional[HaplotypingMetrics] = None  # ONT only


class ReferenceInfo(BaseModel):
    """Information about the reference sequence."""

    name: str
    path: str
    length: int
    segment_count: int = 1
    segments: list[str] = Field(default_factory=list)


class PrimerInfo(BaseModel):
    """Information about primer scheme."""

    provided: bool
    path: Optional[str] = None
    amplicon_count: Optional[int] = None


class RunParameters(BaseModel):
    """Key parameters used for the run."""

    min_depth_coverage: int
    min_consensus_freq: float
    min_variant_frequency: float
    downsample_to: int = 0


class RunMetadata(BaseModel):
    """Metadata about the pipeline run."""

    platform: Platform
    reference: ReferenceInfo
    primers: PrimerInfo
    parameters: RunParameters
    start_time: Optional[datetime] = None
    end_time: Optional[datetime] = None
    duration_seconds: Optional[int] = None
    nextflow_version: Optional[str] = None
    oneroof_version: Optional[str] = None


class Summary(BaseModel):
    """Aggregate summary across all samples."""

    sample_count: int
    samples_pass: int
    samples_warn: int
    samples_fail: int
    mean_coverage_depth: float
    mean_genome_coverage: float
    total_variants_called: int
    unique_variants: int


class OneRoofReport(BaseModel):
    """Top-level report structure."""

    schema_version: str = "1.0.0"
    generated_at: datetime
    run_metadata: RunMetadata
    summary: Summary
    samples: dict[str, SampleMetrics]

    # Optional aggregated data (included at "standard" and "full" levels)
    coverage_by_position: Optional[dict] = None
    variant_matrix: Optional[dict] = None
