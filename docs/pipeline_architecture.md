# OneRoof Pipeline Architecture
OneRoof Development Team
2025-12-16

- [Main Workflow Entry Point](#main-workflow-entry-point)
  - [main.nf](#mainnf)
- [Platform-Specific Workflows](#platform-specific-workflows)
  - [NANOPORE Workflow (workflows/nanopore.nf)](#sec-nanopore)
  - [ILLUMINA Workflow (workflows/illumina.nf)](#sec-illumina)
- [Sub-workflows and Dependencies](#sub-workflows-and-dependencies)
  - [Data Gathering Workflows](#data-gathering-workflows)
  - [Processing Workflows](#processing-workflows)
  - [Optional Feature Workflows](#optional-feature-workflows)
- [Key Modules/Processes](#key-modulesprocesses)
  - [Critical Processes for Testing](#critical-processes-for-testing)
- [Critical Testing Paths](#critical-testing-paths)
  - [Minimal Test Path (No Primers)](#minimal-test-path-no-primers)
  - [Full Test Path (With Primers)](#full-test-path-with-primers)
  - [Key Input Requirements](#key-input-requirements)
- [Output Structure](#output-structure)
  - [Nanopore Output Tree](#nanopore-output-tree)
  - [Illumina Output Tree](#illumina-output-tree)
- [Configuration and Parameters](#configuration-and-parameters)
  - [Platform-Specific Defaults](#platform-specific-defaults)
  - [Resource Management](#resource-management)
  - [Key Process Labels](#key-process-labels)
- [Error Handling and Retry Strategy](#error-handling-and-retry-strategy)
- [Testing Considerations](#testing-considerations)
  - [Critical Validation Points](#critical-validation-points)
  - [Edge Cases to Test](#edge-cases-to-test)
  - [Integration Test Scenarios](#integration-test-scenarios)

This document provides a comprehensive map of the OneRoof Nextflow pipeline structure, including workflow dependencies, data flow, and critical testing points.

## Main Workflow Entry Point

### main.nf

- **Purpose**: Central orchestrator that routes to platform-specific workflows
- **Key Functions**:
  - Platform detection (Nanopore vs Illumina) based on input parameters
  - Input channel initialization for all required files
  - Workflow selection and invocation
  - Email notification on completion

**Input Channels**:

| Channel | Parameter | Description |
|---------|-----------|-------------|
| `ch_primer_bed` | `--primer_bed` | Optional primer BED file |
| `ch_refseq` | `--refseq` | Required reference FASTA |
| `ch_ref_gbk` | `--ref_gbk` | Optional GenBank file for annotation |
| `ch_contam_fasta` | `--contam_fasta` | Optional contamination sequences for dehosting |
| `ch_metagenomics_ref` | `--meta_ref` | Optional metagenomics reference (.syldb or FASTA) |
| `ch_snpeff_config` | `--snpEff_config` | Optional SnpEff configuration |
| `ch_primer_tsv` | `--primer_tsv` | Optional primer TSV file (alternative to BED) |
| `ch_sylph_tax_db` | `--sylph_tax_db` | Optional Sylph taxonomy identifier (e.g., "GTDB_r220") |
| `ch_meta_ref_link` | `--sylph_db_link` | Optional URL to download Sylph database |

## Platform-Specific Workflows

### NANOPORE Workflow (workflows/nanopore.nf)

**Workflow DAG**:

``` mermaid
graph TD
    A[GATHER_NANOPORE] --> B[PRIMER_HANDLING]
    A --> C[ALIGNMENT]
    B --> C
    B --> G[METAGENOMICS]
    A --> G
    C --> D[CONSENSUS]
    C --> E[VARIANTS]
    C --> F[HAPLOTYPING]
    D --> H[PHYLO]
    D --> I[SLACK_ALERT]
    E --> I
    C --> I

    style B fill:#f9f,stroke:#333,stroke-width:2px,stroke-dasharray: 5 5
    style F fill:#f9f,stroke:#333,stroke-width:2px,stroke-dasharray: 5 5
    style G fill:#f9f,stroke:#333,stroke-width:2px,stroke-dasharray: 5 5
    style H fill:#f9f,stroke:#333,stroke-width:2px,stroke-dasharray: 5 5
```

> [!NOTE]
>
> Dashed boxes indicate optional workflow components

**Key Parameters**:

- `platform = "ont"`
- `min_variant_frequency = 0.2`
- `min_qual = 10`
- `max_mismatch = 2` (higher tolerance for Nanopore error rate)

### ILLUMINA Workflow (workflows/illumina.nf)

**Workflow DAG**:

``` mermaid
graph TD
    A[GATHER_ILLUMINA] --> B[ILLUMINA_CORRECTION]
    B --> C[PRIMER_HANDLING]
    B --> D[ALIGNMENT]
    C --> D
    B --> G[METAGENOMICS]
    C --> G
    D --> E[CONSENSUS]
    D --> F[VARIANTS]
    E --> H[PHYLO]
    D --> I[SLACK_ALERT]
    E --> I
    F --> I

    style C fill:#f9f,stroke:#333,stroke-width:2px,stroke-dasharray: 5 5
    style G fill:#f9f,stroke:#333,stroke-width:2px,stroke-dasharray: 5 5
    style H fill:#f9f,stroke:#333,stroke-width:2px,stroke-dasharray: 5 5
```

**Key Parameters**:

- `platform = "illumina"`
- `min_variant_frequency = 0.05`
- `min_qual = 20`
- `max_mismatch = 0` (exact matching for Illumina's low error rate)

## Sub-workflows and Dependencies

### Data Gathering Workflows

#### GATHER_NANOPORE (subworkflows/gather_nanopore.nf)

**Purpose**: Handle multiple Nanopore input formats

**Input Options**:

1.  Remote POD5 monitoring (`remote_pod5_location`)
2.  Local POD5 directory (`pod5_dir`)
3.  Pre-called staging directory (`precalled_staging`)
4.  Pre-processed data directory (`prepped_data`)

**Process Flow**:

``` mermaid
graph LR
    A[POD5 Input] --> B[DOWNLOAD_MODELS]
    B --> C[BASECALL]
    C --> D[MERGE_BAMS]
    D --> E[DEMULTIPLEX]

    F[Pre-called Input] --> G[VALIDATE_NANOPORE]
    E --> G
    G --> H[FILTER_WITH_CHOPPER]
    H --> I[DECONTAMINATE]
    I --> J[COMPRESS_TO_SORTED_FASTA]
    J --> K[FAIDX]
    K --> L[EARLY_RASUSA_READ_DOWNSAMPLING]

    style I fill:#f9f,stroke:#333,stroke-width:2px,stroke-dasharray: 5 5
```

#### GATHER_ILLUMINA (subworkflows/gather_illumina.nf)

**Purpose**: Process paired-end Illumina FASTQ files

**Process Flow**:

``` mermaid
graph LR
    A[Paired FASTQs] --> B[VALIDATE_ILLUMINA]
    B --> C[MERGE_READ_PAIRS]
```

### Processing Workflows

#### ILLUMINA_CORRECTION (subworkflows/illumina_correction.nf)

**Purpose**: Quality control and decontamination for Illumina reads

**Process Flow**:

``` mermaid
graph TD
    A[CORRECT_WITH_FASTP] --> B[DECONTAMINATE]
    B --> C[FASTQC]
    C --> D[MULTIQC]
    B --> E[COMPRESS_TO_SORTED_FASTA]
    E --> F[FAIDX]
    F --> G[EARLY_RASUSA_READ_DOWNSAMPLING]

    style B fill:#f9f,stroke:#333,stroke-width:2px,stroke-dasharray: 5 5
```

#### PRIMER_HANDLING (subworkflows/primer_handling.nf)

**Purpose**: Validate primers and extract complete amplicons using the high-performance Rust-based amplicon finder

**Input Options**:

1.  Primer BED file (`--primer_bed`)
2.  Primer TSV file (`--primer_tsv`)

**Process Flow**:

``` mermaid
graph TD
    A[VALIDATE_PRIMER_BED] --> B[RESPLICE_PRIMERS]
    B --> C[SPLIT_PRIMER_COMBOS]
    C --> D[GET_PRIMER_SEQS]
    D --> E[GET_PRIMER_PATTERNS]
    D --> F[CREATE_PRIMER_TSV]
    F --> G[COLLECT_PRIMER_TSV]
    E --> H[FIND_AND_TRIM_AMPLICONS]
    H --> I[FAIDX]
    I --> J[READ_DOWNSAMPLING]
    J --> K[AMPLICON_STATS]
    J --> L[MERGE_BY_SAMPLE]
    K --> M[CREATE_AMPLICON_TSV]
```

**Key Component**: `FIND_AND_TRIM_AMPLICONS` uses a high-performance Rust script (`bin/find_and_trim_amplicons.rs`) that:
- Performs fuzzy primer matching with configurable mismatch tolerance
- Supports windowed search for performance optimization on long reads
- Handles both forward and reverse orientations
- Outputs trimmed amplicons in FASTA format

#### ALIGNMENT (subworkflows/alignment.nf)

**Purpose**: Map reads to reference and generate coverage statistics

**Process Flow**:

``` mermaid
graph TD
    A[ALIGN_WITH_PRESET] --> B[CONVERT_AND_SORT]
    B --> C[RASUSA_ALN_DOWNSAMPLING]
    C --> D[SORT_BAM]
    D --> E[INDEX]
    E --> F[MOSDEPTH]
    F --> G[PLOT_COVERAGE]
    G --> H[COVERAGE_SUMMARY]
```

#### VARIANTS (subworkflows/variant_calling.nf)

**Purpose**: Call and annotate variants

**Process Flow**:

``` mermaid
graph TD
    A[CALL_VARIANTS] --> B[CONVERT_TO_VCF]
    B --> C[ANNOTATE_VCF]
    C --> D[EXTRACT_FIELDS]
    D --> E[MERGE_VCF_FILES]
```

#### CONSENSUS (subworkflows/consensus_calling.nf)

**Purpose**: Generate consensus sequences

**Process Flow**:

``` mermaid
graph LR
    A[CALL_CONSENSUS] --> B[CONCAT]
```

### Optional Feature Workflows

#### PHYLO (subworkflows/phylo.nf)

**Purpose**: Phylogenetic analysis using Nextclade

**Condition**: Requires `--nextclade_dataset` parameter

**Process Flow**:

``` mermaid
graph LR
    A[CHECK_DATASET] --> B[DOWNLOAD_DATASET]
    B --> C[RUN_NEXTCLADE]
```

#### METAGENOMICS (subworkflows/metagenomics.nf)

**Purpose**: Metagenomic profiling using Sylph with optional taxonomic annotation via sylph-tax

**Condition**: Requires `--meta_ref` or `--sylph_db_link` parameter

**Input Options**:
1. Pre-built `.syldb` file (`--meta_ref`)
2. Custom FASTA to be sketched (`--meta_ref`)
3. URL to download `.syldb` (`--sylph_db_link`)

**Process Flow**:

``` mermaid
graph TD
    A[DOWNLOAD_SYLPH_DB] --> D[PROFILE_SAMPLES]
    B[SKETCH_DATABASE] --> D
    C[SKETCH_SAMPLES] --> D
    D --> E[ADD_TAXONOMY]
    E --> F[MERGE_TAXONOMY]

    style A fill:#f9f,stroke:#333,stroke-width:2px,stroke-dasharray: 5 5
    style B fill:#f9f,stroke:#333,stroke-width:2px,stroke-dasharray: 5 5
    style E fill:#f9f,stroke:#333,stroke-width:2px,stroke-dasharray: 5 5
    style F fill:#f9f,stroke:#333,stroke-width:2px,stroke-dasharray: 5 5
```

**Key Features**:
- Automatic database type detection (pre-built vs custom FASTA)
- Optional taxonomy overlay using sylph-tax with `--no-config` for container compatibility
- Supports taxonomy identifiers like `GTDB_r220`, `IMGVR_4.1`, etc.

#### HAPLOTYPING (subworkflows/haplotyping.nf)

**Purpose**: Viral haplotype reconstruction using devider (Nanopore only)

**Condition**: Number of reference sequences equals number of amplicons

**Process Flow**:

``` mermaid
graph TD
    A[SPLIT_SEGMENTS] --> B[COMPRESS_AND_INDEX_VCF]
    A --> C[PHASE_READS_WITH_DEVIDER]
    B --> C
```

**Key Parameters**:
- `--devider_preset`: Phasing preset (`nanopore-r9`, `nanopore-r10`, `hi-fi`, `old-long-reads`)
- `--min_haplo_reads`: Minimum read support for haplotype reporting

#### DECONTAMINATE (subworkflows/decontaminate.nf)

**Purpose**: Remove host/contaminant reads from samples

**Condition**: Requires `--contam_fasta` parameter

**Process Flow**:

``` mermaid
graph LR
    A[Input Reads] --> B[DEACON]
    B --> C[Filtered Reads]
```

## Key Modules/Processes

### Critical Processes for Testing

<div class="panel-tabset">

#### dorado.nf

- `DOWNLOAD_MODELS`: Model caching
- `BASECALL`: GPU-based basecalling
- `DEMULTIPLEX`: Barcode demultiplexing

#### find_and_trim_amplicons.nf

- `FIND_AND_TRIM_AMPLICONS`: High-performance Rust-based primer finding and amplicon trimming

#### minimap2.nf

- `ALIGN_WITH_PRESET`: Platform-specific alignment

#### ivar.nf

- `CALL_VARIANTS`: Variant detection
- `CALL_CONSENSUS`: Consensus generation
- `CONVERT_TO_VCF`: Format conversion

#### samtools.nf

- `CONVERT_AND_SORT`: BAM processing
- `INDEX`: BAM indexing
- `SPLIT_SEGMENTS`: Split BAM by reference segment (for haplotyping)

#### sylph.nf

- `DOWNLOAD_SYLPH_DB`: Download pre-built database from URL
- `SKETCH_DATABASE`: Sketch custom FASTA into sylph database
- `SKETCH_SAMPLES`: Sketch sample reads
- `PROFILE_SAMPLES`: Metagenomic profiling
- `DOWNLOAD_TAXONOMY`: Download sylph-tax metadata
- `ADD_TAXONOMY`: Add taxonomic annotations
- `MERGE_TAXONOMY`: Merge taxonomy profiles across samples

#### devider.nf

- `PHASE_READS_WITH_DEVIDER`: Haplotype phasing

#### validate.nf

- `VALIDATE_NANOPORE`: Input validation
- `VALIDATE_ILLUMINA`: Paired-end validation
- `VALIDATE_PRIMER_BED`: Primer validation

</div>

## Critical Testing Paths

### Minimal Test Path (No Primers)

1.  **Nanopore**: POD5/FASTQ → Basecall → Align → Consensus/Variants
2.  **Illumina**: Paired FASTQs → Merge → Correct → Align → Consensus/Variants

### Full Test Path (With Primers)

1.  Input validation
2.  Primer handling and amplicon extraction (Rust-based)
3.  Alignment and coverage analysis
4.  Variant calling and annotation
5.  Consensus generation
6.  Optional: Phylogenetics, metagenomics, haplotyping

### Key Input Requirements

> [!IMPORTANT]
>
> ### Minimal Requirements
>
> - Reference FASTA (`--refseq`)
> - Sequencing data:
>   - Nanopore: POD5 files + kit name OR pre-called BAM/FASTQ
>   - Illumina: Paired-end FASTQ directory

> [!TIP]
>
> ### Full Feature Requirements
>
> - Primer BED file (`--primer_bed`) or TSV (`--primer_tsv`)
> - Reference GenBank (`--ref_gbk`) for annotation
> - SnpEff config for variant annotation
> - Contamination FASTA (`--contam_fasta`) for decontamination
> - Metagenomics database (`--meta_ref` or `--sylph_db_link`) for classification
> - Taxonomy identifier (`--sylph_tax_db`) for taxonomic annotation
> - Nextclade dataset (`--nextclade_dataset`) for phylogenetics

## Output Structure

### Nanopore Output Tree

``` bash
nanopore/
├── 01_basecalled_demuxed/
│   ├── bams/
│   └── fastqs/
├── 02_primer_handling/
│   ├── 01_respliced_primers/
│   ├── 02_complete_amplicons/
│   └── 03_merged_by_sample/
├── 03_alignments/
│   ├── 01_genomecov/
│   └── 02_coverage_plots/
├── 04_consensus_seqs/
├── 05_variants/
│   ├── 01_ivar_tables/
│   ├── 02_annotated_vcfs/
│   ├── 03_variant_tsv/
│   └── 04_haplotypes/
├── 06_QC/
├── 07_phylo/
│   └── 01_nextclade/
├── 08_metagenomics/
└── haplotyping/
```

### Illumina Output Tree

``` bash
illumina/
├── 01_merged_reads/
├── 02_primer_handling/
│   ├── 01_respliced_primers/
│   ├── 02_complete_amplicons/
│   └── 03_merged_by_sample/
├── 03_alignments/
│   ├── 01_genomecov/
│   └── 02_coverage_plots/
├── 04_consensus_seqs/
├── 05_variants/
│   ├── 01_ivar_tables/
│   ├── 02_annotated_vcfs/
│   └── 03_variant_tsv/
├── 06_QC/
├── 07_phylo/
│   └── 01_nextclade/
└── 08_metagenomics/
```

## Configuration and Parameters

### Platform-Specific Defaults

| Parameter               | Nanopore  | Illumina |
|-------------------------|-----------|----------|
| `min_variant_frequency` | 0.2       | 0.05     |
| `min_qual`              | 10        | 20       |
| `max_mismatch`          | 2         | 0        |
| Alignment preset        | `map-ont` | `sr`     |

### Resource Management

| Parameter | Description |
|-----------|-------------|
| `pod5_batch_size` | Controls GPU memory usage during basecalling |
| `downsample_to` | Coverage depth limiting (post-alignment) |
| `early_downsample_to` | Coverage depth limiting (pre-alignment) |
| `basecall_max` | Parallel basecalling instances |
| `low_memory` | Resource-constrained mode |
| `forward_window` | Limit forward primer search to first N bases |
| `reverse_window` | Limit reverse primer search to last N bases |

### Key Process Labels

- `big_mem`: Memory-intensive processes (variant calling, consensus)
- GPU requirements: Dorado basecalling

## Error Handling and Retry Strategy

Most processes implement:

``` groovy
errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
maxRetries 2
```

This provides resilience against transient failures while preventing infinite loops.

## Testing Considerations

### Critical Validation Points

1.  Input file validation (exists, correct format)
2.  Primer validation (coordinates, sequences)
3.  Read count filtering (empty file handling)
4.  Platform-specific parameter application
5.  Optional workflow branching

### Edge Cases to Test

- Empty input files
- No reads passing filters
- Missing optional inputs
- Primer mismatches beyond tolerance
- Low coverage regions
- Multiple reference sequences
- Remote file watching timeout
- Pre-built vs custom metagenomics databases

### Integration Test Scenarios

1.  **Minimal run**: reference + reads only
2.  **Full featured run**: all optional inputs
3.  **Real-time processing**: file watching
4.  **Multi-sample processing**
5.  **Platform switching**: same data, different platforms
6.  **Metagenomics with taxonomy**: sylph + sylph-tax integration
7.  **Haplotyping**: multi-segment reference with primers
