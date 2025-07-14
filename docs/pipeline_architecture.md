# OneRoof Pipeline Architecture
OneRoof Development Team
2025-07-14

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

- `ch_primer_bed`: Optional primer BED file
- `ch_refseq`: Required reference FASTA
- `ch_ref_gbk`: Optional GenBank file for annotation
- `ch_contam_fasta`: Optional contamination sequences
- `ch_metagenomics_ref`: Optional metagenomics reference
- `ch_snpeff_config`: Optional SnpEff configuration
- `ch_primer_tsv`: Optional primer TSV file
- `ch_sylph_tax_db`: Optional Sylph taxonomy database

## Platform-Specific Workflows

### NANOPORE Workflow (workflows/nanopore.nf)

**Workflow DAG**:

``` mermaid
graph TD
    A[GATHER_NANOPORE] --> B[PRIMER_HANDLING]
    B --> C[ALIGNMENT]
    C --> D[CONSENSUS]
    C --> E[VARIANTS]
    C --> F[HAPLOTYPING]
    C --> G[METAGENOMICS]
    D --> H[PHYLO]
    E --> I[SLACK_ALERT]

    style B fill:#f9f,stroke:#333,stroke-width:2px,stroke-dasharray: 5 5
```

> [!NOTE]
>
> Dashed boxes indicate optional workflow components

**Key Parameters**:

- `platform = "ont"`
- `min_variant_frequency = 0.2`
- `min_qual = 10`

### ILLUMINA Workflow (workflows/illumina.nf)

**Workflow DAG**:

``` mermaid
graph TD
    A[GATHER_ILLUMINA] --> B[ILLUMINA_CORRECTION]
    B --> C[PRIMER_HANDLING]
    C --> D[ALIGNMENT]
    D --> E[CONSENSUS]
    D --> F[VARIANTS]
    D --> G[PHYLO]
    D --> H[METAGENOMICS]
    E --> I[SLACK_ALERT]
    F --> I

    style C fill:#f9f,stroke:#333,stroke-width:2px,stroke-dasharray: 5 5
```

**Key Parameters**:

- `platform = "illumina"`
- `min_variant_frequency = 0.05`
- `min_qual = 20`

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
    H --> I[COMPRESS_TO_SORTED_FASTA]
    I --> J[FAIDX]
    J --> K[EARLY_RASUSA_READ_DOWNSAMPLING]
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

**Purpose**: Validate primers and extract complete amplicons

**Input Options**:

1.  Primer BED file
2.  Primer TSV file

**Process Flow**:

``` mermaid
graph TD
    A[ORIENT_READS] --> B[GET_PRIMER_PATTERNS]
    B --> C[FIND_COMPLETE_AMPLICONS]
    B --> D[TRIM_ENDS_TO_PRIMERS]
    D --> E[PER_AMPLICON_FILTERS]
    E --> F[MERGE_BY_SAMPLE]
```

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

**Process Flow**:

``` mermaid
graph LR
    A[CHECK_DATASET] --> B[DOWNLOAD_DATASET]
    B --> C[RUN_NEXTCLADE]
```

#### METAGENOMICS (subworkflows/metagenomics.nf)

**Purpose**: Metagenomic classification using Sylph

**Process Flow**:

``` mermaid
graph TD
    A[SKETCH_DATABASE_KMERS] --> C[CLASSIFY_SAMPLE]
    B[SKETCH_SAMPLE_KMERS] --> C
    C --> D[OVERLAY_TAXONOMY]
    D --> E[MERGE_TAXONOMY]

    style D fill:#f9f,stroke:#333,stroke-width:2px,stroke-dasharray: 5 5
    style E fill:#f9f,stroke:#333,stroke-width:2px,stroke-dasharray: 5 5
```

#### HAPLOTYPING (subworkflows/haplotyping.nf)

**Purpose**: Viral haplotype reconstruction (Nanopore only)

**Condition**: Number of reference sequences equals number of amplicons

## Key Modules/Processes

### Critical Processes for Testing

<div class="panel-tabset">

#### dorado.nf

- `DOWNLOAD_MODELS`: Model caching
- `BASECALL`: GPU-based basecalling
- `DEMULTIPLEX`: Barcode demultiplexing

#### minimap2.nf

- `ALIGN_WITH_PRESET`: Platform-specific alignment

#### ivar.nf

- `CALL_VARIANTS`: Variant detection
- `CALL_CONSENSUS`: Consensus generation
- `CONVERT_TO_VCF`: Format conversion

#### samtools.nf

- `CONVERT_AND_SORT`: BAM processing
- `INDEX`: BAM indexing

#### validate.nf

- `VALIDATE_NANOPORE`: Input validation
- `VALIDATE_ILLUMINA`: Paired-end validation
- `VALIDATE_PRIMER_BED`: Primer validation

</div>

## Critical Testing Paths

### Minimal Test Path (No Primers)

1.  **Nanopore**: POD5/FASTQ → Basecall → Align → Consensus/Variants
2.  **Illumina**: Paired FASTQs → Merge → Align → Consensus/Variants

### Full Test Path (With Primers)

1.  Input validation
2.  Primer handling and amplicon extraction
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
> - Contamination FASTA for decontamination
> - Metagenomics database for classification

## Output Structure

### Nanopore Output Tree

``` bash
nanopore/
├── 01_basecalled_demuxed/
├── 02_primer_handling/
├── 03_alignments/
├── 04_consensus_seqs/
├── 05_variants/
├── 06_QC/
├── 07_phylo/
├── metagenomics/
└── haplotyping/
```

### Illumina Output Tree

``` bash
illumina/
├── 01_merged_reads/
├── 02_primer_handling/
├── 03_alignments/
├── 04_consensus_seqs/
├── 05_variants/
├── 06_QC/
├── 07_phylo/
└── metagenomics/
```

## Configuration and Parameters

### Platform-Specific Defaults

| Parameter               | Nanopore  | Illumina |
|-------------------------|-----------|----------|
| `min_variant_frequency` | 0.2       | 0.05     |
| `min_qual`              | 10        | 20       |
| Alignment preset        | `map-ont` | `sr`     |

### Resource Management

- `pod5_batch_size`: Controls GPU memory usage
- `downsample_to`: Coverage depth limiting
- `basecall_max`: Parallel basecalling instances
- `low_memory`: Resource-constrained mode

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
- Primer mismatches
- Low coverage regions
- Multiple reference sequences
- Remote file watching timeout

### Integration Test Scenarios

1.  **Minimal run**: reference + reads only
2.  **Full featured run**: all optional inputs
3.  **Real-time processing**: file watching
4.  **Multi-sample processing**
5.  **Platform switching**: same data, different platforms
