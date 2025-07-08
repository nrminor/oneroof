# Test Data for OneRoof Pipeline

This directory contains minimal test data files for testing the OneRoof pipeline.

## Files

### Reference Files
- `test_reference.fasta` - A synthetic 1500bp reference genome
- `test_reference.gbk` - GenBank format with 4 annotated genes

### Primer Files
- `test_primers.bed` - BED file defining 4 amplicons with left/right primers

### Sequencing Data
- `test_illumina_R1.fastq` / `test_illumina_R2.fastq` - Paired-end Illumina reads (10 reads each)
- `test_nanopore.fastq` - Nanopore long reads (10 reads)

### Additional Files
- `test_pod5/` - Directory placeholder for POD5 test files
- `sample_barcode_mapping.csv` - Sample barcode mapping for multiplexed samples

## Test Data Properties

1. **Reference Genome**: 1500bp synthetic sequence with repeating patterns
   - Position 1-360: ATGGCAACTAT repeat (gene1)
   - Position 361-720: CGACGCTCGACGCT repeat (gene2)
   - Position 721-1080: TTAGGATTCC repeat (gene3)
   - Position 1081-1440: AACCGGTTAA repeat (gene4)

2. **Amplicons**: 4 overlapping amplicons covering the genome
   - Amplicon 1: 50-370 (320bp)
   - Amplicon 2: 300-670 (370bp)
   - Amplicon 3: 600-1020 (420bp)
   - Amplicon 4: 950-1370 (420bp)

3. **Variants**: Several reads contain mutations for testing variant calling
   - Read 2 (Illumina): A->G mutation at position ~100
   - Read 6 (Illumina): A->C mutation at position ~800
   - Read 2 (Nanopore): A->G mutation at position ~280
   - Read 6 (Nanopore): A->G mutation at position ~1200

4. **Quality**: All reads have high quality scores (I = Q40) except one low-quality Nanopore read

## Usage

These test files can be used to verify basic pipeline functionality:

```bash
# Test with Illumina data
nextflow run . \
  --illumina_fastq_dir tests/data/ \
  --primer_bed tests/data/test_primers.bed \
  --refseq tests/data/test_reference.fasta \
  --ref_gbk tests/data/test_reference.gbk

# Test with Nanopore data
nextflow run . \
  --nanopore_fastq tests/data/test_nanopore.fastq \
  --primer_bed tests/data/test_primers.bed \
  --refseq tests/data/test_reference.fasta \
  --ref_gbk tests/data/test_reference.gbk
```
