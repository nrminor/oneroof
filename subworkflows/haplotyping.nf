include { COMPRESS_AND_INDEX_VCF } from "../modules/bcftools"
include { PHASE_READS_WITH_DEVIDER } from "../modules/devider"

workflow HAPLOTYPING {
    take:
    ch_bams      // tuple(sample_id, bam, bai)
    ch_vcfs      // tuple(sample_id, vcf)
    ch_ref       // path(reference.fasta)

    main:

    // Compress and index VCFs for devider
    COMPRESS_AND_INDEX_VCF(ch_vcfs)

    // Join BAMs with VCFs by sample_id, then add reference
    // devider handles multi-contig BAMs natively, processing each contig independently
    ch_devider_input = ch_bams
        .combine(COMPRESS_AND_INDEX_VCF.out, by: 0)  // join by sample_id
        .combine(ch_ref)                              // add reference

    PHASE_READS_WITH_DEVIDER(ch_devider_input)
}

