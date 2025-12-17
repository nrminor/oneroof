include { SPLIT_SEGMENTS } from "../modules/samtools"
include { COMPRESS_AND_INDEX_VCF } from "../modules/bcftools"
include { PHASE_READS_WITH_DEVIDER } from "../modules/devider"

workflow HAPLOTYPING {
    take:
    ch_bams      // tuple(sample_id, bam, bai)
    ch_vcfs      // tuple(sample_id, vcf)
    ch_ref       // path(reference.fasta)

    main:

    // Split BAMs by segment (already indexes them)
    SPLIT_SEGMENTS(ch_bams)

    // Flatten to per-segment tuples
    ch_split_bams = SPLIT_SEGMENTS.out
        .flatMap { sample_id, bams, bais ->
            def bam_list = bams instanceof List ? bams : [bams]
            def bai_list = bais instanceof List ? bais : [bais]
            
            bam_list.indices.collect { i ->
                def bam = bam_list[i]
                def bai = bai_list[i]
                def segment = bam.baseName.replace("${sample_id}_", "")
                tuple(sample_id, segment, bam, bai)
            }
        }

    // Compress and index VCFs for devider
    COMPRESS_AND_INDEX_VCF(ch_vcfs)

    // Join BAMs with VCFs by sample_id, then add reference
    ch_devider_input = ch_split_bams
        .combine(COMPRESS_AND_INDEX_VCF.out, by: 0)  // join by sample_id
        .combine(ch_ref)                              // add reference
        .map { sample_id, segment, bam, bai, vcf_gz, tbi, ref ->
            tuple("${sample_id}_${segment}", bam, bai, vcf_gz, tbi, ref)
        }

    PHASE_READS_WITH_DEVIDER(ch_devider_input)
}
