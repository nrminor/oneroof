include { SPLIT_SEGMENTS ; FASTQ_CONVERSION } from "../modules/samtools"
include { IDENTIFY_HAPLOTYPES      } from "../modules/vsearch"
include { PHASE_READS_WITH_DEVIDER } from "../modules/devider"

workflow HAPLOTYPING {
    take:
    ch_bams
    ch_vcfs
    ch_ref

    main:
    SPLIT_SEGMENTS(
        ch_bams
    )

    PHASE_READS_WITH_DEVIDER(
        SPLIT_SEGMENTS.out.join(ch_vcfs).combine(ch_ref)
    )

    FASTQ_CONVERSION(
        SPLIT_SEGMENTS.out.flatten().map { bam ->
            tuple(file(bam).getSimpleName().split("_")[0..-2].join('_'), bam)
        }
    )

    IDENTIFY_HAPLOTYPES(
        FASTQ_CONVERSION.out
    )
}
