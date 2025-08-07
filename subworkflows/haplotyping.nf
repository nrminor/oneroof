include { SPLIT_SEGMENTS ; FASTQ_CONVERSION } from "../modules/samtools"
include { IDENTIFY_HAPLOTYPES      } from "../modules/vsearch"
include { PHASE_READS_WITH_DEVIDER } from "../modules/devider"
// include { PROCESS_VCF   } from "../modules/bcftools"
// include { PROCESS_REF   } from "../modules/bcftools"

workflow HAPLOTYPING {
    take:
    ch_bams
    _ch_vcfs
    _ch_ref

    main:
    SPLIT_SEGMENTS(
        ch_bams
    )

    FASTQ_CONVERSION(
    SPLIT_SEGMENTS.out.flatten()
        .filter { ~/.bam$/ } // skip .bai files
        .map{bam -> tuple(file(bam).getSimpleName(), file(bam))}
    )

    IDENTIFY_HAPLOTYPES(
        FASTQ_CONVERSION.out
    )
}
