include { CORRECT_WITH_FASTP } from "../modules/fastp"
include { FASTQ_TO_FASTA } from "../modules/seqkit"

workflow ILLUMINA_CORRECTION {
    take:
    ch_uncorrected_reads

    main:
    CORRECT_WITH_FASTP(
        ch_uncorrected_reads
    )

    FASTQ_TO_FASTA(
        CORRECT_WITH_FASTP.out
    )

    emit:
    FASTQ_TO_FASTA.out
}
