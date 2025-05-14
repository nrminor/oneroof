include { CORRECT_WITH_FASTP } from "../modules/fastp"
include { COMPRESS_TO_SORTED_FASTA } from "../modules/seqkit"

workflow ILLUMINA_CORRECTION {
    take:
    ch_uncorrected_reads

    main:
    CORRECT_WITH_FASTP(
        ch_uncorrected_reads
    )

    COMPRESS_TO_SORTED_FASTA(
        CORRECT_WITH_FASTP.out
    )

    emit:
    COMPRESS_TO_SORTED_FASTA.out
}
