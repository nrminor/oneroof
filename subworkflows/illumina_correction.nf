include { CORRECT_WITH_FASTP } from "../modules/fastp"

workflow ILLUMINA_CORRECTION {
    take:
    ch_uncorrected_reads

    main:
    CORRECT_WITH_FASTP(
        ch_uncorrected_reads
    )

    // SORT_AND_COMPRESS(
    //     QUALITY_TRIM.out
    // )

    emit:
    CORRECT_WITH_FASTP.out
}
