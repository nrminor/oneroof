include { INDEX_CONTAMINANTS } from "../modules/deacon"
include { DECON } from "../modules/deacon"

workflow DECONTAMINATE {
    take:
    contam_fasta
    ch_fastqs
    

    main:

    INDEX_CONTAMINANTS (
        contam_fasta
    )

    DECON (
        ch_fastqs.combine(INDEX_CONTAMINANTS.out)
    )
    

    emit:
    DECON.out
}