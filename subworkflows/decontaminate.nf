include { INDEX_CONTAMINANTS } from "../modules/deacon"
include { RUN_DEACON         } from "../modules/deacon"

workflow DECONTAMINATE {
    take:
    contam_fasta
    ch_fastqs

    main:

    INDEX_CONTAMINANTS(contam_fasta)

    RUN_DEACON(ch_fastqs.combine(INDEX_CONTAMINANTS.out))

    emit:
    RUN_DEACON.out
}
