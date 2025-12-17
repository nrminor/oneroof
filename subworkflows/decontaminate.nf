include { INDEX_CONTAMINANTS ; RUN_DEACON ; GET_INDEX } from "../modules/deacon"

workflow DECONTAMINATE {
    take:
    ch_fastqs
    ch_contam_fasta

    main:

    INDEX_CONTAMINANTS(ch_contam_fasta)

    RUN_DEACON(ch_fastqs.combine(INDEX_CONTAMINANTS.out))

    GET_INDEX()

    RUN_DEACON(ch_fastqs.combine(GET_INDEX.out))

    emit:
    RUN_DEACON.out
}
