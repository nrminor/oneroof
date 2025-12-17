include { INDEX_CONTAMINANTS } from "../modules/deacon"
include { RUN_DEACON         } from "../modules/deacon"

workflow DECONTAMINATE {
    take:
    ch_fastqs

    main:

    INDEX_CONTAMINANTS(contam_fasta)

    RUN_DEACON(ch_fastqs.combine(INDEX_CONTAMINANTS.out))

    GET_INDEX()

    DECON (
        ch_fastqs.combine(GET_INDEX.out)
    )

    }
    
    emit:
    RUN_DEACON.out
}
