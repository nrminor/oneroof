include { INDEX_CONTAMINANTS } from "../modules/deacon"
include { DECON } from "../modules/deacon"
include { GET_INDEX } from "../modules/deacon"

workflow DECONTAMINATE {
    take:
    ch_fastqs
    contam_fasta
    

    main:

    if(contam_fasta){

    INDEX_CONTAMINANTS (
        contam_fasta.view()
    )

    DECON (
        ch_fastqs.combine(INDEX_CONTAMINANTS.out)
    )
    } else {

    GET_INDEX()

    DECON (
        ch_fastqs.combine(GET_INDEX.out)
    )

    }
    
    emit:
    DECON.out
}