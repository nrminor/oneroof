#!/usr/bin/env nextflow

include { DOWNLOAD_DATASET; NEXTCLADE_RUN } from "../modules/nextclade"

workflow PHYLO {

    take:
        ch_consensus
    
    main:
        DOWNLOAD_DATASET ( )

        NEXTCLADE_RUN (
            ch_consensus,
            DOWNLOAD_DATASET.out
        )

}