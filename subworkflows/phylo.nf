#!/usr/bin/env nextflow

include { DOWNLOAD_DATASET; RUN_NEXTCLADE } from "../modules/nextclade"

workflow PHYLO {

    take:
        ch_consensus
    
    main:
        DOWNLOAD_DATASET ( )

        RUN_NEXTCLADE (
            ch_consensus,
            DOWNLOAD_DATASET.out
        )

}