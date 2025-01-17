#!/usr/bin/env nextflow

include { CHECK_DATASET ; DOWNLOAD_DATASET ; RUN_NEXTCLADE } from "../modules/nextclade"

workflow PHYLO {
    take:
    ch_consensus

    main:
    CHECK_DATASET()

    DOWNLOAD_DATASET(
        CHECK_DATASET.out
    )

    RUN_NEXTCLADE(
        ch_consensus,
        DOWNLOAD_DATASET.out,
    )
}
