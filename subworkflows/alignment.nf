#!/usr/bin/env nextflow

include { ALIGN_WITH_PRESET } from "../modules/minimap2"
include { CONVERT_AND_SORT; INDEX } from "../modules/samtools"

workflow ALIGNMENT {

    /* */

    take:
        ch_amplicons
        ch_refseq
        ch_platform

    main:
        ALIGN_WITH_PRESET (
            ch_amplicons,
            ch_refseq,
            ch_platform
        )

        CONVERT_AND_SORT (
            ALIGN_WITH_PRESET.out
        )

        INDEX (
            CONVERT_AND_SORT.out
        )

    emit:
        INDEX.out

}