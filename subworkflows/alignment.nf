#!/usr/bin/env nextflow

include { ALIGN_WITH_MINIMAP } from "../modules/minimap2"
include { TRANSFORM_WITH_SAMTOOLS } from "../modules/samtools"

workflow ALIGNMENT {

    /* */

    take:
        ch_amplicons,
        ch_refseq

    main:
        ALIGN_WITH_MINIMAP (
            ch_amplicons,
            ch_refseq
        )

        TRANSFORM_WITH_SAMTOOLS (
            ALIGN_WITH_MINIMAP.out
        )

    emit:
        TRANSFORM_WITH_SAMTOOLS.out

}