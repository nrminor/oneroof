#!/usr/bin/env nextflow

include { ALIGN_WITH_MINIMAP } from "../modules/minimap2"
include { CONVERT_AND_SORT; INDEX } from "../modules/samtools"

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

        CONVERT_AND_SORT (
            ALIGN_WITH_MINIMAP.out
        )

        INDEX (
            CONVERT_AND_SORT.out
        )



    emit:
        INDEX.out

}