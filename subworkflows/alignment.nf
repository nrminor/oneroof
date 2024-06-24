#!/usr/bin/env nextflow

include { ALIGN_WITH_PRESET } from "../modules/minimap2"
include { CONVERT_AND_SORT; INDEX } from "../modules/samtools"
include { MOSDEPTH } from "../modules/mosdepth"
include { PLOT_COVERAGE } from "../modules/plot_coverage"

workflow ALIGNMENT {

    /* */

    take:
        ch_amplicons
        ch_refseq

    main:
        ALIGN_WITH_PRESET (
            ch_amplicons,
            ch_refseq
        )

        CONVERT_AND_SORT (
            ALIGN_WITH_PRESET.out
        )

        INDEX (
            CONVERT_AND_SORT.out
        )

        MOSDEPTH (
            INDEX.out
        )

        PLOT_COVERAGE (
            MOSDEPTH.out
        )

    emit:
        INDEX.out

}
