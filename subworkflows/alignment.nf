#!/usr/bin/env nextflow

include { ALIGN_WITH_PRESET } from "../modules/minimap2"
include { CONVERT_AND_SORT; SORT_BAM; INDEX } from "../modules/samtools"
include { RASUSA_ALN } from "../modules/rasusa"
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

        RASUSA_ALN (
            CONVERT_AND_SORT.out
        )

        SORT_BAM (
            RASUSA_ALN.out
        )

        INDEX (
            SORT_BAM.out
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
