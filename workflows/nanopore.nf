#!/usr/bin/env nextflow

include { GATHER_DATA } from "../subworkflows/gather_data"
// include { ERROR_CORRECTION } from "../subworkflows/error_correction"
include { PRIMER_HANDLING } from "../subworkflows/primer_handling"
include { ALIGNMENT } from "../subworkflows/alignment"
include { QUALITY_CONTROL } from "../subworkflows/quality_control"
include { CONSENSUS } from "../subworkflows/consensus_calling"
include { VARIANTS } from "../subworkflows/variant_calling"
include { PHYLO } from "../subworkflows/phylo"

includeConfig '../config/nanopore.config'

workflow NANOPORE {

    /* */

    take:
        ch_primer_bed
        ch_refseq
        ch_refgbk
        ch_snpeff_config

    main:
        assert params.kit : "Please provide the Nanopore barcoding kit used."

        GATHER_DATA ( )

        // ERROR_CORRECTION (
        //     GATHER_DATA.out
        // )

        PRIMER_HANDLING (
            GATHER_DATA.out,
            ch_primer_bed,
            ch_refseq
        )

        ALIGNMENT (
            PRIMER_HANDLING.out,
            ch_refseq
        )

        QUALITY_CONTROL (
            PRIMER_HANDLING.out,
            ALIGNMENT.out
        )

        CONSENSUS (
            ALIGNMENT.out
        )

        VARIANTS (
            ALIGNMENT.out,
            ch_refseq,
            ch_refgbk,
            ch_snpeff_config
        )

        PHYLO (
            CONSENSUS.out
        )

}