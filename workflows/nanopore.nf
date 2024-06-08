#!/usr/bin/env nextflow

include { GATHER_NANOPORE } from "../subworkflows/gather_nanopore"
// include { ERROR_CORRECTION } from "../subworkflows/error_correction"
include { PRIMER_HANDLING } from "../subworkflows/primer_handling"
include { ALIGNMENT } from "../subworkflows/alignment"
include { QUALITY_CONTROL } from "../subworkflows/quality_control"
include { CONSENSUS } from "../subworkflows/consensus_calling"
include { VARIANTS } from "../subworkflows/variant_calling"
include { PHYLO } from "../subworkflows/phylo"

workflow NANOPORE {

    /* */

    take:
        ch_primer_bed
        ch_refseq
        ch_refgbk
        ch_snpeff_config

    main:
        assert params.platform == "ont"
        assert params.kit : "Please provide the Nanopore barcoding kit used."

        GATHER_NANOPORE ( )

        // ERROR_CORRECTION (
        //     GATHER_NANOPORE.out
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