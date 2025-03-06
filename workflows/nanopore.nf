#!/usr/bin/env nextflow

include { GATHER_NANOPORE } from "../subworkflows/gather_nanopore"
include { PRIMER_HANDLING } from "../subworkflows/primer_handling"
include { ALIGNMENT } from "../subworkflows/alignment"
include { HAPLOTYPING } from "../subworkflows/haplotyping"
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

        GATHER_NANOPORE ( )

        if ( params.primer_bed && params.primer_bed != "" ) {

            PRIMER_HANDLING (
                GATHER_NANOPORE.out,
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

        } else {

            ALIGNMENT (
                GATHER_NANOPORE.out,
                ch_refseq
            )

            QUALITY_CONTROL (
                GATHER_NANOPORE.out,
                ALIGNMENT.out
            )

        }

        CONSENSUS (
            ALIGNMENT.out
        )

        VARIANTS (
            ALIGNMENT.out,
            ch_refseq,
            ch_refgbk,
            ch_snpeff_config
        )

        if ( params.primer_bed && Utils.countFastaHeaders(params.refseq) == Utils.countAmplicons(params.primer_bed) ) {

            HAPLOTYPING (
                ALIGNMENT.out
            )

        }

        PHYLO (
            CONSENSUS.out
        )

}
