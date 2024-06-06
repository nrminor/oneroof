#!/usr/bin/env nextflow

include { GATHER_DATA } from "../subworkflows/gather_data"
include { ALIGNMENT } from "../subworkflows/alignment"
include { CONSENSUS } from "../subworkflows/consensus_calling"
include { VARIANTS } from "../subworkflows/variant_calling"
include { REPORTING } from "../subworkflows/reporting"

workflow NANOPORE {

    /* */

    take:
        ch_primer_bed
        ch_refseq

    main:
        GATHER_DATA ( )

        PRIMER_HANDLING (
            DORADO.out,
            ch_primer_bed,
            ch_refseq
        )

        ALIGNMENT (
            PRIMER_HANDLING.out,
            ch_refseq
        )

        CONSENSUS (
            ALIGNMENT.out.groupTuple( by:  ),
            ch_refseq
        )

        VARIANTS (
            ALIGNMENT.out,
            ch_refseq,
            "ont"
        )

        // if ( params.reporting ) {
            
        //     REPORTING (
        //         ALIGNMENT.out,
        //         CONSENSUS.out,
        //         VARIANTS.out
        //     )

        // }

}