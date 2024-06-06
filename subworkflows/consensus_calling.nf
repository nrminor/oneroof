#!/usr/bin/env nextflow

include { CALL_CONSENSUS } from "../modules/samtools"

workflow CONSENSUS {

    /* */

    take:
        ch_aligned_amplicons

    main:
        CALL_CONSENSUS (
            ch_aligned_amplicons
        )

    emit:
        CALL_CONSENSUS.out

}