#!/usr/bin/env nextflow

include { CALL_CONSENSUS } from "../modules/samtools"

workflow CONSENSUS {

    /* */

    take:
        ch_amplicons
        ch_refseq

    main:
        CALL_CONSENSUS (
            ch_amplicons,
            ch_refseq
        )

    emit:
        CALL_CONSENSUS.out

}