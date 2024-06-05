#!/usr/bin/env nextflow

workflow CONSENSUS {

    /* */

    take:
        ch_amplicons
        ch_refseq

    main:
        MEDAKA_CONSENSUS (
            ch_amplicons,
            ch_refseq
        )

    emit:
        MEDAKA_CONSENSUS.out

}