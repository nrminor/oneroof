#!/usr/bin/env nextflow

workflow VARIANTS {

    /* */

    take:
        ch_amplicons
        ch_refseq
        ch_platform

    main:
        MEDAKA_VARIANTS (
            ch_amplicons,
            ch_refseq
        )
    
    emit:
        MEDAKA_VARIANTS.out
}