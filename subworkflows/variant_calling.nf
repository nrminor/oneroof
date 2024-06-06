#!/usr/bin/env nextflow

include { GENERATE_MPILEUP } from "../modules/samtools"
include { CALL_VARIANTS; CONVERT_TO_VCF } from "../modules/ivar"
include { BUILD_DB; ANNOTATE_VCF } from "../modules/snpeff"

workflow VARIANTS {

    /* */

    take:
        ch_amplicons
        ch_refseq
        ch_platform

    main:
        GENERATE_MPILEUP (
            ch_amplicons,
            ch_refseq
        )

        CALL_VARIANTS (
            GENERATE_MPILEUP.out
        )

        CONVERT_TO_VCF (
            CALL_VARIANTS
        )

        BUILD_DB ( )

        ANNOTATE_VCF (
            BUILD_DB.out,
            ANNOTATE_VCF.out
        )

    // emit:
    //     MEDAKA_VARIANTS.out

}