#!/usr/bin/env nextflow

include { GATHER_ILLUMINA } from "../subworkflows/gather_illumina"

workflow ILLUMINA {

    take:
        ch_primer_bed
        ch_refseq
    
    main:
        assert params.platform == "illumina"

        GATHER_ILLUMINA ( )

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