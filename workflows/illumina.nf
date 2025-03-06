#!/usr/bin/env nextflow

include { GATHER_ILLUMINA } from "../subworkflows/gather_illumina"
include { ILLUMINA_CORRECTION } from "../subworkflows/illumina_correction"
include { PRIMER_HANDLING } from "../subworkflows/primer_handling"
include { ALIGNMENT } from "../subworkflows/alignment"
include { QUALITY_CONTROL } from "../subworkflows/quality_control"
include { CONSENSUS } from "../subworkflows/consensus_calling"
include { VARIANTS } from "../subworkflows/variant_calling"
include { PHYLO } from "../subworkflows/phylo"

workflow ILLUMINA {

    take:
        ch_primer_bed
        ch_refseq
        ch_ref_gbk
        ch_snpeff_config

    main:
        assert params.platform == "illumina"
        assert params.illumina_fastq_dir != "" :
        "Please double check that a directory of Illumina FASTQs or Nanopore POD5s is provided."
        assert file( params.illumina_fastq_dir ).isDirectory() :
        "The provided Illumina FASTQ directory ${params.illumina_fastq_dir} does not exist."

        GATHER_ILLUMINA ( )

        ILLUMINA_CORRECTION (
            GATHER_ILLUMINA.out
        )

        if ( params.primer_bed ) {

            PRIMER_HANDLING (
                ILLUMINA_CORRECTION.out,
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
                ILLUMINA_CORRECTION.out,
                ch_refseq
            )

            QUALITY_CONTROL (
                ILLUMINA_CORRECTION.out,
                ALIGNMENT.out
            )

        }

        CONSENSUS (
            ALIGNMENT.out
        )

        VARIANTS (
            ALIGNMENT.out,
            ch_refseq,
            ch_ref_gbk,
            ch_snpeff_config
        )

        PHYLO (
            CONSENSUS.out
        )

}
