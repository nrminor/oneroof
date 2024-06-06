#!/usr/bin/env nextflow

include { GATHER_DATA } from "../subworkflows/gather_data"
include { ALIGNMENT } from "../subworkflows/alignment"
include { CONSENSUS } from "../subworkflows/consensus_calling"
include { VARIANTS } from "../subworkflows/variant_calling"
// include { REPORTING } from "../subworkflows/reporting"

// parameters derived from user-supplied parameters
// -----------------------------------------------------------------------------

// basecalling and demultiplexing results
params.basecall_results = params.results + "/01_basecalled_demuxed"
params.basecall_bams = params.basecall_results + "bams"
params.basecall_fastqs = params.basecall_results + "fastqs"

// primer handling results
params.primer_handling = params.results + "/02_primer_handling"
params.respliced = params.primer_handling + "/01_respliced_primers"
params.complete_amplicons = params.primer_handling + "/02_complete_amplicons"
params.merged_by_sample = params.primer_handling + "/03_merged_by_sampe"

// alignment results
params.alignment = params.results + "/03_alignment_results"

// consensus results
params.consensus = params.results + "/04_consensus_seqs"

// variant results
params.variants = params.results + "/05_variants"
params.ivar = params.variants + "/01_ivar_tables"
params.vcf = params.variants + "/02_annotated_vcfs"

// -----------------------------------------------------------------------------

workflow NANOPORE {

    /* */

    take:
        ch_primer_bed
        ch_refseq
        ch_refgbk
        ch_snpeff_config

    main:
        assert params.kit : "Please provide the Nanopore barcoding kit used."

        GATHER_DATA ( )

        PRIMER_HANDLING (
            GATHER_DATA.out,
            ch_primer_bed,
            ch_refseq
        )

        ALIGNMENT (
            PRIMER_HANDLING.out,
            ch_refseq,
            "ont"
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

        // if ( params.reporting ) {
            
        //     REPORTING (
        //         ALIGNMENT.out,
        //         CONSENSUS.out,
        //         VARIANTS.out
        //     )

        // }

}