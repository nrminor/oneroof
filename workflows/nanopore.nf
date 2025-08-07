#!/usr/bin/env nextflow

include { GATHER_NANOPORE } from "../subworkflows/gather_nanopore"
include { PRIMER_HANDLING } from "../subworkflows/primer_handling"
include { ALIGNMENT } from "../subworkflows/alignment"
include { HAPLOTYPING } from "../subworkflows/haplotyping"
include { QUALITY_CONTROL } from "../subworkflows/quality_control"
include { CONSENSUS } from "../subworkflows/consensus_calling"
include { VARIANTS } from "../subworkflows/variant_calling"
include { METAGENOMICS } from "../subworkflows/metagenomics"
include { PHYLO } from "../subworkflows/phylo"
include { SLACK_ALERT } from "../subworkflows/slack_alert"

workflow NANOPORE {

    /* */

    take:
        ch_primer_bed
        ch_refseq
        ch_refgbk
        _ch_contam_fasta
        ch_snpeff_config
        ch_metagenome_ref
        ch_primer_tsv
        ch_sylph_tax_db
        ch_sylph_db_link

    main:
        assert params.platform == "ont"

        GATHER_NANOPORE ( )

        if ( params.primer_bed && params.primer_bed != "" || params.primer_tsv && params.primer_tsv != "" ) {

            PRIMER_HANDLING (
                GATHER_NANOPORE.out,
                ch_primer_bed,
                ch_refseq,
                ch_primer_tsv
            )

            METAGENOMICS(
                ch_metagenome_ref,
                ch_sylph_tax_db,
                ch_sylph_db_link,
                PRIMER_HANDLING.out
            )

            alignment_outputs = ALIGNMENT (
                PRIMER_HANDLING.out,
                ch_refseq
            )

        } else {

            METAGENOMICS(
                ch_metagenome_ref,
                ch_sylph_tax_db,
                ch_sylph_db_link,
                GATHER_NANOPORE.out
            )

            alignment_outputs = ALIGNMENT (
                GATHER_NANOPORE.out,
                ch_refseq
            )

        }

        CONSENSUS (
            alignment_outputs.index
        )

        variant_outputs = VARIANTS (
            alignment_outputs.index,
            ch_refseq,
            ch_refgbk,
            ch_snpeff_config
        )

        if ( params.primer_bed && Utils.countFastaHeaders(params.refseq) == Utils.countAmplicons(params.primer_bed) ) {

            HAPLOTYPING (
                alignment_outputs.index,
                variant_outputs.annotate,
                ch_refseq
            )

        }

        PHYLO (
            CONSENSUS.out
        )

        SLACK_ALERT(
            alignment_outputs.coverage_summary,
            CONSENSUS.out.collect(),
            variant_outputs.merge_vcf_files
        )


}
