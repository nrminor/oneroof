#!/usr/bin/env nextflow

include { RESPLICE_PRIMERS } from "../modules/resplice_primers"
include { SPLIT_PRIMER_COMBOS } from "../modules/split_primer_combos"
include { GET_PRIMER_SEQS } from "../modules/get_primer_seqs"
include { FASTQ_CONVERSION } from "../modules/samtools"
include { ORIENT_READS } from "../modules/vsearch"
include { FIND_COMPLETE_AMPLICONS; MERGE_BY_SAMPLE } from "../modules/find_amplicons"

workflow PRIMER_HANDLING {

    /* */

    take:
        ch_basecalls
        ch_primer_bed
        ch_refseq

    main:
        RESPLICE_PRIMERS (
            ch_primer_bed
        )

        SPLIT_PRIMER_COMBOS (
            RESPLICE_PRIMERS.out
        )

        GET_PRIMER_SEQS (
            SPLIT_PRIMER_COMBOS.out.flatten(),
            ch_refseq
        )

        FASTQ_CONVERSION (
            ch_basecalls
        )

        ORIENT_READS (
            FASTQ_CONVERSION.out,
            ch_refseq
        )

        FIND_COMPLETE_AMPLICONS (
            ORIENT_READS.out
                .map { barcode, fastq -> fastq },
            GET_PRIMER_SEQS.out
        )

        MERGE_BY_SAMPLE (
            FIND_COMPLETE_AMPLICONS.out.groupTuple( by: 0 )
        )
    
    emit:
        MERGE_BY_SAMPLE.out
}
