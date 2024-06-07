#!/usr/bin/env nextflow

include { RESPLICE_PRIMERS } from "../modules/resplice_primers"
include { SPLIT_PRIMER_COMBOS } from "../modules/split_primer_combos"
include { GET_PRIMER_SEQS } from "../modules/bedtools"
include { TRIM_ENDS_TO_PRIMERS } from "../modules/cutadapt"
include { FASTQ_CONVERSION } from "../modules/samtools"
include { ORIENT_READS } from "../modules/vsearch"
include { FIND_COMPLETE_AMPLICONS; MERGE_BY_SAMPLE; DOWNSAMPLE_READS } from "../modules/seqkit"

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

        GET_PRIMER_PATTERNS (
            GET_PRIMER_SEQS.out
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
            GET_PRIMER_PATTERNS.out
        )

        TRIM_ENDS_TO_PRIMERS (
            FIND_COMPLETE_AMPLICONS.out
        )

        AMPLICON_STATS (
            TRIM_ENDS_TO_PRIMERS.out.groupTuple( by: 0 )
        )

        MERGE_BY_SAMPLE (
            TRIM_ENDS_TO_PRIMERS.out.groupTuple( by: 0 )
        )

        DOWNSAMPLE_READS (
            MERGE_BY_SAMPLE.out
        )
    
    emit:
        DOWNSAMPLE_READS.out
}
