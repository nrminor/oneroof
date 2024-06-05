#!/usr/bin/env nextflow

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

        FIND_TEMPLATE_AMPLICONS (
            FASTQ_CONVERSION.out,
            GET_PRIMER_SEQS.out
        )

        FIND_COMPLEMENT_AMPLICONS (
            FASTQ_CONVERSION.out,
            GET_PRIMER_SEQS.out
        )

        MERGE_PER_SAMPLE (
            
        )
    
    emit:
        MERGE_AMPLICONS.out
}
