#!/usr/bin/env nextflow

include { SPLIT_SEGMENTS; FASTQ_CONVERSION } from "../modules/samtools"
include { IDENTIFY_HAPLOTYPES } from "../modules/vsearch"

workflow HAPLOTYPING {

    take:
        ch_bams
        ch_refseq

    main:
        COUNT_SEGMENTS (
            ch_bams
        )

        SPLIT_SEGMENTS (
            COUNT_SEGMENTS.out
        )

        FASTQ_CONVERSION (
            SPLIT_SEGMENTS.out
                .flatten()
                .map { bam -> 
                    tuple( file(bam).getSimpleName().split("_")[0],  bam ) 
                }
        )

        IDENTIFY_HAPLOTYPES (
            FASTQ_CONVERSION.out
        )

}
