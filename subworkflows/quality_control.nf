#!/usr/bin/env nextflow

include { FASTQC } from "../modules/fastqc"
include { CRAMINO } from "../modules/cramino"
include { MULTIQC } from "../modules/multiqc"

workflow QUALITY_CONTROL {

    take:
        ch_reads
        ch_alignments

    main:
        FASTQC (
            ch_reads
        )

        // CRAMINO (
        //     ch_alignments
        // )

        MULTIQC (
            FASTQC.out.zip.collect()
        )

}
