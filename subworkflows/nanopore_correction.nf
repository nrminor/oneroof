#!/usr/bin/env nextflow

include { CHOPPER } from "../modules/chopper"
include { DECHAT } from "../modules/dechat"
include { RACON } from "../modules/racon"

workflow NANOPORE_CORRECTION {

    take:
        ch_uncorrected_reads

    main:
        CHOPPER (
            ch_uncorrected_reads
        )

        DECHAT (
            CHOPPER.out
        )

        RACON (
            DECHAT.out
        )

    emit:
        RACON.out

}
