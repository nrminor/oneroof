#!/usr/bin/env nextflow

include { CALL_CONSENSUS } from "../modules/ivar"
include { CONCAT } from "../modules/concat_consensus"

workflow CONSENSUS {

    /* */

    take:
        ch_mpileups

    main:
        CALL_CONSENSUS (
            ch_mpileups
        )

        CONCAT (
            CALL_CONSENSUS.out
                .map { id, fasta -> fasta }
                .collect()
        )

    emit:
        CONCAT.out

}