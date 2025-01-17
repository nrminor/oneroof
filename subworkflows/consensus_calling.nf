#!/usr/bin/env nextflow

include { CALL_CONSENSUS } from "../modules/ivar"
include { CONCAT         } from "../modules/concat_consensus"

workflow CONSENSUS {
    take:
    ch_aligned_amplicons

    main:
    CALL_CONSENSUS(
        ch_aligned_amplicons
    )

    CONCAT(
        CALL_CONSENSUS.out.map { _id, fasta -> fasta }.collect()
    )

    emit:
    CONCAT.out
}
