#!/usr/bin/env nextflow

include { CALL_CONSENSUS } from "../modules/ivar"
include { CONCAT         } from "../modules/concat_consensus"
include { EXTRACT_CONSENSUS_METRICS } from "../modules/reporting"

workflow CONSENSUS {
    take:
    ch_aligned_amplicons

    main:
    CALL_CONSENSUS(
        ch_aligned_amplicons
    )

    EXTRACT_CONSENSUS_METRICS(
        CALL_CONSENSUS.out
    )

    CONCAT(
        CALL_CONSENSUS.out.map { _id, fasta -> fasta }.collect()
    )

    emit:
    consensus_seqs    = CALL_CONSENSUS.out
    concat            = CONCAT.out
    consensus_metrics = EXTRACT_CONSENSUS_METRICS.out
}
