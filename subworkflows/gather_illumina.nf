#!/usr/bin/env nextflow

include { PUBLISH_COMMAND   } from "../modules/reporting"
include { VALIDATE_ILLUMINA } from "../modules/validate"
include { BBMERGE ; CLUMP_READS } from "../modules/bbmap"

workflow GATHER_ILLUMINA {
    main:
    ch_prepped = Channel.fromFilePairs("${params.illumina_fastq_dir}/*{R1,R2}*.fastq.gz", flat: true, maxDepth: 1)


    PUBLISH_COMMAND()

    VALIDATE_ILLUMINA(
        ch_prepped
    )

    BBMERGE(
        VALIDATE_ILLUMINA.out.map { id, reads1, reads2, _report -> tuple(id, file(reads1), file(reads2)) }
    )

    CLUMP_READS(
        BBMERGE.out
    )

    emit:
    CLUMP_READS.out
}
