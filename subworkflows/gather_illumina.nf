#!/usr/bin/env nextflow

include { PUBLISH_COMMAND   } from "../modules/reporting"
include { VALIDATE_ILLUMINA } from "../modules/validate"
include { MERGE_READ_PAIRS  } from "../modules/vsearch.nf"
include { FAIDX             } from "../modules/samtools"
include { EARLY_RASUSA_READ_DOWNSAMPLING } from "../modules/rasusa"



workflow GATHER_ILLUMINA {
    main:
    ch_prepped = Channel.fromFilePairs("${params.illumina_fastq_dir}/*{R1,R2}*.fastq.gz", flat: true, maxDepth: 1)


    PUBLISH_COMMAND()

    VALIDATE_ILLUMINA(
        ch_prepped
    )

    MERGE_READ_PAIRS(
        VALIDATE_ILLUMINA.out.map { id, reads1, reads2, _report -> tuple(id, file(reads1), file(reads2)) }
    )

    FAIDX(
        MERGE_READ_PAIRS.out
            .map { id, fastq -> tuple(id, fastq, file(fastq).countFasta()) }
            .filter { it[2] > 0 }
            .map { id, fastq, _read_count -> tuple(id, file(fastq)) }
    )

    EARLY_RASUSA_READ_DOWNSAMPLING(
        FAIDX.out
    )

    emit:
    EARLY_RASUSA_READ_DOWNSAMPLING.out
}
