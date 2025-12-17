#!/usr/bin/env nextflow

include { PUBLISH_COMMAND   } from "../modules/reporting"
include { VALIDATE_ILLUMINA } from "../modules/validate"
include { MERGE_READ_PAIRS  } from "../modules/vsearch.nf"
include { DEDUP } from "../modules/deduplication.nf"

workflow GATHER_ILLUMINA {

    main:
    ch_prepped = channel.fromFilePairs("${params.illumina_fastq_dir}/*{R1,R2}*.fastq.gz", flat: true, maxDepth: 1)

    PUBLISH_COMMAND()

    VALIDATE_ILLUMINA(
        ch_prepped
    )
    

     if(params.dedup){

        DEDUP(VALIDATE_ILLUMINA.out)

        MERGE_READ_PAIRS(
        DEDUP.out 
    )
    } else {

    MERGE_READ_PAIRS(
        VALIDATE_ILLUMINA.out.map { id, reads1, reads2, _report -> tuple(id, file(reads1), file(reads2)) }
    )
    }


    emit:
     MERGE_READ_PAIRS.out
    
}
