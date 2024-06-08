#!/usr/bin/env nextflow

// include { BBMERGE } from "../modules/bbmap"

workflow GATHER_ILLUMINA {

    /* */

    main:
        ch_prepped = Channel
            .fromFilePairs ( "${params.illumina_fastq_dir}/*{1,2}.fastq.gz", flat: true, maxDepth: 1 )
        
    //     BBMERGE (
    //         ch_prepped
    //     )
    
    // emit:
    //     BBMERGE.out

}