#!/usr/bin/env nextflow

// include { BBMERGE } from "../modules/bbmap"

workflow GATHER_ILLUMINA {

    /* */

    main:
        assert params.illumina_fastq_dir != "" && file(params.illumina_fastq_dir).isDirectory() :
        "No local Illumina paired-end FASTQ directories were provided."

        ch_prepped = Channel
            .fromFilePairs ( "${params.illumina_fastq_dir}/*{1,2}.fastq.gz", flat: true, maxDepth: 1 )
        
    //     BBMERGE (
    //         ch_prepped
    //     )
    
    // emit:
    //     BBMERGE.out

}