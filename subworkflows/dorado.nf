#!/usr/bin/env nextflow

include { WATCH_FOR_POD5S } from "../modules/pod5_watcher"
include { DORADO_BASECALL } from "../modules/basecall"
include { DORADO_DEMUX } from "../modules/demux"

workflow DORADO {

    /* */

    main:
        if ( params.remote_pod5_location ) {

            ch_remote_pod5s = Channel
                .from ( params.remote_pod5_location )

            WATCH_FOR_POD5S (
                ch_remote_pod5s
            )

            ch_staged_pod5s = Channel
                .watchPath ( "${params.pod5_staging}/*.bam", 'create' )
                .subscribe { println "A new POD5 has arrived: $it" }
                // .until( ??? ) TODO

            DORADO_BASECALL (
                ch_staged_pod5s
            )

        } else {

            assert file(params.pod5_dir).isDirectory() : 
            "Please double check that the provided POD5 directory exists: $params.pod5_dir"

            ch_pod5_dir = Channel
                .fromPath ( "${params.pod5_dir}/*.pod5" )

            DORADO_BASECALL (
                ch_pod5_dir
            )

        }

        DORADO_DEMUX (
            DORADO_BASECALL.out.collect()
        )
    
    emit:
        DORADO_DEMUX.out.flatten()
}