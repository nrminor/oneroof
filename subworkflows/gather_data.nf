#!/usr/bin/env nextflow

include { WATCH_FOR_POD5S } from "../modules/pod5_watcher"
include { DORADO_BASECALL } from "../modules/dorado"
include { DORADO_DEMUX } from "../modules/dorado"

workflow gather_data {

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
                .map { pod5 -> tuple( file(pod5).getParent(), file(pod5) ) }
                // .until( ??? ) TODO

            BASECALL (
                ch_staged_pod5s
            )

        } else {

            assert file(params.pod5_dir).isDirectory() : 
            "Please double check that the provided POD5 directory exists: $params.pod5_dir"

            ch_pod5_dir = Channel
                .fromPath ( "${params.pod5_dir}/*.pod5" )
                .map { pod5 -> tuple( file(pod5).getParent(), file(pod5) ) }

            BASECALL (
                ch_pod5_dir
            )

        }

        DEMULTIPLEX (
            BASECALL.out.groupTuple( by: 0 )
        )
    
    emit:
        DEMULTIPLEX.out.flatten()
}