#!/usr/bin/env nextflow

include { WATCH_FOR_POD5S } from "../modules/pod5_watcher"
include { DOWNLOAD_MODELS; BASECALL; DEMULTIPLEX } from "../modules/dorado"

workflow GATHER_DATA {

    /* */

    main:
        DOWNLOAD_MODELS ( )

        if ( params.remote_pod5_location ) {

            error(
                message = "Watching a remote directory for pod5's is not yet supported. Please transfer the directory of pod5's yourself and supply it with the `--pod5_dir` command line argument."
            )

            ch_remote_pod5s = Channel
                .from ( params.remote_pod5_location )

            WATCH_FOR_POD5S (
                ch_remote_pod5s
            )

            ch_staged_pod5s = Channel
                .watchPath ( "${params.pod5_staging}/*.bam", 'create' )
                // .until( ??? ) TODO

            BASECALL (
                DOWNLOAD_MODELS.out,
                ch_staged_pod5s
            )

        } else {

            assert file(params.pod5_dir).isDirectory() : 
            "Please double check that the provided POD5 directory exists: $params.pod5_dir"

            ch_pod5_dir = Channel
                .fromPath ( "${params.pod5_dir}/*.pod5" )
                .take ( params.pod5_batch_size )

            BASECALL (
                DOWNLOAD_MODELS.out,
                ch_pod5_dir
            )

        }

        DEMULTIPLEX (
            BASECALL.out.collect()
        )
    
    emit:
        DEMULTIPLEX.out
            .flatten()
            .map { demux_bam -> tuple( 
                        file(demux_bam).getSimpleName().replace("${params.kit}_", ""),
                        file(demux_bam)
                    ) 
            }
}