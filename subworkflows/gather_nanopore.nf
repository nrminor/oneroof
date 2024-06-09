#!/usr/bin/env nextflow

include { WATCH_FOR_POD5S } from "../modules/pod5_watcher"
include { DOWNLOAD_MODELS; BASECALL; DEMULTIPLEX } from "../modules/dorado"
include { MERGE_BARCODES } from "../modules/samtools"
include { VALIDATE_NANOPORE } from "../modules/validate"

workflow GATHER_NANOPORE {

    /* */

    main:

        // Invoke Dorado basecaller and demultiplexer if pod5's are provided
        // ------------------------------------------------------------------ //

        if ( params.remote_pod5_location != "" || params.pod5_dir != "" ) {

            DOWNLOAD_MODELS ( )

            if ( params.remote_pod5_location != "" ) {

                error(
                    message = "Watching a remote directory for pod5's is not yet supported. Please transfer the directory of pod5's yourself and supply it with the `--pod5_dir` command line argument."
                )

                ch_remote_pod5s = Channel
                    .from ( params.remote_pod5_location )

                WATCH_FOR_POD5S (
                    ch_remote_pod5s
                )

                ch_staged_pod5s = Channel
                    .watchPath ( "${params.pod5_staging}/*.pod5", 'create' )
                    // .until( ??? ) TODO

                BASECALL (
                    DOWNLOAD_MODELS.out.collect(),
                    ch_staged_pod5s
                )

            } else {

                assert file(params.pod5_dir).isDirectory() :
                "Please double check that the provided POD5 directory exists: ${params.pod5_dir}"

                ch_pod5_dir = Channel
                    .fromPath ( "${params.pod5_dir}/*.pod5" )
                    // .take ( params.pod5_batch_size )

                BASECALL (
                    DOWNLOAD_MODELS.out.collect(),
                    ch_pod5_dir
                )

            }

            DEMULTIPLEX (
                BASECALL.out
            )

            MERGE_BARCODES (
                DEMULTIPLEX.out
                    .flatten()
                    .map { demux_bam ->
                            tuple(
                                file(demux_bam).getSimpleName().replace("${params.kit}_", ""),
                                file(demux_bam)
                            )
                    }
                    .groupTuple( by: 0 )
            )

            VALIDATE_NANOPORE (
                MERGE_BARCODES.output
            )

        // ------------------------------------------------------------------ //


        // Otherwise, work with provided FASTQ or BAM files
        // ------------------------------------------------------------------ //

        } else {

            assert params.prepped_data != "" && file(params.prepped_data).isDirectory() :
            """
            No local or remote pod5 directories were supplied. Please supply a directory of
            pre-called and pre-demultiplexed BAM or FASTQ files to run the Nanopore pipeline.
            """

            ch_prepped = Channel
                .fromPath( "${params.prepped_data}/*.fastq.gz" )
                .map { fastq -> tuple( file(fastq).getSimpleName(), file(fastq) ) }

            VALIDATE_NANOPORE (
                ch_prepped
            )

        }

        // ------------------------------------------------------------------ //

    emit:
        VALIDATE_NANOPORE.out
            .map { label, seq_file, status -> tuple( label, file(seq_file) ) }

}
