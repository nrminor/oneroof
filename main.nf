#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { NANOPORE } from "../workflows/nanopore"
include { ILLUMINA } from "../workflows/illumina"

workflow {

    assert params.primer_bed
    assert file(params.primer_bed).isfile()
    assert params.refseq
    assert file(params.refseq).isfile()

    ch_primer_bed = Channel
        .fromPath( params.primer_bed )

    ch_refseq = Channel
        .fromPath( params.refseq )

    if ( params.remote_pod5_location || params.pod5_dir ) {

        NANOPORE (
            ch_primer_bed,
            ch_refseq
        )

    } else (

        assert params.illumina_fastq_dir : 
        "Please double check that a directory of Illumina FASTQs or Nanopore POD5s is provided." 
        assert file( params.illumina_fastq_dir ).isDirectory() : 
        "The provided Illumina FASTQ directory '$params.illumina_fastq_dir' does not exist."

        ILLUMINA (
            ch_primer_bed,
            ch_refseq
        )
    )

}
