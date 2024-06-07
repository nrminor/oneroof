#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { NANOPORE } from "$projectDir/workflows/nanopore"
// include { ILLUMINA } from "$projectDir/workflows/illumina"

workflow {

    // make sure required primer bed is provided and exists
    assert params.primer_bed != "" 
    : "Please provide a primer bed file with the parameter `primer_bed`."
    assert file(params.primer_bed).isFile() 
    : "Please double check that the primer bed file provided with the parameter `primer_bed` exists."
    
    // make sure provided refseq is provided and exists
    assert params.refseq != "" 
    : "Please provide a reference FASTA file with the parameter `refseq`."
    assert file(params.refseq).isFile() 
    : "Please double check that the reference FASTA file provided with the parameter `refseq` exists."
    
    // make sure required reference genbank is provided and exists
    assert params.ref_gbk != "" 
    : "Please provide a reference Genbank file with the parameter `ref_gbk`."
    assert file(params.ref_gbk).isFile() 
    : "Please double check that the reference Genbank file provided with the parameter `ref_gbk` exists."
    
    // make sure required snpeff config is provided and exists
    assert params.snpEff_config != "" 
    : "Please provide a snpEff config file with the parameter `snpEff_config`."
    assert file(params.snpEff_config).isFile() 
    : "Please double check that the snpEff config file provided with the parameter `snpEff_config` exists."

    // initialize input channels
    ch_primer_bed = Channel
        .fromPath( params.primer_bed )

    ch_refseq = Channel
        .fromPath( params.refseq )
    
    ch_ref_gbk = Channel
        .fromPath( params.ref_gbk )
    
    ch_snpeff_config = Channel
        .fromPath( params.snpEff_config )

    if ( params.platform == "ont" ) {

        NANOPORE (
            ch_primer_bed,
            ch_refseq,
            ch_ref_gbk,
            ch_snpeff_config
        )

    } // else {

    //     error(message = "Execution with Illumina reads is not yet supported.")

    //     assert params.illumina_fastq_dir != "" : "Please double check that a directory of Illumina FASTQs or Nanopore POD5s is provided." 
    //     assert file( params.illumina_fastq_dir ).isDirectory() : "The provided Illumina FASTQ directory ${params.illumina_fastq_dir} does not exist."

    //     ILLUMINA (
    //         ch_primer_bed,
    //         ch_refseq
    //     )
    // }

}
