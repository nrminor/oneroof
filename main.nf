#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { NANOPORE } from "$projectDir/workflows/nanopore"
include { ILLUMINA } from "$projectDir/workflows/illumina"

workflow {

    // the sequencing platform used
    params.platform = params.illumina_fastq_dir == "" ? "ont" : "illumina"

    // Checking for required files
    // ---------------------------------------------------------------------- //

    // make sure required primer bed is provided and exists
    assert params.primer_bed != "" :
    "Please provide a primer bed file with the parameter `primer_bed`."
    assert file(params.primer_bed).isFile() :
    "Please double check that the primer bed file provided with the parameter `primer_bed` exists."
    
    // make sure provided refseq is provided and exists
    assert params.refseq != "" :
    "Please provide a reference FASTA file with the parameter `refseq`."
    assert file(params.refseq).isFile() :
    "Please double check that the reference FASTA file provided with the parameter `refseq` exists."
    
    // make sure required reference genbank is provided and exists
    assert params.ref_gbk != "" :
    "Please provide a reference Genbank file with the parameter `ref_gbk`."
    assert file(params.ref_gbk).isFile() :
    "Please double check that the reference Genbank file provided with the parameter `ref_gbk` exists."
    
    // make sure required snpeff config is provided and exists
    assert params.snpEff_config != "" :
    "Please provide a snpEff config file with the parameter `snpEff_config`."
    assert file(params.snpEff_config).isFile() :
    "Please double check that the snpEff config file provided with the parameter `snpEff_config` exists."

    // ---------------------------------------------------------------------- //

    // initialize input channels
    ch_primer_bed = Channel
        .fromPath( params.primer_bed )

    ch_refseq = Channel
        .fromPath( params.refseq )
    
    ch_ref_gbk = Channel
        .fromPath( params.ref_gbk )
    
    ch_snpeff_config = Channel
        .fromPath( params.snpEff_config )


    // decide whether to run the ont or the illumina workflow
    if ( params.platform == "ont" ) {

        NANOPORE (
            ch_primer_bed,
            ch_refseq,
            ch_ref_gbk,
            ch_snpeff_config
        )

    }  else if ( params.platform == "illumina" ) {

        error "Execution with Illumina reads is not yet supported."

        ILLUMINA (
            ch_primer_bed,
            ch_refseq
        )

    } else {

        error """
        Unrecognized platform provided with ${params.platform}. This pipeline only supports
        Nanopore data with the keyword 'ont' and Illumina data with the keyword 'illumina'.
        """

    }

}
