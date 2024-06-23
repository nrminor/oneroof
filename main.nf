#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { NANOPORE } from "$projectDir/workflows/nanopore"
include { ILLUMINA } from "$projectDir/workflows/illumina"

log.info    """
                                                                  .8888b
                                                                  88   "
            .d8888b. 88d888b. .d8888b. 88d888b. .d8888b. .d8888b. 88aaa
            88'  `88 88'  `88 88ooood8 88'  `88 88'  `88 88'  `88 88
            88.  .88 88    88 88.  ... 88       88.  .88 88.  .88 88
            `88888P' dP    dP `88888P' dP       `88888P' `88888P' dP
            ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

            oneroof: Base-, Variant-, and Consensus-calling under One Proverbial Roof
            =========================================================================
            `oneroof` is a bespoke bioinformatic pipeline that can handle basecalling
            Oxford Nanopore POD5 basecalling, Illumina paired-end read-merging, read
            alignment and variant-calling, variant-effect annotation, consensus
            sequence calling, quality reporting, and phylogenetic tree-building, all
            under "one roof."
            (version 0.1.0)
            =========================================================================

            Workflow settings:
            -------------------------------------------------------------------------
            Launch directory            : ${workflow.launchDir}
            Workflow files directory    : ${workflow.projectDir}
            Run start time              : ${workflow.start}
            Run command                 : ${workflow.commandLine}
            Profile                     : ${workflow.profile}
            Nextflow version            : ${nextflow.version}
            Cleanup mode                : ${params.cleanup}

            User-provided settings:
            -------------------------------------------------------------------------
            Sequencing Platform         : ${params.platform}
            Basecalling model           : ${params.model}
            Nanopore barcoding kit      : ${params.kit}
            Permitted primer mismatches : ${params.max_mismatch}
            Desired coverage            : ${params.downsample_to}X
            Minimum coverage            : ${params.min_depth_coverage}X
            NextClade Dataset           : ${params.nextclade_dataset}

            User-provided inputs and outputs:
            -------------------------------------------------------------------------
            Primer BED file             : ${params.primer_bed}
            Reference sequence FASTA    : ${params.refseq}
            Reference sequence GBK      : ${params.ref_gbk}
            File watcher config         : ${params.file_watcher_config}
            Remote POD5 directory       : ${params.remote_pod5_location}
            Local POD5 directory        : ${params.pod5_dir}
            Pre-basecalled directory    : ${params.precalled_staging ?: params.prepped_data }
            Results directory           : ${params.results}

            """
            .stripIndent()

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
