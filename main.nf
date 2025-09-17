#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { NANOPORE } from "./workflows/nanopore"
include { ILLUMINA } from "./workflows/illumina"

// the sequencing platform used
params.platform = params.illumina_fastq_dir == "" ? "ont" : "illumina"

workflow {

    frontMatter = """
                                                                .8888b
                                                                88   "
        .d8888b. 88d888b. .d8888b. 88d888b. .d8888b. .d8888b. 88aaa
        88'  `88 88'  `88 88ooood8 88'  `88 88'  `88 88'  `88 88
        88.  .88 88    88 88.  ... 88       88.  .88 88.  .88 88
        `88888P' dP    dP `88888P' dP       `88888P' `88888P' dP
        ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

        oneroof: Base-, Variant-, and Consensus-calling under One Proverbial Roof
        =========================================================================
        `oneroof` is a bespoke bioinformatic pipeline that can handle Oxford
        Nanopore POD5 basecalling, Illumina paired-end read-merging, read
        alignment and variant-calling, variant-effect annotation, consensus
        sequence calling, quality reporting, and phylogenetic tree-building, all
        under "one roof."
        (version 0.1.0)
        =========================================================================
        """
        .stripIndent()

    if (params.help) {
        Utils.helpMessage(workflow, log, frontMatter)
        exit 0
    }

    log.info frontMatter
    Utils.workflowDisplay(params, workflow, log, nextflow)

    // Checking for required files
    // ---------------------------------------------------------------------- //

    // make sure a reference sequence FASTA, the minimum pipeline dependency, is
    // provided and exists
    assert params.refseq :
    "Please provide a reference FASTA file with the parameter `refseq`."
    assert file(params.refseq).isFile() :
    "Please double check that the reference FASTA file provided with the parameter `refseq` exists."

    // ---------------------------------------------------------------------- //

    // initialize input channels
    ch_primer_bed = params.primer_bed ?
        Channel.fromPath( params.primer_bed ) :
        Channel.empty()

    ch_refseq = Channel
        .fromPath( params.refseq )

    ch_ref_gbk = params.ref_gbk ?
        Channel.fromPath( params.ref_gbk ) :
        Channel.empty()

    ch_contam_fasta = params.contam_fasta && file(params.contam_fasta).isFile()
        ? Channel.fromPath( params.contam_fasta )
        : Channel.empty()

    ch_metagenomics_ref = params.meta_ref
        ? file(params.meta_ref).isFile()
        ? Channel.fromPath( params.meta_ref )
        : Channel.from( params.meta_ref )
        : Channel.empty()

    ch_snpeff_config = params.snpEff_config ?
        Channel.fromPath( params.snpEff_config ) :
        Channel.empty()

    ch_primer_tsv = params.primer_tsv ?
        Channel.fromPath ( params.primer_tsv ) :
        Channel.empty()

    ch_sylph_tax_db = params.sylph_tax_db ?
        Channel.fromPath ( params.sylph_tax_db ) :
        Channel.empty()

    ch_meta_ref_link = params.sylph_db_link ?
        Channel.value ( params.sylph_db_link ) :
        Channel.empty()

    // decide whether to run the ont or the illumina workflowls

    if ( params.platform == "ont" ) {

        NANOPORE (
            ch_primer_bed,
            ch_refseq,
            ch_ref_gbk,
            ch_contam_fasta,
            ch_snpeff_config,
            ch_metagenomics_ref,
            ch_primer_tsv,
            ch_sylph_tax_db,
            ch_meta_ref_link
        )

    }  else if ( params.platform == "illumina" ) {

        ILLUMINA (
            ch_primer_bed,
            ch_refseq,
            ch_ref_gbk,
            ch_contam_fasta,
            ch_snpeff_config,
            ch_metagenomics_ref,
            ch_primer_tsv,
            ch_sylph_tax_db,
            ch_meta_ref_link
        )

    } else {

        error """
        Unrecognized platform provided with ${params.platform}. This pipeline only supports
        Nanopore data with the keyword 'ont' and Illumina data with the keyword 'illumina'.
        """

    }

    if ( params.email ) {
        workflow.onComplete {

            def msg = """\
                Oneroof has finished running with the following settings:

                ${Utils.workflowDisplay(params, workflow, log, nextflow)}
                """
                .stripIndent()

            sendMail(to: params.email, subject: 'Oneroof Execution Report', body: msg)
        }
    }

}
