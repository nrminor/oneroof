#!/usr/bin/env nextflow

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


def helpMessage() {

    log.info """
    $frontMatter

    Usage:

    The typical command for running `oneroof` is as follows:
    nextflow run . --prepped_data inputs/ --primer_bed primers.bed --refseq GENOME.fasta --ref_gbk GENOME.gbk -profile containerless

    Mandatory arguments:
    --refseq                       The reference sequence to be used for mapping in FASTA format.

    Optional arguments:
    --primer_bed                   A bed file of primer coordinates relative to the reference provided with the parameters 'refseq' and 'ref_gbk'.
    --sample_lookup                Lookup JSON file where the key is the barcode or sequencer ID and the value is the desired sample ID.
    --ref_gbk                      The reference sequence to be used for variant annotation in Genbank format.
    --fwd_suffix                   Suffix in the primer bed file denoting whether a primer is forward. Default: '_LEFT'
    --rev_suffix                   Suffix in the primer bed file denoting whether a primer is reverse. Default: '_RIGHT'
    --remote_pod5_location         A remote location to use with an SSH client to watch for pod5 files in real-time as they are generated. Default: ''
    --file_watcher_config          Configuration file for remote file monitoring. Default: ''
    --pod5_staging                 Directory where pod5 files are cached as they arrive from the remote location. Default: '$launchDir/pod5_cache'
    --pod5_dir                     Directory where pod5 files are manually transferred if no remote pod5 location is given. Default: ''
    --precalled_staging            Directory to watch for Nanopore FASTQs or BAMs as they become available. Default: ''
    --prepped_data                 Location of prepped data if pod5 files are already basecalled and demultiplexed. Default: ''
    --illumina_fastq_dir           Location of paired-end Illumina FASTQ files to be processed. Default: ''
    --model                        Nanopore basecalling model to apply to the provided pod5 data. Default: 'sup@latest'
    --model_cache                  Directory to cache basecalling models locally. Default: '$launchDir/work/basecalling_models'
    --kit                          Nanopore barcoding kit used to prepare sequencing libraries. Default: null
    --pod5_batch_size              How many pod5 files to basecall at once. Default: null
    --basecall_max                 Number of parallel instances of the basecaller to run at once. Default: 1
    --max_len                      Maximum acceptable read length. Default: 12345678910
    --min_len                      Minimum acceptable read length. Default: 1
    --min_qual                     Minimum acceptable average quality for a given read. Default: 20
    --secondary                    Enable secondary alignments for each amplicon. Default: null
    --max_mismatch                 Maximum number of mismatches allowed when finding primers. Default: 0
    --downsample_to                Desired coverage to downsample to. Default: 0 (no downsampling)
    --min_consensus_freq           Minimum frequency of a variant base to be included in a consensus sequence. Default: 0.5
    --min_haplo_reads              Minimum read support to report an amplicon-haplotype. Default: 2
    --snpeff_cache                 Directory to cache a custom snpEff database. Default: '$launchDir/work/snpEff_cache'
    --min_depth_coverage           Minimum depth of coverage. Default: 20
    --nextclade_dataset            Nextclade dataset location. Default: null
    --nextclade_cache              Directory to cache Nextclade datasets. Default: '$launchDir/work/nextclade_datasets'
    --results                      Where to place the results. Default: '$launchDir/results'
    --cleanup                      Whether to clean up the work directory after a successful run. Default: null

    Advanced parameters:
    --snpEff_config                Configuration file for snpEff. Default: '$projectDir/conf/snpeff.config'
    """
    .stripIndent()
}

if (params.help) {
    helpMessage()
    exit 0
}

nextflow.enable.dsl = 2

include { NANOPORE } from "$projectDir/workflows/nanopore"
include { ILLUMINA } from "$projectDir/workflows/illumina"

log.info    """
            $frontMatter

            Workflow settings:
            -------------------------------------------------------------------------
            Launch directory            : ${workflow.launchDir}
            Workflow files directory    : ${workflow.projectDir}
            Run start time              : ${workflow.start}
            Run command                 : ${workflow.commandLine}
            Profile                     : ${workflow.profile}
            Nextflow version            : ${nextflow.version}
            Cleanup mode                : ${params.cleanup ?: ""}

            User-provided settings:
            -------------------------------------------------------------------------
            Sequencing Platform         : ${params.platform}
            Basecalling model           : ${params.model ?: ""}
            Nanopore barcoding kit      : ${params.kit ?: ""}
            Minimum read length         : ${params.min_len}
            Maximum read length         : ${params.max_len}
            Minimum avg. read quality   : ${params.min_qual}
            Permitted primer mismatches : ${params.max_mismatch}
            Desired coverage            : ${params.downsample_to != 0 ? params.downsample_to + "X" : "No downsampling"}
            Secondary alignments        : ${params.secondary ? "on" : "off"}
            Minimum coverage            : ${params.min_depth_coverage}X
            Minimum consensus frequency : ${params.min_consensus_freq}
            NextClade Dataset           : ${params.nextclade_dataset ?: ""}

            User-provided inputs and outputs:
            -------------------------------------------------------------------------
            Primer BED file             : ${params.primer_bed ?: ""}
            Reference sequence FASTA    : ${params.refseq}
            Reference sequence GBK      : ${params.ref_gbk ?: ""}
            File watcher config         : ${params.file_watcher_config ?: ""}
            Remote POD5 directory       : ${params.remote_pod5_location ?: ""}
            Local POD5 directory        : ${params.pod5_dir ?: ""}
            Pre-basecalled directory    : ${params.precalled_staging ?: params.prepped_data}
            Sample ID lookup            : ${params.sample_lookup ?: ""}
            Results directory           : ${params.results}

            """
            .stripIndent()

workflow {

    // the sequencing platform used
    params.platform = params.illumina_fastq_dir == "" ? "ont" : "illumina"

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

    ch_snpeff_config = params.snpEff_config ?
        Channel.fromPath( params.snpEff_config ) :
        Channel.empty()

    ch_sample_lookup = params.sample_lookup ?
        Channel.fromPath( params.sample_lookup ) :
        Channel.empty()

    // decide whether to run the ont or the illumina workflow
    if ( params.platform == "ont" ) {

        NANOPORE (
            ch_primer_bed,
            ch_refseq,
            ch_ref_gbk,
            ch_snpeff_config,
            ch_sample_lookup
        )

    }  else if ( params.platform == "illumina" ) {

        ILLUMINA (
            ch_primer_bed,
            ch_refseq,
            ch_ref_gbk,
            ch_snpeff_config,
            ch_sample_lookup
        )

    } else {

        error """
        Unrecognized platform provided with ${params.platform}. This pipeline only supports
        Nanopore data with the keyword 'ont' and Illumina data with the keyword 'illumina'.
        """

    }

}
