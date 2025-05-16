//
// This file holds a couple Groovy functions for counting amplicons and reference contigs
//

import java.nio.file.Files
import java.nio.file.Paths
import java.util.stream.Collectors

class Utils {

    public static Integer countFastaHeaders(String filePath) {
        def headerCount = 0
        def file = new File(filePath)

        if (!file.exists()) {
            throw new FileNotFoundException("The file $filePath does not exist.")
        }

        file.eachLine { line ->
            if (line.trim().startsWith('>')) {
                headerCount++
            }
        }

        return headerCount
    }

    public static Integer countAmplicons(String filePath) {
        // Read all lines from the file
        def lines = Files.readAllLines(Paths.get(filePath))

        // Extract the fourth column, replace hyphens, and remove _LEFT and _RIGHT
        def processedItems = lines.collect { line ->
            def columns = line.split('\t')
            if (columns.size() >= 4) {
                def item = columns[3]
                item = item.replace('-', '_')
                item = item.replace('_LEFT', '').replace('_RIGHT', '')
                return item
            }
            return null
        }.findAll { it != null }

        // Convert to a set to find unique items
        def uniqueItems = processedItems.toSet()

        // Return the number of unique items
        return uniqueItems.size()

    }

    public static void helpMessage(workflow, log, frontMatter) {

        log.info frontMatter
        log.info """
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
        --pod5_staging                 Directory where pod5 files are cached as they arrive from the remote location. Default: '${workflow.launchDir}/pod5_cache'
        --pod5_dir                     Directory where pod5 files are manually transferred if no remote pod5 location is given. Default: ''
        --precalled_staging            Directory to watch for Nanopore FASTQs or BAMs as they become available. Default: ''
        --prepped_data                 Location of prepped data if pod5 files are already basecalled and demultiplexed. Default: ''
        --illumina_fastq_dir           Location of paired-end Illumina FASTQ files to be processed. Default: ''
        --model                        Nanopore basecalling model to apply to the provided pod5 data. Default: 'sup@latest'
        --model_cache                  Directory to cache basecalling models locally. Default: '${workflow.launchDir}/work/basecalling_models'
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
        --snpeff_cache                 Directory to cache a custom snpEff database. Default: '${workflow.launchDir}/work/snpEff_cache'
        --min_depth_coverage           Minimum depth of coverage. Default: 10
        --nextclade_dataset            Nextclade dataset location. Default: null
        --nextclade_cache              Directory to cache Nextclade datasets. Default: '${workflow.launchDir}/work/nextclade_datasets'
        --results                      Where to place the results. Default: '${workflow.launchDir}/results'
        --email                        Email address to notify on completion. Multiple comma dilimeted emails may be provided.
        --cleanup                      Whether to clean up the work directory after a successful run. Default: null

        Advanced parameters:
        --snpEff_config                Configuration file for snpEff. Default: '${workflow.projectDir}/conf/snpeff.config'
        """
        .stripIndent()
    }

    public static void workflowDisplay(params, workflow, log, nextflow) {
        log.info """
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
                Illumina FASTQs directory   : ${params.illumina_fastq_dir ?: ""}
                Sample ID lookup            : ${params.sample_lookup ?: ""}
                Results directory           : ${params.results}
                Email Address(es)           : ${params.email ?: ""}

                """
                .stripIndent()
    }

    public static String reverseComplement(String seq) {
        def complementMap = [
            'A': 'T', 'T': 'A',
            'C': 'G', 'G': 'C',
            'a': 't', 't': 'a',
            'c': 'g', 'g': 'c',
            'N': 'N', 'n': 'n',
            'U': 'A', 'u': 'a',
        ]
        return seq.reverse().collect { base ->
            complementMap.get(base, 'N')  // default to 'N' if unknown base
        }.join()
    }

}
