process REPORT_REFERENCES {

    /* */

    publishDir "${params.results}/reference_assets", mode: "copy", overwrite: true

	errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
	maxRetries 1

    cpus 3

    input:
    path files

    shell:
    '''
    if [ ! -d !{params.results}/reference_assets ]; then
        mkdir -p !{params.results}/reference_assets
    fi

    find . -name '*.fa*' -o -name '*.bed' -o -name '*.g*' > ref_files.txt

    for file in "$(cat ref_files.txt)"; do
        cp `realpath ${file}` !{params.results}/reference_assets/
    done
    '''
}

process PUBLISH_COMMAND {

    /* */

    publishDir params.run_command, mode: "copy", overwrite: true

	errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
	maxRetries 1

    cpus 1

    output:
    path "run_command.sh"

    script:
    command = workflow.commandLine
    """
    echo "${command}" > run_command.sh
    """

}

process EXTRACT_COVERAGE_METRICS {

    /*
     * Extract coverage metrics from bedtools genomecov output.
     * Produces a JSON file with coverage statistics for downstream reporting.
     */

    tag "${sample_id}"

    errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
    maxRetries 1

    input:
    tuple val(sample_id), path(coverage_bed)

    output:
    tuple val(sample_id), path("${sample_id}_coverage_metrics.json")

    script:
    """
    extract_metrics.py coverage \\
        --sample-id ${sample_id} \\
        --bed ${coverage_bed} \\
        --output ${sample_id}_coverage_metrics.json
    """
}

process EXTRACT_ALIGNMENT_METRICS {

    /*
     * Extract alignment metrics from BAM files.
     * Produces a JSON file with read counts and mapping statistics.
     */

    tag "${sample_id}"

    errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
    maxRetries 1

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}_alignment_metrics.json")

    script:
    """
    extract_metrics.py alignment \\
        --sample-id ${sample_id} \\
        --bam ${bam} \\
        --output ${sample_id}_alignment_metrics.json
    """
}

process EXTRACT_VARIANT_METRICS {

    /*
     * Extract variant metrics from SnpSift effects TSV.
     * Produces a JSON file with variant counts and effect statistics.
     */

    tag "${sample_id}"

    errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
    maxRetries 1

    input:
    tuple val(sample_id), path(vcf), path(effects_tsv)

    output:
    tuple val(sample_id), path("${sample_id}_variant_metrics.json")

    script:
    """
    extract_metrics.py variants \\
        --sample-id ${sample_id} \\
        --effects-tsv ${effects_tsv} \\
        --output ${sample_id}_variant_metrics.json
    """
}

process EXTRACT_CONSENSUS_METRICS {

    /*
     * Extract consensus metrics from FASTA files.
     * Produces a JSON file with sequence quality statistics.
     */

    tag "${sample_id}"

    errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
    maxRetries 1

    input:
    tuple val(sample_id), path(fasta)

    output:
    tuple val(sample_id), path("${sample_id}_consensus_metrics.json")

    script:
    """
    extract_metrics.py consensus \\
        --sample-id ${sample_id} \\
        --fasta ${fasta} \\
        --output ${sample_id}_consensus_metrics.json
    """
}

process EXTRACT_METAGENOMICS_METRICS {

    /*
     * Extract metagenomics metrics from Sylph profile output.
     * Produces a JSON file with taxonomic profiling results.
     */

    tag "${sample_id}"

    errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
    maxRetries 1

    input:
    tuple val(sample_id), path(profile_tsv)

    output:
    tuple val(sample_id), path("${sample_id}_metagenomics_metrics.json")

    script:
    """
    extract_metrics.py metagenomics \\
        --sample-id ${sample_id} \\
        --profile-tsv ${profile_tsv} \\
        --output ${sample_id}_metagenomics_metrics.json
    """
}

process EXTRACT_HAPLOTYPING_METRICS {

    /*
     * Extract haplotyping metrics from Devider output directory.
     * Produces a JSON file with haplotype phasing results.
     * Nanopore only.
     */

    tag "${devider_dir.baseName}"

    errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
    maxRetries 1

    input:
    path devider_dir

    output:
    tuple val(sample_id), path("${sample_id}_haplotyping_metrics.json")

    script:
    // Extract sample_id from directory name (e.g., "sample1_devider" -> "sample1")
    sample_id = devider_dir.baseName.replaceAll(/_devider$/, '')
    """
    extract_metrics.py haplotyping \\
        --sample-id ${sample_id} \\
        --devider-dir ${devider_dir} \\
        --output ${sample_id}_haplotyping_metrics.json
    """
}

process ASSEMBLE_REPORT {

    /*
     * Assemble the final OneRoof report from collected metrics.
     * Produces JSON reports, MultiQC custom content files, and visualizations.
     *
     * The amplicon_summary input is optional - pass a placeholder file named
     * 'NO_AMPLICON_SUMMARY' when primers are not provided.
     */

    publishDir "${params.report_dir}", mode: 'copy', overwrite: true

    errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
    maxRetries 1

    input:
    path metrics_files, stageAs: "metrics/*"
    val platform
    val reference_name
    val qc_thresholds
    path multiqc_config_template
    path amplicon_summary

    output:
    path "oneroof_report.json", emit: full_report
    path "oneroof_report_summary.json", emit: summary_report
    path "multiqc/*.tsv", emit: multiqc_tsv
    path "multiqc/multiqc_config.yaml", emit: multiqc_config
    path "visualizations/*", emit: visualizations, optional: true

    script:
    def thresholds_json = groovy.json.JsonOutput.toJson(qc_thresholds)
    def config_arg = multiqc_config_template.name != 'NO_CONFIG' ? "--config-template ${multiqc_config_template}" : ""
    def amplicon_arg = amplicon_summary.name != 'NO_AMPLICON_SUMMARY' ? "--amplicon-summary ${amplicon_summary}" : ""
    """
    mkdir -p multiqc visualizations

    assemble_report.py assemble \\
        --metrics-dir metrics \\
        --output-dir . \\
        --multiqc-dir multiqc \\
        --viz-dir visualizations \\
        --platform ${platform} \\
        --reference-name "${reference_name}" \\
        --qc-thresholds '${thresholds_json}' \\
        ${config_arg} \\
        ${amplicon_arg}
    """
}
