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

process ASSEMBLE_REPORT {

    /*
     * Assemble the final OneRoof report from collected metrics.
     * Produces JSON reports, MultiQC custom content files, and visualizations.
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

    output:
    path "oneroof_report.json", emit: full_report
    path "oneroof_report_summary.json", emit: summary_report
    path "multiqc/*", emit: multiqc_files
    path "visualizations/*", emit: visualizations, optional: true

    script:
    def thresholds_json = groovy.json.JsonOutput.toJson(qc_thresholds)
    def config_arg = multiqc_config_template.name != 'NO_CONFIG' ? "--config-template ${multiqc_config_template}" : ""
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
        ${config_arg}
    """
}
