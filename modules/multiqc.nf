/*
 * MULTIQC - Aggregate QC reports into a single HTML report
 *
 * This process aggregates FastQC reports and optionally OneRoof custom
 * content files into a unified MultiQC report.
 *
 * Inputs:
 *   - fastqc_files: FastQC output files (required)
 *   - oneroof_multiqc_files: OneRoof custom content TSV files (optional, use [] for none)
 *   - multiqc_config: MultiQC config YAML (optional, use file('NO_CONFIG') for none)
 *
 * Example usage:
 *   // With OneRoof files and config:
 *   MULTIQC(fastqc_ch.collect(), oneroof_ch.collect(), config_ch)
 *
 *   // Without OneRoof files (backward compatible):
 *   MULTIQC(fastqc_ch.collect(), Channel.of([]), Channel.of(file('NO_CONFIG')))
 */
process MULTIQC {

    publishDir params.qc, mode: 'copy', overwrite: true

    input:
    path fastqc_files, stageAs: "fastqc/*"
    path oneroof_multiqc_files, stageAs: "oneroof/*"
    path multiqc_config

    output:
    path "oneroof_multiqc_report.html", emit: report
    path "oneroof_multiqc_report_data", emit: data, optional: true

    script:
    def config_arg = multiqc_config.name != 'NO_CONFIG' ? "--config ${multiqc_config}" : ""
    """
    multiqc . ${config_arg} --filename oneroof_multiqc_report
    """

}