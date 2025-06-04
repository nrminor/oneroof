process CALL_SLACK_ALERT{

    script:
    """
     slack_alerts.py \
        --exp_num ${workflow.launchDir} \
        --input_tsv "${params.results}/${params.platform}/03_alignments/coverage_summary.tsv" \
        --depth ${params.min_depth_coverage}
    """

}