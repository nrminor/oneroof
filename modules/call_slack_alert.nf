process CALL_SLACK_ALERT{

    input: 
        val(alignment)

    script:
    """
     slack_alerts.py \
        --exp_num ${workflow.launchDir} \
        --input_tsv_dir ${params.results} \
        --platform ${params.platform} \
        --depth ${params.min_depth_coverage}
    """

}