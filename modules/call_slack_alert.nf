process CALL_SLACK_ALERT{

    input: 
    val alignment

    script:
    // TODO: Replace experiment number flag with `--run_label` and make it a string, and then
    // use groovy/nextflow to parse that out here
    """
     slack_alerts.py \
    --exp_num ${workflow.launchDir} \
    --input_tsv_dir ${params.results} \
    --platform ${params.platform} \
    --depth ${params.min_depth_coverage}
    """

}
