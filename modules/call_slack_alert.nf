process CALL_SLACK_ALERT{

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    input: 
        val alignment
        val consensus
        val variants


    script:
    // TODO: Replace experiment number flag with `--run_label` and make it a string, and then
    // use groovy/nextflow to parse that out here
    def run = workflow.launchDir.getName()
    """
     slack_alerts.py \
    --run_label ${run} \
    --input_tsv_dir ${alignment} \
    --depth ${params.min_depth_coverage}
    """

}
