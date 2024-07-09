process PLOT_COVERAGE {

    /* */

    tag "${sample_id}"
    publishDir "${params.cov_plots}/${sample_id}", mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    tuple val(sample_id), path(mosdepth_files)

    output:
    path "${sample_id}.*.pdf", emit: plots
    path "${sample_id}*.tsv", emit: passing_cov

    script:
    """
    plot_coverage.py \
    --label ${sample_id} \
    --input ${sample_id}.per-base.bed \
    --depth ${params.min_depth_coverage}
    """

}
