process PLOT_COVERAGE {

    /* */

    tag "${sample_id}"
    publishDir params.cov_plots, mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    tuple val(sample_id), path(mosdepth_files)

    output:
    path "${sample_id}.coverage.pdf"

    script:
    """
    cat ${sample_id}.per-base.bed.gz | gzip -d -c > ${sample_id}.per-base.bed.gz && \
    plot_coverage.py \
    --label ${sample_id} \
    --input ${sample_id}.per-base.bed
    """

}
