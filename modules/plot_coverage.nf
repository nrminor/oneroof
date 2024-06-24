process PLOT_COVERAGE {

    /* */

    tag "${sample_id}"
    publishDir params.cramino, mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    tuple val(sample_id), path(mosdepth_files)

    output:
    path "*"

    script:
    """
    plot_coverage.py --label ${sample_id}
    """

}
