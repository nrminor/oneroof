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

process MULTI_SAMPLE_PLOT {

    /* */

    publishDir params.cov_plots, mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    path "inputs/*"

    output:
    path "*.pdf"

    script:
    if (params.log)
        """
        multisample_plot.py \
        --input_dir inputs \
        --log
        """
    else
        """
        multisample_plot.py \
        --input_dir inputs \
        """

}

process COVERAGE_SUMMARY {

    /* */

    publishDir params.alignment, mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    path collected_passing_cov

    output:
    path "coverage_summary.tsv", emit: coverage_summary

    script:
    """
    awk -F"\t" 'BEGIN {print "sample id\tproportion ≥ specified depth"} \$1!="sample id"{print \$1 "\t" \$4;}' *.tsv > coverage_summary.tsv
    """

}