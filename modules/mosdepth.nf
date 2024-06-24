process MOSDEPTH {

    /* */

    tag "${sample_id}"
    publishDir params.mosdepth, mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    cpus 3

    input:
    tuple val(sample_id), path(bam), path(index)

    output:
    tuple val(sample_id), path("${sample_id}*")

    script:
    """
    mosdepth \
    --threads ${task.cpus} \
    ${sample_id} \
    ${bam} && \
    gunzip ${sample_id}.per-base.bed.gz
    """

}
