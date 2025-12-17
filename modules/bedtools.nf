process BEDTOOLS_GENOMECOV {

    /* */

    tag "${sample_id}"
    publishDir params.genomecov, mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    tuple val(sample_id), path(bam), path(index)

    output:
    tuple val(sample_id), path("${sample_id}.per-base.bed")

    script:
    """
    bedtools genomecov -bga -ibam ${bam} > ${sample_id}.per-base.bed
    """

}
