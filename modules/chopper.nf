process FILTER_WITH_CHOPPER {

    /* */

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	cpus 4

    input:
    tuple val(label), path(fastq)

    output:
    tuple val(label), path("${label}.filtered.fastq.gz")

    script:
    """
    gunzip -c ${fastq} \
    | chopper \
    --maxlength ${params.max_len} \
    --minlength ${params.min_len} \
    --quality ${params.min_qual} \
    --threads ${task.cpus} \
    | gzip -c > ${label}.filtered.fastq.gz
    """
}
