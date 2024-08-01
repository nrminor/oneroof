process FILTER_WITH_CHOPPER {

    /* */

    tag "${label}"
	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	cpus 4

    input:
    tuple val(label), path(fastq)

    output:
    tuple val(label), path("${new_id}.filtered.fastq.gz")

    script:
    new_id = file(fastq).getName().replace(".fastq.gz", "")
    """
    gunzip -c ${fastq} \
    | chopper \
    --maxlength ${params.max_len} \
    --minlength ${params.min_len} \
    --quality ${params.min_qual} \
    --threads ${task.cpus} \
    | bgzip -c > ${new_id}.filtered.fastq.gz
    """
}
