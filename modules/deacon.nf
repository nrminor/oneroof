process INDEX_CONTAMINANTS {

    /* */

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	cpus 1

    input:
    // tuple path(contam_fasta), val(sample_id)
    tuple path(contam_fasta)

    output:
    path "*.idx"

    script:

    def dbName = file(contam_fasta).getSimpleName()
    """
    deacon index build ${contam_fasta} > ${dbName}.idx
    """

}

process RUN_DEACON {

    /* */

	tag "${barcode}"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	cpus 1

    input:
	tuple val(barcode), path(reads), path(contam_index)

    output:
	tuple val(barcode), path("${barcode}.decontam.fastq.gz")

    script:
    def dbName = file(contam_index).getSimpleName()
    """
    zcat ${reads} \
    | deacon filter ${contam_index} --log ${dbName}.log.json \
    | pigz > ${barcode}.decontam.fastq.gz
    """

}
