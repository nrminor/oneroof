process FIND_COMPLETE_AMPLICONS {
    /* */

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	cpus 3

    input:
	tuple path(reads), path(patterns)

    output:
    tuple val(barcode), path(patterns), path("${barcode}_amplicons.fasta.gz")

    script:
	barcode = file(reads).getSimpleName()
    """
    grepq \
    --read-gzip \
    --write-gzip \
    --fasta \
    --fast \
    
    """
}
