process CONCAT {

    /* */

    publishDir params.consensus, mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    path consensus_fastas

    output:
    path "*.fasta"

    script:
    """
    concat_consensus.py
    """
}