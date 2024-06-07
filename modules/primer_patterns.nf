process GET_PRIMER_PATTERNS {

    /* */

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	input:
	path primer_fasta

	output:
	path "${primer_combo}.txt"

	script:
	primer_combo = file(primer_fasta.toString()).getSimpleName()
	"""
	make_primer_patterns.py ${primer_combo}.fasta ${primer_combo}
	"""

}