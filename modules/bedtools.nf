process GET_PRIMER_SEQS {

    /* */

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	input:
	path bed
	each path(refseq)

	output:
	path "${primer_combo}.fasta"

	script:
	primer_combo = file(bed.toString()).getSimpleName()
	"""
	bedtools getfasta -fi ${refseq} -bed ${bed} > ${primer_combo}.fasta
	
	# pull out the correct regex patterns for the amplicon
	make_primer_patterns.py ${primer_combo}.fasta ${primer_combo}
	"""

}