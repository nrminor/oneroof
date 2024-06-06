process SPLIT_PRIMER_COMBOS {

    /*
    */

	errorStrategy { task.attempt < 3 ? 'retry' : params.errorMode }
	maxRetries 2

	input:
	path all_primer_combos

	output:
	path "*.bed"

	script:
	"""
	split_primer_combos.py -i ${all_primer_combos}
	"""

}