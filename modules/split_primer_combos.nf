process SPLIT_PRIMER_COMBOS {

    /*
    */

	tag "${params.desired_amplicon}"
    label "general"

	errorStrategy { task.attempt < 3 ? 'retry' : params.errorMode }
	maxRetries 2

	input:
	path all_primer_combos

	output:
	path "${params.desired_amplicon}*.bed"

	script:
	"""
	split_primer_combos.py -i ${all_primer_combos}
	"""

}