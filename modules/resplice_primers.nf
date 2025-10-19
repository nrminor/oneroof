process RESPLICE_PRIMERS {

    /*
    */

    publishDir params.respliced, mode: 'copy', overwrite: true

    errorStrategy 'finish'

    input:
    path bed_file

    output:
    path "respliced.bed"

    script:
    """
    resplice_primers.py \
    --input_bed ${bed_file} \
    --fwd_suffix ${params.fwd_suffix} \
    --rev_suffix ${params.rev_suffix} \
    -vv
    """

}

process SPLIT_PRIMER_COMBOS {

    /*
    */

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
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
