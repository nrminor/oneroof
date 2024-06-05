process RESPLICE_PRIMERS {

    /*
    */

	// tag "${params.desired_amplicon}"
    // label "general"
	// publishDir params.results, mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : params.errorMode }
	maxRetries 2

    input:
    path bed_file

    output:
    path "respliced.bed"

    script:
    """
	resplice_primers.py -i ${bed_file}
    """

}