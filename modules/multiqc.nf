process MULTIQC {

    publishDir params.qc, mode: 'copy', overwrite: true

    input:
	path fastqc_files, stageAs: "?/*"

	output:
	path "*.html"

    script:
    """
    multiqc .
    """

}