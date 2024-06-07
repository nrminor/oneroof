process MULTIQC {

    publishDir params.qc, mode: 'copy', mode: true

    input:
	path fastqc_files, stageAs: "?/*"

	output:
	path "*.html"

    script:
    """
    multiqc .
    """

}