process MULTIQC {

    publishDir params.qc, mode: 'copy', mode: true

    input:
	path fastqc_files

	output:
	path "*.html"

    script:
    """
    multiqc ${fastq_files}
    """

}