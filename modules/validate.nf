process VALIDATE_NANOPORE {

    /* */

    tag "${label}"

	errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
	maxRetries 1

    cpus 3

    input:
    tuple val(label), path(seq_file)

    output:
    tuple val(label), path("${label}.validated.fastq.gz"), val("success!")

    script:
    if ( file(seq_file).getName().endsWith(".fastq.gz") )
        """
        seqkit seq --validate-seq ${seq_file} -o ${label}.validated.fastq.gz
        """
    else if ( file(seq_file).getName().endsWith(".bam") )
        """
        samtools quickcheck -v -u ${seq_file} && \
        samtools fastq ${seq_file} | bgzip -o ${label}.validated.fastq.gz
        """
    else
        error "Unrecognized file format provided."

}

process VALIDATE_ILLUMINA {

    /* */

    tag "${label}"

	errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
	maxRetries 1

    cpus 3

    input:
    tuple val(label), path(reads1), path(reads2)

    output:
    tuple val(label), path(reads1), path(reads2), path("${barcode}.report.txt")

    script:
    """
    seqfu check --deep --verbose --thousands \
    ${reads1} ${reads2} \
    > ${barcode}.report.txt
    """

}
