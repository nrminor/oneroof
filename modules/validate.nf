process VALIDATE_NANOPORE {

    /* */

    tag "${label}"

	errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
	maxRetries 1

    cpus 3

    input:
    tuple val(label), path(seq_file)

    output:
    tuple val(label), path("${label}.validated.bam"), val("success!")

    script:
    if ( file(seq_file).getName().endsWith(".fastq.gz") )
        """
        seqkit seq -validate-seq ${seq_file} &
        samtools import -0 ${seq_file} -o ${label}.validated.bam
        wait
        """
    else if ( file(seq_file).getName().endsWith(".bam") )
        """
        samtools quickcheck -v -u ${seq_file} && \
        cp `realpath ${seq_file}` ${label}.validated.bam
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
