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
        # fq lint \
        # --lint-mode panic \
        # --single-read-validation-level high \
        # ${seq_file} && ln -s \$(realpath ${seq_file}) ${label}.validated.fastq.gz

        cat ${seq_file} | seqkit seq --validate-seq --out-file ${label}.validated.fastq.gz
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

    // errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
    // maxRetries 1

    cpus 3

    input:
    tuple val(label), path(reads1), path(reads2)

    output:
    tuple val(label), path(reads1), path(reads2), path("${label}.report.txt")

    script:
    """
    # fq lint \\
    # --lint-mode panic \\
    # --single-read-validation-level high \\
    # --paired-read-validation-level high \\
    # ${reads1} ${reads2} \\
    # > ${label}.report.txt

    reformat.sh in=${reads1} in2=${reads2} vpair 2> ${label}.report.txt
    """

}
