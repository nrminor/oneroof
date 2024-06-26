process TRIM_ENDS_TO_PRIMERS {

    /* */

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	cpus 3

    input:
    tuple val(barcode), path(patterns_file), path(untrimmed)

    output:
    tuple val(barcode), path("${barcode}*.trimmed.amplicons.fastq.gz")

    script:
    amplicon = file(patterns_file).getSimpleName()
    """
    FORWARD_PATTERN=\$(head -n 1 ${patterns_file})
    REVERSE_PATTERN=\$(tail -n 1 ${patterns_file})
    cutadapt \
    -j ${task.cpus} \
    -a \$FORWARD_PATTERN \
    ${untrimmed} \
    -o tmp.fq && \
    cutadapt \
    -j ${task.cpus} \
    -g \$REVERSE_PATTERN \
    tmp.fq \
    -o ${barcode}.${amplicon}.trimmed.fastq.gz && \
    rm tmp.fq
    """

}
