process TRIM_ENDS_TO_PRIMERS {

    /* */

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	cpus 3

    input:
    tuple val(barcode), path(patterns_file), path(untrimmed)

    output:
    tuple val(barcode), path("${barcode}.trimmed.amplicons.fastq.gz")

    script:
    patterns = file(patterns_file).readLines()
    assert patterns.size() == 2 : "Too many patterns provided in ${patterns_file}"
    forward_pattern = patterns[0]
    reverse_pattern = patterns[1]
    """
    cutadapt \
    -j ${task.cpus} \
    -a ${reverse_pattern} \
    ${untrimmed} \
    -o tmp.fq && \
    cutadapt \
    -j ${task.cpus} \
    -g ${forward_pattern} \
    tmp.fq \
    -o ${barcode}.trimmed.amplicons.fastq.gz && \
    rm tmp.fq
    """

}