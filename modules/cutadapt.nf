process TRIM_ENDS_TO_PRIMERS {

    /* */

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    cpus 3

    input:
    tuple val(barcode), path(patterns_file), path(untrimmed)

    output:
    tuple val(barcode), path("${barcode}*.trimmed.fastq.gz")

    script:
    amplicon = file(patterns_file).getSimpleName()
    """
    FORWARD_PATTERN=\$(head -n 1 ${patterns_file})
    REVERSE_PATTERN=\$(tail -n 1 ${patterns_file})

    echo "Removing the forward primer: \$FORWARD_PATTERN"
    cutadapt \
    ${untrimmed} \
    -j ${task.cpus} \
    --front \$FORWARD_PATTERN \
    --revcomp \
    --match-read-wildcards \
    --errors ${params.max_mismatch} \
    --untrimmed-output ${barcode}_${amplicon}_fwd.untrimmed.fastq.gz \
    -o tmp.fq

    echo "Removing the reverse primer: \$REVERSE_PATTERN"
    cutadapt \
    tmp.fq \
    -j ${task.cpus} \
    --adapter \$REVERSE_PATTERN \
    --revcomp \
    --match-read-wildcards \
    --errors ${params.max_mismatch} \
    --minimum-length ${params.min_len} \
    --maximum-length ${params.max_len} \
    --untrimmed-output ${barcode}_${amplicon}_rev.untrimmed.fastq.gz \
    -o ${barcode}.${amplicon}.trimmed.fastq.gz && \
    rm tmp.fq
    """

}
