process FIND_COMPLETE_AMPLICONS {

    /* */

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	cpus 3

    input:
	tuple path(reads), path(patterns)

    output:
    tuple val(barcode), path(patterns), path("${barcode}_amplicons.fastq.gz")

    script:
	barcode = file(reads).getSimpleName()
    """
	cat ${reads} | \
    seqkit grep \
	--threads ${task.cpus} \
	--max-mismatch ${params.max_mismatch} \
	--by-seq \
	--pattern-file ${patterns} \
	-o ${barcode}_amplicons.fastq.gz
    """

}

process AMPLICON_STATS {

    /* */

	tag "${barcode}"
	publishDir params.complete_amplicons, mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	cpus 3

	input:
	tuple val(barcode), path("amplicons/*")

    output:
    path "${barcode}.per_amplicon_stats.tsv"

    script:
    """
    seqkit stats \
	--threads ${task.cpus} \
	--all --basename --tabular \
	amplicons/*.fastq.gz > ${barcode}.per_amplicon_stats.tsv
    """

}

process MERGE_BY_SAMPLE {

    /* */

	tag "${barcode}"
	publishDir params.complete_amplicons, mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	cpus 3

	input:
	tuple val(barcode), path("fastqs/*")

	output:
	tuple val(barcode), path("${barcode}.amplicons.fastq.gz")

	script:
	"""
	seqkit scat \
	--find-only \
	--threads ${task.cpus} \
	fastqs/ \
	| bgzip -o ${barcode}.amplicons.fastq.gz
	"""
}

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

    seqkit amplicon \
    -f -r 1:-1 \
    --forward \$FORWARD_PATTERN \
    --reverse \$REVERSE_PATTERN \
    --max-mismatch ${params.max_mismatch} \
    --strict-mode \
    --threads ${task.cpus} \
    --out-file ${barcode}.${amplicon}.trimmed.fastq.gz \
    ${untrimmed}
    """

}
