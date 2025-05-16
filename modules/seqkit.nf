process COMPRESS_TO_SORTED_FASTA {

    tag "${barcode}"

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    input:
    tuple val(barcode), path(fastq_reads)

    output:
    tuple val(barcode), path("${barcode}.fasta.gz")

    script:
    """
    seqkit fq2fa ${fastq_reads} \
    | seqkit seq --only-id \
    | seqkit sort --two-pass -o "${barcode}.fasta.gz"
    """
}

process FIND_COMPLETE_AMPLICONS {

    /* */

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	cpus 3

    input:
	tuple path(reads), path(patterns)

    output:
    tuple val(barcode), path(patterns), path("${barcode}_amplicons.fasta.gz")

    script:
	barcode = file(reads).getSimpleName()
    """
	cat ${reads} | \
    seqkit grep \
	--threads ${task.cpus} \
	--max-mismatch ${params.max_mismatch} \
	--by-seq \
	--pattern-file ${patterns} \
	-o ${barcode}_amplicons.fasta.gz
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
    tuple val(barcode), path("${barcode}*.trimmed.fasta.gz")

    script:
    amplicon = file(patterns_file).getSimpleName()
    """
    FORWARD_PATTERN=\$(head -n 1 ${patterns_file})
    REVERSE_PATTERN=\$(tail -n 1 ${patterns_file})
    FORWARD_LENGTH=\${#FORWARD_PATTERN}
    REVERSE_LENGTH=\${#REVERSE_PATTERN}

    seqkit amplicon \
    --region \${FORWARD_LENGTH}:-\${REVERSE_LENGTH} \
    --forward \$FORWARD_PATTERN \
    --reverse \$REVERSE_PATTERN \
    --max-mismatch ${params.max_mismatch} \
    --strict-mode \
    --threads ${task.cpus} \
    --out-file ${barcode}.${amplicon}.trimmed.fasta.gz \
    ${untrimmed}
    """

}

process PER_AMPLICON_FILTERS {

    /* */

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    cpus 4

    input:
    tuple val(label), path(fasta)

    output:
    tuple val(label), path("${new_id}.filtered.fasta.gz")

    script:
    new_id = file(fasta).getName().replace(".fasta.gz", "")
    """
    seqkit seq \
    --max-len ${params.max_len} \
    --min-len ${params.min_len} \
    --min-qual ${params.min_qual} \
    --threads ${task.cpus} \
    -o ${new_id}.filtered.fasta.gz
    ${fasta}
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
	tuple val(barcode), path("fastas/*")

	output:
	tuple val(barcode), path("${barcode}.amplicons.fasta.gz")

	script:
	"""
	seqkit scat \
	--find-only \
	--threads ${task.cpus} \
	fastas/ \
	| bgzip -o ${barcode}.amplicons.fasta.gz
	"""
}
