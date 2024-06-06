process CONVERT_AND_SORT {

    tag "${barcode}"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    tuple val(barcode), path(sam)

    output:
    tuple val(barcode), path("${barcode}.bam")

    script:
    """
    samtools view -bS ${sam} \
    | samtools sort -o ${barcode}.bam
    """

}

process INDEX {

    tag "${barcode}"
    publishDir params.alignment, mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    tuple val(barcode), path(bam)

    output:
    tuple val(barcode), path(alignment), path("${alignment}.bai")

    script:
    """
    samtools index ${bam}
    """

}

process CALL_CONSENSUS {

    tag "${barcode}"
    publishDir params.consensus, mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    tuple val(barcode), path(bam), path(bai)

    output:
    tuple val(barcode), path("${barcode}.consensus.fasta")

    script:
    """
    samtools consensus \
    -m simple \
    -c ${params.min_variant_frequency} \
    -d ${params.min_depth_coverage} \
    ${bam} \
    > ${barcode}.consensus.fasta
    """

}

process GENERATE_MPILEUP {

    tag "${barcode}"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    tuple val(barcode), path(bam)

    output:
    tuple val(barcode), path("${barcode}.mpileup")

    script:
    """
    samtools mpileup -aa -A -Q 0 -d 0 ${bam} > ${barcode}.mpileup
    """

}

process FASTQ_CONVERSION {

    tag "${barcode}"
    publishDir params.basecall_fastqs

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    tuple val(barcode), path(bam)

    output:
    tuple val(barcode), path("${barcode}.fastq.gz")

    script:
    """
    samtools fastq ${bam} | gzip -c > ${barcode}.fastq.gz
    """

}

