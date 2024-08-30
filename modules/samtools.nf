process MERGE_BARCODES {

    /* */

    tag "${barcode}"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    tuple val(barcode), path("to_merge/???.bam")

    output:
    tuple val(barcode), path("${barcode}.bam")

    script:
    """
    samtools merge -o ${barcode}.bam to_merge/*.bam
    """
}

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
    | samtools sort -M -o ${barcode}.bam
    """

}

process SORT_BAM {

    tag "${barcode}"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    tuple val(barcode), path(bam)

    output:
    tuple val(barcode), path("${barcode}.sorted.bam")

    script:
    """
    cat ${bam} \
    | samtools sort -M -o ${barcode}.sorted.bam
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
    tuple val(barcode), path(bam), path("${barcode}*.bam.bai")

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
    -A \
    -c ${params.min_consensus_freq} \
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
    tuple val(barcode), path(bam), path(bai)

    output:
    tuple val(barcode), path("${barcode}.mpileup")

    script:
    """
    samtools mpileup -aa -A -Q 0 -d 0 ${bam} > ${barcode}.mpileup
    """

}

process FASTQ_CONVERSION {

    tag "${barcode}"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    tuple val(barcode), path(bam)

    output:
    tuple val(barcode), path("${barcode}.fastq.gz")

    script:
    """
    samtools fastq ${bam} | bgzip -o ${barcode}.fastq.gz
    """

}

process FAIDX {

    tag "${barcode}"
    publishDir params.basecall_fastqs, mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    tuple val(barcode), path(fastq)

    output:
    tuple val(barcode), path(fastq), path("*.fai")

    script:
    index_name = file(fastq).getName()
    """
    samtools faidx --fastq ${fastq}
    """

}

process SPLIT_SEGMENTS {

    tag "${sample_id}"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    cpus 4

    input:
    tuple val(sample_id), path(bam), path(index)

    output:
    path "${sample_id}_*.bam"

    shell:
    '''
    # Get list of contigs
    contigs=$(samtools idxstats "!{bam}" | cut -f 1 | grep -v '*')

    # Iterate through contigs and create separate BAM files
    for contig in $contigs
    do
        output_bam="!{sample_id}_${contig}.bam"
        (
            samtools view -b "!{bam}" "$contig" > "$output_bam"
            echo "Created $output_bam"
        ) &
    done

    # Wait for all remaining background jobs to finish
    wait

    echo "Splitting the !{sample_id} BAM is complete."
    '''

}
