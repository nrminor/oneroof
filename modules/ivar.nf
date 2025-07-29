process CALL_VARIANTS {

    tag "${barcode}"
    label "big_mem"
    publishDir params.ivar, mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    tuple val(barcode), path(bam), path(bai)
    each path(refseq)

    output:
    tuple val(barcode), path("${barcode}.tsv")

    script:
    """
    samtools mpileup -aa -A -Q 0 -d 0 ${bam} \
    | ivar variants \
    -t ${params.min_variant_frequency} \
    -m ${params.min_depth_coverage} \
    -p ${barcode} \
    -r ${refseq}
    """

}

process CALL_CONSENSUS {

    tag "${barcode}"
    label "big_mem"
    publishDir params.consensus, mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    tuple val(barcode), path(bam), path(bai)

    output:
    tuple val(barcode), path("${barcode}.consensus.fa*")

    script:
    """
    samtools mpileup -aa -A -Q 0 -d 0 ${bam} \
    | ivar consensus \
    -p ${barcode}.consensus \
    -t ${params.min_consensus_freq} \
    -m ${params.min_depth_coverage} \
    -q 0 \
    -n N
    """
}

process CONVERT_TO_VCF {

    tag "${barcode}"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    tuple val(barcode), path(ivar_table)

    output:
    tuple val(barcode), path("${barcode}.vcf")

    script:
    """
    ivar_variants_to_vcf.py convert ${ivar_table} ${barcode}.vcf
    """

}
