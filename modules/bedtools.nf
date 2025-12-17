process GET_PRIMER_SEQS {

    /* */

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	input:
	path bed
	each path(refseq)

	output:
	path "${primer_combo}.fasta"

	script:
	primer_combo = file(bed.toString()).getSimpleName()
	"""
	bedtools getfasta -fi ${refseq} -bed ${bed} > ${primer_combo}.fasta
	"""

}

process BEDTOOLS_GENOMECOV {

    /* */

    tag "${sample_id}"
    publishDir params.genomecov, mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    tuple val(sample_id), path(bam), path(index)

    output:
    tuple val(sample_id), path("${sample_id}.per-base.bed")

    script:
    """
    bedtools genomecov -bga -ibam ${bam} > ${sample_id}.per-base.bed
    """

}
