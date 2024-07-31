process ORIENT_READS {

    tag "${barcode}"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    tuple val(barcode), path(reads)
    each path(refseq)

    output:
    tuple val(barcode), path("${barcode}.oriented.fastq")

    script:
    """
    vsearch --orient ${reads} \
    --db ${refseq} \
    --fastqout ${barcode}.oriented.fastq
    """

}

process IDENTIFY_HAPLOTYPES {

	/*
    */

	tag "${sample_id}"
	publishDir "${params.haplo}/${sample_id}", mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : params.errorMode }
	maxRetries 2

	cpus 1

	input:
	tuple val(sample_id), path(reads)

	output:
	tuple val(sample_id), path("${sample_id}_deduped.fasta"), emit: deduped_fasta
	path "${sample_id}_haplotype_metadata.tsv", emit: metadata

	script:
	"""
	vsearch \
	--fastx_uniques ${reads} \
	--fastaout ${sample_id}_deduped.fasta \
	--sizeout \
	--minuniquesize ${params.min_reads} \
	--tabbedout tmp.tsv \
	--strand both

	csvtk add-header -t \
	--names orig_label,clust_label,clust_index,seq_index_in_clust,clust_abundance,first_seq_label \
	tmp.tsv \
	| csvtk filter -t --filter "clust_abundance>=${params.min_reads}" \
	> ${sample_id}_haplotype_metadata.tsv

	rm tmp.tsv
	"""

}
