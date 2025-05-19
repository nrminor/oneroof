process MERGE_READ_PAIRS {

    /* */

	tag "${sample_id}"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	cpus 4

	input:
	tuple val(sample_id), path(reads1), path(reads2)

	output:
	tuple val(sample_id), path("${sample_id}.merged.fastq.gz")

	script:
	"""
	vsearch \
	--fastq_mergepairs ${reads1} \
	--reverse ${reads2} \
	--fastqout_notmerged_fwd ${sample_id}.unmerged_fwd.fastq.gz \
	--fastqout_notmerged_rev ${sample_id}.unmerged_rev.fastq.gz \
	--eetabbedout ${sample_id}.merge_stats.tsv \
	--log ${sample_id}.merging.log \
	--threads ${task.cpus} \
	--eeout \
	--fastqout - \
	| gzip -c - > ${sample_id}.merged.fastq.gz
	"""

}

process ORIENT_READS {

    tag "${barcode}"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    tuple val(barcode), path(reads)
    each path(refseq)

    output:
    tuple val(barcode), path("${barcode}.oriented.fasta.gz")

    script:
    """
    vsearch \
	--orient ${reads} \
    --db ${refseq} \
    --fastaout - \
	| gzip -c > ${barcode}.oriented.fasta.gz
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
	tuple val(sample_id), path("${sample_id}_haplotypes.fasta"), emit: deduped_fasta
	path "${sample_id}_haplotype_metadata.tsv", emit: metadata

	script:
	"""
	vsearch \
	--fastx_uniques ${reads} \
	--fastaout ${sample_id}_haplotypes.fasta \
	--sizeout \
	--minuniquesize ${params.min_haplo_reads} \
	--tabbedout tmp.tsv \
	--strand both

	csvtk add-header -t \
	--names orig_label,clust_label,clust_index,seq_index_in_clust,clust_abundance,first_seq_label \
	tmp.tsv \
	| csvtk filter -t --filter "clust_abundance>=${params.min_haplo_reads}" \
	> ${sample_id}_haplotype_metadata.tsv

	rm tmp.tsv
	"""

}
