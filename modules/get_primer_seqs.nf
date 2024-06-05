process GET_PRIMER_SEQS {

    /*
    */

	tag "${params.desired_amplicon}"
    label "general"

	errorStrategy { task.attempt < 3 ? 'retry' : params.errorMode }
	maxRetries 2

	input:
	path bed
	each path(refseq)

	output:
	path "patterns/"

	script:
	primer_combo = file(bed.toString()).getSimpleName()
	"""
	bedtools getfasta -fi ${refseq} -bed ${bed} > primer_seqs.fasta
	seqkit seq --complement --validate-seq --seq-type DNA primer_seqs.fasta | \
	grep -v "^>" > ${primer_combo}_comp_patterns.txt
	cat primer_seqs.fasta | grep -v "^>" > ${primer_combo}_patterns.txt
	mkdir patterns
	mv ${primer_combo}_comp_patterns.txt ${primer_combo}_patterns.txt patterns/
	"""
}