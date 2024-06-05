process FIND_COMPLETE_AMPLICONS {

    /*
    */

    tag "${sample_id}"
    // label "general"
	// publishDir params.amplicon_reads, mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : params.errorMode }
	maxRetries 2

	cpus params.max_cpus

    input:
	tuple val(sample_id), path(reads)
    each path(search_patterns)

    output:
    tuple val(sample_id), env(primer_combo), path("${sample_id}_*_amplicons.fastq.gz")

    script:
    """
	primer_combo=`ls patterns/ | head -n 1 | sed 's/_comp_patterns.txt//'`

	mkdir primer_orients

	# find primers in the "template" orientation
	cat ${reads} | \
    seqkit grep \
	--threads ${task.cpus} \
	--max-mismatch 3 \
	--by-seq \
	--delete-matched \
	--pattern-file ${search_patterns}/\${primer_combo}_patterns.txt \
	-o primer_orients/${sample_id}_\${primer_combo}_amplicons_temp.fastq.gz

	# find primers in the complement orientation
	cat ${reads} | \
    seqkit grep \
	--threads ${task.cpus} \
	--max-mismatch 3 \
	--by-seq \
	--delete-matched \
	--pattern-file ${search_patterns}/\${primer_combo}_comp_patterns.txt \
	-o primer_orients/${sample_id}_\${primer_combo}_amplicons_comp.fastq.gz

	# combine them into one fastq
	seqkit scat -j ${task.cpus} --find-only primer_orients/ -o ${sample_id}_\${primer_combo}_amplicons.fastq.gz
    """

}