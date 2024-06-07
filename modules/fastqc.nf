process FASTQC {

    /* */

	tag "${sample_id}"

	errorStrategy { task.attempt < 3 ? 'retry' : params.errorMode }
	maxRetries 2

	cpus 1

	input:
	path reads

	output:
    path "*.html", emit: html
    path "*.zip", emit: zip

	script:
	"""
	fastqc \
	--threads ${task.cpus} \
	--memory ${fastqc_memory} \
	${reads}
	"""

}

process FASTQC_RS {

    /* */

	tag "${sample_id}"

	errorStrategy { task.attempt < 3 ? 'retry' : params.errorMode }
	maxRetries 2

	cpus 1

	input:
	path reads

	output:
	path "${sample_id}_qc.html", emit: html
	path "${sample_id}/", emit: multiqc_data

	script:
    label = file(reads.toString()).getSimpleName()
	"""
	fqc -q ${reads} -s . > ${label}_qc.html
	mkdir ${sample_id}
	mv fastqc_data.txt ${sample_id}/fastqc_data.txt
	"""

}