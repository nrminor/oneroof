process FASTQC {

    /* */

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
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
	${reads}
	"""

}

process FASTQC_RS {

    /* */

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	cpus 1

	input:
	path reads

	output:
	path "${label}_qc.html", emit: html
	path "${label}/", emit: multiqc_data

	script:
    label = file(reads.toString()).getSimpleName()
	"""
	fqc -q ${reads} -s . > ${label}_qc.html
	mkdir ${label}
	mv fastqc_data.txt ${label}/fastqc_data.txt
	"""

}