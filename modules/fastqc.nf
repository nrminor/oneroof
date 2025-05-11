process FASTQC {

    /* */

	tag "${barcode}"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	cpus 1

	input:
	tuple val(barcode), path(reads)

	output:
    path "*.html", emit: html
    path "*.zip", emit: zip
    tuple val(barcode), path(reads), emit: fastq

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
