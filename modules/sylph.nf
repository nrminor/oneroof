process SKETCH_DATABASE_KMERS {

	/* */

	// storeDir params.db_cache

	errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
	maxRetries 1

    cpus 6

	input:
	path input_db

	output:
	path "*.syldb"

	when:
    params.sylph_db

	script:
	"""
	sylph sketch -t ${task.cpus} -k 31 -i -c 200 -g ${input_db} -o ${input_db}.syldb
	"""
}

process SKETCH_SAMPLE_KMERS {

	/*
	-c is the subsampling rate, also referred to as the compression parameter.
	A higher -c is faster but less sensitive at low coverage 
	*/

	tag "${sample_id}"

	errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
	maxRetries 1

	cpus 4

	input:
	tuple val(sample_id), path(reads)

	output:
	tuple val(sample_id), path("${sample_id}*.sylsp")

	when:
    params.sylph_db

	script:
	"""
	sylph sketch -t ${task.cpus} -k ${params.k} -c 100 -r ${reads} -o ${sample_id}
	"""
}

process CLASSIFY_SAMPLE {

	/* */

	tag "${sample_id}"

	errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
	maxRetries 1

	input:
	tuple val(sample_id), path(sample_sketches), path(syldb)

	output:
	tuple val(sample_id), path("${sample_id}*.tsv")

	script:
	"""
	sylph profile \
	-t ${task.cpus} --minimum-ani 90 --estimate-unknown -M 3 --read-seq-id 0.80 \
	${sample_sketches} ${syldb} > ${sample_id}_sylph_results.tsv
	"""
}

process SYLPH_TAX_DOWNLOAD {

	errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
	maxRetries 1

	output:
	val "ready"

	script:
	"""
	if [ ! -d ${params.sylph_tax_dir} ]; then
		mkdir -p ${params.sylph_tax_dir}
	fi
	sylph-tax download --download-to ${params.sylph_tax_dir}
	"""
}

process OVERLAY_TAXONOMY {

	/* */

	tag "${sample_id}"

	errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
	maxRetries 1

	input:
	tuple val(sample_id), path(tsv_dir), val(ready)

	script:
	"""
	sylph-tax taxprof sylph_results/*.tsv -t ${params.sylph_taxonomy} -o ${sample_id}
	"""
}
