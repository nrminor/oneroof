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
    params.meta_ref

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
    params.meta_ref

	script:
	"""
	sylph sketch -t ${task.cpus} -k ${params.k} -c 100 -r ${reads} -o ${sample_id}
	"""
}

// process CLASSIFY_SAMPLE {

// 	publishDir "${params.results}/sylph_results", mode: 'copy'

//     tag "${sample_id}"

//     errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
//     maxRetries 1

//     input:
//     tuple val(sample_id), path(sample_sketches), path(syldb)

//     output:
//     tuple val(sample_id), path("${sample_id}_sylph_results.tsv")

//    script:
// 	"""
// 	mkdir -p sylph_results
// 	sylph profile \
// 	-t ${task.cpus} --minimum-ani 90 --estimate-unknown -M 3 --read-seq-id 0.80 \
// 	${syldb} ${sample_sketches} -o ${sample_id}_sylph_results.tsv
// 	"""
// }
process CLASSIFY_SAMPLE {

    publishDir params.metagenomics, mode: 'copy'

    tag "${sample_id}"

    errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
    maxRetries 1

    input:
    tuple val(sample_id), path(sample_sketches), path(syldb)

    output:
    tuple val(sample_id), path("sylph_results/${sample_id}_sylph_results.tsv")

    script:
    """
	mkdir -p sylph_results
	sylph profile \
	-t ${task.cpus} --minimum-ani 90 --estimate-unknown -M 3 --read-seq-id 0.80 \
	${sample_sketches} ${syldb} > sylph_results/${sample_id}_sylph_results.tsv
    """
}



process SYLPH_TAX_DOWNLOAD {

	// errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
	// maxRetries 1

	input: 
	tuple val(sample_id), path("sylph_results/${sample_id}_sylph_results.tsv")

	output:
	val "ready"

	script:
	results = params.results
	"""
	cd ${results}
	mkdir -p sylph_tax_databases
	sylph-tax download --download-to sylph_tax_databases
	"""
}

process OVERLAY_TAXONOMY {
    publishDir params.metagenomics, mode: 'copy'
    tag "${sample_id}"
    
    input:
    tuple val(sample_id), path(tsv_path)
    path(input_db)
    val "ready"
    
    output:
    tuple val(sample_id), path("sylph_tax_results/${sample_id}*.sylphmpa")
    
    script:
    """
    mkdir -p sylph_tax_results
    sylph-tax taxprof ${tsv_path} -t ${input_db} -o "sylph_tax_results/${sample_id}"
    """
}



process MERGE_TAXONOMY {

    publishDir params.metagenomics, mode: 'copy'

    tag "merge"

    errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
    maxRetries 1

    input:
    tuple val(sample_id), path(input_dir)

    script:
    """
	sylph-tax merge ${input_dir} --column relative_abundance -o merged_taxonomy.tsv

    """
}


