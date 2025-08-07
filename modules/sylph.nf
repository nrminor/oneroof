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
	tuple val(sample_id), path("${reads}.sylsp")

	when:
    params.meta_ref

	script:
	"""
	sylph sketch -t ${task.cpus} -k ${params.k} -c 100 -r ${reads} -o ${sample_id}
	"""
}

process CLASSIFY_SAMPLE {

    publishDir params.metagenomics, mode: 'copy'

    tag "${sample_id}"

    //errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
    maxRetries 1

    input:
    tuple val(sample_id), path(sample_sketches), path(syldb)

    output:
    // tuple val(sample_id), path("sylph_results/${sample_id}_sylph_results.tsv")
	tuple val(sample_id), path("${sample_id}_sylph_results.tsv")

    script:
    """

    sylph profile \\
        -t ${task.cpus} --minimum-ani 90 --estimate-unknown -M 3 -c 10 --read-seq-id 0.80 \\
        ${sample_sketches} ${syldb} > ${sample_id}_sylph_results.tsv
    """
}



process SYLPH_TAX_DOWNLOAD {

	// errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
	// maxRetries 1
    input: 
	val(reads)

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
    tuple val(sample_id), path(tsv_path), val ("ready"), path(input_db)
    
    output:
    tuple val(sample_id), path("${sample_id}*.sylphmpa")
    
    script:
    """
	echo ${tsv_path}
    mkdir -p sylph_tax_results

	conflict_file=\$(find sylph_tax_results -type f -name '${sample_id}*.sylphmpa')
	if [[ -n "\$conflict_file" ]]; then
	echo "Removing conflicting file: \$conflict_file"
	rm -f "\$conflict_file"
	fi

	sylph-tax taxprof ${tsv_path} -t ${input_db} --add-folder-information -o ${sample_id}
    """
}



process MERGE_TAXONOMY {

    publishDir params.metagenomics, mode: 'copy'

    tag "merge"

    errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
    maxRetries 1

    input:
    tuple val(sample_id), path(input_dir)

	output:
	path "merged_taxonomy.tsv"

    script:
    """
	sylph-tax merge ${input_dir} --column relative_abundance -o merged_taxonomy.tsv

    """
}

process DOWNLOAD_DB_LINK {

	tag "download"

	input: 
	val db_link

	output:

	path "*.syldb"

	script: 
	// this should just download into cwd 
	"""
	wget ${db_link} -O database.syldb
	"""

}


