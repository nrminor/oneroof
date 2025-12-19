/*
 * Sylph metagenomic profiling processes
 *
 * Sylph is an ultrafast tool for metagenomic profiling that uses k-mer sketching.
 * Sylph-tax adds taxonomic annotations to sylph output.
 *
 * Workflow:
 *   1. Sketch database (FASTA → .syldb) OR use pre-built .syldb
 *   2. Sketch sample reads (FASTQ → .sylsp)
 *   3. Profile samples against database (→ TSV)
 *   4. Optionally add taxonomy with sylph-tax (→ .sylphmpa)
 *   5. Optionally merge taxonomy profiles (→ merged TSV)
 */


// Download a pre-built sylph database from a URL
process DOWNLOAD_SYLPH_DB {

    tag "${db_url}"

    // Use publishDir instead of storeDir to avoid cross-run cache issues
    publishDir "${params.results}/sylph_cache", mode: 'copy', overwrite: true

    errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
    maxRetries 1

    input:
    val db_url

    output:
    path "*.syldb", emit: syldb

    script:
    // Extract filename from URL for meaningful output name
    def db_name = db_url.tokenize('/')[-1]
    """
    wget -q "${db_url}" -O ${db_name}
    """
}


// Sketch a custom FASTA into a sylph database
process SKETCH_DATABASE {

    tag "${fasta.baseName}"

    storeDir "${params.results}/sylph_cache"

    errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
    maxRetries 1

    cpus 6

    input:
    path fasta

    output:
    path "*.syldb", emit: syldb

    script:
    """
    sylph sketch \\
        -t ${task.cpus} \\
        -k ${params.k} \\
        -c 200 \\
        ${fasta} \\
        -o ${fasta.baseName}
    """
}


// Sketch sample reads into a sylph sample sketch
process SKETCH_SAMPLES {

    tag "${sample_id}"

    errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
    maxRetries 1

    cpus 4

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}.sylsp"), emit: sketch, optional: true

    script:
    """
    sylph sketch \\
        -t ${task.cpus} \\
        -k ${params.k} \\
        -c 100 \\
        -r ${reads} \\
        -o ${sample_id}
    """
}


// Profile samples against a sylph database
process PROFILE_SAMPLES {

    tag "${sample_id}"

    publishDir params.metagenomics, mode: 'copy'

    errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
    maxRetries 1

    cpus 4

    input:
    tuple val(sample_id), path(sample_sketch), path(syldb)

    output:
    tuple val(sample_id), path("${sample_id}_sylph.tsv"), emit: profile, optional: true

    script:
    """
    sylph profile \\
        -t ${task.cpus} \\
        --minimum-ani 90 \\
        --estimate-unknown \\
        ${sample_sketch} ${syldb} \\
        > ${sample_id}_sylph.tsv
    """
}


// Download sylph-tax taxonomy metadata (one-time, cached)
process DOWNLOAD_TAXONOMY {

    storeDir "${params.results}/sylph_cache"

    errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
    maxRetries 1

    output:
    path "sylph_taxonomy", emit: taxonomy_dir

    script:
    """
    mkdir -p sylph_taxonomy
    sylph-tax --no-config download --download-to sylph_taxonomy
    """
}


// Add taxonomy annotations to sylph profile results
process ADD_TAXONOMY {

    tag "${sample_id}"

    publishDir params.metagenomics, mode: 'copy'

    errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
    maxRetries 1

    input:
    tuple val(sample_id), path(profile_tsv)
    path taxonomy_dir
    val db_identifier  // e.g., "GTDB_r220", "IMGVR_4.1"

    output:
    tuple val(sample_id), path("${sample_id}*.sylphmpa"), emit: taxprofile

    script:
    """
    sylph-tax --no-config --taxonomy-dir ${taxonomy_dir} \\
        taxprof ${profile_tsv} \\
        -t ${db_identifier} \\
        -o ${sample_id}
    """
}


// Merge all taxonomy profiles into a single table
process MERGE_TAXONOMY {

    publishDir params.metagenomics, mode: 'copy'

    errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
    maxRetries 1

    input:
    path taxprofiles  // collected .sylphmpa files

    output:
    path "merged_taxonomy.tsv", emit: merged

    script:
    """
    sylph-tax merge ${taxprofiles} --column relative_abundance -o merged_taxonomy.tsv
    """
}
