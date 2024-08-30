process REPORT_REFERENCES {

    /* */

    publishDir "${params.results}/reference_assets", mode: "copy", overwrite: true

	errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
	maxRetries 1

    cpus 3

    input:
    path files

    shell:
    '''
    if [ ! -d !{params.results}/reference_assets ]; then
        mkdir -p !{params.results}/reference_assets
    fi

    REF_LIST=$(find . -type f -name '*.fasta' -o -name '*.bed' -o -name '*.gbk' -o -name '*.gb' -o -name '*.fa')

    for file in "${REF_LIST}"; do
        cp `realpath ${file}` !{params.results}/reference_assets/
    done
    '''
}