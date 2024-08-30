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
    REF_LIST=$(find . -type f -name '*.fasta' -o -name '*.bed' -o -name '*.gbk' -o -name '*.gb' -o -name '*.fa')
    if [ ! -d !{params.results}/reference_assets ]; then
        mkdir -p !{params.results}/reference_assets
    fi

    for file in "${REF_LIST}"; do
        cp ${file} !{params.results}/reference_assets/ &
    done

    wait
    '''
}