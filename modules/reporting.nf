process REPORT_REFERENCES {

    /* */

    publishDir "${params.results}/reference_assets", mode: "copy", overwrite: true

	errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
	maxRetries 1

    cpus 3

    input:
    path files

    output:
    path "*"

    shell:
    '''
    REF_LIST=$(find . -type f -name '*.fasta' -o -name '*.bed' -o -name '*.gbk' -o -name '*.gb' -o -name '*.fa')

    for file in "${{REF_LIST}}"; do
        cp ${file} ./ &
    done
    
    wait
    '''
}