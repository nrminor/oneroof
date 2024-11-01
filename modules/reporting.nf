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

    find . -name '*.fa*' -o -name '*.bed' -o -name '*.gb*' > ref_files.txt

    for file in "$(cat ref_files.txt)"; do
        cp `realpath ${file}` !{params.results}/reference_assets/
    done
    '''
}

process PUBLISH_COMMAND {

    /* */

    publishDir params.results, mode: "copy", overwrite: true

	errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
	maxRetries 1

    cpus 1

    output:
    path "run_command.sh"

    script:
    command = workflow.commandLine
    """
    echo "${command}" > run_command.sh
    """

}
