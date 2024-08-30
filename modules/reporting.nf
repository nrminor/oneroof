process REPORT_REFERENCES {

    /* */

    publishDir "${params.results}/reference_assets", mode: "copy", overwrite: true

	errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
	maxRetries 1

    cpus 3

    output:
    path "*"

    shell:
    '''
    if [ -n !{params.ref_fasta} ] && [ -e `realpath !{params.ref_fasta}` ]; then
        cp `realpath !{ref_fasta}` . &
    fi
    if [ -n !{params.ref_gbk} ] && [ -e `realpath !{params.ref_gbk}` ]; then
        cp `realpath !{params.ref_gbk}` . &
    fi
    if [ -n !{params.primer_bed} ] && [ -e `realpath !{params.primer_bed}` ]; then
        cp `realpath !{params.primer_bed}` . &
    fi
    wait
    '''
}
