process DOWNLOAD_MODELS {

    storeDir params.model_cache

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    output:
    path "*"

    script:
    """
    dorado download --verbose
    """

}

process BASECALL {

    maxForks params.basecall_max

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    path models
    path "pod5s/???.pod5"

    output:
    path "basecalled.bam"

    script:
    """
    dorado basecaller \
    ${params.model} pod5s/ \
    --kit-name ${params.kit} \
    > basecalled.bam
    """

}

process DEMULTIPLEX {

    /* */

    publishDir params.basecall_bams, mode: 'copy', overwrite: true
    maxForks params.basecall_max

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    path bam

    output:
    path "demux/*barcode*"

    script:
    """
    dorado demux ${bam} --no-classify --output-dir demux
    """

}
