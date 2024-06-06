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
    each path(models)
    path pod5_dir

    output:
    val "basecalled.bam"

    script:
    """
    dorado basecaller \
    ${params.model} ${pod5_dir} \
    --kit-name ${params.kit} \
    > "basecalled.bam"
    """

}

process DEMULTIPLEX {

    publishDir params.basecall_bams, mode: 'copy', overwrite: true
    maxForks params.basecall_max

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    path "bams/???.bam"

    output:
    path "demux/*barcode*"

    script:
    """
    dorado demux \
    bams/ \
    --kit-name ${params.kit} \
    --output-dir demux/
    """

}