process DOWNLOAD_MODELS {

    storeDir params.model_cache

    output:
    path "*"

    script:
    println "Downloading basecaller models."
    """
    dorado download --verbose
    """

}

process BASECALL {

    tag "${sample_id}"
    maxForks params.basecall_max

    input:
    each path(models)
    tuple val(sample_id), path("pod5s/???.pod5")

    output:
    tuple val(sample_id), path("")

    script:
    """
    dorado basecaller \
    ${params.model} pod5s/ \
    --kit-name ${params.kit} \
    > "${sample_id}.bam"
    """

}

process DEMULTIPLEX {

    publishDir params.basecall_bams, mode: 'copy', overwrite: true
    maxForks params.basecall_max

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