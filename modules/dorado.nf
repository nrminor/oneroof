process BASECALL (

    tag 
    maxForks params.basecall_max

    input:
    tuple val(sample_id), path(pod5)

    output:
    tuple val(sample_id), path("")

    script:
    """
    dorado basecall
    """

)

process DEMULTIPLEX (

    tag 
    maxForks params.basecall_max

    input:
    tuple val(sample_id), path(pod5)

    output:
    tuple val(sample_id), path("")

    script:
    """
    dorado demux
    """

)