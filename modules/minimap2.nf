process ALIGN_WITH_PRESET {

    tag "${barcode}, ${params.platform}"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    tuple val(barcode), path(reads)
    each path(refseq)

    output:
    tuple val(barcode), path("${barcode}.sam")

    script:
    preset = params.platform == "ont" ? "map-ont" : "sr"
    secondary = params.secondary ? "no" : "yes"
    """
    minimap2 \
    --secondary=${secondary} \
    -ax ${preset} \
    ${refseq} ${reads} \
    | samtools view -h -e '(rlen)>=${params.min_len}' \
    | samtools view -h -e '(rlen)<=${params.max_len}' > ${barcode}.sam
    """

}
