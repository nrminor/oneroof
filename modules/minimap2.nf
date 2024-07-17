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
    secondary = params.secondary ? "no" : "ues"
    """
    minimap2 \
    --secondary=${secondary} \
    -ax ${preset} \
    ${refseq} ${reads} > ${barcode}.sam
    """

}
