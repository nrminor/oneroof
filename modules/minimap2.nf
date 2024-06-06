process ALIGN_WITH_PRESET {

    tag "${barcode}, ${platform}"

    input:
    tuple val(barcode), path(fastq)
    each path(refseq)
    each path(platform)

    output:
    tuple val(barcode), path("${barcode}.sam")

    script:
    String preset = platform == "ont" ? "map-ont" : "sr"
    """
    minimap2 -ax ${preset} ${refseq} <(zcat ${fastq}) > ${barcode}.sam
    """

}

