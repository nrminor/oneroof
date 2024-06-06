process ALIGN_WITH_PRESET {

    script:
    preset = platform == "ont" ? "map-ont" : "sr"
    """
    minimap2 -ax ${preset} ${refseq} <(zcat ${fastq}) > output.sam
    """

}

