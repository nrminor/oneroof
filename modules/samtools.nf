process CONVERT_AND_SORT {

    script:
    """
    samtools view -bS ${sam} \
    | samtools sort -o {output.bam}
    """

}

process INDEX {

    output:
    tuple path(alignment), path("${alignment}.bai")

    script:
    """
    samtools index ${alignment}
    """

}

process CALL_CONSENSUS {

    script:
    """
    samtools consensus -m simple -c {params.min_variant_frequency} -d {params.min_depth_coverage} {input.bam} > {output.consensus}
    """

}

process GENERATE_MPILEUP {

    script:
    """
    samtools mpileup -aa -A -Q 0 -d 0 {input.bam} > ???
    """

}

