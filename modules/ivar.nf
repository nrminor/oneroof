process CALL_VARIANTS {

    script:
    """
    cat ${mpileup} \
    | ivar variants -t {params.min_variant_frequency} -m {params.min_depth_coverage} -p temp/{wildcards.sample} -r {input.reference}
    cp temp/{wildcards.sample}.tsv {output.tsv}
    """

}

process CONVERT_TO_VCF {

    script:
    """
    ivar_variants_to_vcf.py temp/{wildcards.sample}.tsv {output.vcf}
    """

}