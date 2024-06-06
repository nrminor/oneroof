process CALL_VARIANTS {

    tag "${barcode}"
    publishDir params.ivar, mode: 'copy', overwrite: true

    input:
    tuple val(barcode), path("${barcode}.mpileup")
    each path(refseq)

    output:
    tuple val(barcode), path("${barcode}.tsv")

    script:
    """
    cat ${mpileup} \
    | ivar variants \
    -t ${params.min_variant_frequency} \
    -m ${params.min_depth_coverage} \
    -p ${barcode} \
    -r ${refseq}
    """

}

process CONVERT_TO_VCF {

    tag "${barcode}"

    input:
    tuple val(barcode), path(ivar_table)

    output:
    tuple val(barcode), path("${barcode}.vcf")

    script:
    """
    ivar_variants_to_vcf.py ${ivar_table} ${barcode}.vcf
    """

}