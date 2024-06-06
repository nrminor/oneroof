process BUILD_DB {

    cache params.snpeff_cache

    input:
    path refseq
    path genbank
    path snpeff_config

    output:
    path "config"

    script:
    String config_dir = "config/genome/ref_genome/"
    """
    # 
    mkdir -p ${config_dir}

    # 
    cp ${genbank} ${config_dir}/genes.gbk

    # 
    snpEff build \
    -c ${snpeff_config} \
    -dataDir genome/ \
    -genbank \
    -v ref_genome
    """

}

process ANNOTATE_VCF {

    tag "${barcode}"
    publishDir params.variants, mode: 'copy', overwrite: true

    input:
    each path(snpeff_db)
    each path(snpeff_config)
    tuple val(barcode), path(vcf)

    output:

    tuple val(barcode), path("${barcode}.annotated.vcf")

    script:
    """
    snpEff \
    -c ${snpeff_config} \
    -dataDir genome/ \
    -v ref_genome \
    ${vcf} \
    > ${barcode}.annotated.vcf
    """

}