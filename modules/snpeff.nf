process BUILD_DB {

    storeDir params.snpeff_cache

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    path refseq
    path genbank
    path snpeff_config

    output:
    path "config"

    script:
    config_dir = "config/genome/ref_genome/"
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
    publishDir params.vcf, mode: 'copy', overwrite: true

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