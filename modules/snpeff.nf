process BUILD_DB {

    storeDir params.snpeff_cache

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    path refseq
    path genbank
    path snpeff_config

    output:
    path "genome"

    script:
    config_dir = "genome/ref_genome"
    """
    # 
    mkdir -p ${config_dir}/

    # 
    cp ${genbank} ${config_dir}/genes.gbk

    #
    cp ${snpeff_config} local.config

    #
    snpEff build -c local.config -dataDir genome/ -genbank -v ref_genome
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
    cp ${snpeff_config} local.config && \
    snpEff \
    -c local.config \
    -dataDir genome/ \
    -v ref_genome \
    ${vcf} \
    > ${barcode}.annotated.vcf
    """

}

process EXTRACT_FIELDS {

    tag "${barcode}"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    tuple val(sample_id), path(vcf)

    output:
    tuple val(sample_id), path(vcf), path("${sample_id}_variant_effects.tsv")

    script:
    """
    snpsift extractFields \
    ${vcf} \
    CHROM REF POS ALT AF AC DP MQ ANN[*].GENE ANN[*].GENEID ANN[*].EFFECT \
    ANN[*].HGVS_P ANN[*].CDNA_POS ANN[*].CDNA_LEN ANN[*].CDS_POS ANN[*].AA_POS \
    > ${sample_id}_variant_effects.tsv
    """
}
