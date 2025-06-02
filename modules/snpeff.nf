process BUILD_DB {

    // storeDir params.snpeff_cache

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

    if(genbank.name.endsWith(".gbk"))
        """
        #
        mkdir -p ${config_dir}

        #
        cp ${genbank} ${config_dir}/genes.gbk

        #
        cp ${snpeff_config} local.config

        #
        snpEff build -c local.config -dataDir genome/ -genbank -v ref_genome
        """

    else if (genbank.name.endsWith(".gff") | genbank.name.endsWith(".gff3") )
        """
        #
        mkdir -p ${config_dir}

        #
        cp ${genbank} ${config_dir}/genes.gff
        cp ${refseq} ${config_dir}/sequences.fa
   

        #
        cp ${snpeff_config} local.config

        #
        snpEff build -c local.config -dataDir genome/ -gff3 -v -noCheckCds -noCheckProtein ref_genome 

        """

}

process ANNOTATE_VCF {

    tag "${barcode}"
    publishDir params.vcf, mode: 'copy', overwrite: true

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    input:
    each path(snpeff_db)
    each path(snpeff_config)
    tuple val(barcode), path(vcf)

    output:
    tuple val(barcode), path("${barcode}.annotated.vcf")

    script:
    """
    cp ${snpeff_config} local.config && \
    snpEff -Xmx3g \
    -c local.config \
    -dataDir genome/ \
    -v ref_genome \
    ${vcf} \
    > ${barcode}.annotated.vcf
    """

}

process EXTRACT_FIELDS {

    tag "${sample_id}"

    publishDir params.variant_tsv, mode: 'copy', overwrite: true, pattern: "*.tsv"

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    input:
    tuple val(sample_id), path(vcf)

    output:
    tuple val(sample_id), path(vcf), path("${sample_id}_variant_effects.tsv")

    script:
    """
    SnpSift -Xmx3g extractFields \
    -s "," \
    ${vcf} \
    CHROM REF POS ALT AF AC DP GEN[0].REF_DP GEN[0].ALT_DP GEN[0].ALT_FREQ MQ ANN[0].GENE ANN[0].EFFECT ANN[0].HGVS_P ANN[0].CDS_POS ANN[0].AA_POS \
    > ${sample_id}_variant_effects.tsv
    """
}
