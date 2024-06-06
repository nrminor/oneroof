process BUILD_DB {

    script:
    """
    mkdir -p config/genome/ref_genome
    cp {input.gbk} {output.snpeff_gbk}
    snpEff build -c {input.config} -dataDir genome/ -genbank -v ref_genome
    """

}

process ANNOTATE_VCF {

    script:
    """
    snpEff -c {input.config} -dataDir genome/ -v ref_genome {input.vcf} > {output.vcf}
    """

}