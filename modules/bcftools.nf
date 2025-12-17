process MERGE_VCF_FILES {

    /* */

    publishDir params.vcf, mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
	maxRetries 1

    input:
    path vcf_files

    output:
    tuple path("all_sample_vcf.gz"), path("all_sample_vcf.gz.tbi")

    script:
    """
    bcftools merge ${vcf_files} | bgzip > all_sample_vcf.gz
    tabix -p vcf all_sample_vcf.gz
    """

}

process COMPRESS_AND_INDEX_VCF {

    tag "${sample_id}"

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    input:
    tuple val(sample_id), path(vcf)

    output: 
    tuple val(sample_id), path("${vcf}.gz"), path("${vcf}.gz.tbi")

    script:
    """
    bgzip -c "${vcf}" > "${vcf}.gz"
    tabix -p vcf "${vcf}.gz"
    """
}