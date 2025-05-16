process CORRECT_WITH_FASTP {

    tag "${sample_id}"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}.corrected.fastq.gz")

    script:
    """
    fastp -i `realpath ${reads}` -o "${sample_id}.corrected.fastq.gz" \
    --qualified_quality_phred 20 \
    --unqualified_percent_limit 30 \
    --length_required 50 \
    --trim_poly_g \
    --trim_poly_x
    """

}
