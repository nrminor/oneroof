
process DEDUP {

    input:
    tuple val(sample_id), path(reads1), path(reads2), path("${sample_id}.report.txt")

    output:
    tuple val(sample_id), path("${sample_id}.dedup.R1.fastq.gz"), path("${sample_id}.dedup.R2.fastq.gz")

    script:
    """
    fastp -D \
      -i ${reads1} -I ${reads2} \
      -o ${sample_id}.dedup.R1.fastq.gz \
      -O ${sample_id}.dedup.R2.fastq.gz 
    """
}
