process PHASE_READS_WITH_DEVIDER {

    /* */

    tag "${sample_id}"
    publishDir params.haplotyping, mode: 'copy'

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    cpus 4

    input:
    tuple val(sample_id), path(bam), path(bai), path(vcf), path(tbi), path (fasta_ref) 

    output:
    path "${sample_id}_devider"

    script:
    """
    devider \
    --bam-file ${bam} \
    --vcf-file ${vcf} \
    --reference-fasta ${fasta_ref} \
    -o ${sample_id}_devider \
    -t ${task.cpus} \
    --preset ${params.devider_preset}
    """
}
