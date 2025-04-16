process PHASE_READS_WITH_DEVIDER {

    /* */

    tag "${sample_id}"
    // publishDir

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    cpus 4

    input:
    tuple val(sample_id), path(reads), path(vcf), path(fasta_ref)

    output:
    path "${sample_id}_devider"

    script:
    """
    devider
    -b ${reads} \
    -v ${vcf} \
    -r ${fasta_ref} \
    -o ${sample_id}_devider \
    -t ${task.cpus} \
    --preset ${params.devider_preset}
    """
}
