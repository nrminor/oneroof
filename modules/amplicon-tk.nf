process AMPLICON_TK_TRIM {

    /* */

    tag "${}"
    publishDir , mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    tuple val(sample_id), val(fastq)
    path primer_bed
    path ref

    output:
    path "${sample_id}.amplicons.fastq.gz"

    script:
    """
    amplicon-tk trim \
    --input-file ${fastq} \
    --bed-file ${primer_bed} \
    --fasta-ref ${ref} \
    --left-suffix "_LEFT" \
    --right-suffix "_RIGHT" \
    --min-freq 0.0 \
    --output ${sample_id}.amplicons.fastq.gz
    """

}

// process AMPLICON_TK_SORT {}

// process AMPLICON_TK_CONSENSUS {}
