process RASUSA_READS {

	/* */

	tag "${barcode}"
	publishDir params.complete_amplicons, mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	cpus 3

	input:
	tuple val(barcode), path(amplicons), path(faidx)

	output:
	tuple val(barcode), path("${barcode}.*.fastq.gz")

    script:
    if ( params.downsample_to == 0 )
        """
        cp ${amplicons} ${amplicons}.no_downsampling.fastq.gz
        """
    else
        """
        rasusa reads \
        --coverage ${params.downsample_to} \
        --genome-size ${faidx} \
        --seed 14 \
        --output-type g \
        --output ${barcode}.${params.downsample_to}x.fastq.gz \
        ${amplicons}
        """

}

process RASUSA_ALN {

	/* */

	tag "${barcode}"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	cpus 3

	input:
	tuple val(barcode), path(bam)

	output:
	tuple val(barcode), path("${barcode}.*.bam")
    script:
    if ( params.downsample_to == 0 )
        """
        cp ${bam} ${barcode}.no_downsampling.bam
        """
    else
        """
        samtools index ${bam} && \
        rasusa aln \
        --coverage ${params.downsample_to} \
        --seed 14 \
        --output ${barcode}.${params.downsample_to}x.bam \
        ${bam}
        """

}
