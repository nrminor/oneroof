process RASUSA_READ_DOWNSAMPLING {

	/* */

	tag "${barcode}"
	publishDir params.complete_amplicons, mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	cpus 3

	input:
	tuple val(barcode), path(amplicons), path(faidx)

	output:
	tuple val(barcode), path("${barcode}*.fastq.gz")

    script:
    basename = file(amplicons).getName().replace(".fastq.gz", "")
    if ( params.downsample_to == 0 )
        """
        cp ${amplicons} ${basename}.no_downsampling.fastq.gz
        """
    else
        """
        rasusa reads \
        --coverage ${params.downsample_to} \
        --genome-size ${faidx} \
        --seed 14 \
        --output-type g \
        --output ${basename}.${params.downsample_to}x.fastq.gz \
        ${amplicons}
        """

}

process RASUSA_ALN_DOWNSAMPLING {

	/* */

	tag "${barcode}"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	cpus 3

	input:
	tuple val(barcode), path(bam)

	output:
	tuple val(barcode), path("${barcode}*.bam")

    script:
    basename = file(bam).getName().replace(".bam", "")
    if ( params.downsample_to == 0 )
        """
        cp ${bam} ${basename}.no_downsampling.bam
        """
    else
        """
        samtools index ${bam} && \
        rasusa aln \
        --coverage ${params.downsample_to} \
        --seed 14 \
        --output ${basename}.${params.downsample_to}x.bam \
        ${bam}
        """

}
