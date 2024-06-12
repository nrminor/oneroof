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
	tuple val(barcode), path("${barcode}.${params.downsample_to}x.fastq.gz")

    script:
    """
    rasusa reads \
    --coverage ${params.downsample_to} \
    --genome-size ${faidx} \
    --seed 14 \
    --output-type b \
    --output ${barcode}.${params.downsample_to}x.fastq.gz
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
	tuple val(barcode), path("${barcode}.${params.downsample_to}x.bam")

    script:
    """
    rasusa aln \
    --coverage ${params.downsample_to} \
    --seed 14 \
    --output ${barcode}.${params.downsample_to}x.bam
    """

}
