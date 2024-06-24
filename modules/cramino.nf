process CRAMINO {

    /* */

    tag "${label}"
    publishDir params.cramino, mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    path bam
    each path(refseq)

    output:
    path "${label}.cramino.arrow"

    script:
    label = file(bam.toString()).getSimpleName()
    """
    cramino \
    --threads ${task.cpus} \
    --reference ${refseq} \
    --hist \
    --checksum \
    --arrow ${label}.cramino.arrow \
    --phased \
    --spliced \
    --ubam \
    ${bam}
    """

}
