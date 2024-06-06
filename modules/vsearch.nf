process ORIENT_READS {

    tag "${barcode}"

    input:
    tuple val(barcode), path(reads)
    each path(refseq)

    output:
    tuple val(barcode), path("${barcode}.oriented.fastq")

    script:
    """
    vsearch --orient ${reads} \
    --db ${refseq} \
    --fastqout ${barcode}.oriented.fastq
    """

}