include { DECONTAMINATE } from "../subworkflows/decontaminate"
include { FASTQC  } from "../modules/fastqc"

workflow QUALITY_CONTROL {
    /*
     * Run quality control on reads: optional decontamination and FastQC.
     * 
     * Note: MULTIQC has been moved to the main workflow level to allow
     * combining FastQC outputs with OneRoof metrics from later pipeline stages.
     */

    take:
    ch_reads
    ch_contam_fasta

    main:
    if ( params.contam_fasta && file(params.contam_fasta).isFile() ) {
        DECONTAMINATE(
            ch_reads,
            ch_contam_fasta
        )

        FASTQC(
            DECONTAMINATE.out
        )
    } else {
        FASTQC(
            ch_reads
        )

    }

    emit:
    fastq = FASTQC.out.fastq
    fastqc_zip = FASTQC.out.zip
}
