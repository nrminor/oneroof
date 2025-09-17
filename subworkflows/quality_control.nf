include { DECONTAMINATE } from "../subworkflows/decontaminate"
include { FASTQC  } from "../modules/fastqc"
include { MULTIQC } from "../modules/multiqc"

workflow QUALITY_CONTROL {
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

    MULTIQC(
        FASTQC.out.zip.collect()
    )

    emit:
    FASTQC.out.fastq
}
