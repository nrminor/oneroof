include { FASTQC  } from "../modules/fastqc"
include { INDEX_CONTAMINANTS; DECONTAMINATE } from "../modules/deacon"
include { MULTIQC } from "../modules/multiqc"

workflow QUALITY_CONTROL {
    take:
    ch_reads
    ch_contam_fasta

    main:
    if ( params.contam_fasta && file(params.contam_fasta).isFile() ) {
        INDEX_CONTAMINANTS(ch_contam_fasta)

        DECONTAMINATE(
            ch_reads.combine(INDEX_CONTAMINANTS.out)
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
