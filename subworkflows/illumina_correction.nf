include { CORRECT_WITH_FASTP } from "../modules/fastp"
include { FASTQC  } from "../modules/fastqc"
include { INDEX_CONTAMINANTS; DECONTAMINATE } from "../modules/deacon"
include { MULTIQC } from "../modules/multiqc"
include { COMPRESS_TO_SORTED_FASTA } from "../modules/seqkit"
include { FAIDX             } from "../modules/samtools"
include { EARLY_RASUSA_READ_DOWNSAMPLING } from "../modules/rasusa"

workflow ILLUMINA_CORRECTION {
    take:
    ch_uncorrected_reads
    ch_contam_fasta

    main:
    CORRECT_WITH_FASTP(
        ch_uncorrected_reads
    )

    if ( params.contam_fasta && file(params.contam_fasta).isFile() ) {
        INDEX_CONTAMINANTS(ch_contam_fasta)

        DECONTAMINATE(
            CORRECT_WITH_FASTP.out.combine(INDEX_CONTAMINANTS.out)
        )

        FASTQC(
            DECONTAMINATE.out
        )
    } else {
        FASTQC(
            CORRECT_WITH_FASTP.out
        )
       
    }

    MULTIQC(
        FASTQC.out.zip.collect()
    )

    COMPRESS_TO_SORTED_FASTA(
        CORRECT_WITH_FASTP.out
    )

    FAIDX(
        COMPRESS_TO_SORTED_FASTA.out
            .map { id, fastq -> tuple(id, fastq, file(fastq).countFasta()) }
            .filter { it[2] > 0 }
            .map { id, fastq, _read_count -> tuple(id, file(fastq)) }
    )

    EARLY_RASUSA_READ_DOWNSAMPLING(
        FAIDX.out
    )

    emit:
    EARLY_RASUSA_READ_DOWNSAMPLING.out
}
