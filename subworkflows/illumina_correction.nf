include { CORRECT_WITH_FASTP } from "../modules/fastp"
include { QUALITY_CONTROL } from "../subworkflows/quality_control"
include { COMPRESS_TO_SORTED_FASTA } from "../modules/seqkit"
include { FAIDX             } from "../modules/samtools"
include { EARLY_READ_DOWNSAMPLING } from "../modules/rasusa"

workflow ILLUMINA_CORRECTION {
    take:
    ch_uncorrected_reads
    ch_contam_fasta

    main:
    CORRECT_WITH_FASTP(
        ch_uncorrected_reads
    )

    QUALITY_CONTROL(
        CORRECT_WITH_FASTP.out,
        ch_contam_fasta
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

    EARLY_READ_DOWNSAMPLING(
        FAIDX.out
    )

    emit:
    EARLY_READ_DOWNSAMPLING.out
}
