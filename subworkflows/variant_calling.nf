#!/usr/bin/env nextflow

include { REPORT_REFERENCES } from "../modules/reporting"
include { GENERATE_MPILEUP  } from "../modules/samtools"
include { CALL_VARIANTS ; CONVERT_TO_VCF } from "../modules/ivar"
include {
    BUILD_DB ;
    ANNOTATE_VCF ;
    EXTRACT_FIELDS
} from "../modules/snpeff"
include { MERGE_VCF_FILES   } from "../modules/bcftools"

workflow VARIANTS {
    take:
    ch_amplicons
    ch_refseq
    ch_genbank
    ch_snpeff_config

    main:
    REPORT_REFERENCES(
        ch_genbank
    )

    CALL_VARIANTS(
        ch_amplicons,
        ch_refseq,
    )

    CONVERT_TO_VCF(
        CALL_VARIANTS.out
    )

    BUILD_DB(
        ch_refseq,
        ch_genbank,
        ch_snpeff_config,
    )

    ANNOTATE_VCF(
        BUILD_DB.out,
        ch_snpeff_config,
        CONVERT_TO_VCF.out,
    )

    EXTRACT_FIELDS(
        ANNOTATE_VCF.out
    )

    MERGE_VCF_FILES(
        ANNOTATE_VCF.out.map { _label, vcf -> vcf }.collect()
    )

    emit:
    ANNOTATE_VCF.out
}
