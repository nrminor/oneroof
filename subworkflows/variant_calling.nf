include { REPORT_REFERENCES ; EXTRACT_VARIANT_METRICS } from "../modules/reporting"
include { GENERATE_MPILEUP  } from "../modules/samtools"
include { CALL_VARIANTS ; CONVERT_TO_VCF } from "../modules/ivar"
include {
    BUILD_DB ;
    ANNOTATE_VCF ;
    EXTRACT_FIELDS
} from "../modules/snpeff"
include { MERGE_VCF_FILES   } from "../modules/bcftools"
include { COLLECT_FULL_VARIANT_TABLE } from "../modules/variant_table"

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
        CALL_VARIANTS.out.filter { _id, tsv -> file(tsv).countLines() > 1 },
        ch_refseq
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

    EXTRACT_VARIANT_METRICS(
        EXTRACT_FIELDS.out
    )

    MERGE_VCF_FILES(
        ANNOTATE_VCF.out.map { _label, vcf -> vcf }.collect()
    )

    COLLECT_FULL_VARIANT_TABLE(
        EXTRACT_FIELDS.out.map { _sample_id, _vcf, tsv -> tsv }.collect()
    )

    emit:
    vcf = CONVERT_TO_VCF.out  // unannotated VCFs for haplotyping
    annotate = ANNOTATE_VCF.out
    merge_vcf_files = MERGE_VCF_FILES.out.collect()
    full_variant_table = COLLECT_FULL_VARIANT_TABLE.out
    variant_metrics = EXTRACT_VARIANT_METRICS.out
}
