include { REPORT_REFERENCES ; EXTRACT_COVERAGE_METRICS ; EXTRACT_ALIGNMENT_METRICS } from "../modules/reporting"
include { ALIGN_WITH_PRESET      } from "../modules/minimap2"
include { CONVERT_AND_SORT ; SORT_BAM ; INDEX } from "../modules/samtools"
include { ALIGNMENT_DOWNSAMPLING } from "../modules/rasusa"
include { BEDTOOLS_GENOMECOV     } from "../modules/bedtools"
include { PLOT_COVERAGE ; MULTI_SAMPLE_PLOT ; COVERAGE_SUMMARY } from "../modules/plot_coverage"



workflow ALIGNMENT {
    take:
    ch_amplicons
    ch_refseq

    main:
    REPORT_REFERENCES(
        ch_refseq
    )

    ALIGN_WITH_PRESET(
        ch_amplicons,
        ch_refseq,
    )

    CONVERT_AND_SORT(
        ALIGN_WITH_PRESET.out
    )

    ALIGNMENT_DOWNSAMPLING(
        CONVERT_AND_SORT.out
    )

    SORT_BAM(
        ALIGNMENT_DOWNSAMPLING.out
    )

    INDEX(
        SORT_BAM.out
    )

    BEDTOOLS_GENOMECOV(
        INDEX.out
    )

    EXTRACT_COVERAGE_METRICS(
        BEDTOOLS_GENOMECOV.out
    )

    EXTRACT_ALIGNMENT_METRICS(
        INDEX.out
    )

    PLOT_COVERAGE(
        BEDTOOLS_GENOMECOV.out
    )

    MULTI_SAMPLE_PLOT(
        BEDTOOLS_GENOMECOV.out.map { _sample_id, bed_file -> bed_file }.collect()
    )

    COVERAGE_SUMMARY(
        PLOT_COVERAGE.out.passing_cov.collect()
    )

    emit:
    index              = INDEX.out
    coverage_summary   = COVERAGE_SUMMARY.out
    coverage_metrics   = EXTRACT_COVERAGE_METRICS.out
    alignment_metrics  = EXTRACT_ALIGNMENT_METRICS.out
}
