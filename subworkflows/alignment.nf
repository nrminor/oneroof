#!/usr/bin/env nextflow

include { REPORT_REFERENCES       } from "../modules/reporting"
include { ALIGN_WITH_PRESET       } from "../modules/minimap2"
include { CONVERT_AND_SORT ; SORT_BAM ; INDEX } from "../modules/samtools"
include { RASUSA_ALN_DOWNSAMPLING } from "../modules/rasusa"
include { MOSDEPTH                } from "../modules/mosdepth"
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

    RASUSA_ALN_DOWNSAMPLING(
        CONVERT_AND_SORT.out
    )

    SORT_BAM(
        RASUSA_ALN_DOWNSAMPLING.out
    )

    INDEX(
        SORT_BAM.out
    )

    MOSDEPTH(
        INDEX.out
    )

    PLOT_COVERAGE(
        MOSDEPTH.out
    )

    MULTI_SAMPLE_PLOT(
        MOSDEPTH.out.map { _sample_id, files -> files }.flatten().filter { mosdepth_file -> mosdepth_file.toString().endsWith(".per-base.bed") }.collect(),
    )
    
    COVERAGE_SUMMARY(
        PLOT_COVERAGE.out.passing_cov.collect()
    )

    emit:
    index = INDEX.out
    coverage_summary = COVERAGE_SUMMARY.out
}
