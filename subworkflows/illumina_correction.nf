#!/usr/bin/env nextfllow

include {
    FIND_ADAPTER_SEQS;
    TRIM_ADAPTERS;
    REMOVE_OPTICAL_DUPLICATES;
    REMOVE_LOW_QUALITY_REGIONS;
    REMOVE_ARTIFACTS;
    ERROR_CORRECT_PHASE_ONE;
    ERROR_CORRECT_PHASE_TWO;
    ERROR_CORRECT_PHASE_THREE;
    QUALITY_TRIM;
    CLUMP_READS;
} from "../modules/bbmap"

workflow ILLUMINA_CORRECTION {

    take:
        ch_uncorrected_reads

    main:
        FIND_ADAPTER_SEQS (
            ch_uncorrected_reads
        )

        TRIM_ADAPTERS (
            FIND_ADAPTER_SEQS.out
        )

        REMOVE_OPTICAL_DUPLICATES (
            TRIM_ADAPTERS.out
        )

        REMOVE_LOW_QUALITY_REGIONS (
            REMOVE_OPTICAL_DUPLICATES.out
        )

        REMOVE_ARTIFACTS (
            REMOVE_LOW_QUALITY_REGIONS.out
        )

        ERROR_CORRECT_PHASE_ONE (
            REMOVE_ARTIFACTS.out
        )

        ERROR_CORRECT_PHASE_TWO (
            ERROR_CORRECT_PHASE_ONE.out
        )

        ERROR_CORRECT_PHASE_THREE (
            ERROR_CORRECT_PHASE_TWO.out
        )

        QUALITY_TRIM (
            ERROR_CORRECT_PHASE_THREE.out
        )

        CLUMP_READS (
            QUALITY_TRIM.out
        )

    emit:
        CLUMP_READS.out

}
