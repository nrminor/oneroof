#!/usr/bin/env nextflow

include { GATHER_NANOPORE } from "../subworkflows/gather_nanopore"
include { PRIMER_HANDLING } from "../subworkflows/primer_handling"
include { ALIGNMENT } from "../subworkflows/alignment"
include { HAPLOTYPING } from "../subworkflows/haplotying"
include { QUALITY_CONTROL } from "../subworkflows/quality_control"
include { CONSENSUS } from "../subworkflows/consensus_calling"
include { VARIANTS } from "../subworkflows/variant_calling"
include { PHYLO } from "../subworkflows/phylo"

workflow NANOPORE {

    /* */

    take:
        ch_primer_bed
        ch_refseq
        ch_refgbk
        ch_snpeff_config

    main:
        assert params.platform == "ont"

        GATHER_NANOPORE ( )

        if ( params.primer_bed && params.primer_bed != "" ) {

            PRIMER_HANDLING (
                GATHER_NANOPORE.out,
                ch_primer_bed,
                ch_refseq
            )

            ALIGNMENT (
                PRIMER_HANDLING.out,
                ch_refseq
            )

            if ( countFastaHeaders(params.refseq) == countAmplicons(params.primer_bed) ) {

                HAPLOTYPING (
                    ALIGNMENT.out,
                    ch_refseq
                )

            }

            QUALITY_CONTROL (
                PRIMER_HANDLING.out,
                ALIGNMENT.out
            )

        } else {

            ALIGNMENT (
                GATHER_NANOPORE.out,
                ch_refseq
            )

            QUALITY_CONTROL (
                GATHER_NANOPORE.out,
                ALIGNMENT.out
            )

        }

        CONSENSUS (
            ALIGNMENT.out
        )

        VARIANTS (
            ALIGNMENT.out,
            ch_refseq,
            ch_refgbk,
            ch_snpeff_config
        )

        PHYLO (
            CONSENSUS.out
        )

}


// Groovy code
import java.nio.file.Files
import java.nio.file.Paths
import java.util.stream.Collectors

def countFastaHeaders(String filePath) {
    def headerCount = 0
    def file = new File(filePath)

    if (!file.exists()) {
        throw new FileNotFoundException("The file $filePath does not exist.")
    }

    file.eachLine { line ->
        if (line.trim().startsWith('>')) {
            headerCount++
        }
    }

    return headerCount
}

def countAmplicons(String filePath) {
    // Read all lines from the file
    def lines = Files.readAllLines(Paths.get(filePath))

    // Extract the fourth column, replace hyphens, and remove _LEFT and _RIGHT
    def processedItems = lines.collect { line ->
        def columns = line.split('\t')
        if (columns.size() >= 4) {
            def item = columns[3]
            item = item.replace('-', '_')
            item = item.replace('_LEFT', '').replace('_RIGHT', '')
            return item
        }
        return null
    }.findAll { it != null }

    // Convert to a set to find unique items
    def uniqueItems = processedItems.toSet()

    // Return the number of unique items
    return uniqueItems.size()
}
