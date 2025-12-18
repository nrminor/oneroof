include { GATHER_NANOPORE } from "../subworkflows/gather_nanopore"
include { PRIMER_HANDLING } from "../subworkflows/primer_handling"
include { ALIGNMENT } from "../subworkflows/alignment"
include { HAPLOTYPING } from "../subworkflows/haplotyping"
include { CONSENSUS } from "../subworkflows/consensus_calling"
include { VARIANTS } from "../subworkflows/variant_calling"
include { METAGENOMICS } from "../subworkflows/metagenomics"
include { PHYLO } from "../subworkflows/phylo"
include { SLACK_ALERT } from "../subworkflows/slack_alert"

workflow NANOPORE {

    /* */

    take:
        ch_primer_bed
        ch_refseq
        ch_refgbk
        ch_contam_fasta
        ch_snpeff_config
        ch_metagenome_ref
        ch_primer_tsv
        ch_sylph_tax_db
        ch_meta_ref_link

    main:
        assert params.platform == "ont"

        GATHER_NANOPORE ( ch_contam_fasta )

        if ( params.primer_bed || params.primer_tsv ) {

            PRIMER_HANDLING(
                GATHER_NANOPORE.out,
                ch_primer_bed,
                ch_refseq,
                ch_primer_tsv
            )

            METAGENOMICS(
                PRIMER_HANDLING.out.merged_reads,
                ch_metagenome_ref,
                ch_meta_ref_link,
                ch_sylph_tax_db
            )

            alignment_outputs = ALIGNMENT (
                PRIMER_HANDLING.out.merged_reads,
                ch_refseq
            )

        } else {

            METAGENOMICS(
                GATHER_NANOPORE.out,
                ch_metagenome_ref,
                ch_meta_ref_link,
                ch_sylph_tax_db
            )

            alignment_outputs = ALIGNMENT (
                GATHER_NANOPORE.out,
                ch_refseq
            )

        }

        CONSENSUS (
            alignment_outputs.index
        )

        variant_outputs = VARIANTS (
            alignment_outputs.index,
            ch_refseq,
            ch_refgbk,
            ch_snpeff_config
        )

        if ( params.devider_preset ) {

            HAPLOTYPING (
                alignment_outputs.index,
                variant_outputs.vcf,
                ch_refseq
            )

        }

        PHYLO (
            CONSENSUS.out
        )

        SLACK_ALERT(
            alignment_outputs.coverage_summary,
            CONSENSUS.out.collect(),
            variant_outputs.merge_vcf_files
        )

        // Assemble OneRoof report from collected metrics
        if ( params.generate_json_report ) {

            // Build QC thresholds map from params
            qc_thresholds = [
                coverage_pass: params.qc_coverage_pass,
                coverage_warn: params.qc_coverage_warn,
                completeness_pass: params.qc_completeness_pass,
                completeness_warn: params.qc_completeness_warn
            ]

            // Get reference name from refseq filename
            reference_name = file(params.refseq).baseName

            // Collect all coverage metrics JSON files
            ch_metrics = alignment_outputs.coverage_metrics
                .map { sample_id, metrics_json -> metrics_json }
                .collect()

            // MultiQC config template (use default if not specified)
            ch_multiqc_template = params.multiqc_config_template
                ? Channel.fromPath(params.multiqc_config_template)
                : Channel.fromPath("${projectDir}/assets/multiqc_config.yaml")

            ASSEMBLE_REPORT(
                ch_metrics,
                params.platform,
                reference_name,
                qc_thresholds,
                ch_multiqc_template.first()
            )

        }

}
