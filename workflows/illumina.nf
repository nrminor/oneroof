include { GATHER_ILLUMINA } from "../subworkflows/gather_illumina"
include { ILLUMINA_CORRECTION } from "../subworkflows/illumina_correction"
include { PRIMER_HANDLING } from "../subworkflows/primer_handling"
include { ALIGNMENT } from "../subworkflows/alignment"
include { CONSENSUS } from "../subworkflows/consensus_calling"
include { VARIANTS } from "../subworkflows/variant_calling"
include { METAGENOMICS } from "../subworkflows/metagenomics"
include { PHYLO } from "../subworkflows/phylo"
include { SLACK_ALERT } from "../subworkflows/slack_alert"
include { ASSEMBLE_REPORT } from "../modules/reporting"

workflow ILLUMINA {

    take:
        ch_primer_bed
        ch_refseq
        ch_ref_gbk
        ch_contam_fasta
        ch_snpeff_config
        ch_metagenome_ref
        ch_primer_tsv
        ch_sylph_tax_db
        ch_meta_ref_link

    main:
        assert params.platform == "illumina"
        assert params.illumina_fastq_dir != "" :
        "Please double check that a directory of Illumina FASTQs or Nanopore POD5s is provided."
        assert file( params.illumina_fastq_dir ).isDirectory() :
        "The provided Illumina FASTQ directory ${params.illumina_fastq_dir} does not exist."

        GATHER_ILLUMINA ( )

        ILLUMINA_CORRECTION (
            GATHER_ILLUMINA.out,
            ch_contam_fasta
        )

        if ( params.primer_bed || params.primer_tsv ) {

            PRIMER_HANDLING (
                ILLUMINA_CORRECTION.out,
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
                ILLUMINA_CORRECTION.out,
                ch_metagenome_ref,
                ch_meta_ref_link,
                ch_sylph_tax_db
            )

            alignment_outputs = ALIGNMENT (
                ILLUMINA_CORRECTION.out,
                ch_refseq
            )

        }


        CONSENSUS (
            alignment_outputs.index
        )

        variant_outputs = VARIANTS (
            alignment_outputs.index,
            ch_refseq,
            ch_ref_gbk,

            ch_snpeff_config
        )

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
