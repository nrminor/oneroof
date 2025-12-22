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
include { MULTIQC } from "../modules/multiqc"
include { RESOLVE_REFSEQ_LOCAL ; RESOLVE_REFSEQ_ACCESSION } from "../modules/reference"
include { RESOLVE_REF_GBK_LOCAL ; RESOLVE_REF_GBK_ACCESSION } from "../modules/reference"

workflow ILLUMINA {

    take:
        ch_primer_bed
        ch_contam_fasta
        ch_snpeff_config
        ch_metagenome_ref
        ch_primer_tsv
        ch_sylph_tax_db
        ch_meta_ref_link

    main:
        assert params.platform == "illumina"

        // Resolve reference FASTA (local file or NCBI accession)
        if (file(params.refseq).exists()) {
            RESOLVE_REFSEQ_LOCAL(Channel.fromPath(params.refseq))
            ch_refseq = RESOLVE_REFSEQ_LOCAL.out.refseq
        } else {
            RESOLVE_REFSEQ_ACCESSION(Channel.value(params.refseq))
            ch_refseq = RESOLVE_REFSEQ_ACCESSION.out.refseq
        }

        // Resolve reference GenBank (local file or NCBI accession)
        if (params.ref_gbk) {
            if (file(params.ref_gbk).exists()) {
                RESOLVE_REF_GBK_LOCAL(Channel.fromPath(params.ref_gbk))
                ch_ref_gbk = RESOLVE_REF_GBK_LOCAL.out.ref_gbk
            } else {
                RESOLVE_REF_GBK_ACCESSION(Channel.value(params.ref_gbk))
                ch_ref_gbk = RESOLVE_REF_GBK_ACCESSION.out.ref_gbk
            }
        } else {
            ch_ref_gbk = Channel.empty()
        }
        assert params.illumina_fastq_dir != "" :
        "Please double check that a directory of Illumina FASTQs or Nanopore POD5s is provided."
        assert file( params.illumina_fastq_dir ).isDirectory() :
        "The provided Illumina FASTQ directory ${params.illumina_fastq_dir} does not exist."

        GATHER_ILLUMINA ( )

        correction_outputs = ILLUMINA_CORRECTION (
            GATHER_ILLUMINA.out,
            ch_contam_fasta
        )

        if ( params.primer_bed || params.primer_tsv ) {

            PRIMER_HANDLING (
                correction_outputs.reads,
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

            ch_amplicon_summary = PRIMER_HANDLING.out.amplicon_summary

        } else {

            METAGENOMICS(
                correction_outputs.reads,
                ch_metagenome_ref,
                ch_meta_ref_link,
                ch_sylph_tax_db
            )

            alignment_outputs = ALIGNMENT (
                correction_outputs.reads,
                ch_refseq
            )

            ch_amplicon_summary = Channel.empty()

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
            CONSENSUS.out.concat
        )


        SLACK_ALERT(
            alignment_outputs.coverage_summary,
            CONSENSUS.out.concat.collect(),
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

            // Collect all metrics JSON files from all subworkflows
            // Each metrics channel emits tuple(sample_id, metrics_json)
            // We extract just the JSON file and collect all into a single channel
            ch_coverage_metrics = alignment_outputs.coverage_metrics
                .map { sample_id, metrics_json -> metrics_json }

            ch_alignment_metrics = alignment_outputs.alignment_metrics
                .map { sample_id, metrics_json -> metrics_json }

            ch_variant_metrics = variant_outputs.variant_metrics
                .map { sample_id, metrics_json -> metrics_json }

            ch_consensus_metrics = CONSENSUS.out.consensus_metrics
                .map { sample_id, metrics_json -> metrics_json }

            ch_metagenomics_metrics = METAGENOMICS.out.metagenomics_metrics
                .map { sample_id, metrics_json -> metrics_json }

            // Combine all metrics channels (no haplotyping for Illumina)
            ch_metrics = ch_coverage_metrics
                .mix(ch_alignment_metrics)
                .mix(ch_variant_metrics)
                .mix(ch_consensus_metrics)
                .mix(ch_metagenomics_metrics)
                .collect()

            // MultiQC config template (use default if not specified)
            multiqc_template = params.multiqc_config_template
                ? file(params.multiqc_config_template)
                : file("${projectDir}/conf/multiqc_config.yaml")

            // Amplicon summary (placeholder if primers not provided)
            ch_amplicon_summary_for_report = ch_amplicon_summary
                .ifEmpty(file("NO_AMPLICON_SUMMARY"))

            ASSEMBLE_REPORT(
                ch_metrics,
                params.platform,
                reference_name,
                qc_thresholds,
                multiqc_template,
                ch_amplicon_summary_for_report
            )

            // Run MultiQC with FastQC outputs and OneRoof custom content
            if ( params.generate_multiqc ) {
                MULTIQC(
                    correction_outputs.fastqc_zip.collect(),
                    ASSEMBLE_REPORT.out.multiqc_tsv.collect(),
                    ASSEMBLE_REPORT.out.multiqc_config
                )
            }

        }

}
