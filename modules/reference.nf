/*
 * Reference sequence resolution processes
 *
 * These processes resolve reference inputs that can be either:
 * - Local file paths (validated and normalized via *_LOCAL processes)
 * - NCBI accessions (fetched via Entrez API via *_ACCESSION processes)
 *
 * Sequences are normalized to contain only valid IUPAC nucleotide characters,
 * with invalid characters replaced by 'N'.
 *
 * The workflow must determine which process to call based on whether the input
 * is a local file path or an NCBI accession. This separation ensures proper
 * Nextflow file staging for local files.
 *
 * Usage in workflows:
 *   if (file(params.refseq).exists()) {
 *       RESOLVE_REFSEQ_LOCAL(Channel.fromPath(params.refseq))
 *       ch_refseq = RESOLVE_REFSEQ_LOCAL.out.refseq
 *   } else {
 *       RESOLVE_REFSEQ_ACCESSION(Channel.value(params.refseq))
 *       ch_refseq = RESOLVE_REFSEQ_ACCESSION.out.refseq
 *   }
 */


/*
 * Resolve a FASTA reference from a local file path.
 *
 * The file is properly staged by Nextflow, validated, and normalized.
 * This process is REQUIRED - if it fails, the pipeline should stop.
 */
process RESOLVE_REFSEQ_LOCAL {

    tag "${local_ref.name}"

    publishDir "${params.results}/reference_assets", mode: 'copy', overwrite: true

    errorStrategy { task.attempt < 3 ? 'retry' : 'finish' }
    maxRetries 2

    input:
    path local_ref

    output:
    path "*.fasta", emit: refseq

    script:
    def output_name = local_ref.baseName + '.resolved.fasta'
    """
    fetch_reference.py fasta "${local_ref}" --output "${output_name}"
    """
}


/*
 * Resolve a FASTA reference from an NCBI accession.
 *
 * The accession is fetched via Entrez API, validated, and normalized.
 * This process is REQUIRED - if it fails, the pipeline should stop.
 */
process RESOLVE_REFSEQ_ACCESSION {

    tag "${accession}"

    publishDir "${params.results}/reference_assets", mode: 'copy', overwrite: true

    errorStrategy { task.attempt < 3 ? 'retry' : 'finish' }
    maxRetries 2

    input:
    val accession

    output:
    path "*.fasta", emit: refseq

    script:
    def output_name = accession.replaceAll(/[^A-Za-z0-9._-]/, '_') + '.fasta'
    """
    fetch_reference.py fasta "${accession}" --output "${output_name}"
    """
}


/*
 * Resolve a GenBank reference from a local file path.
 *
 * The file is properly staged by Nextflow, validated, and normalized.
 * This process is OPTIONAL - if it fails, the pipeline continues without
 * variant annotation.
 */
process RESOLVE_REF_GBK_LOCAL {

    tag "${local_ref.name}"

    publishDir "${params.results}/reference_assets", mode: 'copy', overwrite: true

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    input:
    path local_ref

    output:
    path "*.gbk", emit: ref_gbk, optional: true

    script:
    def output_name = local_ref.baseName + '.resolved.gbk'
    """
    fetch_reference.py genbank "${local_ref}" --output "${output_name}"
    """
}


/*
 * Resolve a GenBank reference from an NCBI accession.
 *
 * The accession is fetched via Entrez API, validated, and normalized.
 * This process is OPTIONAL - if it fails, the pipeline continues without
 * variant annotation.
 */
process RESOLVE_REF_GBK_ACCESSION {

    tag "${accession}"

    publishDir "${params.results}/reference_assets", mode: 'copy', overwrite: true

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    input:
    val accession

    output:
    path "*.gbk", emit: ref_gbk, optional: true

    script:
    def output_name = accession.replaceAll(/[^A-Za-z0-9._-]/, '_') + '.gbk'
    """
    fetch_reference.py genbank "${accession}" --output "${output_name}"
    """
}
