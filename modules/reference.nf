/*
 * Reference sequence resolution processes
 *
 * These processes resolve reference inputs that can be either:
 * - Local file paths (validated and normalized)
 * - NCBI accessions (fetched via Entrez API)
 *
 * Sequences are normalized to contain only valid IUPAC nucleotide characters,
 * with invalid characters replaced by 'N'.
 *
 * Usage in main.nf:
 *   RESOLVE_REFSEQ(Channel.value(params.refseq))
 *   RESOLVE_REF_GBK(Channel.value(params.ref_gbk))
 */


/*
 * Resolve a FASTA reference from a local path or NCBI accession.
 *
 * This process is REQUIRED - if it fails, the pipeline should stop.
 * Uses 'finish' error strategy to allow other running processes to complete
 * before terminating.
 */
process RESOLVE_REFSEQ {

    tag "${ref_input}"

    // Use standard Nextflow caching (respects -resume) rather than storeDir
    // to avoid cross-run cache conflicts when switching references
    publishDir "${params.results}/reference_assets", mode: 'copy', overwrite: true

    errorStrategy { task.attempt < 3 ? 'retry' : 'finish' }
    maxRetries 2

    input:
    val ref_input  // Local path string or NCBI accession

    output:
    path "*.fasta", emit: refseq

    script:
    // Determine output filename based on input
    // For accessions: NC_045512.2.fasta
    // For local files: preserve original basename
    def output_name = ref_input.contains('/') || ref_input.contains('\\')
        ? file(ref_input).baseName + '.fasta'
        : ref_input.replaceAll(/[^A-Za-z0-9._-]/, '_') + '.fasta'

    """
    fetch_reference.py fasta "${ref_input}" --output "${output_name}"
    """
}


/*
 * Resolve a GenBank reference from a local path or NCBI accession.
 *
 * This process is OPTIONAL - if it fails, the pipeline continues without
 * variant annotation. Uses 'ignore' error strategy as final fallback.
 */
process RESOLVE_REF_GBK {

    tag "${ref_input}"

    // Use standard Nextflow caching (respects -resume) rather than storeDir
    // to avoid cross-run cache conflicts when switching references
    publishDir "${params.results}/reference_assets", mode: 'copy', overwrite: true

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    input:
    val ref_input  // Local path string or NCBI accession

    output:
    path "*.gbk", emit: ref_gbk, optional: true

    script:
    // Determine output filename based on input
    // For accessions: NC_045512.2.gbk
    // For local files: preserve original basename
    def output_name = ref_input.contains('/') || ref_input.contains('\\')
        ? file(ref_input).baseName + '.gbk'
        : ref_input.replaceAll(/[^A-Za-z0-9._-]/, '_') + '.gbk'

    """
    fetch_reference.py genbank "${ref_input}" --output "${output_name}"
    """
}
