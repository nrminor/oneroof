/*
 * Metagenomics subworkflow
 *
 * Performs metagenomic profiling using Sylph and optionally adds
 * taxonomic annotations using sylph-tax.
 *
 * Input options for the reference database:
 *   - Pre-built .syldb file (params.meta_ref)
 *   - Custom FASTA to be sketched (params.meta_ref)
 *   - URL to download a .syldb (params.sylph_db_link)
 *
 * Taxonomy is added if params.sylph_tax_db is set to a valid
 * taxonomy identifier (e.g., "GTDB_r220", "IMGVR_4.1").
 */

include {
    DOWNLOAD_SYLPH_DB ;
    SKETCH_DATABASE ;
    SKETCH_SAMPLES ;
    PROFILE_SAMPLES ;
    DOWNLOAD_TAXONOMY ;
    ADD_TAXONOMY ;
    MERGE_TAXONOMY
} from "../modules/sylph"

include { EXTRACT_METAGENOMICS_METRICS } from "../modules/reporting"


workflow METAGENOMICS {

    take:
    ch_sample_reads      // tuple(sample_id, reads)
    ch_meta_ref          // path to .syldb or .fasta, or empty channel
    ch_sylph_db_link     // URL string to download .syldb, or empty channel
    ch_sylph_tax_db      // taxonomy identifier string (e.g., "GTDB_r220"), or empty channel

    main:

    // -------------------------------------------------------------------------
    // Step 1: Resolve the sylph database
    // -------------------------------------------------------------------------
    // Three options:
    //   a) Download from URL if sylph_db_link is provided
    //   b) Use pre-built .syldb file directly
    //   c) Sketch custom FASTA into .syldb

    ch_downloaded_db = params.sylph_db_link
        ? DOWNLOAD_SYLPH_DB(ch_sylph_db_link)
        : Channel.empty()

    // Branch the meta_ref channel by file type
    ch_meta_ref
        .branch { ref ->
            prebuilt: ref.name.endsWith('.syldb')
            custom: true  // .fasta, .fa, .fna, etc.
        }
        .set { ch_ref_branched }

    // Sketch custom FASTA references
    ch_sketched_db = SKETCH_DATABASE(ch_ref_branched.custom)

    // Combine all database sources into one channel
    ch_syldb = ch_downloaded_db
        .mix(ch_ref_branched.prebuilt)
        .mix(ch_sketched_db)

    // -------------------------------------------------------------------------
    // Step 2: Sketch sample reads (only if we have a database to profile against)
    // -------------------------------------------------------------------------
    // Combine samples with database first - if ch_syldb is empty, nothing runs
    ch_samples_with_db = ch_sample_reads.combine(ch_syldb)

    SKETCH_SAMPLES(
        ch_samples_with_db.map { sample_id, reads, _syldb -> tuple(sample_id, reads) }
    )

    // -------------------------------------------------------------------------
    // Step 3: Profile samples against database
    // -------------------------------------------------------------------------
    // Rejoin sketches with database for profiling
    ch_sketches_with_db = SKETCH_SAMPLES.out.sketch
        .combine(ch_syldb)

    PROFILE_SAMPLES(ch_sketches_with_db)

    // Filter out empty results (only header line)
    ch_profiles = PROFILE_SAMPLES.out.profile
        .filter { _sample_id, tsv -> tsv.countLines() > 1 }

    // Extract metagenomics metrics from non-empty profiles
    EXTRACT_METAGENOMICS_METRICS(ch_profiles)

    // -------------------------------------------------------------------------
    // Step 4: Add taxonomy (optional)
    // -------------------------------------------------------------------------
    ch_taxprofiles = Channel.empty()
    ch_merged = Channel.empty()

    if (params.sylph_tax_db) {
        // Download taxonomy metadata (cached via storeDir)
        DOWNLOAD_TAXONOMY()

        // Add taxonomy to each profile
        ADD_TAXONOMY(
            ch_profiles,
            DOWNLOAD_TAXONOMY.out.taxonomy_dir,
            ch_sylph_tax_db
        )

        ch_taxprofiles = ADD_TAXONOMY.out.taxprofile

        // -------------------------------------------------------------------------
        // Step 5: Merge all taxonomy profiles
        // -------------------------------------------------------------------------
        MERGE_TAXONOMY(
            ch_taxprofiles
                .map { _sample_id, sylphmpa -> sylphmpa }
                .collect()
        )

        ch_merged = MERGE_TAXONOMY.out.merged
    }

    emit:
    profiles             = ch_profiles
    taxonomy             = ch_taxprofiles
    merged               = ch_merged
    metagenomics_metrics = EXTRACT_METAGENOMICS_METRICS.out
}
