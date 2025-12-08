include {
    SKETCH_DATABASE_KMERS ;
    SKETCH_SAMPLE_KMERS ;
    CLASSIFY_SAMPLE ;
    OVERLAY_TAXONOMY ;
    SYLPH_TAX_DOWNLOAD ;
    MERGE_TAXONOMY ;
    DOWNLOAD_DB_LINK
} from "../modules/sylph"

workflow METAGENOMICS {

    take:
    ch_metagenome_ref
    ch_sylph_tax_db
    ch_meta_ref_link
    ch_sample_reads

    main:

     if (params.sylph_db_link) {
        // This captures the output from the process into a new variable
        ch_downloaded_db = DOWNLOAD_DB_LINK(ch_meta_ref_link)

       ch_prebuilt_ref = ch_downloaded_db.filter { fa ->
        !file(fa).getBaseName().endsWith("syldb") || !file(fa).getBaseName().endsWith("syldb")
    }
     } else {
        ch_prebuilt_ref = ch_metagenome_ref.filter { fa ->
        !file(fa).getBaseName().endsWith("syldb") || !file(fa).getBaseName().endsWith("syldb")
    }
     }

    ch_custom_ref = ch_metagenome_ref.filter { fa ->
        file(fa).getBaseName().endsWith("fasta") || file(fa).getBaseName().endsWith("fasta")
    }

    SKETCH_DATABASE_KMERS(ch_custom_ref)

    ch_sylph_db_queue = SKETCH_DATABASE_KMERS.out.mix(ch_prebuilt_ref)

    SKETCH_SAMPLE_KMERS(ch_sample_reads.combine(ch_sylph_db_queue))
    
    CLASSIFY_SAMPLE(
        SKETCH_SAMPLE_KMERS.out.combine(ch_sylph_db_queue)
    )

    ch_classified = CLASSIFY_SAMPLE.out.filter { _id, sketch -> file(sketch).countLines() > 1 }


    // SYLPH_TAX_DOWNLOAD(ch_classified)

    // OVERLAY_TAXONOMY(
    //     ch_classified.combine(ch_sylph_tax_db).combine(SYLPH_TAX_DOWNLOAD.out)
    // )
    // Channel
    //     .of(1)
    //     .set { sylph_tax_trig_ch }

    // // Run SYLPH_TAX_DOWNLOAD once, driven by the trigger
    // SYLPH_TAX_DOWNLOAD(sylph_tax_trig_ch)

    // // Capture its output channel (assumed to be the taxonomy dir or files)
    // SYLPH_TAX_DOWNLOAD.out .map { tax_output -> tax_output.toString() }   // turn into a simple val
    //     .set { ch_tax_ready }
    //------------------------------------------------------------------
    // Overlay taxonomy per sample, combining:
    //   - ch_classified: (sample_id, tsv_path)
    //   - ch_sylph_taxonomy_dir: shared Sylph taxonomy directory
    //   - ch_sylph_tax_db: DB path
    //------------------------------------------------------------------

    // OVERLAY_TAXONOMY(
    //     ch_classified
    //         .combine(ch_sylph_taxonomy_dir)
    //         .combine(ch_sylph_tax_db)
    // )

    Channel
    .of(1)
    .set { sylph_tax_trig_ch }

SYLPH_TAX_DOWNLOAD(sylph_tax_trig_ch)

// this will now be a channel that emits the *directory* sylph_taxonomy
SYLPH_TAX_DOWNLOAD.out
    .set { ch_sylph_taxonomy_dir }


    // OVERLAY_TAXONOMY(
    //     ch_classified
    //         .combine(ch_tax_ready)
    //         .combine(ch_sylph_tax_db)
    // )
    OVERLAY_TAXONOMY(
    ch_classified          // (sample_id, tsv_path)
        .combine(ch_sylph_taxonomy_dir)  // (taxonomy_dir)
        .combine(ch_sylph_tax_db)        // (input_db)
)


    MERGE_TAXONOMY(
        OVERLAY_TAXONOMY.out
    )

}