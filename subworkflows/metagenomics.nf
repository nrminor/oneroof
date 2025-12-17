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
        ch_metagenome_link = DOWNLOAD_DB_LINK(ch_meta_ref_link)
        ch_prebuilt_ref = ch_metagenome_link.filter { fa ->
            !file(fa).getBaseName().endsWith("fasta") || !file(fa).getBaseName().endsWith("fasta")
        }
    }
    else {
        ch_prebuilt_ref = ch_metagenome_ref.filter { fa ->
            !file(fa).getBaseName().endsWith("fasta") || !file(fa).getBaseName().endsWith("fasta")
        }
    }

    ch_custom_ref = ch_metagenome_ref.filter { fa ->
        file(fa).getBaseName().endsWith("fasta") || file(fa).getBaseName().endsWith("fasta")
    }

    SKETCH_DATABASE_KMERS(ch_custom_ref)

    ch_sylph_db_queue = SKETCH_DATABASE_KMERS.out.mix(ch_prebuilt_ref)

    SKETCH_SAMPLE_KMERS(ch_sample_reads)

    CLASSIFY_SAMPLE(
        SKETCH_SAMPLE_KMERS.out.combine(ch_sylph_db_queue)
    )

    ch_classified = CLASSIFY_SAMPLE.out.filter { _id, sketch -> file(sketch).countLines() > 1 }

    if (params.sylph_tax_db) {

    SYLPH_TAX_DOWNLOAD(sylph_tax_trig_ch)

    ch_sylph_taxonomy_dir = SYLPH_TAX_DOWNLOAD.out

        MERGE_TAXONOMY(
            OVERLAY_TAXONOMY.out
        )
    }
}
