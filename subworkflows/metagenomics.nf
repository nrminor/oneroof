include {
    SKETCH_DATABASE_KMERS ;
    SKETCH_SAMPLE_KMERS ;
    PROFILE_SAMPLE ;
    OVERLAY_TAXONOMY
} from "../modules/sylph"

workflow METAGENOMICS {

    take:
    ch_ref
    ch_sample_reads
    ch_taxonomy_files

    main:
    SKETCH_DATABASE_KMERS(ch_ref)

    SKETCH_SAMPLE_KMERS(ch_sample_reads)

    PROFILE_SAMPLE(
        SKETCH_SAMPLE_KMERS.out.combine(SKETCH_DATABASE_KMERS.out)
    )

    OVERLAY_TAXONOMY(
        PROFILE_SAMPLE.out.combine(ch_taxonomy_files)
    )

    // emit:

}
