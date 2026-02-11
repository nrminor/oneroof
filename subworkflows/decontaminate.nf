include { INDEX_CONTAMINANTS ; RUN_DEACON ; GET_INDEX } from "../modules/deacon"

workflow DECONTAMINATE {
    take:
    ch_fastqs
    ch_contam_fasta

    main:
    def has_local_fasta = params.contam_fasta && file(params.contam_fasta).isFile()

    ch_built_index = has_local_fasta
        ? INDEX_CONTAMINANTS(ch_contam_fasta)
        : Channel.empty()

    ch_fetched_index = (!has_local_fasta && params.contam_link)
        ? GET_INDEX()
        : Channel.empty()

    ch_index = ch_built_index.mix(ch_fetched_index)

    ch_decontaminated = RUN_DEACON(ch_fastqs.combine(ch_index))

    emit:
    ch_decontaminated
}
