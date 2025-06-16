process CHECK_DATASET {

    errorStrategy 'ignore'

    output:
    val "${params.nextclade_dataset}"

    when:
    params.nextclade_dataset

    script:
    '''
    dataset_options=$(nextclade dataset list)

    query_name_line="\"name\"=\"!{params.nextclade_dataset}\""

    if grep -Fq "$query_name_line" <<< "$dataset_options"; then
        echo "$query_name_line found in supported Nextclade datasets. It is safe to proceed."
    else
        echo "$query_name_line not found in the list of supported Nextclade datasets."
    fi
    '''

}

process DOWNLOAD_DATASET {

    storeDir params.nextclade_cache

    input:
    val valid_dataset

    output:
    path dataset_label

    script:
    dataset_label = params.nextclade_dataset.split("/")[-1]
    """
    nextclade dataset get \
    --name ${valid_dataset} \
    --output-dir ${dataset_label}
    """

}

process RUN_NEXTCLADE {

    publishDir params.nextclade, mode: 'copy'

    input:
    path sequences
    each path(nextclade_dataset)

    output:
    path "${label}/"

    script:
    label = file(sequences.toString()).getSimpleName()
    """
    nextclade run \
    --input-dataset ${nextclade_dataset} \
    --output-all=${label}/ \
    ${sequences}
    """

}
