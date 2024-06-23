process RENAME_BY_BARCODE {

    /* */

	errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
	maxRetries 1

    cpus 8

    inputs:
    path files_to_rename
    path samplesheet

    output:
    path "*"

    shell:
    '''
    while IFS=$'\t' read -r barcode sample_id; do
        # Find files in the current directory with the base_before name
        for old_filename in "${barcode}".*; do
            if [ -e "$old_filename" ]; then
                # Extract the file extension
                extension="${old_filename##*.}"

                # Construct the new filename with the same extension
                new_file_with_extension="${sample_id}.${extension}"

                # Rename the file
                mv "$old_filename" "$new_file_with_extension" &

                echo "Renamed $old_filename to $new_file_with_extension"
            else
                echo "File $old_filename not found."
            fi
        done
    done < !{samplesheet}

    wait
    '''
}
