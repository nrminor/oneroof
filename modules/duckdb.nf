process FASTA_TO_PARQUET {

    tag "${sample_id}"

    input:
    tuple val(sample_id), path(amplicons)

    output:
    tuple val(sample_id), path("${sample_id}.parquet")

    script:
    """
    seqkit fx2tab ${amplicons} \
    | duckdb :memory: <<'SQL'
        copy (
            select * from read_csv(
                '/dev/stdin', 
                columns = {'id': 'TEXT', 'seq': 'TEXT', 'qual': 'TEXT'}
            )
        ) 
        to '${sample_id}.parquet'
        with (format 'parquet', compression 'ZSTD');
    SQL
    """
}

process QUERY_FOR_PRIMERS {

    cpus 12

    input:
    tuple val(sample_id), path(parquet), val(amplicon), val(fwd), val(fwd_rc), val(rev), val(rev_rc)

    output:
    tuple val(sample_id), path("*.parquet")

    script:
    """
    duckdb :memory: <<SQL
    
    SET threads TO ${task.cpus};
    
    -- Read input Parquet (assumed to have 'id' and 'seq' columns)
    WITH input AS (
      SELECT * FROM '${parquet}'
    ),

    -- Window sequences for forward primer length
    fwd_substrings AS (
      SELECT
        id,
        i AS pos,
        SUBSTRING(seq, i, LENGTH('${fwd}')) AS window
      FROM input,
           generate_series(1, LENGTH(seq) - LENGTH('${fwd}') + 1) AS t(i)
    ),

    -- Window sequences for reverse primer length
    rev_substrings AS (
      SELECT
        id,
        i AS pos,
        SUBSTRING(seq, i, LENGTH('${rev}')) AS window
      FROM input,
           generate_series(1, LENGTH(seq) - LENGTH('${rev}') + 1) AS t(i)
    ),


    -- Forward primer matches (with position + strand info)
    fwd_matches AS (
      SELECT id, pos, primer, mismatches
          FROM (
          SELECT id, pos, 'FWD' AS primer, mismatches
          FROM (
            SELECT id, pos, hamming(window, '${fwd}') AS mismatches
            FROM fwd_substrings
          ) WHERE mismatches <= ${params.max_mismatches}

          UNION ALL

          SELECT id, pos, 'FWD_RC' AS primer, mismatches
          FROM (
            SELECT id, pos, hamming(window, '${fwd_rc}') AS mismatches
            FROM fwd_substrings
          ) WHERE mismatches <= ${params.max_mismatches}
      )
      QUALIFY ROW_NUMBER() OVER (PARTITION BY id ORDER BY pos) = 1
    ),

    -- Reverse primer matches (with position + strand info)
    rev_matches AS (
      SELECT id, pos, primer, mismatches
      FROM (
          SELECT id, pos, 'REV' AS primer, mismatches
          FROM (
            SELECT id, pos, hamming(window, '${rev}') AS mismatches
            FROM fwd_substrings
          ) WHERE mismatches <= ${params.max_mismatches}

          UNION ALL

          SELECT id, pos, 'REV_RC' AS primer, mismatches
          FROM (
            SELECT id, pos, hamming(window, '${rev_rc}') AS mismatches
            FROM fwd_substrings
          ) WHERE mismatches <= ${params.max_mismatches}
      )
      QUALIFY ROW_NUMBER() OVER (PARTITION BY id ORDER BY pos) = 1
    ),

    -- Only keep rows that matched both forward and reverse
    matches AS (
      SELECT id
      FROM fwd_matches
      INTERSECT
      SELECT id
      FROM rev_matches
    ),

    -- run filtering
    filtered AS (
        SELECT i.*
        FROM input i
        JOIN matches m ON i.id = m.id
    )

    -- put the results of the whole query into a new parquet file
    COPY filtered TO '${sample_id}.${amplicon}.parquet' (FORMAT PARQUET, COMPRESSION 'ZSTD');

    SQL
    """
}

process PRIMER_SEARCH_ENGINE {

    tag "${amplicon}, ${sample_count} samples"

    cpus 12

    input:
    tuple val(sample_count), path(parquet), val(amplicon), val(fwd), val(fwd_rc), val(rev), val(rev_rc)

    output:
    tuple path("*.parquet")

    script:
    """
    duckdb :memory: <<SQL
    
    SET threads TO ${task.cpus};
    
    -- Read input Parquet (assumed to have 'id' and 'seq' columns)
    WITH input AS (
      SELECT * FROM '${parquet}'
    ),

    -- Window sequences for forward primer length
    fwd_substrings AS (
      SELECT
        id,
        i AS pos,
        SUBSTRING(seq, i, LENGTH('${fwd}')) AS window
      FROM input,
           generate_series(1, LENGTH(seq) - LENGTH('${fwd}') + 1) AS t(i)
    ),

    -- Window sequences for reverse primer length
    rev_substrings AS (
      SELECT
        id,
        i AS pos,
        SUBSTRING(seq, i, LENGTH('${rev}')) AS window
      FROM input,
           generate_series(1, LENGTH(seq) - LENGTH('${rev}') + 1) AS t(i)
    ),


    -- Forward primer matches (with position + strand info)
    fwd_matches AS (
      SELECT id, pos, primer, mismatches
          FROM (
          SELECT id, pos, 'FWD' AS primer, mismatches
          FROM (
            SELECT id, pos, hamming(window, '${fwd}') AS mismatches
            FROM fwd_substrings
          ) WHERE mismatches <= ${params.max_mismatches}

          UNION ALL

          SELECT id, pos, 'FWD_RC' AS primer, mismatches
          FROM (
            SELECT id, pos, hamming(window, '${fwd_rc}') AS mismatches
            FROM fwd_substrings
          ) WHERE mismatches <= ${params.max_mismatches}
      )
      QUALIFY ROW_NUMBER() OVER (PARTITION BY id ORDER BY pos) = 1
    ),

    -- Reverse primer matches (with position + strand info)
    rev_matches AS (
      SELECT id, pos, primer, mismatches
      FROM (
          SELECT id, pos, 'REV' AS primer, mismatches
          FROM (
            SELECT id, pos, hamming(window, '${rev}') AS mismatches
            FROM fwd_substrings
          ) WHERE mismatches <= ${params.max_mismatches}

          UNION ALL

          SELECT id, pos, 'REV_RC' AS primer, mismatches
          FROM (
            SELECT id, pos, hamming(window, '${rev_rc}') AS mismatches
            FROM fwd_substrings
          ) WHERE mismatches <= ${params.max_mismatches}
      )
      QUALIFY ROW_NUMBER() OVER (PARTITION BY id ORDER BY pos) = 1
    ),

    -- Only keep rows that matched both forward and reverse
    matches AS (
      SELECT id
      FROM fwd_matches
      INTERSECT
      SELECT id
      FROM rev_matches
    ),

    -- run filtering
    filtered AS (
        SELECT i.*
        FROM input i
        JOIN matches m ON i.id = m.id
    )

    -- put the results of the whole query into a new parquet file
    COPY filtered TO '${amplicon}.multisample.parquet' (FORMAT PARQUET, COMPRESSION 'ZSTD');

    SQL
    """
}

process FILTER_BY_QUERY {

    tag "${sample_id}"

    input:
    tuple val(sample_id), path(parquet)

    script:
    """
    duckdb :memory: <<'SQL' |
        copy (
            select * from ${parquet}
            where length(sequence) >= ${params.min_len}
            and length(sequence) <= ${params.max_len}
        )
        to '${sample_id}.filtered.parquet'
        with (format 'parquet', compression 'ZSTD');
    SQL
    """
}

process SPLIT_HIVE_PARQUET {
    input:
    path hive_dir

    output:
    path "*.parquet"

    script:
    """
    # Find all sample_ids in the hive
    duckdb /dev/null <<SQL
      CREATE TEMP TABLE sample_ids AS
      SELECT DISTINCT sample_id
      FROM '${hive_dir}/*.parquet';
  
      COPY (
        SELECT DISTINCT sample_id
        FROM sample_ids
      ) TO 'sample_list.tsv' (DELIMITER '\t', HEADER FALSE);
    SQL

    # Loop over each sample_id and extract its partition
    while read sample_id; do
      duckdb /dev/null <<SQL
      COPY (
        SELECT *
        FROM '${hive_dir}/*.parquet'
        WHERE sample_id = '\$sample_id'
      ) TO '\$sample_id.parquet' (
        FORMAT PARQUET,
        COMPRESSION ZSTD
      );
      SQL
    done < sample_list.tsv
    """
}

process PARQUET_TO_FASTA {

    tag "${sample_id}"

    input:
    tuple val(sample_id), path(parquet)

    output:
    tuple val(sample_id), path("${sample_id}.fastq.gz")

    script:
    """
    duckdb :memory: <<'SQL' |
        copy (
            select * from '${parquet}'
        ) 
        to '/dev/stdout' 
        with (format csv, header false, delimiter 't');
    SQL
    seqkit tab2fx -o ${sample_id}.fasta.gz
    """
}

