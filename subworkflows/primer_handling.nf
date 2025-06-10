include { VALIDATE_PRIMER_BED      } from "../modules/validate"
include { REPORT_REFERENCES        } from "../modules/reporting"
include { RESPLICE_PRIMERS         } from "../modules/resplice_primers"
include { SPLIT_PRIMER_COMBOS      } from "../modules/split_primer_combos"
include { GET_PRIMER_PATTERNS      } from "../modules/primer_patterns"
include { GET_PRIMER_SEQS          } from "../modules/bedtools"
include { FAIDX                    } from "../modules/samtools"
include { ORIENT_READS             } from "../modules/vsearch"
include { TRIM_ENDS_TO_PRIMERS     } from "../modules/cutadapt"
include {
    FIND_COMPLETE_AMPLICONS ;
    // TRIM_ENDS_TO_PRIMERS ;
    PER_AMPLICON_FILTERS ;
    MERGE_BY_SAMPLE ;
    AMPLICON_STATS
} from "../modules/seqkit"
include { FILTER_WITH_CHOPPER      } from "../modules/chopper"
include { RASUSA_READ_DOWNSAMPLING } from "../modules/rasusa"
include { WRITE_PRIMER_FASTA } from "../modules/write_primer_fasta"
include { CREATE_PRIMER_TSV } from "../modules/output_primer_tsv"
include { COLLECT_PRIMER_TSV } from "../modules/output_primer_tsv"


workflow PRIMER_HANDLING {
    take:
    ch_basecalls
    ch_primer_bed
    ch_refseq
    ch_primer_tsv

    main:

    ORIENT_READS(
        ch_basecalls
            .map { barcode, reads -> tuple(barcode, file(reads), file(reads).countFasta()) }
            .filter { it[2] > 100 }
            .map { barcode, reads, _read_count -> tuple(barcode, file(reads)) },
        ch_refseq,
    )
  
    if (params.primer_bed && params.primer_bed  != "") {
        
        VALIDATE_PRIMER_BED(
            ch_primer_bed
        )

        REPORT_REFERENCES(
            VALIDATE_PRIMER_BED.out
        )

        RESPLICE_PRIMERS(
            VALIDATE_PRIMER_BED.out
        )

        SPLIT_PRIMER_COMBOS(
            RESPLICE_PRIMERS.out
        )

        GET_PRIMER_SEQS(
            SPLIT_PRIMER_COMBOS.out.flatten(),
            ch_refseq,
        )

        CREATE_PRIMER_TSV(
            GET_PRIMER_SEQS.out
        )

        COLLECT_PRIMER_TSV(
            CREATE_PRIMER_TSV.out.primer_pair.collect()
        )

        GET_PRIMER_PATTERNS(
            GET_PRIMER_SEQS.out
        )

    } else if (params.primer_tsv && params.primer_tsv != "") {
        ch_tsv_to_fasta = ch_primer_tsv
            .splitCsv(header:true, sep:'\t', strip: true)
            .filter { row -> !row.startsWith("#")}
            .filter { row ->
                row.amplicon_name && row.fwd_sequence && row.reverse_sequence
            }
            .map { row ->
                def name = row.amplicon_name
                def fwd = row.fwd_sequence
                def rev = row.reverse_sequence

                def header1 = ">${name}:0-${fwd.size() - 1}"
                def header2 = ">${name}:${fwd.size()}-${fwd.size() + rev.size() - 1}"

                def fasta = "${header1}\n${fwd}\n${header2}\n${rev}"

                tuple(name, fasta)
          }
    
        WRITE_PRIMER_FASTA(
            ch_tsv_to_fasta
        )

         GET_PRIMER_PATTERNS(
            WRITE_PRIMER_FASTA.out
        )
    }


    FIND_COMPLETE_AMPLICONS(
        ORIENT_READS.out.map { _barcode, fastq -> fastq }.combine(GET_PRIMER_PATTERNS.out)
    )


    TRIM_ENDS_TO_PRIMERS(
        FIND_COMPLETE_AMPLICONS.out
    )

    PER_AMPLICON_FILTERS(
        TRIM_ENDS_TO_PRIMERS.out
            .map { id, fastq -> tuple(id, fastq, file(fastq).countFasta()) }
            .filter { it[2] > 0 }
            .map { id, fastq, _read_count -> tuple(id, file(fastq)) }
    )

    FAIDX(
        PER_AMPLICON_FILTERS.out
            .map { id, fastq -> tuple(id, fastq, file(fastq).countFasta()) }
            .filter { it[2] > 0 }
            .map { id, fastq, _read_count -> tuple(id, file(fastq)) }
    )

    RASUSA_READ_DOWNSAMPLING(
        FAIDX.out
    )

    AMPLICON_STATS(
        RASUSA_READ_DOWNSAMPLING.out.groupTuple(by: 0)
    )

    MERGE_BY_SAMPLE(
        RASUSA_READ_DOWNSAMPLING.out.groupTuple(by: 0)
    )

    emit:
    MERGE_BY_SAMPLE.out
}

