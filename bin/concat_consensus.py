#!/usr/bin/env python3

import glob

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def main() -> None:
    consensus_files = glob.glob("*.consensus.fasta")

    assert (
        len(consensus_files) > 0
    ), "Please double check that the working directory provided contains some files with the extension '.consensus.fasta'."

    consensus_records = []
    for consensus_file in consensus_files:
        sample_name = consensus_file.split("/")[-1].replace(".consensus.fasta", "")
        concatenated_sequence = ""
        with open(consensus_file, "r", encoding="utf8") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                concatenated_sequence += str(record.seq)

        # Create a SeqRecord with the concatenated sequence
        concatenated_record = SeqRecord(
            Seq(concatenated_sequence), id=sample_name, description=""
        )
        consensus_records.append(concatenated_record)

    # Write the concatenated sequences to the output FASTA file
    with open("all_sample_consensus.fasta", "w", encoding="utf8") as output_handle:
        SeqIO.write(consensus_records, output_handle, "fasta")


if __name__ == "__main__":
    main()
