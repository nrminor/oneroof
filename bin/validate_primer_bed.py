#!/usr/bin/env python3

import sys
from pathlib import Path


def fix_bed_files(input_bed: Path, output_prefix: str = "validated") -> None:
    with open(input_bed, encoding="utf8") as file, open(
        f"{output_prefix}.bed",
        "w",
        encoding="utf8",
    ) as output:
        lines = [row.strip().split("\t") for row in file if row != ""]

        fixed_file = []
        for i, row in enumerate(lines):
            required_col_count = 6
            assert (
                len(row) >= required_col_count
            ), f"BED file row {i} has fewer than the required 6 columns: \n\n {row}"
            start = int(row[1])
            stop = int(row[2])
            sign = row[5]
            if start < stop:
                fixed_file.append("\t".join(row))
                continue
            row[2] = str(start)
            row[1] = str(stop)
            if sign == "+":
                row[5] = "-"
            elif sign == "-":
                row[5] = "+"
            else:
                message = f"unsupported value encountered in the strand column of BED file row {i}: {row}"
                raise ValueError(message)
            fixed_file.append("\t".join(row))
        output.write("\n".join(fixed_file))


def main() -> None:
    arguments = sys.argv
    assert len(arguments) >= 2, f"Inputs were less than the required 2: {arguments}"  # noqa: PLR2004
    input_bed = Path(arguments[1])
    assert input_bed.is_file(), f"Inputed path is not file: {input_bed}"
    output_prefix = arguments[2] if len(arguments) >= 3 else "validated"
    fix_bed_files(input_bed, output_prefix)


if __name__ == "__main__":
    main()
