#!/usr/bin/env python3
# Licensed under the Apache 2.0 license, see LICENSE file for details
"""A tool to flip GenBank format DNA records to make sure genes of interest are on the right strand"""

from argparse import ArgumentParser, FileType
from pathlib import Path
import sys
from typing import Set, TextIO

from Bio import SeqIO


def main():
    parser = ArgumentParser()
    parser.add_argument("id_file", type=FileType("r", encoding="utf-8"),
                        help="File containing the record IDs of the records to flip, one per line")
    parser.add_argument("gbk_list_file", type=FileType("r", encoding="utf-8"),
                        help="File containing the paths and filenames of GenBank files to process")
    parser.add_argument("-s", "--suffix", default="_flipped",
                        help="Suffix to attach to the filename for the new file with flipped records.")
    args = parser.parse_args()

    run(args.id_file, args.gbk_list_file, args.suffix)


def run(id_file: TextIO, gbk_list_file: TextIO, suffix: str):
    """Run the record flipper"""
    records_to_flip: Set[str] = set()
    for line in id_file:
        records_to_flip.add(line.strip())

    for line in gbk_list_file:
        file = Path(line.strip())
        handle_file(file, records_to_flip, suffix)


def handle_file(file: Path, records_to_flip: Set[str], suffix: str):
    """ Handle a single file. """
    records = list(SeqIO.parse(file, "genbank"))

    changed = False
    for i, record in enumerate(records):
        if record.name not in records_to_flip:
            print(record.name, "not in", records_to_flip, file=sys.stderr)
            continue
        changed = True
        flipped = record.reverse_complement(id=True, name=True, description=True,
                                            annotations=True, dbxrefs=True)
        break

    if changed:
        records[i] = flipped

        new_name = file.parent / f"{file.stem}{suffix}.gbk"
        SeqIO.write(records, new_name, "genbank")


if  __name__ == "__main__":
    main()
