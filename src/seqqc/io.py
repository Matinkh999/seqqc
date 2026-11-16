"""Readers for FASTQ and FASTA files.

Status: the FASTA reader is implemented; the FASTQ reader is planned (week 2).

Design notes
------------
* Readers are generators: they yield one record at a time, so memory use does
  not depend on the number of records in the file.
* gzip-compressed files (.gz) are opened transparently.
* A FASTQ record is always 4 lines. Do NOT detect records by the first
  character: a quality line may start with '@' or '+'.
* Malformed input raises FormatError.
"""

import gzip
from dataclasses import dataclass


class FormatError(ValueError):
    """Raised when a FASTQ/FASTA file is malformed."""


@dataclass(frozen=True)
class FastqRecord:
    name: str  # line 1, without the leading '@'
    seq: str   # line 2
    qual: str  # line 4, same length as seq


@dataclass(frozen=True)
class FastaRecord:
    name: str  # the whole header line, without the leading '>'
    seq: str   # all sequence lines joined together


def open_text(path):
    """Open a text file for reading. Files ending in .gz are decompressed on the fly."""
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "rt", encoding="utf-8")


def read_fastq(path):
    """Yield FastqRecord objects from a FASTQ file (.fastq / .fq / .gz)."""
    raise NotImplementedError("planned for week 2")


def read_fasta(path):
    """Yield one FastaRecord per sequence in a FASTA file (.fasta / .fa / .fna / .gz).

    * A sequence may span many lines; the lines are joined.
    * Blank lines are ignored, and Windows line endings (CRLF) are handled.
    * Sequence text before the first '>' header raises FormatError.
    """
    name = None
    chunks = []
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield FastaRecord(name, "".join(chunks))
                name = line[1:].strip()
                chunks = []
            elif name is None:
                raise FormatError(
                    f"line {line_number}: sequence found before the first '>' header"
                )
            else:
                chunks.append(line)
    if name is not None:
        yield FastaRecord(name, "".join(chunks))
