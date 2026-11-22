"""Readers for FASTQ and FASTA files.

Both readers work line by line and are generators: they hold only the current
record in memory, so memory use does not depend on the size of the file.
gzip-compressed files (.gz) are opened transparently.

FASTQ validation (see check_fastq_record): every record must have 4 lines,
a header that starts with '@', a third line that starts with '+', a sequence
made only of A/C/G/T/N, and a quality string with exactly one valid character
per base.
"""

import gzip
import re
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


# ----------------------------------------------------------------------------
# FASTQ
# ----------------------------------------------------------------------------

_VALID_SEQUENCE = re.compile(r"[ACGTNacgtn]*")
_HIGHEST_QUALITY_CHAR = "~"  # ASCII 126, the last printable character


def check_fastq_record(record_number, line_number, header, seq, plus, qual, offset=33):
    """Raise FormatError if the four lines of one FASTQ record are not valid.

    record_number and line_number (the line of the header) are only used in
    the error message. Line endings must already be removed.
    """
    where = f"record {record_number} (line {line_number})"
    if not header.startswith("@"):
        raise FormatError(f"{where}: the header line must start with '@'")
    if not plus.startswith("+"):
        raise FormatError(f"{where}: the third line must start with '+'")
    if len(seq) != len(qual):
        raise FormatError(
            f"{where}: sequence length ({len(seq)}) is not equal to "
            f"quality length ({len(qual)})"
        )
    if not _VALID_SEQUENCE.fullmatch(seq):
        raise FormatError(
            f"{where}: the sequence contains characters other than A, C, G, T, N"
        )
    if qual and (min(qual) < chr(offset) or max(qual) > _HIGHEST_QUALITY_CHAR):
        raise FormatError(
            f"{where}: the quality string contains a character outside the "
            f"valid range for offset {offset}"
        )


def read_fastq(path, offset=33):
    """Yield one validated FastqRecord at a time from a FASTQ file (.fastq / .fq / .gz).

    * The file is read line by line; one record is 4 lines.
    * Blank lines between records and at the end of the file are ignored.
    * A bad record raises FormatError that names the record and its line.
      Records before the bad one have already been returned by then.
    """
    if offset not in (33, 64):
        raise ValueError("offset must be 33 or 64")
    with open_text(path) as handle:
        lines = enumerate(handle, start=1)
        record_number = 0
        for line_number, header in lines:
            if not header.strip():
                continue
            record_number += 1
            try:
                seq = next(lines)[1]
                plus = next(lines)[1]
                qual = next(lines)[1]
            except StopIteration:
                raise FormatError(
                    f"record {record_number} (line {line_number}): "
                    "the file ends in the middle of a record"
                ) from None
            header, seq, plus, qual = (s.rstrip("\r\n") for s in (header, seq, plus, qual))
            check_fastq_record(record_number, line_number, header, seq, plus, qual, offset)
            yield FastqRecord(header[1:], seq, qual)


# ----------------------------------------------------------------------------
# FASTA
# ----------------------------------------------------------------------------

def iter_fasta_lines(path):
    """Yield (True, header_text) for every '>' line and (False, line) for every sequence line.

    This is the only place that knows how a FASTA file is laid out. Blank lines
    are skipped, and sequence text before the first '>' header raises FormatError.
    """
    seen_header = False
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                seen_header = True
                yield True, line[1:].strip()
            elif not seen_header:
                raise FormatError(
                    f"line {line_number}: sequence found before the first '>' header"
                )
            else:
                yield False, line


def read_fasta(path):
    """Yield one FastaRecord per sequence in a FASTA file (.fasta / .fa / .fna / .gz).

    A sequence may span many lines; the lines are joined, so one whole sequence
    is held in memory at a time (use fasta_stats.iter_fasta_stats when only
    counts are needed).
    """
    name = None
    chunks = []
    for is_header, text in iter_fasta_lines(path):
        if is_header:
            if name is not None:
                yield FastaRecord(name, "".join(chunks))
            name, chunks = text, []
        else:
            chunks.append(text)
    if name is not None:
        yield FastaRecord(name, "".join(chunks))
