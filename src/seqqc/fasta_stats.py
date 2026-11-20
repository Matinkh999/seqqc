"""Length and GC content of every sequence in a FASTA file.

Definitions
-----------
* length = number of characters in the sequence (N and other letters included)
* GC%    = (G + C) / (A + C + G + T) x 100

N and other ambiguity codes are counted in the length, but not in the GC%
denominator. Lower-case letters are treated like upper-case ones.

The file is processed line by line and the sequence itself is never stored:
only a few running counters per record are kept, so memory use does not depend
on the size of the file (or of a single sequence).

Run it:
    python -m seqqc.fasta_stats file.fasta     (package installed)
    python scripts/fasta_gc.py file.fasta      (no installation needed)
"""

import argparse
import sys
from dataclasses import dataclass

from .io import FormatError, iter_fasta_lines

NAME_WIDTH = 16


@dataclass(frozen=True)
class SeqStats:
    name: str
    length: int  # all characters of the sequence
    gc: int      # number of G and C
    acgt: int    # number of A, C, G and T (the denominator of GC%)

    @property
    def gc_percent(self) -> float:
        """GC percentage; 0.0 when the sequence has no A/C/G/T at all."""
        return 100 * self.gc / self.acgt if self.acgt else 0.0


def sequence_stats(name, seq) -> SeqStats:
    """Length and base counts of one sequence that is already in memory."""
    seq = seq.upper()
    g, c = seq.count("G"), seq.count("C")
    a, t = seq.count("A"), seq.count("T")
    return SeqStats(name, len(seq), g + c, g + c + a + t)


def iter_fasta_stats(path):
    """Yield the SeqStats of each sequence, one at a time, reading line by line."""
    name = None
    length = g = c = a = t = 0
    for is_header, text in iter_fasta_lines(path):
        if is_header:
            if name is not None:
                yield SeqStats(name, length, g + c, g + c + a + t)
            name = text
            length = g = c = a = t = 0
        else:
            line = text.upper()
            length += len(line)
            g += line.count("G")
            c += line.count("C")
            a += line.count("A")
            t += line.count("T")
    if name is not None:
        yield SeqStats(name, length, g + c, g + c + a + t)


def fasta_stats(path) -> list[SeqStats]:
    """All SeqStats of a file as a list (convenient for small files and tests)."""
    return list(iter_fasta_stats(path))


def total_stats(stats) -> SeqStats:
    """Sum of several sequences, returned as one SeqStats named 'TOTAL'."""
    stats = list(stats)
    return SeqStats(
        "TOTAL",
        sum(s.length for s in stats),
        sum(s.gc for s in stats),
        sum(s.acgt for s in stats),
    )


def table_header() -> str:
    return f"{'sequence':<{NAME_WIDTH}}  {'length':>12}  {'GC%':>7}"


def table_row(stats) -> str:
    label = (stats.name.split() or ["(unnamed)"])[0]
    return f"{label:<{NAME_WIDTH}}  {stats.length:>12,}  {stats.gc_percent:>7.2f}"


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog="fasta_gc",
        description="Print the length and GC% of every sequence in a FASTA file.",
    )
    parser.add_argument("fasta", help="FASTA file (.fasta, .fa, .fna; .gz is fine)")
    args = parser.parse_args(argv)

    count = length = gc = acgt = 0
    try:
        for stats in iter_fasta_stats(args.fasta):
            if count == 0:
                print(table_header())
            print(table_row(stats))
            count += 1
            length += stats.length
            gc += stats.gc
            acgt += stats.acgt
    except (OSError, FormatError, UnicodeDecodeError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    if count == 0:
        print(f"error: no sequences found in {args.fasta}", file=sys.stderr)
        return 1

    print(table_row(SeqStats("TOTAL", length, gc, acgt)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
