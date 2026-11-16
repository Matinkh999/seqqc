"""Length and GC content of every sequence in a FASTA file.

Definitions
-----------
* length = number of characters in the sequence (N and other letters included)
* GC%    = (G + C) / (A + C + G + T) x 100

N and other ambiguity codes are counted in the length, but not in the GC%
denominator. Lower-case letters are treated like upper-case ones.

Run it:
    python -m seqqc.fasta_stats file.fasta     (package installed)
    python scripts/fasta_gc.py file.fasta      (no installation needed)
"""

import argparse
import sys
from dataclasses import dataclass

from .io import FormatError, read_fasta


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
    """Length and base counts of one sequence."""
    seq = seq.upper()
    g, c = seq.count("G"), seq.count("C")
    a, t = seq.count("A"), seq.count("T")
    return SeqStats(name, len(seq), g + c, g + c + a + t)


def fasta_stats(path) -> list[SeqStats]:
    """Statistics of every sequence in a FASTA file."""
    return [sequence_stats(record.name, record.seq) for record in read_fasta(path)]


def total_stats(stats) -> SeqStats:
    """Sum of several sequences, returned as one SeqStats named 'TOTAL'."""
    stats = list(stats)
    return SeqStats(
        "TOTAL",
        sum(s.length for s in stats),
        sum(s.gc for s in stats),
        sum(s.acgt for s in stats),
    )


def format_table(stats) -> str:
    """Text table: one row per sequence plus a TOTAL row."""
    rows = list(stats) + [total_stats(stats)]
    labels = [(s.name.split() or ["(unnamed)"])[0] for s in rows]
    width = max(len(label) for label in ["sequence", *labels])
    lines = [f"{'sequence':<{width}}  {'length':>12}  {'GC%':>7}"]
    for label, s in zip(labels, rows):
        lines.append(f"{label:<{width}}  {s.length:>12,}  {s.gc_percent:>7.2f}")
    return "\n".join(lines)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog="fasta_gc",
        description="Print the length and GC% of every sequence in a FASTA file.",
    )
    parser.add_argument("fasta", help="FASTA file (.fasta, .fa, .fna; .gz is fine)")
    args = parser.parse_args(argv)

    try:
        stats = fasta_stats(args.fasta)
    except (OSError, FormatError, UnicodeDecodeError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    if not stats:
        print(f"error: no sequences found in {args.fasta}", file=sys.stderr)
        return 1

    print(format_table(stats))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
