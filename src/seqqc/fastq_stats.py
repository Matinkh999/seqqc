"""Per-read mean Phred score and file-level Q20 / Q30 of a FASTQ file.

Definitions
-----------
* Phred score of a base: Q = (ASCII code of its quality character) - offset
  (offset 33 for modern Illumina data, 64 for old data)
* mean Phred of a read:  the arithmetic mean of the Q values of its bases
* Q20 / Q30 of a file:   the percentage of ALL bases in the file whose Q is at
  least 20 / at least 30 (a base with Q = 20 counts as Q20)

Everything is computed in a single pass over the file. Only a handful of
running counters are kept, so memory use does not depend on the number of reads.

Run it:
    python -m seqqc.fastq_stats reads.fastq     (package installed)
    python scripts/fastq_qc.py reads.fastq      (no installation needed)
"""

import argparse
import csv
import sys
from functools import lru_cache
from pathlib import Path

from .io import FormatError, read_fastq


def phred_scores(qual, offset=33):
    """Phred score of every base of a quality string (the plain, easy-to-read version)."""
    return [ord(ch) - offset for ch in qual]


def mean_phred(qual, offset=33):
    """Mean Phred score of one read (0.0 for an empty read).

    sum(bytes) adds up the ASCII codes of all characters at C speed; subtracting
    offset once per base turns that sum into the sum of Phred scores.
    """
    if not qual:
        return 0.0
    return sum(qual.encode("ascii")) / len(qual) - offset


@lru_cache(maxsize=None)
def _threshold_table(cutoff):
    """256-entry table: 1 for byte values >= cutoff, 0 for the others."""
    return bytes(1 if value >= cutoff else 0 for value in range(256))


def count_at_least(qual, q, offset=33):
    """Number of bases in a quality string whose Phred score is >= q.

    bytes.translate replaces every character by 1 (good enough) or 0, and
    count(1) counts the ones; both steps run at C speed.
    """
    return qual.encode("ascii").translate(_threshold_table(offset + q)).count(1)


class QualityAccumulator:
    """Running totals for a FASTQ file; memory use does not depend on the number of reads."""

    def __init__(self, offset=33):
        self.offset = offset
        self.reads = 0
        self.bases = 0
        self.q20_bases = 0
        self.q30_bases = 0
        self.min_length = None
        self.max_length = 0
        self.min_read_mean = None
        self.max_read_mean = None
        self._phred_sum = 0         # sum of the Phred scores of all bases
        self._read_mean_sum = 0.0   # sum of the per-read means
        self._reads_with_bases = 0
        self._table20 = _threshold_table(offset + 20)
        self._table30 = _threshold_table(offset + 30)

    def add(self, qual):
        """Add one read to the totals and return its mean Phred score."""
        n = len(qual)
        data = qual.encode("ascii")
        phred_sum = sum(data) - self.offset * n
        mean = phred_sum / n if n else 0.0

        self.reads += 1
        self.bases += n
        self._phred_sum += phred_sum
        self.q20_bases += data.translate(self._table20).count(1)
        self.q30_bases += data.translate(self._table30).count(1)
        self.min_length = n if self.min_length is None else min(self.min_length, n)
        self.max_length = max(self.max_length, n)
        if n:  # empty reads have no mean, so they are left out of the per-read statistics
            self._reads_with_bases += 1
            self._read_mean_sum += mean
            if self.min_read_mean is None or mean < self.min_read_mean:
                self.min_read_mean = mean
            if self.max_read_mean is None or mean > self.max_read_mean:
                self.max_read_mean = mean
        return mean

    @property
    def mean_length(self):
        return self.bases / self.reads if self.reads else 0.0

    @property
    def mean_phred(self):
        """Mean Phred score over all bases of the file."""
        return self._phred_sum / self.bases if self.bases else 0.0

    @property
    def mean_read_phred(self):
        """Mean of the per-read mean Phred scores."""
        if not self._reads_with_bases:
            return 0.0
        return self._read_mean_sum / self._reads_with_bases

    @property
    def q20_percent(self):
        """Percentage of all bases with Phred score >= 20."""
        return 100 * self.q20_bases / self.bases if self.bases else 0.0

    @property
    def q30_percent(self):
        """Percentage of all bases with Phred score >= 30."""
        return 100 * self.q30_bases / self.bases if self.bases else 0.0


def summarize_fastq(path, offset=33):
    """Read a FASTQ file once and return its QualityAccumulator."""
    acc = QualityAccumulator(offset)
    for record in read_fastq(path, offset):
        acc.add(record.qual)
    return acc


def iter_read_means(path, offset=33):
    """Yield (read_id, length, mean_phred) for every read, one at a time."""
    for record in read_fastq(path, offset):
        yield record.name, len(record.seq), mean_phred(record.qual, offset)


def _number(value):
    return "n/a" if value is None else f"{value:.2f}"


def format_summary(acc, path, first_reads=()):
    """Text report of a finished QualityAccumulator."""
    lines = [
        f"file                : {Path(path).name}",
        f"reads               : {acc.reads:,}",
        f"bases               : {acc.bases:,}",
        f"read length         : mean {acc.mean_length:.1f} (min {acc.min_length}, max {acc.max_length})",
        f"mean Phred (bases)  : {acc.mean_phred:.2f}",
        f"mean Phred per read : min {_number(acc.min_read_mean)} / mean {acc.mean_read_phred:.2f}"
        f" / max {_number(acc.max_read_mean)}",
        f"Q20 (Q >= 20)       : {acc.q20_percent:.2f} %",
        f"Q30 (Q >= 30)       : {acc.q30_percent:.2f} %",
    ]
    if first_reads:
        lines += ["", f"first {len(first_reads)} reads", f"{'read_id':<16}{'length':>8}{'mean_Phred':>12}"]
        lines += [f"{rid:<16}{length:>8}{mean:>12.2f}" for rid, length, mean in first_reads]
    return "\n".join(lines)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog="fastq_qc",
        description="Per-read mean Phred score and Q20/Q30 percentages of a FASTQ file.",
    )
    parser.add_argument("fastq", help="FASTQ file (.fastq, .fq; .gz is fine)")
    parser.add_argument("--offset", type=int, choices=(33, 64), default=33,
                        help="Phred offset of the quality characters (default: 33)")
    parser.add_argument("--per-read", metavar="CSV",
                        help="also write the mean Phred score of every read to this CSV file")
    parser.add_argument("--show", type=int, default=5, metavar="N",
                        help="print the mean Phred score of the first N reads (default: 5)")
    args = parser.parse_args(argv)

    acc = QualityAccumulator(args.offset)
    first_reads = []
    csv_file = writer = None
    try:
        if args.per_read:
            csv_file = open(args.per_read, "w", newline="", encoding="utf-8")
            writer = csv.writer(csv_file)
            writer.writerow(["read_id", "length", "mean_phred"])
        for record in read_fastq(args.fastq, args.offset):
            mean = acc.add(record.qual)
            read_id = (record.name.split() or ["(unnamed)"])[0]
            if len(first_reads) < args.show:
                first_reads.append((read_id, len(record.seq), mean))
            if writer:
                writer.writerow([read_id, len(record.seq), f"{mean:.4f}"])
    except (OSError, FormatError, UnicodeDecodeError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    finally:
        if csv_file:
            csv_file.close()

    if acc.reads == 0:
        print(f"error: no reads found in {args.fastq}", file=sys.stderr)
        return 1
    print(format_summary(acc, args.fastq, first_reads))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
