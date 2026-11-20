"""Streaming tests: memory must not grow with the size of the file."""

import inspect
import random
import tracemalloc

from seqqc.fasta_stats import iter_fasta_stats
from seqqc.fastq_stats import summarize_fastq


def make_fastq(path, reads, length=100, seed=1):
    """Write a FASTQ file read by read; nothing is kept in memory."""
    rng = random.Random(seed)
    pool = [
        ("".join(rng.choices("ACGT", k=length)), "".join(rng.choices("IHG5+", k=length)))
        for _ in range(500)
    ]
    with open(path, "w", newline="\n") as handle:
        for i in range(reads):
            seq, qual = pool[i % 500]
            handle.write(f"@read{i}\n{seq}\n+\n{qual}\n")


def peak_memory(function, *args):
    """Largest amount of memory (bytes) that Python allocated while function ran."""
    tracemalloc.start()
    try:
        function(*args)
        _, peak = tracemalloc.get_traced_memory()
    finally:
        tracemalloc.stop()
    return peak


def consume(iterator):
    for _ in iterator:
        pass


def test_fastq_summary_uses_little_memory(tmp_path):
    path = tmp_path / "big.fastq"
    make_fastq(path, reads=40_000)                       # about 9 MB on disk
    assert path.stat().st_size > 8_000_000
    assert peak_memory(summarize_fastq, path) < 1_000_000   # less than 1 MB


def test_fastq_memory_does_not_grow_with_the_file_size(tmp_path):
    small, large = tmp_path / "small.fastq", tmp_path / "large.fastq"
    make_fastq(small, reads=5_000)
    make_fastq(large, reads=40_000)                      # 8 times bigger
    assert peak_memory(summarize_fastq, large) < peak_memory(summarize_fastq, small) + 100_000


def test_fasta_stats_do_not_keep_the_sequence_in_memory(tmp_path):
    path = tmp_path / "chromosome.fasta"
    with open(path, "w", newline="\n") as handle:
        handle.write(">chr1 one very long sequence\n")
        for _ in range(50_000):                          # 50,000 lines x 80 bases = 4 million bases
            handle.write("ACGT" * 20 + "\n")
    stats = next(iter_fasta_stats(path))
    assert stats.length == 4_000_000
    assert stats.gc_percent == 50.0
    assert peak_memory(lambda p: consume(iter_fasta_stats(p)), path) < 1_000_000


def test_fasta_stats_are_produced_lazily(tmp_path):
    path = tmp_path / "two.fasta"
    path.write_bytes(b">a\nACGT\n>b\nGGCC\n")
    assert inspect.isgenerator(iter_fasta_stats(path))
