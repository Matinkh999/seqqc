import gzip
from pathlib import Path

import pytest

from seqqc.fasta_stats import fasta_stats, main, sequence_stats, total_stats
from seqqc.io import FormatError, read_fasta

DATA = Path(__file__).parent / "data"


def test_read_fasta_joins_multiline_sequences():
    records = list(read_fasta(DATA / "tiny.fasta"))
    assert [r.name for r in records] == ["seq1 demo", "seq2 no_gc"]
    assert records[0].seq == "GATTACAGATTAGGCCTTAA"
    assert records[1].seq == "ATATATAT"


def test_stats_of_the_sample_file():
    seq1, seq2 = fasta_stats(DATA / "tiny.fasta")
    assert (seq1.length, seq1.gc, seq1.acgt) == (20, 7, 20)
    assert seq1.gc_percent == pytest.approx(35.0)
    assert (seq2.length, seq2.gc_percent) == (8, 0.0)


def test_total_stats():
    total = total_stats(fasta_stats(DATA / "tiny.fasta"))
    assert total.length == 28
    assert total.gc_percent == pytest.approx(25.0)  # 7 of 28 bases


def test_lowercase_blank_lines_and_windows_line_endings(tmp_path):
    path = tmp_path / "messy.fasta"
    path.write_bytes(b">s1 lower\r\nacgtNNNN\r\n\r\ngg\r\n>s2\r\nCCCC\r\n")
    s1, s2 = fasta_stats(path)
    assert (s1.length, s1.gc, s1.acgt) == (10, 4, 6)
    assert s1.gc_percent == pytest.approx(66.6667, abs=1e-3)
    assert s2.gc_percent == 100.0


def test_n_counts_in_length_but_not_in_gc_percent():
    stats = sequence_stats("x", "ACGTNNNN")
    assert stats.length == 8
    assert stats.acgt == 4
    assert stats.gc_percent == 50.0


def test_gzip_file(tmp_path):
    path = tmp_path / "tiny.fasta.gz"
    with gzip.open(path, "wt") as handle:
        handle.write((DATA / "tiny.fasta").read_text())
    assert [s.length for s in fasta_stats(path)] == [20, 8]


def test_header_without_sequence(tmp_path):
    path = tmp_path / "empty.fasta"
    path.write_text(">only_header\n>next\nAC\n")
    empty, nxt = fasta_stats(path)
    assert empty.length == 0
    assert empty.gc_percent == 0.0
    assert nxt.length == 2


def test_sequence_before_first_header_is_an_error(tmp_path):
    path = tmp_path / "bad.fasta"
    path.write_text("ACGT\n>seq1\nAC\n")
    with pytest.raises(FormatError):
        list(read_fasta(path))


def test_cli_prints_a_table(capsys):
    assert main([str(DATA / "tiny.fasta")]) == 0
    out = capsys.readouterr().out
    assert "seq1" in out and "35.00" in out
    assert "TOTAL" in out and "25.00" in out


def test_cli_reports_a_missing_file(capsys):
    assert main(["no_such_file.fasta"]) == 1
    assert "no_such_file.fasta" in capsys.readouterr().err
