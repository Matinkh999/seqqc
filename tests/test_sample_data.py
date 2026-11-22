"""Week-1 exercise: check the sample files and the formulas by hand.

These tests do not use the seqqc package yet. They only verify the format
rules and the hand-calculated values explained in docs/formulas.md.
"""

from pathlib import Path

import pytest

DATA = Path(__file__).parent / "data"


def fastq_records(path):
    lines = path.read_text().splitlines()
    assert len(lines) % 4 == 0, "FASTQ must have exactly 4 lines per record"
    return [tuple(lines[i:i + 4]) for i in range(0, len(lines), 4)]


def phred(qual, offset=33):
    return [ord(c) - offset for c in qual]


def test_fastq_structure():
    records = fastq_records(DATA / "tiny.fastq")
    assert len(records) == 5
    for header, seq, plus, qual in records:
        assert header.startswith("@")
        assert plus.startswith("+")
        assert len(seq) == len(qual)


def test_quality_line_may_start_with_at_sign():
    # read5: the quality line starts with '@', so records must be split by
    # "4 lines at a time", not by looking at the first character.
    _, _, _, qual = fastq_records(DATA / "tiny.fastq")[4]
    assert qual.startswith("@")


def test_phred_characters():
    # Q = ASCII code - 33
    assert phred("!+5?I") == [0, 10, 20, 30, 40]


def test_hand_calculated_read1():
    _, seq, _, qual = fastq_records(DATA / "tiny.fastq")[0]
    gc = (seq.count("G") + seq.count("C")) / len(seq) * 100
    q = phred(qual)
    assert gc == 25.0
    assert sum(q) / len(q) == pytest.approx(395 / 12)  # about 32.92
    assert sum(1 for x in q if x >= 30) == 10          # 10 of 12 bases are Q30+


def test_fasta_structure_and_gc():
    records, name = {}, None
    for line in (DATA / "tiny.fasta").read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            records[name] = ""
        else:
            records[name] += line
    assert set(records) == {"seq1", "seq2"}
    gc = {k: (v.count("G") + v.count("C")) / len(v) * 100 for k, v in records.items()}
    assert gc["seq1"] == 35.0   # 7 of 20 bases
    assert gc["seq2"] == 0.0
