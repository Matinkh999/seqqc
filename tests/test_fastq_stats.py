import inspect
from pathlib import Path

import pytest

from seqqc.fastq_stats import (
    QualityAccumulator,
    count_at_least,
    iter_read_means,
    main,
    mean_phred,
    phred_scores,
    summarize_fastq,
)

DATA = Path(__file__).parent / "data"
TINY = DATA / "tiny.fastq"
SAMPLE = DATA / "sample_dna_reads.fastq"

# Hand-calculated values for tests/data/tiny.fastq (5 reads of 12 bases, see docs/formulas.md):
#   Phred sums per read: 395, 480, 400, 304, 402  -> 1981 in total over 60 bases
#   bases with Q >= 20: 10+12+10+8+10 = 50;  bases with Q >= 30: 10+12+10+6+10 = 48


# ---- per-base and per-read functions ------------------------------------------

def test_phred_scores():
    assert phred_scores("!+5?I") == [0, 10, 20, 30, 40]


def test_phred_scores_with_offset_64():
    assert phred_scores("@Th", offset=64) == [0, 20, 40]


def test_mean_phred_of_one_read_by_hand():
    assert mean_phred("IIIIIHHHGG#!") == pytest.approx(395 / 12)


def test_mean_phred_with_offset_64():
    assert mean_phred("hhhh", offset=64) == 40.0


def test_mean_phred_of_an_empty_read_is_zero():
    assert mean_phred("") == 0.0


def test_count_at_least_q20_and_q30():
    assert count_at_least("IIIIIHHHGG#!", 30) == 10
    assert count_at_least("IIIIII55++##", 20) == 8
    assert count_at_least("IIIIII55++##", 30) == 6


def test_a_base_with_exactly_q20_counts_as_q20():
    assert count_at_least("5555555555", 20) == 10
    assert count_at_least("5555555555", 30) == 0


# ---- the accumulator -----------------------------------------------------------

def test_add_returns_the_mean_of_that_read():
    acc = QualityAccumulator()
    assert acc.add("IIII") == 40.0
    assert acc.add("5555") == 20.0


def test_summary_of_the_small_sample_file():
    acc = summarize_fastq(TINY)
    assert acc.reads == 5
    assert acc.bases == 60
    assert (acc.min_length, acc.max_length) == (12, 12)
    assert acc.mean_phred == pytest.approx(1981 / 60)
    assert acc.q20_percent == pytest.approx(100 * 50 / 60)
    assert acc.q30_percent == pytest.approx(100 * 48 / 60)


def test_summary_of_the_course_sample_file():
    acc = summarize_fastq(SAMPLE)
    assert (acc.reads, acc.bases) == (10, 500)
    assert acc.mean_phred == pytest.approx(35.7)
    assert acc.mean_read_phred == pytest.approx(35.7)
    assert (acc.min_read_mean, acc.max_read_mean) == (20.0, 40.0)
    assert acc.q20_percent == 100.0      # every base has Q >= 20 (some have exactly 20)
    assert acc.q30_percent == 80.0       # the two reads made of '5' (Q20) are below Q30


def test_empty_reads_are_counted_but_do_not_disturb_the_means():
    acc = QualityAccumulator()
    acc.add("IIII")
    acc.add("")
    assert acc.reads == 2
    assert acc.bases == 4
    assert acc.min_length == 0
    assert acc.mean_read_phred == 40.0
    assert acc.q30_percent == 100.0


def test_accumulator_without_reads_reports_zeros():
    acc = QualityAccumulator()
    assert (acc.mean_phred, acc.q20_percent, acc.q30_percent) == (0.0, 0.0, 0.0)


# ---- mean Phred of every read ----------------------------------------------------

def test_mean_phred_of_every_read_of_the_course_sample():
    means = [mean for _, _, mean in iter_read_means(SAMPLE)]
    assert means == [40, 39, 20, 40, 40, 39, 40, 20, 40, 39]


def test_mean_phred_of_every_read_of_the_small_sample():
    means = [mean for _, _, mean in iter_read_means(TINY)]
    expected = [395 / 12, 40.0, 400 / 12, 304 / 12, 402 / 12]
    assert means == pytest.approx(expected)


def test_per_read_means_are_produced_lazily():
    assert inspect.isgenerator(iter_read_means(TINY))


# ---- command line ---------------------------------------------------------------

def test_cli_prints_the_summary(capsys):
    assert main([str(SAMPLE)]) == 0
    out = capsys.readouterr().out
    assert "reads               : 10" in out
    assert "mean Phred (bases)  : 35.70" in out
    assert "Q20 (Q >= 20)       : 100.00 %" in out
    assert "Q30 (Q >= 30)       : 80.00 %" in out
    assert "READ_003" in out


def test_cli_writes_the_mean_of_every_read_to_csv(tmp_path):
    out_file = tmp_path / "means.csv"
    assert main([str(SAMPLE), "--per-read", str(out_file), "--show", "0"]) == 0
    lines = out_file.read_text().splitlines()
    assert lines[0] == "read_id,length,mean_phred"
    assert len(lines) == 11
    assert lines[1] == "READ_001,50,40.0000"
    assert lines[3] == "READ_003,50,20.0000"


def test_cli_reports_a_malformed_file(tmp_path, capsys):
    path = tmp_path / "bad.fastq"
    path.write_bytes(b"@r1\nACGT\n+\nIII\n")
    assert main([str(path)]) == 1
    err = capsys.readouterr().err
    assert "record 1" in err and "is not equal to" in err


def test_cli_reports_a_missing_file(capsys):
    assert main(["no_such_file.fastq"]) == 1
    assert "no_such_file.fastq" in capsys.readouterr().err


def test_cli_rejects_a_file_without_reads(tmp_path, capsys):
    path = tmp_path / "empty.fastq"
    path.write_bytes(b"")
    assert main([str(path)]) == 1
    assert "no reads" in capsys.readouterr().err
