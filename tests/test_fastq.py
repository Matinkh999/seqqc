import gzip
import inspect
import re
from pathlib import Path

import pytest

from seqqc.io import FastqRecord, FormatError, check_fastq_record, read_fastq

DATA = Path(__file__).parent / "data"

GOOD = "@r1\nACGT\n+\nIIII\n"


def write(tmp_path, text, name="reads.fastq"):
    path = tmp_path / name
    path.write_bytes(text.encode("ascii"))
    return path


# ---- reading ---------------------------------------------------------------

def test_reads_every_record_of_the_small_sample_file():
    records = list(read_fastq(DATA / "tiny.fastq"))
    assert len(records) == 5
    assert records[0] == FastqRecord("read1 demo", "GATTACAGATTA", "IIIIIHHHGG#!")


def test_reads_the_course_sample_file():
    records = list(read_fastq(DATA / "sample_dna_reads.fastq"))
    assert len(records) == 10
    assert records[0].name == "READ_001"
    assert all(len(r.seq) == len(r.qual) == 50 for r in records)


def test_reader_is_a_generator():
    assert inspect.isgenerator(read_fastq(DATA / "tiny.fastq"))


def test_records_are_validated_one_at_a_time(tmp_path):
    path = write(tmp_path, GOOD + "@r2\nACGT\n+\nIII\n")  # record 2: quality too short
    records = read_fastq(path)
    assert next(records).name == "r1"      # record 1 is returned first ...
    with pytest.raises(FormatError):       # ... the error appears only when record 2 is reached
        next(records)


def test_quality_line_that_starts_with_at_sign_is_not_a_new_record():
    records = list(read_fastq(DATA / "tiny.fastq"))
    assert records[4].name == "read5 qual_starts_with_at"
    assert records[4].qual == "@@IIIIIIII++"


def test_blank_lines_between_records_and_at_the_end_are_ignored(tmp_path):
    path = write(tmp_path, GOOD + "\n\n" + GOOD.replace("r1", "r2") + "\n\n")
    assert [r.name for r in read_fastq(path)] == ["r1", "r2"]


def test_windows_line_endings(tmp_path):
    path = write(tmp_path, GOOD.replace("\n", "\r\n"))
    (record,) = read_fastq(path)
    assert (record.seq, record.qual) == ("ACGT", "IIII")


def test_gzip_file(tmp_path):
    path = tmp_path / "reads.fastq.gz"
    with gzip.open(path, "wt") as handle:
        handle.write(GOOD)
    assert [r.seq for r in read_fastq(path)] == ["ACGT"]


def test_empty_file_has_no_records(tmp_path):
    assert list(read_fastq(write(tmp_path, ""))) == []


def test_empty_read_is_allowed(tmp_path):
    (record,) = read_fastq(write(tmp_path, "@empty\n\n+\n\n"))
    assert (record.seq, record.qual) == ("", "")


# ---- validation ------------------------------------------------------------

@pytest.mark.parametrize(
    "text, expected",
    [
        ("r1\nACGT\n+\nIIII\n", "must start with '@'"),
        ("@r1\nACGT\nIIII\nIIII\n", "must start with '+'"),
        ("@r1\nACGT\n+\nIII\n", "is not equal to"),
        ("@r1\nACXT\n+\nIIII\n", "characters other than"),
        ("@r1\nACGT\n+\nII I\n", "outside the valid range"),
        ("@r1\nACGT\n+\n", "ends in the middle"),
    ],
)
def test_invalid_records_raise_format_error(tmp_path, text, expected):
    with pytest.raises(FormatError, match=re.escape(expected)):
        list(read_fastq(write(tmp_path, text)))


def test_error_message_names_the_record_and_the_line(tmp_path):
    text = GOOD + GOOD.replace("r1", "r2") + "@r3\nACGT\n+\nIII\n"
    with pytest.raises(FormatError) as error:
        list(read_fastq(write(tmp_path, text)))
    assert "record 3" in str(error.value)
    assert "line 9" in str(error.value)


def test_check_fastq_record_accepts_a_valid_record():
    check_fastq_record(1, 1, "@r1", "ACGT", "+", "IIII")  # no error


def test_check_fastq_record_rejects_different_lengths():
    with pytest.raises(FormatError, match="is not equal to"):
        check_fastq_record(1, 1, "@r1", "ACGT", "+", "III")


def test_lower_case_bases_are_valid():
    check_fastq_record(1, 1, "@r1", "acgtn", "+", "IIIII")


def test_quality_range_depends_on_the_offset():
    sample = DATA / "sample_dna_reads.fastq"          # its qualities include '5' (ASCII 53)
    assert len(list(read_fastq(sample, offset=33))) == 10
    with pytest.raises(FormatError):                  # '5' is below '@' (ASCII 64)
        list(read_fastq(sample, offset=64))


def test_offset_must_be_33_or_64():
    with pytest.raises(ValueError):
        list(read_fastq(DATA / "tiny.fastq", offset=50))
