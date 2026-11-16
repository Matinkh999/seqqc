# seqqc

A small, extensible Python tool for quality control of DNA sequencing data (FASTQ / FASTA).

> Status: week 1 of 8 - project skeleton, FASTA reader, and a first length/GC% report.

## Planned features

- Streaming FASTQ/FASTA reader (gzip supported); memory use does not grow with file size
- Metrics: GC content, mean Phred quality (per position and per read), nucleotide distribution, Q20/Q30
- Self-contained HTML report (plots embedded) and a JSON file with all values
- Plugin mechanism: add a new metric as one small Python file, without touching the core
- Validation against FastQC on public datasets

## Development setup

```bash
python -m venv .venv
source .venv/bin/activate          # Windows (CMD): .venv\Scripts\activate
pip install -e ".[dev]"
pytest
seqqc --version
```

## Length and GC% of a FASTA file

```bash
python scripts/fasta_gc.py tests/data/tiny.fasta
```

```
sequence        length      GC%
seq1                20    35.00
seq2                 8     0.00
TOTAL               28    25.00
```

GC% = (G + C) / (A + C + G + T) x 100. N and other ambiguity codes count in the
length but not in the GC% denominator. Gzip-compressed files (.gz) work too.

## Layout

```
src/seqqc/          the package (io, fasta_stats, metrics, report, plugins, cli)
tests/              unit tests; tests/data holds tiny sample files
docs/formulas.md    definition of every metric (Persian)
scripts/            helper scripts (fasta_gc.py; later: comparison with FastQC)
```

## Roadmap

1. Design, formulas, repository (done)
2. FASTQ reader and basic statistics (the FASTA reader is already done)
3. Phred metrics
4. GC and nucleotide metrics, FASTA mode, command-line interface
5. Plots and HTML/JSON report
6. Plugin mechanism and two example plugins
7. Comparison with FastQC on three public datasets
8. Documentation, final report, buffer
