"""Command-line interface (full implementation planned for week 4)."""

import argparse

from . import __version__


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="seqqc",
        description="Quality control for DNA sequencing data (FASTQ/FASTA).",
    )
    parser.add_argument("--version", action="version", version=f"seqqc {__version__}")
    return parser


def main(argv=None) -> int:
    build_parser().parse_args(argv)
    print("seqqc: nothing implemented yet (project skeleton, week 1).")
    return 0
