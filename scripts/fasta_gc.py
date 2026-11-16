"""Print the length and GC% of every sequence in a FASTA file.

Works without installing the package:

    python scripts/fasta_gc.py tests/data/tiny.fasta
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from seqqc.fasta_stats import main  # noqa: E402

if __name__ == "__main__":
    raise SystemExit(main())
