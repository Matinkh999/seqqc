"""Per-read mean Phred score and Q20/Q30 percentages of a FASTQ file.

Works without installing the package:

    python scripts/fastq_qc.py tests/data/sample_dna_reads.fastq
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from seqqc.fastq_stats import main  # noqa: E402

if __name__ == "__main__":
    raise SystemExit(main())
