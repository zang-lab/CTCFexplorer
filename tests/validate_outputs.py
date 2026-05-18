from __future__ import annotations

import csv
import sys
from pathlib import Path


def require(path: Path) -> Path:
    if not path.exists():
        raise SystemExit(f"Missing expected output: {path}")
    return path


def validate_bed(path: Path) -> None:
    with path.open() as handle:
        rows = [line.strip().split("	") for line in handle if line.strip()]
    if not rows:
        raise SystemExit(f"BED file is empty: {path}")
    for idx, row in enumerate(rows, start=1):
        if len(row) < 3:
            raise SystemExit(f"BED row {idx} in {path} has fewer than 3 columns")
        try:
            start = int(row[1])
            end = int(row[2])
        except ValueError as exc:
            raise SystemExit(f"BED row {idx} in {path} has non-integer coordinates") from exc
        if start >= end:
            raise SystemExit(f"BED row {idx} in {path} has start >= end")


def validate_tsv(path: Path, min_rows: int = 1) -> None:
    with path.open() as handle:
        reader = csv.reader(handle, delimiter="	")
        rows = list(reader)
    if len(rows) <= min_rows:
        raise SystemExit(f"TSV file has too few rows: {path}")


def main(outdir: str) -> None:
    root = Path(outdir)
    require(root / "01_fastq")
    require(root / "02_fastqc")
    validate_bed(require(root / "06_union" / "union_sites.bed"))
    validate_tsv(require(root / "07_occupancy" / "union_occupancy_matrix.tsv"))
    validate_tsv(require(root / "07_occupancy" / "union_occupancy_scores.tsv"))
    validate_tsv(require(root / "08_signal" / "union_RPKM.tsv"))
    validate_bed(require(root / "09_final" / "high_confidence.bed"))
    validate_bed(require(root / "09_final" / "constitutive.bed"))
    validate_tsv(require(root / "09_final" / "union_occupancy_high_confidence.tsv"))
    validate_tsv(require(root / "09_final" / "union_RPKM_high_confidence.tsv"))


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("Usage: validate_outputs.py <nextflow-results-dir>")
    main(sys.argv[1])
