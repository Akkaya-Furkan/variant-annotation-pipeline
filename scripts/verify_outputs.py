#!/usr/bin/env python3
"""Regenerate all five outputs offline and compare SHA-256 to committed files."""
import hashlib
from pathlib import Path
import tempfile

from variant_annotation.pipeline import run

ROOT = Path(__file__).resolve().parents[1]
OUTPUTS = ("annotated.tsv", "filtered.tsv", "summary.json", "report.md", "summary.png")


def main():
    with tempfile.TemporaryDirectory() as directory:
        regenerated = Path(directory)
        run(ROOT / "data/sites.vcf", ROOT / "data/clinvar_source.vcf", regenerated)
        for name in OUTPUTS:
            expected = (ROOT / "results" / name).read_bytes()
            actual = (regenerated / name).read_bytes()
            if hashlib.sha256(expected).digest() != hashlib.sha256(actual).digest():
                raise SystemExit(f"Reproducibility mismatch: {name}")
            print(f"MATCH {name} {hashlib.sha256(actual).hexdigest()}")
    print(f"Verified {len(OUTPUTS)}/{len(OUTPUTS)} reproducible outputs")


if __name__ == "__main__":
    main()
