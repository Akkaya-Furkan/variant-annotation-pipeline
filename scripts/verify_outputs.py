#!/usr/bin/env python3
"""Verify PNG decoded pixels and all other output bytes against committed files."""
import hashlib
from pathlib import Path
import tempfile

from PIL import Image

from variant_annotation.pipeline import run

ROOT = Path(__file__).resolve().parents[1]
OUTPUTS = ("annotated.tsv", "filtered.tsv", "summary.json", "report.md", "summary.png")


def output_signature(path):
    if path.name == "summary.png":
        # Preserve mode and dimensions: no conversion, resampling or tolerance.
        with Image.open(path) as image:
            return image.mode, image.size, hashlib.sha256(image.tobytes()).hexdigest()
    return (hashlib.sha256(path.read_bytes()).hexdigest(),)


def main():
    with tempfile.TemporaryDirectory() as directory:
        regenerated = Path(directory)
        run(ROOT / "data/sites.vcf", ROOT / "data/clinvar_source.vcf", regenerated)
        for name in OUTPUTS:
            expected = output_signature(ROOT / "results" / name)
            actual = output_signature(regenerated / name)
            if expected != actual:
                raise SystemExit(f"Reproducibility mismatch: {name}")
            detail = ""
            if name == "summary.png":
                mode, (width, height), _ = actual
                detail = f" (decoded pixels, {mode}, {width}x{height})"
            print(f"MATCH {name} {actual[-1]}{detail}")
    print(f"Verified {len(OUTPUTS)}/{len(OUTPUTS)} reproducible outputs")


if __name__ == "__main__":
    main()
