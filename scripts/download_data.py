#!/usr/bin/env python3
"""Stream a fixed ClinVar prefix; no full database download is needed."""
import argparse
from datetime import datetime, timezone
import gzip
import hashlib
import io
import json
from pathlib import Path
from urllib.request import Request, urlopen

SOURCE = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/2026/clinvar_20261004.vcf.gz"


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def download(destination, count=300):
    if count < 1:
        raise ValueError("count must be positive")
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    # Build in memory first: download/parse failure cannot overwrite an existing snapshot.
    headers, records = [], []
    request = Request(SOURCE, headers={"User-Agent": "variant-annotation-pipeline/0.1"})
    with urlopen(request, timeout=60) as response:
        last_modified = response.headers.get("Last-Modified")
        with gzip.GzipFile(fileobj=response) as compressed:
            with io.TextIOWrapper(compressed, encoding="utf-8") as lines:
                for line in lines:
                    if line.startswith("#"):
                        headers.append(line)
                    elif line.strip():
                        if len(line.rstrip().split("\t")) < 8:
                            raise ValueError("Malformed upstream VCF")
                        records.append(line)
                        if len(records) == count:
                            break
    if len(records) != count or not any(h.startswith("#CHROM\t") for h in headers):
        raise ValueError("Upstream file lacks requested records or VCF header")
    source_text = "".join(headers + records)
    sites = ["##fileformat=VCFv4.2\n", "##reference=GRCh38\n",
             "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"]
    for line in records:
        fields = line.rstrip("\n").split("\t")[:8]
        fields[2], fields[7] = ".", "."
        sites.append("\t".join(fields) + "\n")
    payloads = {"clinvar_source.vcf": source_text, "sites.vcf": "".join(sites)}
    manifest = {
        "source_url": SOURCE, "release_date": "2026-10-04", "assembly": "GRCh38",
        "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
        "http_last_modified": last_modified, "source_records": count,
        "selection": f"First {count} records in source order; no clinical selection",
        "source_file_date": next((h.strip().split("=", 1)[1] for h in headers
                                  if h.startswith("##fileDate=")), None),
        "sha256": {name: hashlib.sha256(body.encode()).hexdigest()
                   for name, body in payloads.items()},
        "checksum_scope": "Uncompressed extracted files; full remote gzip not downloaded or verified",
    }
    for name, body in payloads.items():
        (destination / name).write_text(body, encoding="utf-8")
    (destination / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("data"))
    parser.add_argument("--count", type=int, default=300)
    args = parser.parse_args()
    print(json.dumps(download(args.output, args.count), indent=2))
