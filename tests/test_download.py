"""Exercise downloader without network, including preservation on failure."""
import gzip
import io
import json

import pytest
from scripts import download_data

HEADER = "##fileformat=VCFv4.2\n##reference=GRCh38\n##fileDate=2026-10-04\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"


class Response(io.BytesIO):
    headers = {"Last-Modified": "test"}


def test_streamed_prefix_and_stripping(tmp_path, monkeypatch):
    body = HEADER + "1\t1\t42\tA\tG\t.\t.\tGENEINFO=GENE:1\n" + "1\t2\t43\tA\tC\t.\t.\t.\n"
    monkeypatch.setattr(download_data, "urlopen", lambda *a, **kw: Response(gzip.compress(body.encode())))
    manifest = download_data.download(tmp_path, count=1)
    assert manifest["source_records"] == 1
    assert "GENEINFO=GENE:1" in (tmp_path / "clinvar_source.vcf").read_text()
    assert (tmp_path / "sites.vcf").read_text().endswith("1\t1\t.\tA\tG\t.\t.\t.\n")
    assert "1\t2\t" not in (tmp_path / "clinvar_source.vcf").read_text()
    assert json.loads((tmp_path / "manifest.json").read_text()) == manifest


@pytest.mark.parametrize("body", [b"not gzip", gzip.compress(HEADER.encode()),
                                  gzip.compress((HEADER + "malformed\n").encode())])
def test_download_failure_preserves_existing_files(tmp_path, monkeypatch, body):
    before = {name: f"previous {name}".encode() for name in
              ["clinvar_source.vcf", "sites.vcf", "manifest.json"]}
    for name, value in before.items():
        (tmp_path / name).write_bytes(value)
    monkeypatch.setattr(download_data, "urlopen", lambda *a, **kw: Response(body))
    with pytest.raises((ValueError, OSError)):
        download_data.download(tmp_path, count=1)
    assert {p.name: p.read_bytes() for p in tmp_path.iterdir()} == before


def test_invalid_count_does_not_create_directory(tmp_path):
    output = tmp_path / "new"
    with pytest.raises(ValueError):
        download_data.download(output, count=0)
    assert not output.exists()
