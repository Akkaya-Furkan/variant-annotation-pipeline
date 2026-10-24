"""Synthetic edge cases and independent public-snapshot acceptance checks."""
import csv
import hashlib
import json
from pathlib import Path

import pytest
from variant_annotation.pipeline import annotate, read_vcf, run

ROOT = Path(__file__).resolve().parents[1]
HEADER = "##fileformat=VCFv4.2\n##reference=GRCh38\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"


def vcf(tmp_path, rows, name="test.vcf"):
    path = tmp_path / name
    path.write_text(HEADER + "\n".join(rows) + "\n")
    return path


def record(pos=10, ref="A", alt="G", info=".", quality="."):
    return f"1\t{pos}\t42\t{ref}\t{alt}\t.\t{quality}\t{info}"


@pytest.mark.parametrize(("ref", "alt", "kind", "substitution"), [
    ("A", "G", "SNV", "transition"), ("C", "T", "SNV", "transition"),
    ("A", "C", "SNV", "transversion"), ("A", "AT", "insertion", "NA"),
    ("AT", "A", "deletion", "NA"), ("AC", "GT", "MNV", "NA"),
])
def test_variant_types(tmp_path, ref, alt, kind, substitution):
    variants = read_vcf(vcf(tmp_path, [record(ref=ref, alt=alt)]))
    row = annotate(variants, variants)[0]
    assert (row["variant_class"], row["substitution"]) == (kind, substitution)


def test_missing_and_unmatched(tmp_path):
    variants = read_vcf(vcf(tmp_path, [record(), record(pos=11)]))
    rows = annotate(variants, variants[:1])
    assert rows[0]["matched"] and rows[0]["gene"] == "NA"
    assert not rows[1]["matched"] and rows[1]["clinvar_id"] == "NA"
    assert not any(row["selected"] for row in rows)


@pytest.mark.parametrize(("sig", "status", "quality", "passes"), [
    ("Pathogenic", "reviewed_by_expert_panel", "PASS", (True, True, True)),
    ("Pathogenic/Likely_pathogenic", "criteria_provided,_multiple_submitters,_no_conflicts", ".", (True, True, True)),
    ("Likely_pathogenic", "criteria_provided,_single_submitter", "PASS", (True, True, False)),
    ("Pathogenic", "practice_guideline", "LowQual", (False, False, False)),
    ("Conflicting_classifications_of_pathogenicity", "reviewed_by_expert_panel", ".", (True, False, False)),
    ("Pathogenic|Benign", "reviewed_by_expert_panel", ".", (True, False, False)),
    ("Pathogenic", "unknown", ".", (True, True, False)),
    ("Uncertain_significance", "practice_guideline", ".", (True, False, False)),
])
def test_cumulative_filters(tmp_path, sig, status, quality, passes):
    records = read_vcf(vcf(tmp_path, [record(info=f"CLNSIG={sig};CLNREVSTAT={status}", quality=quality)]))
    row = annotate(records, records)[0]
    assert tuple(row[k] for k in ["pass_quality", "pass_pathogenic", "pass_review"]) == passes
    assert row["selected"] == passes[-1]


@pytest.mark.parametrize(("ref", "alt"), [("A", "G,T"), ("A", "<DEL>"), ("A", "*"), ("N", "A"), ("A", "A"), ("A", ".")])
def test_unsupported_alleles(tmp_path, ref, alt):
    with pytest.raises(ValueError, match="Unsupported alleles"):
        read_vcf(vcf(tmp_path, [record(ref=ref, alt=alt)]))


@pytest.mark.parametrize("body", [
    "1\t10\t.\tA\tG\t.\t.\t.\n",
    HEADER + "1\t10\tA\n",
    HEADER + record(pos=0),
    HEADER + record(pos="bad"),
    HEADER + record() + "\n" + record(),
    HEADER.replace("GRCh38", "GRCh37") + record(),
    HEADER.replace("GRCh38", "GRCh380") + record(),
    HEADER + record(info="CLNSIG"),
    HEADER + record(info="CLNSIG=Benign;CLNSIG=Pathogenic"),
    HEADER,
])
def test_malformed_vcf(tmp_path, body):
    path = tmp_path / "invalid.vcf"
    path.write_text(body)
    with pytest.raises(ValueError):
        read_vcf(path)


def test_exact_key_no_chr_alias_or_indel_normalization(tmp_path):
    source = read_vcf(vcf(tmp_path, [record(ref="AT", alt="A")]))
    different = read_vcf(vcf(tmp_path, [record(ref="AT", alt="T")], "different.vcf"))
    assert not annotate(different, source)[0]["matched"]
    path = vcf(tmp_path, [record().replace("1\t", "chr1\t", 1)], "alias.vcf")
    assert not annotate(read_vcf(path), source)[0]["matched"]


def test_end_to_end_positive_and_negative(tmp_path):
    info = "GENEINFO=GENE:123;MC=SO:0001583|missense_variant;CLNSIG=Pathogenic;CLNREVSTAT=reviewed_by_expert_panel"
    resource = vcf(tmp_path, [record(info=info), record(pos=11)], "resource.vcf")
    inputs = vcf(tmp_path, [record(), record(pos=11), record(pos=12)], "input.vcf")
    output = tmp_path / "output"
    summary = run(inputs, resource, output)
    assert summary["filter_stages"] == {"input": 3, "quality_eligible": 3,
                                        "pathogenic_aggregate": 1, "review_at_least_two_stars": 1}
    assert summary["matched_records"] == 2
    assert summary["unmatched_records"] == 1
    assert summary["missing_gene"] == 2
    with (output / "filtered.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 1 and rows[0]["gene"] == "GENE:123"
    assert (output / "summary.png").read_bytes().startswith(b"\x89PNG\r\n\x1a\n")


def test_snapshot_integrity_and_stripped_input():
    manifest = json.loads((ROOT / "data/manifest.json").read_text())
    for filename, checksum in manifest["sha256"].items():
        assert hashlib.sha256((ROOT / "data" / filename).read_bytes()).hexdigest() == checksum
    source = read_vcf(ROOT / "data/clinvar_source.vcf")
    sites = read_vcf(ROOT / "data/sites.vcf")
    assert len(source) == len(sites) == manifest["source_records"] == 300
    assert len({r.key for r in sites}) == 300
    assert [r.key for r in source] == [r.key for r in sites]
    assert all(not r.info and r.identifier == "." for r in sites)


def test_committed_summary_counts_and_tables():
    summary = json.loads((ROOT / "results/summary.json").read_text())
    with (ROOT / "results/annotated.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == summary["input_records"] == 300
    assert sum(r["matched"] == "True" for r in rows) == summary["matched_records"] == 300
    assert sum(summary["variant_classes"].values()) == 300
    assert sum(summary["clinical_significance"].values()) == 300
    counts = list(summary["filter_stages"][key] for key in
                  ["input", "quality_eligible", "pathogenic_aggregate", "review_at_least_two_stars"])
    assert counts == sorted(counts, reverse=True)
    with (ROOT / "results/filtered.tsv").open() as handle:
        selected = list(csv.DictReader(handle, delimiter="\t"))
    assert len(selected) == sum(r["selected"] == "True" for r in rows) == counts[-1]


def test_invalid_input_creates_no_partial_outputs(tmp_path):
    invalid = vcf(tmp_path, [record(alt="G,T")])
    output = tmp_path / "output"
    with pytest.raises(ValueError):
        run(invalid, ROOT / "data/clinvar_source.vcf", output)
    assert not output.exists()


def test_percent_encoded_info(tmp_path):
    records = read_vcf(vcf(tmp_path, [record(info="GENEINFO=A%3BB:1;SOMATIC")]))
    assert records[0].info["GENEINFO"] == "A;B:1"
    assert records[0].info["SOMATIC"] is True


def test_explicit_missing_info_values(tmp_path):
    records = read_vcf(vcf(tmp_path, [record(info="GENEINFO=.;MC=;CLNSIG=.;CLNREVSTAT=.")]))
    row = annotate(records, records)[0]
    assert all(row[name] == "NA" for name in ["gene", "consequence", "clinical_significance", "review_status"])
    assert not row["selected"]
