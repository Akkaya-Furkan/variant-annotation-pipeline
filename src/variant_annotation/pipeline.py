"""Offline exact-allele annotation, explicit filters and reproducible reporting."""
import argparse
from collections import Counter
import csv
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import re
from urllib.parse import unquote

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


@dataclass(frozen=True)
class Variant:
    chrom: str
    pos: int
    identifier: str
    ref: str
    alt: str
    filter: str
    info: dict

    @property
    def key(self):
        return self.chrom, self.pos, self.ref, self.alt


def read_vcf(path):
    """Require declared GRCh38 and simple biallelic, literal A/C/G/T alleles.

    Fail on unsupported alleles instead of assigning record-level annotations to
    individual ALTs. No chromosome alias conversion or indel normalization occurs.
    """
    records, seen = [], set()
    header_seen, assembly_seen = False, False
    with Path(path).open(encoding="utf-8") as handle:
        for number, line in enumerate(handle, 1):
            line = line.rstrip("\r\n")
            if line == "##reference=GRCh38":
                assembly_seen = True
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM\t"):
                if line.split("\t")[:8] != ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]:
                    raise ValueError(f"Invalid VCF header at line {number}")
                header_seen = True
                continue
            if not line:
                continue
            if line.startswith("#") or not header_seen:
                raise ValueError(f"Missing or invalid VCF header at line {number}")
            fields = line.split("\t")
            if len(fields) < 8:
                raise ValueError(f"Expected at least eight VCF columns at line {number}")
            chrom, pos, identifier, ref, alt, _, quality_filter, raw_info = fields[:8]
            try:
                pos = int(pos)
            except ValueError as error:
                raise ValueError(f"Invalid position at line {number}") from error
            if pos < 1 or not chrom:
                raise ValueError(f"Invalid coordinate at line {number}")
            if not re.fullmatch("[ACGT]+", ref) or not re.fullmatch("[ACGT]+", alt) or ref == alt:
                raise ValueError(f"Unsupported alleles at line {number}: {ref}>{alt}")
            info = {}
            if raw_info != ".":
                for item in raw_info.split(";"):
                    key, separator, value = item.partition("=")
                    if not key or key in info:
                        raise ValueError(f"Empty or duplicate INFO key at line {number}")
                    if key in {"GENEINFO", "MC", "CLNSIG", "CLNREVSTAT"} and not separator:
                        raise ValueError(f"Annotation INFO value required at line {number}")
                    info[key] = unquote(value) if separator else True
            record = Variant(chrom, pos, identifier, ref, alt, quality_filter, info)
            if record.key in seen:
                raise ValueError(f"Duplicate coordinate/allele key at line {number}")
            seen.add(record.key)
            records.append(record)
    if not header_seen or not assembly_seen:
        raise ValueError("VCF must declare a column header and GRCh38 reference")
    if not records:
        raise ValueError("VCF contains no records")
    return records


def variant_class(ref, alt):
    if len(ref) == len(alt) == 1:
        return "SNV"
    if len(ref) == len(alt):
        return "MNV"
    return "insertion" if len(alt) > len(ref) else "deletion"


def review_stars(status):
    # Exact aggregate statuses from ClinVar; unknown status is not evidence.
    return {
        "practice_guideline": 4,
        "reviewed_by_expert_panel": 3,
        "criteria_provided,_multiple_submitters,_no_conflicts": 2,
        "criteria_provided,_single_submitter": 1,
        "criteria_provided,_conflicting_classifications": 1,
        "criteria_provided,_conflicting_interpretations": 1,
        "no_assertion_criteria_provided": 0,
        "no_assertion_provided": 0,
        "no_classification_provided": 0,
    }.get(status, 0)


def pathogenic(significance):
    # Do not use a substring match: 'Conflicting_...Pathogenic' must never pass.
    terms = set(re.split(r"[/|,]", significance))
    return bool(terms) and terms <= {"Pathogenic", "Likely_pathogenic"}


FIELDS = ["chrom", "pos", "ref", "alt", "variant_class", "substitution", "matched",
          "clinvar_id", "gene", "consequence", "clinical_significance", "review_status",
          "review_stars", "source_filter", "pass_quality", "pass_pathogenic",
          "pass_review", "selected"]


def annotate(variants, resource):
    lookup = {record.key: record for record in resource}
    rows = []
    for record in variants:
        match = lookup.get(record.key)
        info = {key: value for key, value in (match.info if match else {}).items()
                if value not in {"", "."}}
        kind = variant_class(record.ref, record.alt)
        substitution = "NA"
        if kind == "SNV":
            substitution = "transition" if {record.ref, record.alt} in ({"A", "G"}, {"C", "T"}) else "transversion"
        significance = info.get("CLNSIG", "NA")
        status = info.get("CLNREVSTAT", "NA")
        stars = review_stars(status)
        # A dot means unfiltered, not a successful upstream quality test.
        quality = record.filter in {"PASS", "."}
        sig_pass = quality and match is not None and pathogenic(significance)
        review_pass = sig_pass and stars >= 2
        rows.append(dict(zip(FIELDS, [
            record.chrom, record.pos, record.ref, record.alt, kind, substitution,
            match is not None, match.identifier if match else "NA",
            info.get("GENEINFO", "NA"), info.get("MC", "NA"), significance, status,
            stars, record.filter, quality, sig_pass, review_pass, review_pass,
        ], strict=True)))
    return rows


def write_table(path, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def run(input_path, resource_path, output):
    variants, resource = read_vcf(input_path), read_vcf(resource_path)
    rows = annotate(variants, resource)
    selected = [row for row in rows if row["selected"]]
    stages = {"input": len(rows), "quality_eligible": sum(r["pass_quality"] for r in rows),
              "pathogenic_aggregate": sum(r["pass_pathogenic"] for r in rows),
              "review_at_least_two_stars": len(selected)}
    summary = {
        "assembly": "GRCh38", "input_records": len(rows),
        "unique_keys": len({v.key for v in variants}), "resource_records": len(resource),
        "matched_records": sum(r["matched"] for r in rows),
        "unmatched_records": sum(not r["matched"] for r in rows),
        "missing_gene": sum(r["gene"] == "NA" for r in rows),
        "missing_consequence": sum(r["consequence"] == "NA" for r in rows),
        "filter_stages": stages,
        "variant_classes": dict(sorted(Counter(r["variant_class"] for r in rows).items())),
        "clinical_significance": dict(sorted(Counter(r["clinical_significance"] for r in rows).items())),
        "review_stars": dict(sorted(Counter(str(r["review_stars"]) for r in rows).items())),
        "regions": {chrom: {"min_pos": min(v.pos for v in variants if v.chrom == chrom),
                            "max_pos": max(v.pos for v in variants if v.chrom == chrom)}
                    for chrom in sorted({v.chrom for v in variants})},
        "input_sha256": hashlib.sha256(Path(input_path).read_bytes()).hexdigest(),
        "resource_sha256": hashlib.sha256(Path(resource_path).read_bytes()).hexdigest(),
    }
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    write_table(output / "annotated.tsv", rows)
    write_table(output / "filtered.tsv", selected)
    (output / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), layout="constrained")
    for ax, counts, title in [(axes[0], summary["variant_classes"], "Input variant classes"),
                               (axes[1], stages, "Cumulative demonstration filters")]:
        labels = list(counts)
        values = list(counts.values())
        bars = ax.bar(range(len(labels)), values, color="#237c8a")
        ax.set_xticks(range(len(labels)), [label.replace("_", "\n") for label in labels], fontsize=8)
        ax.set_ylabel("Variants")
        ax.set_title(title)
        ax.set_ylim(0, max(values, default=0) * 1.2 + 1)
        ax.bar_label(bars, padding=3)
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("ClinVar GRCh38 prefix — exact-allele lookup, no clinical interpretation")
    fig.savefig(output / "summary.png", dpi=150, metadata={"Software": "variant-annotation-pipeline"})
    plt.close(fig)
    report = ["# Annotation report", "", "Exact GRCh38 allele lookup against the bundled ClinVar snapshot.", "",
              f"Input: {len(rows)} records; unique keys: {summary['unique_keys']}; "
              f"matched: {summary['matched_records']}; unmatched: {summary['unmatched_records']}.", "",
              "| Cumulative stage | Variants |", "| --- | ---: |"]
    report.extend(f"| {name} | {value} |" for name, value in stages.items())
    report.extend(["", f"Missing gene: {summary['missing_gene']}; missing consequence: {summary['missing_consequence']}.",
                   "", "![Variant classes and filters](summary.png)", "",
                   "FILTER '.' is unfiltered and is admitted by this demonstration rule.",
                   "Unknown review status receives zero stars. Clinical interpretation is not done.",
                   "The input was stripped from this same source: match coverage is a round-trip check, not accuracy."])
    (output / "report.md").write_text("\n".join(report) + "\n")
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=Path("data/sites.vcf"))
    parser.add_argument("--resource", type=Path, default=Path("data/clinvar_source.vcf"))
    parser.add_argument("--output", type=Path, default=Path("results"))
    args = parser.parse_args()
    print(json.dumps(run(args.input, args.resource, args.output), indent=2))


if __name__ == "__main__":
    main()
