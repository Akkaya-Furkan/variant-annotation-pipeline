# One-day scope (2026-10-24)

## Question
Can a small public VCF be reproducibly annotated against a lightweight local
ClinVar snapshot, filtered in explicit layers, and summarized without a large
VEP/SnpEff installation?

## Implementation boundary
Use the first 300 records in the archived 2026-10-04 ClinVar GRCh38 VCF,
in source order, without selecting on clinical significance. Stream only the
required prefix of the compressed download. Keep the original annotated records
as a local lookup resource and generate an INFO-free sites VCF from them. Join
on exact GRCh38 chromosome, position, REF, ALT; copy gene, molecular consequence,
clinical significance and review status. Derive variant class and SNV substitution
class. This deliberately measures known-variant lookup, not independent biological
prediction. Reject multiallelic and symbolic records rather than misassign annotations.

Filter cumulatively: supported literal alleles and FILTER PASS or '.' (unfiltered);
then an explicitly pathogenic/likely pathogenic aggregate without conflict;
then a ClinVar review status equivalent to at least two stars. Missing or unknown
review status fails the last layer. These are demonstration rules, not clinical advice.

## Acceptance criteria (fixed before implementation)
- Commit a public-data snapshot with exactly 300 source records and 300 stripped
  input records; record source URL, assembly, release date and SHA-256 checksums.
- Produce exactly 300 annotated rows and 300 unique coordinate/allele keys,
  each with a match to the local source; preserve missing annotations explicitly.
- Produce an annotated TSV, filtered TSV, JSON summary, Markdown report and PNG
  summary figure (5 outputs); report all filter-stage counts, including zero.
- Offline reruns produce identical tables, summary, report and figure hashes.
- At least 12 pytest cases pass, including missing annotations, unmatched variants,
  indels, conflict/review rules, unsupported alleles, malformed VCF and full fixture
  integration; CI runs tests and verifies committed outputs without a network download.
- Provide a pinned Python environment, data download script, usage documentation
  and GitHub Actions workflow. Verify pushed main against local HEAD and read CI.

## Explicit exclusions
Whole-genome analysis, novel-variant consequence prediction, VEP/SnpEff databases,
reference-based normalization, population frequency annotation, genotypes, disease
interpretation and clinical validation are **not done**. Results from this positional
prefix are not representative of all ClinVar or of a human genome. The source VCF
already contains annotations; stripping and rejoining them is a transparent technical
round trip, not a benchmark of annotation accuracy or clinical sensitivity.
