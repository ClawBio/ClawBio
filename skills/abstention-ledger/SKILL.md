---
name: abstention-ledger
description: >-
  Parse legacy SnpEff EFF annotations per transcript and emit an auditable
  abstention ledger for a family segregation table: which records can be reviewed,
  which are withheld, and the checked evidence behind every refusal.
license: MIT
metadata:
  version: 0.1.0
  author: ClawBio Berlin hackathon
  domain: genomics
  inputs:
    - name: input_file
      type: file
      format:
        - tsv
      description: Segregation table with per-role genotype columns and a legacy EFF column
      required: true
    - name: vcf
      type: file
      format:
        - vcf
      description: Matching VCF; enables the sample-to-role resolution check
      required: false
    - name: evidence
      type: file
      format:
        - json
      description: Cached build-matched annotation from fetch_evidence.py
      required: false
  outputs:
    - name: report
      type: file
      format: md
      description: Ledger, reproduced segregation, and cohort abstentions
    - name: result
      type: file
      format: json
      description: Machine-readable verdicts with per-gate evidence
  dependencies:
    python: ">=3.10"
  tags:
    - abstention
    - segregation
    - snpeff
    - legacy-eff
    - provenance
    - family
    - transcript-selection
  demo_data:
    - path: examples/demo_segregation.tsv
      description: Synthetic 8-record pedigree table exercising every Tier-1 gate
  endpoints:
    cli: python skills/abstention-ledger/abstention_ledger.py --input {input_file} --output {output_dir}
  openclaw:
    requires:
      bins:
        - python3
    always: false
    emoji: "\U0001F6D1"
    homepage: https://github.com/ClawBio/ClawBio
    os:
      - darwin
      - linux
    install: []
    trigger_keywords:
      - abstention
      - abstention list
      - legacy EFF
      - SnpEff EFF
      - segregation
      - parent of origin
      - what can this data not support
      - withheld variants
      - transcript selection
---

# Abstention Ledger

## Trigger

**Fire when:**

- A variant table must yield a review list *and* an explicit account of what it
  cannot support.
- The input carries a legacy SnpEff `EFF` column. This is the only skill in the
  library that parses it; every other skill reads `MC`, `Consequence` or `ANN`.
- A family segregation table needs its parent-of-origin labels re-derived rather
  than trusted.
- Someone asks "why is this variant *not* on the list?", "what can this data not
  support?", or "which of these were you not entitled to rank?".
- Sample-to-role assignment is uncertain or contradicted by its documentation.
- A reviewer needs a transcript-level view because `max(impact)` is hiding a
  disagreement.

**Do NOT fire when:**

- The question is what a variant *means clinically*. Route to
  `clinical-variant-reporter` for ACMG/AMP classification.
- The input is a modern `ANN=`/VEP-annotated VCF with no `EFF` column and no
  pedigree — `vcf-annotator` or `variant-annotation` fit better.
- Copy-number or structural variants are the subject. Route to
  `cnv-acmg-classifier`.
- A phenotype and HPO terms exist and phenotype-driven ranking is wanted. This
  skill cannot do it, and no skill in the library currently can.
- The user wants a ranked shortlist *without* the withheld set. Refusing to
  separate the two is the whole point; do not fire and then suppress half the
  output.

## Why This Exists

Two reasons, and the first is a gap the library documents about itself.

**1. Nothing here parses legacy `EFF`.** `rare-high-impact-variants` reads its
consequence from `MC`, `Consequence` or `ANN` (see
`rare_high_impact_variants.py`, the `info.get` chain). Classic SnpEff writes
`EFF`. So the skill whose report explains the population-frequency blind spot
receives an empty consequence for every record in that format, which is why the
hackathon brief warns: *"Do not feed the historical quartet directly to
rare-high-impact-variants: it does not parse the legacy EFF field."*
`legacy_eff.py` closes that gap.

**2. Parsing `EFF` properly turns out to be a gate.** `EFF` carries one
annotation per transcript. Pipelines routinely reduce that to `max(impact)`,
which converts a transcript-selection decision into an apparent biological fact.
A variant can be `START_LOST(HIGH)` on a minor transcript and
`MISSENSE(MODERATE)` on two canonical ones — same variant, same gene, same
record. Keeping every annotation is what makes the disagreement visible.

The output is therefore two artefacts. A review list, and a ledger of everything
withheld with the check, the triggering value, and its source. A reader can
disagree with a check. They cannot disagree with a conclusion that shows no
working.

## Core Capabilities

- Per-transcript legacy `EFF` parsing, with comma handling that survives commas
  inside codon and amino-acid fields.
- Sample-to-role resolution decided **from genotypes**, by testing all 24
  assignments of four samples to four roles and keeping those that reproduce
  every genotype in the table.
- Parent-of-origin re-derivation, reported alongside what a swapped role
  assignment would have produced.
- Cohort abstentions: dataset-level facts that bound every per-record claim.
- Per-record gates, each attaching the value that fired it.
- Optional build-matched evidence layer via the Ensembl GRCh37 REST endpoint,
  kept in separate columns from any supplied labels.

## Scope

**In scope.** Deciding what a variant table can and cannot support, and showing
the working.

**Out of scope.** Clinical interpretation, phenotype-driven prioritisation
(no HPO ingestion exists here), phasing, structural variants, and any statement
about an individual's health. A `REVIEWABLE` verdict means the record survived
the checks this skill could run — nothing more.

## Input Formats

Tab-separated, with `CHROM POS ID REF ALT`, one `<ROLE>_GT_DP_GQ` column per
family member (`SON`, `FATHER`, `SISTER`, `MOTHER`), an optional supplied
parent-of-origin column, and an `EFF` column in classic SnpEff format.

## Workflow

1. Load the table; compute the input checksum.
2. Resolve sample-to-role mapping against the VCF, if supplied.
3. Re-derive parent-of-origin; report agreement with any supplied labels, and the
   cost of the alternative role assignment.
4. Parse `EFF` per transcript for every record.
5. Apply cohort abstentions, then per-record gates.
6. Join the cached evidence layer if present; otherwise mark records
   `NO_EVIDENCE_LAYER` rather than assuming a value.
7. Write `report.md`, `result.json`, `tables/`, `reproducibility/`.

## CLI Reference

```bash
# Synthetic demo, no network
python skills/abstention-ledger/abstention_ledger.py --demo --output /tmp/al

# Real table, with the sample-map check enabled
python skills/abstention-ledger/abstention_ledger.py \
    --input data/challenge1/challenge1-b37-segregation.tsv \
    --vcf   data/challenge1/challenge1-b37-segregation.vcf.gz \
    --output out/run1

# Fetch the build-matched evidence layer once, then reuse the cache
python skills/abstention-ledger/fetch_evidence.py \
    --input data/challenge1/challenge1-b37-segregation.tsv \
    --output out/vep_grch37_cache.json
```

## Demo

`--demo` runs `examples/demo_segregation.tsv`: eight synthetic records, no
network. Expected outcome — 3 withheld (one transcript-severity disagreement, one
olfactory-receptor family locus, one mucin family locus with a second
disagreement), 1 further withheld on position (extended MHC), 4 reviewable.

## Reason Codes

### Cohort abstentions — properties of the dataset

| Code | Blocks |
|---|---|
| `NO_PHENOTYPE` | Any statement tying a variant to this person's clinical picture |
| `UNPHASED` | Any statement about two variants sharing a copy |
| `SELECTION_BIAS` | Any count or burden presented as describing the family |
| `HISTORICAL_ANNOTATION` | Any claim that supplied effect labels are current |
| `NON_PROBAND_DISCLOSURE` | Treating the output as being about the proband alone |

### Per-record gates

| Code | Fires when |
|---|---|
| `TRANSCRIPT_ARTIFACT` | Impact disagrees across transcripts of the same gene |
| `LOW_COMPLEXITY_LOCUS` | Gene family or position where short-read calls are unreliable |
| `NO_FREQUENCY_RECORD` | Build-matched query returned no frequency — a fact about the database, not the variant |
| `FREQUENCY_DOCUMENTED_COMMON` | Documented frequency at or above 1% |
| `ANNOTATION_SUPERSEDED` | Current annotation does not reproduce the supplied impact |
| `NO_EVIDENCE_LAYER` | No build-matched annotation was available |

## Domain Decisions

- **Transcript disagreement is scoped to a single gene.** A HIGH call in gene A
  beside a MODIFIER call in overlapping gene B is two genes, not a disagreement.
- **`NON_PROBAND_DISCLOSURE` is a cohort abstention, not a per-record gate.** In
  a one-carrier-parent file every record discloses a parent by construction, so as
  a veto it fires on 100% of records and discriminates nothing. It was written
  that way first here, and it made the review list trivially empty.
- **Absence of a frequency record is never treated as evidence of anything.** It
  is reported as a property of the reference database.
- **The evidence layer is build-matched on purpose.** The data is GRCh37; gnomAD
  v4 is GRCh38-native, so a liftover would sit in the path and a failed liftover
  is indistinguishable from a missing record. Querying the GRCh37 endpoint removes
  that ambiguity — at the cost of an older, exome-only frequency set, which the
  report states rather than hides.
- **`REVIEWABLE` is a claim about our evidence, not about the variant.**

## Example Output

Real excerpt from `report.md` on the four-person teaching pedigree:

```markdown
# Abstention Ledger

Records: **68**
Assembly: GRCh37/b37 (contigs without chr prefix)

**61 of 68 records are withheld from the review list.** 7 survived every check.

## 1. Which sample is which person

Resolved from genotypes: **1 of 24** possible assignments of sample IDs to family
roles reproduces every genotype in the table.

| Role assignment          | paternal | maternal | ambiguous | no carrier parent | disagreements |
|--------------------------|----------|----------|-----------|-------------------|---------------|
| as labelled              | 30       | 38       | 0         | 0                 | **0 / 68**    |
| sister and mother swapped| 11       | 25       | 19        | 13                | **32 / 68**   |

## 4. Per-record gates

| Reason code                   | Records |
|-------------------------------|---------|
| `FREQUENCY_DOCUMENTED_COMMON` | 54      |
| `TRANSCRIPT_ARTIFACT`         | 22      |
| `LOW_COMPLEXITY_LOCUS`        | 19      |
| `ANNOTATION_SUPERSEDED`       | 4       |
| `NO_FREQUENCY_RECORD`         | 2       |

### Checks that ran and found nothing

- **ACMG secondary findings:** all 68 records screened against 81 genes — **0 hits**.
```

And one ledger row, showing the shape of a refusal:

```
variant_id  22:32875190 rs11107 G>A
genes       FBXO7
verdict     WITHHELD
codes       TRANSCRIPT_ARTIFACT
evidence    supplied max impact HIGH is not reproduced across transcripts of the
            same gene — FBXO7: NM_012179.3=MODERATE(NON_SYNONYMOUS_CODING);
            NM_001033024.1=MODERATE(NON_SYNONYMOUS_CODING);
            NM_001257990.1=HIGH(START_LOST)
source      legacy EFF field of the input record, parsed per transcript
```

## Safety Rules

Per ClawBio safety rule 2, every generated report carries the standard
disclaimer verbatim: *"ClawBio is a research and educational tool. It is not a
medical device and does not provide clinical diagnoses. Consult a healthcare
professional before making any medical decisions."* It is asserted by
`tests/test_abstention_ledger.py`, not left to reviewer discipline.

- Never assert that a variant causes, or is likely to cause, disease.
- Never describe a frequency as low without a documented value.
- Never present an inferred parent-of-origin label as molecular phase.
- Never let a failed network call become an empty result: fetch failures are
  recorded with their reason, and downstream records carry `NO_EVIDENCE_LAYER`.
- Never let an unrun check look like a passed one. If the secondary-findings list
  fails to load, records carry `SF_LIST_UNAVAILABLE`; if no evidence layer was
  fetched, the report says `not checked` rather than `no record found`.
- Demo data is synthetic. Do not ship real individual genotypes in this
  directory.

## Agent Boundary

The skill decides what is *supportable*, never what is *true of a person*. It
returns verdicts and evidence; it does not counsel, does not rank by predicted
severity, and does not produce a clinical report. Where a reviewer would need
phenotype, consent scope, or phase to proceed, the skill names the missing input
instead of substituting a default.

## Output Structure

```
output/
├── report.md
├── result.json
├── tables/
│   ├── segregation.tsv
│   └── abstention_ledger.tsv
└── reproducibility/
    ├── commands.sh
    └── checksums.sha256
```

## Dependencies

Python standard library only. No pip install, no compiled dependency, no API key.
The evidence layer uses one batched HTTP POST and caches the result to disk.

## Gotchas

- The GRCh37 Ensembl endpoint answers **503** for parameters it does not support
  (`mane`, `af_1kg`, `af_gnomade`, `af_gnomadg`). It looks like an outage and is
  not one. `PARAMS` in `fetch_evidence.py` is verified parameter by parameter.
- b37 contigs have no `chr` prefix. Annotators defaulting to UCSC naming return
  nothing rather than erroring.
- A supplied parent-of-origin column is a label, not a measurement. Re-derive it.

## Chaining Partners

| Skill | Direction | How |
|---|---|---|
| `rare-high-impact-variants` | downstream | `eff_to_info.py` translates legacy `EFF` into the `MC`/`GENEINFO` keys it reads, so it works on legacy-annotated input instead of silently returning zero. `prove_gap.py` measures the before and after. |
| `clinical-variant-reporter` | upstream | We import its `ACMG_SF_V32_GENES` and `is_secondary_finding_gene` rather than keeping a second copy of the list. Route to it for actual ACMG/AMP classification. |
| `cnv-acmg-classifier` | sibling | Handles what this skill structurally cannot see. An SNV-only input should carry an explicit note that copy-number evidence was never examined. |
| `nfcore-sarek-wrapper` | upstream | Produces the VCF; annotate, then feed the segregation table here. |
| `rare-disease-rnaseq` | sibling | Expression outliers are the tie-breaker this skill lacks. Where a record is withheld for want of evidence, RNA is one way to get some. |
| `lit-synthesizer`, `clinical-trial-finder` | downstream | Only for records that reached `REVIEWABLE`. Running literature search over a withheld variant manufactures the significance the ledger just declined to assign. |

## Maintenance

**Review cadence:** whenever the ACMG SF list version changes, or the Ensembl
GRCh37 REST parameter surface changes.

**Staleness signals:**

- `ACMG_SF_V32_GENES` is ACMG SF **v3.2** (81 genes). The current statement is
  **v3.3** (84 genes; added `ABCD1`, `CYP27A1`, `PLN`). When the upstream skill
  updates, this skill inherits it — screening against a superseded list is exactly
  the class of error this skill exists to name, so the version is printed in every
  report rather than assumed.
- `LOW_COMPLEXITY_FAMILIES` is a curated prefix list, not a coordinate track. If
  reviewers start disagreeing with its calls, replace it with a real segmental
  duplication / repeat-masker overlap.
- `fetch_evidence.PARAMS` is verified parameter by parameter against the live
  endpoint. If a fetch starts returning 503, suspect a parameter before an outage.

**Deprecation criteria:** retire the `legacy_eff.py` half if upstream skills gain
native `EFF` support. The gates stay useful regardless of annotation format, so
they should be split out rather than retired with it.

## Citations

- SnpEff classic `EFF` field specification — SnpEff documentation.
- McLaren W, et al. The Ensembl Variant Effect Predictor. *Genome Biol.* 2016;17:122. doi:10.1186/s13059-016-0974-4
- Lee K, Abul-Husn NS, et al. ACMG SF v3.3 list for reporting of secondary findings in clinical exome and genome sequencing. *Genet Med.* 2025;27(8):101454. doi:10.1016/j.gim.2025.101454
- Corpas M, et al. Crowdsourced direct-to-consumer genomic analysis of a family quartet. *BMC Genomics.* 2015;16:910. doi:10.1186/s12864-015-1973-7
