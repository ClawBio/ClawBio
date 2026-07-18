---
name: population-equity-auditor
description: >-
  Audit a cohort's small-variant (SNV) calls for reference-panel referral inflation: re-classify each
  variant under gnomAD-blind vs population-aware allele frequency, quantify how many actionable ACMG
  calls are reference-panel artifacts, and show which safeguard (ClinGen PVS1 gene-mechanism gate or
  population frequency) corrects each. Ancestry-stratified, assay-scope-aware, with a
  never-suppress-a-true-pathogenic safety invariant.
version: 0.1.0
author: Manuel Corpas
domain: genomics
license: MIT

inputs:
  - name: input_file
    type: file
    format: [tsv, csv, vcf]
    description: Annotated cohort table (consequence, CADD, gnomAD AF, gnomAD-AMR AF, per-subpopulation allele counts)
    required: true

outputs:
  - name: report
    type: file
    format: md
    description: Equity-audit report with the referral-inflation table and safety section
  - name: result
    type: file
    format: json
    description: Machine-readable audit with per-variant 2x2 classification grid and content hash

dependencies:
  python: ">=3.11"
  packages: []

tags: [acmg, equity, ancestry, population-frequency, referral, gnomad, reference-panel-bias, variant-classification]

demo_data:
  - path: demo_input.txt
    description: Real Peru high-impact-LoF cohort (150 WGS SNV genomes, 7 subpopulations; Borda/Guio/O'Connor)

endpoints:
  cli: python skills/population-equity-auditor/population_equity_auditor.py --input {input_file} --output {output_dir}

metadata:
  openclaw:
    requires:
      bins:
        - python3
      env: []
      config: []
    always: false
    emoji: ⚖️
    homepage: https://github.com/ClawBio/ClawBio
    os: [macos, linux]
    install: []
    trigger_keywords:
      - reference-panel bias
      - population-aware allele frequency
      - gnomAD-blind
      - referral inflation
      - VUS inflation
      - ancestry equity audit
      - false pathogenic
      - underrepresented population
      - PM2 over-call
      - population equity auditor
---

# ⚖️ Population Equity Auditor

You are **Population Equity Auditor**, a specialised ClawBio agent. Your role is to measure how much of a
cohort's ACMG *actionable* burden is an artifact of using a Euro-biased population reference (gnomAD) to
classify variants from an under-represented population, and to show that the artifact is removed by two
independent safeguards without ever suppressing a true pathogenic variant.

## Trigger

**Fire this skill when the user says any of:**
- "audit this cohort for reference-panel / gnomAD bias"
- "how many of these actionable calls are artifacts of gnomAD-blindness?"
- "re-classify under population-aware allele frequency"
- "is our VUS / referral rate inflated for this ancestry?"
- "run population-equity-auditor"

**Do NOT fire when:**
- The user wants plain ACMG classification of one VCF (use `clinical-variant-reporter`)
- The user wants CNV dosage classification (use `cnv-acmg-classifier` — a different variant class)
- The user wants diversity/FST summary statistics (use `equity-scorer`)
- The user wants GWAS risk by ancestry (use `ancestry-risk-profiler`)

## Why This Exists

- **Without it**: A hospital deploying agentic ACMG automation cannot tell whether an "actionable"
  finding in an under-represented patient is real or an artifact of the reference panel. gnomAD
  under-samples non-European ancestries, so a variant that is common (and benign) in the patient's own
  population can look "absent / ultra-rare" (PM2, moderate pathogenic) and, combined with an ungated PVS1
  loss-of-function call, be reported as **Likely Pathogenic**. The bundled synthetic demo reproduces exactly
  this artifact; on real diverse cohorts a substantial fraction of high-impact LoF SNVs are affected.
- **With it**: Every variant is classified under a 2x2 grid — {gnomAD-global, gnomAD-AMR, population} x
  {naive, ClinGen-PVS1-gated} — and the report states how many actionable calls are reference-panel
  artifacts, which safeguard corrects each, the referral rate under each configuration, and any variant
  where population frequency would *mask* a genuine pathogenic (flagged, never auto-benigned).
- **Why ClawBio**: This is the trust-and-equity layer for agentic genomics. It reuses the shipped
  `clinical-variant-reporter` ACMG engine and ClawBio's ClinGen/gnomAD reference data — no ungrounded
  scoring — and turns the 2025 diversity-gap thesis into a per-cohort, reproducible, auditable measurement.

## Core Capabilities

1. **2x2 configuration audit**: classify each variant under three frequency sources (gnomAD-global,
   gnomAD-AMR, cohort population) x two pathogenic-evidence safeguards (naive; ClinGen PVS1 gene-mechanism
   gated) using the ACMG/AMP 2015 combining rules.
2. **Reference-panel referral inflation**: count false *actionable* (P/LP) calls under naive + gnomAD-blind
   that are not actionable under the trustworthy (gated + population-aware) config, and attribute each
   correction to the PVS1 gate, population frequency, or both.
3. **Ancestry stratification**: per-subpopulation allele frequency with Wilson 95% confidence intervals.
4. **Assay-scope comparability**: each cohort declares assay/caller/variant-classes/build; cross-cohort
   comparisons are restricted to the shared variant subset (e.g. SNV-only) and the restriction is reported.
5. **Safety invariant**: population frequency must never silently downgrade a variant that is actionable
   from gene-mechanism-gated evidence (a possible real founder pathogenic common in the cohort); such
   variants are flagged for expert review.

## Input Formats

The skill is **cohort-agnostic**. Provide a TSV plus a metadata sidecar:

| File | Contents |
|------|----------|
| `<cohort>.tsv` | `LOCATION REF ALT SYMBOL CADD_PHRED Consequence CLIN_SIG <one ALT-allele-count column per subpopulation> GNOMAD_AF GNOMAD_AF_AMR GNOMAD_AN GNOMAD_AN_AMR` |
| `<cohort>.meta.json` | `cohort_id, assay, caller, build, variant_classes, subpopulation_sizes {name: N}, notes` |

No cohort structure is hard-coded — subpopulation sizes, assay, caller, and build all come from the metadata,
so any hospital points the skill at its own annotated cohort. Cohort allele frequency =
sum(subpopulation ALT allele counts) / (2 x total N). See `data/demo_cohort.*` for the schema.

## Data governance & model independence

- **Model-independent**: classification is a **deterministic** rule engine (ACMG/AMP Richards 2015 + ClinGen
  PVS1 gate + Tavtigian points). No LLM, no ML inference — reproducible and auditable, with a content hash on
  every result. The skill is the model-independent *oracle* against which agentic LLMs (closed and open) can
  be benchmarked; it does not itself run a model.
- **Sealed / local-first**: the skill makes **zero network calls**. It operates on a pre-annotated cohort
  table. For a clinical deployment, produce the gnomAD/VEP/CADD annotation **locally** so no cohort-derived
  coordinates or counts ever leave the environment.
- **No patient data bundled or transmitted**: ships only public gene-level reference tables and a synthetic
  demo. See `data/README.md`.

## Workflow

1. **Validate**: detect the annotated-cohort format; confirm required columns; build the `CohortScope`
   descriptor (assay, caller, variant classes, build).
2. **Classify (2x2 grid)**: for each variant, run the ACMG engine five ways — naive+gnomAD-global,
   naive+gnomAD-AMR, naive+population, gated+gnomAD-global, gated+population.
3. **Attribute corrections**: mark false-actionable variants and whether the gate, frequency, or both
   removed the actionable call.
4. **Check safety**: flag any variant where population frequency masks gate-surviving pathogenic evidence.
5. **Report**: write the referral-inflation table, the corrected-variant table, per-subpopulation strata,
   the safety section, a content hash, and a reproducibility bundle.

## CLI Reference

```bash
# Audit an annotated cohort TSV
python skills/population-equity-auditor/population_equity_auditor.py \
  --input <cohort.tsv> --output <report_dir>

# Demo on the bundled real Peru high-impact-LoF cohort
python skills/population-equity-auditor/population_equity_auditor.py --demo --output /tmp/pea_demo

# Via ClawBio runner
python clawbio.py run population-equity-auditor --demo
```

## Demo

```bash
python clawbio.py run population-equity-auditor --demo
```

Expected output: an audit of the bundled synthetic demo cohort. The report shows false-actionable calls
under naive + gnomAD-blind automation (LoF-tolerant genes common in the cohort), each corrected by the
ClinGen PVS1 gate, population-aware frequency, and the PM2-strength update; a control that is not flagged;
and one flagged safety case (a ClinGen-haploinsufficient gene where population frequency would mask a real
pathogenic — surfaced for review, never auto-benigned).

## Algorithm / Methodology

Frequency enters ACMG classification only through the population-AF-driven criteria, so substituting the
frequency source is the entire equity lever:

| Code | Strength | Fires when |
|------|----------|-----------|
| BA1 | Stand-alone benign | AF > 5% |
| BS1 | Strong benign | AF > 1% |
| PM2 | Moderate pathogenic | AF < 1e-4 or absent |
| PVS1 | Very strong pathogenic | LoF consequence — **gated** by ClinGen gene-mechanism in the hardened config |

A LoF SNV that is common in the cohort (BA1/BS1) but absent/ultra-rare in gnomAD (PM2) flips between
Benign and Likely Pathogenic purely on which allele frequency the engine is given. The ClinGen PVS1
gene-mechanism gate (Abou Tayoun 2018) independently withholds PVS1 unless loss-of-function is an
established/inferred disease mechanism for the gene (ClinGen HI=3 or gnomAD LOEUF<0.35 / pLI>0.9).

**Key thresholds**: BA1 AF>0.05; BS1 AF>0.01; PM2 AF<1e-4; PVS1 gate LOEUF<0.35 or pLI>0.9. Combining rules:
ACMG/AMP Richards et al. 2015.

## Example Queries

- "Audit this Peruvian cohort — how many actionable calls are gnomAD-bias artifacts?"
- "Re-classify these LoF variants under the population's own allele frequency."
- "Is our VUS rate inflated for under-represented ancestries, and which safeguard fixes it?"

## Output Structure

```
output_directory/
├── report.md                       # Referral-inflation table, corrected variants, safety section
├── result.json                     # Per-variant 2x2 grid, cohort scope, safety, content hash
├── tables/
│   └── per_variant_audit.tsv       # One row per variant: naive-blind vs trustworthy class + correction
└── reproducibility/
    └── commands.sh                 # Exact command + engine/data provenance + content hash
```

## Dependencies

**Required**:
- Python 3.11+ (standard library only for the auditor)
- The sibling `clinical-variant-reporter` skill (its `acmg_engine.py` is imported for classification)

**Bundled reference data** (public, shipped under `data/`):
- ClinGen haploinsufficiency table (GRCh38) + provenance
- gnomAD per-gene constraint (LOEUF / pLI)

## Gotchas

- **Frequency source is the whole point**: never collapse the three frequency baselines into one. gnomAD
  *has* an AMR panel; the audit reports harm under gnomAD-global **and** gnomAD-AMR so a reviewer cannot
  dismiss the finding as "use the AMR panel".
- **SNV-only cohorts**: a caller run in SNP-only mode (e.g. GATK UnifiedGenotyper) omits frameshift indels,
  a major LoF class. The `CohortScope` records this and comparisons are restricted to SNVs. Do not compare
  an array cohort's referral rate to a WGS cohort's as if equal.
- **Population frequency can mask a real founder pathogenic**: high cohort AF must not auto-benign a
  gene-mechanism-supported pathogenic (BA1 has documented ClinGen exceptions). These are flagged, never
  silently downgraded.

## Safety

- **Local-first**: classification and reference lookups run locally; no patient identifiers leave the machine
- **Disclaimer**: every report includes the ClawBio medical disclaimer
- **Never-suppress-a-true-pathogenic**: population frequency cannot remove gate-surviving pathogenic
  evidence; masked candidates are flagged for expert review
- **Provenance**: every result carries a content hash; every ACMG code traces to its engine, database, and threshold

## Agent Boundary

The agent (LLM) dispatches and explains. The skill (Python) executes the ACMG grid and the safety check.
The agent must NOT override thresholds, re-label a flagged masked variant, or invent associations.

## Integration with Bio Orchestrator

**Trigger conditions**: the orchestrator routes here when the user asks about reference-panel / gnomAD bias,
population-aware re-classification, ancestry-stratified VUS/referral inflation, or false actionable calls in
an under-represented cohort.

**Chaining partners**:
- `variant-annotation` / `vcf-annotator`: upstream — supply the VEP consequence + gnomAD AF + CADD the audit consumes
- `clinical-variant-reporter`: upstream engine — the ACMG classifier this skill re-runs under each configuration
- `equity-scorer`: complementary — cohort-level diversity/FST metrics alongside this per-variant referral audit
- `profile-report`: downstream — flagged variants feed expert review

## Maintenance

- **Review cadence**: re-evaluate when ClinGen HI or gnomAD constraint releases update (refresh `data/`)
- **Staleness signals**: new gnomAD major version (frequency baselines shift), ClinGen dosage curation updates
- **Deprecation**: archive to `skills/_deprecated/` only if superseded by a more comprehensive equity auditor

## Citations

- [Richards et al. (2015)](https://pubmed.ncbi.nlm.nih.gov/25741868/) — ACMG/AMP variant interpretation standards. *Genet Med* 17:405–424
- [Abou Tayoun et al. (2018)](https://pubmed.ncbi.nlm.nih.gov/30192042/) — ClinGen PVS1 refinement. *Hum Mutat* 39:1517–1524
- [Pejaver et al. (2022)](https://pubmed.ncbi.nlm.nih.gov/36413997/) — calibrated in-silico (PP3/BP4) thresholds. *Am J Hum Genet*
- [gnomAD](https://gnomad.broadinstitute.org/) — Genome Aggregation Database (allele frequency + constraint)
- [ClinGen Dosage Sensitivity](https://clinicalgenome.org/) — haploinsufficiency curation
- Borda, Guio, O'Connor et al. (2025), PMID 41031014 — Peruvian population genomics (source of the demo cohort)
