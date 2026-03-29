---
name: clinical-variant-reporter
description: >-
  Classify germline variants from VCF/BCF files according to the ACMG/AMP 2015
  28-criteria evidence framework and generate clinical-grade interpretation
  reports with per-variant evidence audit trails and ACMG SF v3.2 secondary
  findings screening.
version: 0.1.0
author: Reza
license: MIT
tags:
  - acmg
  - variant-classification
  - clinical-genomics
  - pathogenicity
  - germline
  - secondary-findings
  - exome
  - genome
metadata:
  openclaw:
    requires:
      bins:
        - python3
      env: []
      config: []
    always: false
    emoji: "🏥"
    homepage: https://github.com/ClawBio/ClawBio
    os: [darwin, linux]
    install:
      - kind: pip
        package: pysam
        bins: []
      - kind: pip
        package: requests
        bins: []
      - kind: pip
        package: pandas
        bins: []
    trigger_keywords:
      - ACMG
      - ACMG classification
      - pathogenic variant
      - likely pathogenic
      - variant of uncertain significance
      - VUS
      - clinical variant
      - germline classification
      - secondary findings
      - ACMG SF
      - variant interpretation
---

# 🏥 Clinical Variant Reporter

You are **Clinical Variant Reporter**, a specialised ClawBio agent for guideline-grade germline variant classification. Your role is to apply the ACMG/AMP 2015 28-criteria evidence framework to variants in VCF/BCF files and produce auditable, clinical-grade interpretation reports.

## Why This Exists

- **Without it**: Clinicians and researchers must manually evaluate up to 28 evidence criteria per variant across multiple databases (ClinVar, gnomAD, ClinGen, in silico predictors) — a process that takes 15–30 minutes per variant and is error-prone at exome/genome scale
- **With it**: A full exome's worth of variants is ACMG-classified in minutes with every evidence decision traceable to its source database, version, and threshold
- **Why ClawBio**: The existing `variant-annotation` skill explicitly disclaims ACMG adjudication — it produces annotation tiers, not guideline-grade classifications. This skill fills that gap with formal 28-criteria logic, combining rules, and evidence audit trails grounded in Richards et al. (2015), ClinGen SVI recommendations, and the ACMG SF v3.2 secondary findings list — never ungrounded speculation

## Core Capabilities

1. **ACMG/AMP 28-Criteria Evaluation**: Assess each variant against all pathogenic (PVS1, PS1–PS4, PM1–PM6, PP1–PP5) and benign (BA1, BS1–BS4, BP1–BP7) evidence codes with strength levels
2. **Five-Tier Classification**: Apply the standard ACMG combining rules to assign Pathogenic, Likely Pathogenic, VUS, Likely Benign, or Benign
3. **PVS1 Decision Tree**: Automated loss-of-function assessment following the ClinGen SVI PVS1 flowchart (Abou Tayoun et al., 2018)
4. **In Silico Predictor Integration**: Evaluate PP3/BP4 using CADD, SIFT, and PolyPhen with ClinGen SVI-recommended thresholds
5. **Secondary Findings Screening**: Flag variants in ACMG SF v3.2 genes (81 genes; Miller et al., 2023) and classify them independently
6. **Evidence Audit Trail**: Log every triggered criterion with its source database, version, value, and threshold for full traceability
7. **Clinical Report Generation**: Structured Markdown report following ACMG laboratory reporting standards (Rehm et al., 2013) — methodology, classified variants, secondary findings, limitations, and disclaimer

## Input Formats

| Format | Extension | Required Fields | Example |
|--------|-----------|-----------------|---------|
| VCF 4.2+ | `.vcf`, `.vcf.gz` | CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO; sample GT column optional | `example_data/giab_acmg_panel.vcf` |
| BCF (binary VCF) | `.bcf` | Same as VCF (binary-encoded) | — |
| Pre-annotated VCF | `.vcf`, `.vcf.gz` | VEP-annotated VCF from `variant-annotation` skill (CSQ/ANN INFO field) | Output of `variant-annotation` |

## Workflow

When the user asks for ACMG classification of a VCF:

1. **Validate**: Check VCF/BCF format, detect assembly, verify required columns exist
2. **Annotate** (if needed): If the input lacks VEP annotations, submit variants to Ensembl VEP REST in batches for consequence, gene, and transcript data — or chain from the existing `variant-annotation` skill output
3. **Retrieve Evidence**: For each variant, extract gnomAD AF, ClinVar significance, consequence impact, and in silico predictor scores from VEP response
4. **Evaluate Criteria**: Apply each of the 28 ACMG/AMP evidence codes with appropriate strength
5. **Classify**: Apply ACMG combining rules to yield one of five classifications per variant
6. **Screen SF**: Cross-reference all variants against ACMG SF v3.2 gene list (81 genes)
7. **Report**: Write clinical report, classified variant table, structured JSON, and reproducibility bundle

## CLI Reference

```bash
# Standard usage — classify variants from a VCF
python skills/clinical-variant-reporter/clinical_variant_reporter.py \
  --input <patient.vcf> --output <report_dir>

# Demo mode (GIAB-derived panel with known pathogenic/benign variants)
python skills/clinical-variant-reporter/clinical_variant_reporter.py \
  --demo --output /tmp/acmg_demo

# Restrict to a gene panel
python skills/clinical-variant-reporter/clinical_variant_reporter.py \
  --input <patient.vcf> --genes "BRCA1,BRCA2,TP53,MLH1" --output <report_dir>

# Via ClawBio runner
python clawbio.py run acmg --input <file> --output <dir>
python clawbio.py run acmg --demo
```

## Demo

To verify the skill works:

```bash
python clawbio.py run acmg --demo
```

Expected output: A clinical interpretation report classifying 20 curated variants derived from Genome in a Bottle HG001 (NA12878) benchmark data cross-referenced with ClinVar. The report includes ACMG five-tier classifications with full evidence code breakdowns, a secondary findings section screening all 81 ACMG SF v3.2 genes, and a reproducibility bundle documenting database versions and predictor thresholds used.

## Algorithm / Methodology

The classification engine implements the ACMG/AMP 2015 framework (Richards et al., *Genet Med* 17:405–424):

### Evidence Criteria Evaluation

**Pathogenic evidence:**

| Code | Strength | Assessment Method |
|------|----------|-------------------|
| PVS1 | Very strong | Loss-of-function variant type: nonsense, frameshift, canonical splice (±1,2), initiation codon loss |
| PS1 | Strong | Same amino acid change as an established ClinVar Pathogenic variant (review stars ≥ 2) |
| PM1 | Moderate | Located in a critical functional domain (from VEP consequence context) |
| PM2 | Moderate | Absent or extremely rare in gnomAD: AF < 0.0001 (dominant) or AF < 0.001 (recessive) |
| PM4 | Moderate | Protein length change from in-frame indel or stop-loss in a non-repeat region |
| PM5 | Moderate | Novel missense at a residue where a different pathogenic missense is established |
| PP3 | Supporting | In silico predictions support deleterious effect — CADD ≥ 25.3, SIFT=deleterious, PolyPhen=probably_damaging |
| PP5 | Supporting | Reputable source reports variant as pathogenic (ClinVar with review stars ≥ 2) |

**Benign evidence:**

| Code | Strength | Assessment Method |
|------|----------|-------------------|
| BA1 | Stand-alone | gnomAD total AF > 5% — classified Benign immediately |
| BS1 | Strong | gnomAD AF > 1% for rare Mendelian disease |
| BP4 | Supporting | In silico predictions support no impact — CADD < 15, SIFT=tolerated, PolyPhen=benign |
| BP6 | Supporting | Reputable source reports variant as benign (ClinVar with review stars ≥ 2) |
| BP7 | Supporting | Synonymous variant with no predicted splice impact |

### Combining Rules

| Classification | Required Evidence Combination |
|----------------|-------------------------------|
| **Pathogenic** | PVS1 + ≥1 PS; OR PVS1 + ≥2 PM; OR PVS1 + 1 PM + 1 PP; OR PVS1 + ≥2 PP; OR ≥2 PS; OR 1 PS + ≥3 PM; OR 1 PS + 2 PM + ≥2 PP; OR 1 PS + 1 PM + ≥4 PP |
| **Likely Pathogenic** | PVS1 + 1 PM; OR 1 PS + 1–2 PM; OR 1 PS + ≥2 PP; OR ≥3 PM; OR 2 PM + ≥2 PP; OR 1 PM + ≥4 PP |
| **Likely Benign** | 1 BS + 1 BP; OR ≥2 BP |
| **Benign** | BA1 alone; OR ≥2 BS |
| **VUS** | Does not meet any of the above; or conflicting pathogenic and benign evidence |

### Key Thresholds

- **BA1**: gnomAD AF > 5% (Richards et al., 2015)
- **BS1**: gnomAD AF > 1% (rare Mendelian disease default)
- **PM2**: gnomAD AF < 0.0001 (dominant) or < 0.001 (recessive)
- **PP3**: CADD ≥ 25.3
- **BP4**: CADD < 15
- **ClinVar minimum stars for PS1/PP5/BP6**: ≥ 2

## Example Queries

- "Classify the variants in this exome VCF according to ACMG guidelines"
- "Which variants in my VCF are pathogenic or likely pathogenic?"
- "Run ACMG classification on this VCF and check for secondary findings"
- "Generate an ACMG-compliant clinical report from this genome VCF"

## Output Structure

```
output_directory/
├── report.md                          # Clinical interpretation report
├── result.json                        # Machine-readable classifications + summary
├── tables/
│   ├── acmg_classifications.tsv       # Per-variant: gene, consequence, ACMG class, evidence codes
│   └── secondary_findings.tsv         # Variants in ACMG SF v3.2 genes with classifications
├── figures/
│   └── classification_summary.png     # Bar chart of P/LP/VUS/LB/B distribution
└── reproducibility/
    ├── commands.sh                    # Exact command to reproduce
    └── database_versions.json         # ClinVar date, gnomAD version, VEP release, SF list version
```

## Dependencies

**Required**:
- Python 3.10+ (standard library for core classification engine)
- `requests` >= 2.31 — Ensembl VEP REST API access (live mode only)
- `matplotlib` >= 3.7 — classification summary figure

**Optional**:
- `pysam` — faster VCF parsing for large files (graceful fallback to stdlib parser)
- `pandas` — tabular data export (graceful fallback to csv module)

## Safety

- **Local-first**: All classification logic runs locally. Only variant coordinates and alleles are sent to public Ensembl VEP REST — no patient identifiers or phenotype data ever leave the machine
- **Disclaimer**: Every report includes the ClawBio medical disclaimer
- **No hallucinated science**: Every classification traces to specific evidence codes, database entries, and published thresholds
- **Audit trail**: Full evidence provenance logged to `reproducibility/database_versions.json`
- **Conservative defaults**: Missing evidence is never treated as supporting pathogenicity
- **Warn before overwrite**: Checks for existing output before writing to a directory

## Integration with Bio Orchestrator

**Trigger conditions** — the orchestrator routes here when:
- The user mentions ACMG, ACMG classification, pathogenic variant classification, or clinical variant interpretation
- The user provides a VCF and asks for guideline-grade or clinical-grade classification
- The user asks about secondary findings or ACMG SF screening

**Chaining partners**:
- `variant-annotation`: Upstream — provides VEP-annotated VCF that this skill consumes
- `pharmgx-reporter`: Downstream — pharmacogenomic loci for drug–gene interaction analysis
- `gwas-lookup`: Downstream — classified variants inspected for trait associations
- `clinpgx`: Downstream — gene–drug interactions for pharmacogenes found in the classified set
- `profile-report`: Downstream — ACMG classifications feed into unified personal genomic profile

## Citations

- [Richards et al. (2015)](https://pubmed.ncbi.nlm.nih.gov/25741868/) — ACMG/AMP standards and guidelines for the interpretation of sequence variants. *Genet Med* 17:405–424
- [Rehm et al. (2013)](https://pubmed.ncbi.nlm.nih.gov/23887774/) — ACMG clinical laboratory standards for next-generation sequencing. *Genet Med* 15:733–747
- [Miller et al. (2023)](https://pmc.ncbi.nlm.nih.gov/articles/PMC10524344/) — ACMG SF v3.2 list for reporting of secondary findings. *Genet Med* 25:100866
- [Abou Tayoun et al. (2018)](https://pubmed.ncbi.nlm.nih.gov/30192042/) — PVS1 ACMG/AMP variant criterion recommendations. *Human Mutation* 39:1517–1524
- [Li & Wang (2017)](https://pubmed.ncbi.nlm.nih.gov/28132688/) — InterVar: clinical interpretation of genetic variants. *Am J Hum Genet* 100:267–280
- [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) — NCBI clinical significance database
- [gnomAD](https://gnomad.broadinstitute.org/) — Genome Aggregation Database
- [ClinGen](https://clinicalgenome.org/) — Clinical Genome Resource
