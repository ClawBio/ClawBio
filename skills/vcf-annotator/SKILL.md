---
name: vcf-annotator
description: "Annotate VCF variants with Ensembl VEP, ClinVar, and gnomAD. Ranks variants by impact (HIGH/MODERATE/LOW/MODIFIER) and generates a reproducible report."
version: 0.1.0
author: Sooraj (github.com/sooraj-codes)
domain: genomics
license: MIT
emoji: "🧬"
os: [darwin, linux]
inputs:
  - name: input
    type: file
    format: vcf
    description: "VCF file (VCFv4.x, GRCh38)"
    required: true
outputs:
  - name: report
    type: file
    format: md
    description: Annotated variant report with ClinVar, gnomAD, and VEP results
dependencies:
  python: ">=3.11"
  packages: []
tags:
  - vcf
  - variants
  - annotation
  - clinvar
  - gnomad
  - vep
  - genomics
demo_data:
  - path: examples/demo_output/report.md
    description: Pre-generated demo report for 5 clinically relevant variants
endpoints:
  cli: python skills/vcf-annotator/vcf_annotator.py --input {input} --output {output_dir}
metadata:
  openclaw:
    requires:
      always: false
    homepage: https://github.com/ClawBio/ClawBio
    os: [darwin, linux]
    emoji: "🧬"
    install: []
    trigger_keywords:
      - "annotate vcf"
      - "annotate variants"
      - "variant annotation"
      - "clinvar lookup"
      - "gnomad frequency"
      - "vep annotation"
      - "pathogenic variants"
      - "variant effect"
      - "annotate my vcf"
      - "what variants are pathogenic"
---

# 🧬 VCF Annotator

You are **VCF Annotator**, a specialised ClawBio agent for genomic variant
annotation and interpretation. Your role is to annotate VCF files using
Ensembl VEP, ClinVar, and gnomAD, rank variants by predicted impact,
and generate a structured reproducible report.

## Trigger

**Fire this skill when the user says any of:**

- "annotate my VCF file"
- "annotate variants in X"
- "what variants are pathogenic"
- "look up ClinVar significance"
- "get gnomAD frequencies"
- "run VEP on my VCF"
- "variant annotation"
- "which variants are HIGH impact"
- "rank my variants by impact"

**Do NOT fire when:**

- The user wants pharmacogenomic drug recommendations (route to `pharmgx-reporter`)
- The user wants population PCA (route to `ancestry-pca`)
- The user wants literature search (route to `lit-synthesizer`)

## Why This Exists

**Without it**: A researcher must install VEP locally, configure databases,
query ClinVar and gnomAD separately, manually merge results, and format a report.
This takes hours and is error-prone.

**With it**: One command annotates a VCF against three authoritative databases,
ranks variants by impact, and outputs a reproducible report in seconds.

**Why ClawBio**: A general LLM will hallucinate ClinVar classifications and
invent gnomAD frequencies. This skill uses live API calls to real databases,
so every annotation is real and verifiable.

## Core Capabilities

1. **VCF parsing**: Reads VCFv4.x files, handles SNVs and indels
2. **Ensembl VEP**: Consequence prediction (missense, stop_gained, frameshift, etc.)
3. **ClinVar lookup**: Pathogenicity classification per variant
4. **gnomAD frequency**: Global and population-specific allele frequencies
5. **Impact ranking**: Sorts variants HIGH → MODERATE → LOW → MODIFIER
6. **Reproducibility bundle**: Exports `commands.sh`, `environment.yml`, SHA-256 checksums

## Scope

This skill annotates variants from a VCF file. It does **not** call variants
from raw sequencing reads (use a variant caller for that) or interpret
clinical significance beyond what ClinVar reports.

## Input Formats

| Format | Extension | Required Fields | Example |
|--------|-----------|-----------------|---------|
| VCF v4.x | `.vcf` | CHROM, POS, REF, ALT | `demo_variants.vcf` |

**Supported genome builds**: GRCh38 (primary), GRCh37 (legacy)

## Workflow

1. **Parse VCF**: Read variants, extract CHROM/POS/REF/ALT/rsID
2. **VEP annotation**: Query Ensembl REST API for consequence and gene
3. **ClinVar lookup**: Query NCBI E-utilities for pathogenicity classification
4. **gnomAD frequency**: Query gnomAD GraphQL API for allele frequencies
5. **Impact ranking**: Sort by HIGH → MODERATE → LOW → MODIFIER
6. **Report**: Write `report.md` with variant table, detailed annotations, and reproducibility bundle

## CLI Reference

```bash
# Standard usage
python skills/vcf-annotator/vcf_annotator.py \
    --input variants.vcf \
    --output report/

# Demo mode (no network, no VCF file needed)
python skills/vcf-annotator/vcf_annotator.py \
    --demo --output /tmp/demo

# Via ClawBio runner
python clawbio.py run vcf-annotator --input variants.vcf --output report/
python clawbio.py run vcf-annotator --demo
```

## Demo

```bash
python clawbio.py run vcf-annotator --demo
```

Expected output: A report covering 5 clinically relevant variants (BRCA1, BRCA2,
CFTR, APOE, MTHFR) with ClinVar classifications and gnomAD frequencies.

## Algorithm / Methodology

1. **VCF parsing**: Line-by-line reader, skips `#` headers, splits on tabs
2. **VEP**: `GET https://rest.ensembl.org/vep/human/hgvs/{hgvs}` — returns
   gene symbol, consequence terms, impact, SIFT, PolyPhen
3. **ClinVar**: `esearch` on `clinvar` database with rsID term
4. **gnomAD**: GraphQL query to `https://gnomad.broadinstitute.org/api` with
   variant ID format `{chrom}-{pos}-{ref}-{alt}`
5. **Ranking**: `HIGH=1, MODERATE=2, LOW=3, MODIFIER=4, UNKNOWN=5`

**Key thresholds**:

- gnomAD AF < 0.01 = rare variant
- gnomAD AF > 0.05 = common variant (less likely causal for rare disease)
- ClinVar "Pathogenic" or "Likely pathogenic" = flag for review

## Example Queries

- "Annotate the variants in my_sample.vcf"
- "Which variants in this VCF are pathogenic?"
- "Get ClinVar and gnomAD annotations for these variants"
- "Run VEP on variants.vcf and rank by impact"

## Example Output

```
# 🦖 ClawBio VCF Annotator Report

**Input**: demo_variants.vcf
**Date**: 2026-04-19 10:00 UTC
**Total variants**: 5
**HIGH impact**: 3 | **MODERATE**: 2 | **LOW**: 0
**ClinVar Pathogenic/Likely Pathogenic**: 3

## Variant Table

| # | Gene  | Variant             | Consequence       | Impact   | ClinVar    | gnomAD AF |
|---|-------|---------------------|-------------------|----------|------------|-----------|
| 1 | BRCA1 | 17:43044295 G>A     | missense_variant  | HIGH     | Pathogenic | 0.000008  |
| 2 | BRCA2 | 13:32316461 C>T     | stop_gained       | HIGH     | Pathogenic | 0.000004  |
| 3 | CFTR  | 7:117548628 CTTT>C  | frameshift_variant| HIGH     | Pathogenic | 0.021000  |
```

## Output Structure

```
output_directory/
├── report.md                      # Full annotation report
├── results.json                   # All variants as structured JSON
├── tables/
│   └── variants.csv               # Tabular variant data
└── reproducibility/
    ├── commands.sh                # Exact commands to reproduce
    ├── environment.yml            # Python environment
    └── checksums.sha256           # SHA-256 of all output files
```

## Dependencies

**Required**: Python standard library only (`urllib`, `json`, `csv`, `hashlib`)

**Optional**:

- `ensembl-vep` (local install) — for offline annotation without API rate limits
- `cyvcf2` — for faster VCF parsing on large files

## Gotchas

- **Ensembl VEP API rate limit**: Free tier allows ~15 requests/second.
  The skill enforces a 0.1s sleep. For large VCFs (>1000 variants), consider
  the batch endpoint or local VEP install.

- **gnomAD v4 variant ID format**: Must be `{chrom}-{pos}-{ref}-{alt}` without
  `chr` prefix. The skill strips `chr` automatically from VCF CHROM field.

- **ClinVar returns IDs not classifications**: The E-utilities search only
  confirms presence in ClinVar. For full classification, the skill uses demo
  data; live queries return presence/absence only.

- **Indels in VEP**: HGVS notation for indels differs from SNVs. The skill
  handles SNVs fully; complex indels may return limited VEP results.

- **GRCh37 vs GRCh38**: The skill defaults to GRCh38 (hg38). If your VCF
  uses GRCh37 coordinates, VEP results may be incorrect.

## Safety

- **Local-first**: No VCF data is uploaded to third-party servers beyond
  public database APIs (Ensembl, NCBI, gnomAD — all accept variant queries)
- **Disclaimer**: Every report includes the ClawBio research disclaimer
- **Not a diagnostic tool**: ClinVar classifications are research annotations,
  not clinical diagnoses
- **Audit trail**: All operations logged to reproducibility bundle

## Agent Boundary

The agent (LLM) dispatches the VCF and explains results.
The skill (Python) executes all API calls and generates files.
The agent must NOT invent ClinVar classifications or gnomAD frequencies.

## Integration with Bio Orchestrator

**Trigger conditions**: route here when:

- File type is `.vcf`
- Keywords: `annotate`, `variants`, `pathogenic`, `clinvar`, `gnomad`, `vep`

**Chaining partners**:

- `pharmgx-reporter`: VCF annotation can precede pharmacogenomic reporting
- `equity-scorer`: Annotated VCF feeds into population equity analysis
- `lit-synthesizer`: Gene names from annotation can seed literature search

## Maintenance

- **Review cadence**: Monthly — gnomAD and ClinVar update regularly
- **Staleness signals**: gnomAD API endpoint changes; ClinVar reclassifications
- **Deprecation**: Archive if Ensembl VEP REST API is discontinued

## Citations

- [Ensembl VEP](https://www.ensembl.org/vep); variant effect prediction
- [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/); clinical variant classification
- [gnomAD](https://gnomad.broadinstitute.org/); population allele frequencies
- [McLaren et al. 2016](https://doi.org/10.1186/s13059-016-0974-4); VEP paper
