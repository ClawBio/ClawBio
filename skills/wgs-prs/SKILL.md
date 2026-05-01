---
name: wgs-prs
description: >-
  End-to-end WGS to polygenic risk score pipeline. Takes paired-end FASTQ files
  (or a pre-existing VCF) through nf-core/sarek for variant calling, applies VCF
  QC (normalisation, hard filtering, Ti/Tv and Het/Hom checks), then computes
  polygenic risk scores via the PGS Catalog. Fills the FASTQ→VCF gap upstream
  of the gwas-prs skill.
version: 0.1.0
author: David de Lorenzo
license: MIT
tags:
  - wgs
  - whole-genome-sequencing
  - polygenic-risk-scores
  - prs
  - sarek
  - nf-core
  - nextflow
  - variant-calling
  - vcf-qc
  - gatk
metadata:
  openclaw:
    requires:
      bins:
        - python3
        - nextflow
      anyBins:
        - docker
        - singularity
      env: []
      config: []
    always: false
    emoji: "🧬"
    homepage: https://github.com/ClawBio/ClawBio
    os: [darwin, linux]
    install:
      - kind: url
        url: https://get.nextflow.io
        bins: [nextflow]
        note: "curl -s https://get.nextflow.io | bash"
      - kind: conda
        package: bcftools
        bins: [bcftools]
        note: "conda install -c bioconda bcftools  (optional — enables full VCF QC)"
    trigger_keywords:
      - WGS
      - whole genome sequencing
      - FASTQ to PRS
      - variant calling
      - nf-core sarek
      - sarek
      - FASTQ VCF
      - WGS polygenic risk
      - germline variant calling
      - GATK HaplotypeCaller
      - WGS pipeline
      - raw sequencing to risk scores
inputs:
  - name: fastq_r1
    type: file
    format: fastq.gz
    description: "Forward reads FASTQ.gz (paired-end WGS)"
    required: false
  - name: fastq_r2
    type: file
    format: fastq.gz
    description: "Reverse reads FASTQ.gz"
    required: false
  - name: input_vcf
    type: file
    format: vcf.gz
    description: "Pre-existing VCF — skips sarek, starts at QC stage"
    required: false
  - name: sample_id
    type: string
    description: "Sample identifier used throughout the pipeline"
    required: false
    default: "SAMPLE"
  - name: sex
    type: string
    description: "Biological sex: XX or XY (affects sex-chromosome calling)"
    required: false
    default: "XX"
outputs:
  - name: bridge_report.md
    description: "Human-readable summary of all pipeline stages"
  - name: bridge_report.json
    description: "Machine-readable stage status and QC metrics"
  - name: vcf_qc/qc_metrics.json
    description: "Ti/Tv ratio, Het/Hom ratio, variant counts, pass/fail"
  - name: vcf_qc/canonical_pass.vcf.gz
    description: "Normalised, filtered canonical VCF ready for PRS scoring"
  - name: prs_output/report.md
    description: "PRS narrative report from gwas-prs"
  - name: prs_output/tables/scores.csv
    description: "Per-trait PRS scores, percentiles, and risk categories"
---

# 🧬 WGS-PRS Pipeline

**Author**: David de Lorenzo (ClawBio Community)
**Requires**: Python 3.9+, nextflow, docker or singularity, bcftools (recommended)

---

You are the **WGS-PRS** skill, an end-to-end pipeline agent for whole-genome sequencing data. Your role is to take a user from raw FASTQ files (or a pre-existing VCF) all the way to polygenic risk scores, with robust QC at every stage.

## Pipeline Stages

1. **Variant calling** — nf-core/sarek (FASTQ → BAM → VCF via GATK HaplotypeCaller)
2. **VCF QC** — bcftools normalisation, hard filtering, Ti/Tv and Het/Hom evaluation
3. **PRS scoring** — ClawBio `gwas-prs` skill (PGS Catalog, 6 curated + 3,000+ live scores)
4. **Aggregated report** — Markdown + JSON summary of all stages

## Entry Points

Users may enter the pipeline at two points:

- **FASTQ entry** (full pipeline): provide `--fastq-r1` and optionally `--fastq-r2`
- **VCF entry** (skip sarek): provide `--input-vcf` with a pre-existing single-sample GRCh38 VCF

## Usage

```bash
# Full pipeline from paired FASTQ
python wgs_prs.py --fastq-r1 sample_R1.fastq.gz --fastq-r2 sample_R2.fastq.gz \
    --sample-id HG001 --output-dir results/

# Start from an existing VCF
python wgs_prs.py --input-vcf sample.vcf.gz --output-dir results/

# Dry run — generate samplesheet and preview commands only
python wgs_prs.py --fastq-r1 sample_R1.fastq.gz --dry-run

# Score a specific trait
python wgs_prs.py --input-vcf sample.vcf.gz --trait "type 2 diabetes"
```

## Key Design Decisions

- **Reference genome**: GRCh38 (GATK.GRCh38 sarek alias). Older GRCh37 VCFs require liftover before PRS scoring.
- **Variant caller**: GATK HaplotypeCaller (default). DeepVariant available via `--tools deepvariant`.
- **VCF QC thresholds**: Ti/Tv 1.8–2.5, Het/Hom 1.0–3.0, QUAL ≥ 30, DP ≥ 10.
- **Fail-fast**: pipeline aborts on QC failure by default. Use `--no-fail-fast` to continue with a warning.
- **Canonical VCF contract**: the handoff point between stages is a normalised, PASS-filtered, single-sample GRCh38 VCF — consistent with what `gwas-prs`, `variant-annotation`, and `pharmgx-reporter` all accept.

## Chaining with other ClawBio Skills

After WGS-PRS completes, the canonical VCF can be passed to:
- `variant-annotation` — Ensembl VEP, ClinVar, gnomAD
- `pharmgx-reporter` — pharmacogenomics from the same VCF
- `claw-ancestry-pca` — ancestry estimation to validate PRS reference population
- `clinical-variant-reporter` — ACMG/AMP pathogenicity classification

## Dependencies

| Tool | Required | Purpose |
|------|----------|---------|
| nextflow | Yes | Executes nf-core/sarek |
| docker or singularity | Yes | Container runtime for sarek |
| bcftools ≥ 1.17 | Recommended | VCF normalisation and stats (falls back to Python if absent) |
| python3 ≥ 3.9 | Yes | Runtime |

## Integration with Bio Orchestrator

This skill is invoked when:
- The user mentions WGS, whole-genome sequencing, FASTQ files, or raw sequencing data
- The user asks to run the full pipeline "from scratch" or "from reads"
- Keywords: WGS, FASTQ, sarek, variant calling, germline variants, raw reads to PRS

It chains downstream to `gwas-prs` automatically. For users who already have a VCF,
the bio-orchestrator should route directly to `gwas-prs` or `variant-annotation` instead.
