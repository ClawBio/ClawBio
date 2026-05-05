---
name: wgs-prs
description: >-
  End-to-end WGS to polygenic risk score pipeline. Takes paired-end FASTQ files
  (or a pre-existing VCF) through nf-core/sarek for variant calling, applies VCF
  QC (normalisation, hard filtering, Ti/Tv and Het/Hom checks), then computes
  polygenic risk scores via the PGS Catalog. Fills the FASTQ to VCF gap upstream
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
        note: "conda install -c bioconda bcftools  (optional, enables full VCF QC)"
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
    description: "Pre-existing VCF, skips sarek, starts at QC stage"
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

## Trigger

**Fire this skill when the user says any of:**
- "run WGS analysis"
- "whole genome sequencing"
- "FASTQ to PRS" / "FASTQ to polygenic risk scores"
- "variant calling from raw reads"
- "nf-core sarek" / "run sarek"
- "germline variant calling"
- "GATK HaplotypeCaller"
- "raw sequencing to risk scores"
- "WGS pipeline" / "WGS polygenic risk"

**Do NOT fire when:**
- The user already has a VCF and wants PRS only: route to `gwas-prs` instead.
- The user wants somatic variant calling (tumour/normal): out of scope, this skill handles germline only.
- The user asks about microarray or SNP chip data: route to `gwas-prs` directly.
- The user wants metagenomics or RNA-seq: wrong pipeline.

## Scope

**One skill, one task.** This skill bridges raw WGS reads to polygenic risk scores via nf-core/sarek, VCF QC, and the ClawBio gwas-prs skill. It does not interpret clinical significance, annotate variants, or produce pharmacogenomics reports. Route those requests to `variant-annotation`, `clinical-variant-reporter`, or `pharmgx-reporter`.

## Pipeline Stages

1. **Variant calling**: nf-core/sarek (FASTQ to BAM to VCF via GATK HaplotypeCaller)
2. **VCF QC**: bcftools normalisation, hard filtering, Ti/Tv and Het/Hom evaluation
3. **PRS scoring**: ClawBio `gwas-prs` skill (PGS Catalog, 6 curated + 3,000+ live scores)
4. **Aggregated report**: Markdown + JSON summary of all stages

## Entry Points

Users may enter the pipeline at two points:

- **FASTQ entry** (full pipeline): provide `--fastq-r1` and optionally `--fastq-r2`
- **VCF entry** (skip sarek): provide `--input-vcf` with a pre-existing single-sample GRCh38 VCF

## Workflow

When the user provides WGS input (FASTQ or VCF):

1. **Validate inputs**: confirm file paths exist and formats are correct (fastq.gz or vcf.gz). Abort with a clear message if required inputs are missing.
2. **Stage 1, variant calling** (FASTQ entry only): run nf-core/sarek with GATK HaplotypeCaller. Generate samplesheet CSV, invoke nextflow, confirm VCF output exists.
3. **Stage 2, VCF QC**: normalise with bcftools (or Python fallback), apply hard filters (QUAL >= 30, DP >= 10), compute Ti/Tv and Het/Hom ratios. Fail fast if thresholds are violated, unless `--no-fail-fast` is set.
4. **Stage 3, PRS scoring**: pass the canonical VCF to `gwas-prs`. Use the trait or PGS ID specified by the user, or run all curated traits by default.
5. **Stage 4, aggregated report**: write `bridge_report.md` and `bridge_report.json` combining stage statuses, QC metrics, and PRS summary.
6. **Surface results**: show the user the report path and key metrics. Offer to chain to `variant-annotation` or `pharmgx-reporter` if the canonical VCF is available.

**Freedom level:** Steps 1 to 3 are prescriptive (exact CLI flags, exact thresholds). Steps 5 to 6 allow interpretive flexibility in the report narrative.

## Usage

```bash
# Full pipeline from paired FASTQ
python wgs_prs.py --fastq-r1 sample_R1.fastq.gz --fastq-r2 sample_R2.fastq.gz \
    --sample-id HG001 --output-dir results/

# Start from an existing VCF
python wgs_prs.py --input-vcf sample.vcf.gz --output-dir results/

# Dry run: generate samplesheet and preview commands only
python wgs_prs.py --fastq-r1 sample_R1.fastq.gz --dry-run

# Score a specific trait
python wgs_prs.py --input-vcf sample.vcf.gz --trait "type 2 diabetes"
```

## Key Design Decisions

- **Reference genome**: GRCh38 (GATK.GRCh38 sarek alias). Older GRCh37 VCFs require liftover before PRS scoring.
- **Variant caller**: GATK HaplotypeCaller (default). DeepVariant available via `--tools deepvariant`.
- **VCF QC thresholds**: Ti/Tv 1.8 to 2.5, Het/Hom 1.0 to 3.0, QUAL >= 30, DP >= 10.
- **Fail-fast**: pipeline aborts on QC failure by default. Use `--no-fail-fast` to continue with a warning.
- **Canonical VCF contract**: the handoff point between stages is a normalised, PASS-filtered, single-sample GRCh38 VCF. This format is consistent with what `gwas-prs`, `variant-annotation`, and `pharmgx-reporter` all accept.

## Example Output

```markdown
# ClawBio WGS-PRS Bridge Report

**Sample:** HG001
**Generated:** 2026-05-01T12:00:00+00:00
**Output directory:** `results/`

## Pipeline Stages

| Stage  | Status     | Duration |
|--------|------------|----------|
| sarek  | success    | 142.3s   |
| vcf_qc | success    | 8.1s     |
| gwas   | success    | 23.5s    |
| report | success    | 0.4s     |

## VCF QC Metrics

**QC Status:** PASS

| Metric            | Value   |
|-------------------|---------|
| Total variants    | 4,821   |
| SNPs              | 4,103   |
| Indels            | 718     |
| Ti/Tv ratio       | 2.12    |
| Het/Hom ratio     | 1.74    |
| Filtered variants | 203     |

## Polygenic Risk Scores

| Trait              | Score  | Percentile | Risk Category |
|--------------------|--------|------------|---------------|
| Type 2 diabetes    | 0.82   | 73rd       | Above average |
| Coronary artery    | 0.61   | 54th       | Average       |

*ClawBio is a research and educational tool. It is not a medical device.*
```

## Chaining with other ClawBio Skills

After WGS-PRS completes, the canonical VCF can be passed to:
- `variant-annotation`: Ensembl VEP, ClinVar, gnomAD
- `pharmgx-reporter`: pharmacogenomics from the same VCF
- `claw-ancestry-pca`: ancestry estimation to validate PRS reference population
- `clinical-variant-reporter`: ACMG/AMP pathogenicity classification

## Dependencies

| Tool | Required | Purpose |
|------|----------|---------|
| nextflow | Yes | Executes nf-core/sarek |
| docker or singularity | Yes | Container runtime for sarek |
| bcftools >= 1.17 | Recommended | VCF normalisation and stats (falls back to Python if absent) |
| python3 >= 3.9 | Yes | Runtime |

## Gotchas

- **Do not skip VCF QC even when the user provides their own VCF.** Users often pass unfiltered or unnormalised VCFs from external pipelines. Always run Stage 2 unless the user explicitly opts out with `--skip-qc`. Skipping QC silently produces unreliable PRS scores.
- **Do not route to this skill when the user already has a VCF and wants PRS only.** The model will be tempted to use wgs-prs because it mentions PRS. If there is no FASTQ and the user has not asked for variant calling, route directly to `gwas-prs` to avoid unnecessary sarek overhead.
- **Do not invent QC thresholds.** Ti/Tv and Het/Hom cut-offs are fixed at 1.8 to 2.5 and 1.0 to 3.0 respectively. Do not adjust these based on the user's wishes or apparent sample quality. If thresholds are debated, surface the metrics and let the user decide whether to continue with `--no-fail-fast`.
- **GRCh37 VCFs will silently produce wrong PRS scores.** The PGS Catalog scores are aligned to GRCh38. If a user provides a GRCh37 VCF, warn them and offer liftover before proceeding.
- **Sarek can take hours on a full genome.** Set expectations with the user before launching Stage 1. For testing, recommend `--dry-run` first.

## Safety

- **Local-first**: all data is processed locally. No reads or variants leave the user's machine.
- **Disclaimer**: every report includes the ClawBio medical disclaimer: *"ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions."*
- **No hallucinated parameters**: all QC thresholds and PGS Catalog identifiers trace to documented sources.
- **Audit trail**: stage durations, commands, and output paths are logged to `bridge_report.json`.

## Agent Boundary

The agent (LLM) dispatches, explains results, and surfaces next steps. The skill (Python) executes all variant calling, QC, and scoring. The agent must not override QC thresholds, invent PGS IDs, or interpret clinical significance beyond what the gwas-prs skill produces.

## Integration with Bio Orchestrator

This skill is invoked when:
- The user mentions WGS, whole-genome sequencing, FASTQ files, or raw sequencing data
- The user asks to run the full pipeline "from scratch" or "from reads"
- Keywords: WGS, FASTQ, sarek, variant calling, germline variants, raw reads to PRS

It chains downstream to `gwas-prs` automatically. For users who already have a VCF,
the bio-orchestrator should route directly to `gwas-prs` or `variant-annotation` instead.
