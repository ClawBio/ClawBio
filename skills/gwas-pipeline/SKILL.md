---
name: gwas-pipeline
description: End-to-end GWAS automation wrapping PLINK2 for genotype QC and REGENIE for two-step whole-genome regression association
  testing. Produces Manhattan plots, QQ plots, clumped lead variants, and structured summary statistics.
license: MIT
metadata:
  version: 0.1.0
  author: Reza
  tags:
  - gwas
  - association-testing
  - regenie
  - plink2
  - population-genetics
  - biobank
  - quality-control
  openclaw:
    requires:
      bins:
      - python3
      - plink2
      - regenie
    always: false
    emoji: 📊
    homepage: https://github.com/ClawBio/ClawBio
    os:
    - darwin
    - linux
    install:
    - kind: conda
      package: plink2
      bins:
      - plink2
    - kind: conda
      package: regenie
      bins:
      - regenie
    trigger_keywords:
    - GWAS
    - genome-wide association
    - PLINK
    - REGENIE
    - association testing
    - Manhattan plot
    - QQ plot
    - polygenic
    - case-control
---

# 📊 GWAS Pipeline

You are **GWAS Pipeline**, a specialised ClawBio agent for genome-wide association studies. Your role is to automate best-practice QC and association testing from genotype files to publication-ready results.

## Why This Exists

- **Without it**: Researchers must orchestrate PLINK2 and REGENIE manually, writing hundreds of lines of bash, managing dozens of parameters, and applying field-standard QC thresholds by hand
- **With it**: A single command runs the full QC cascade, REGENIE two-step regression, and post-GWAS visualisation on any genotype dataset
- **Why ClawBio**: Grounded in Anderson et al. (2010) QC thresholds and Mbatchou et al. (2021) REGENIE methodology — not ad hoc parameter choices. Every command logged for reproducibility

## Core Capabilities

1. **Genotype QC via PLINK2**: Sample/variant missingness, MAF, HWE, LD pruning
2. **REGENIE Step 1**: Whole-genome ridge regression with LOCO predictions
3. **REGENIE Step 2**: Single-variant association (Firth logistic / linear)
4. **Visualisation**: Manhattan plot, QQ plot with lambda GC
5. **Post-GWAS**: Lead variant extraction at genome-wide significance (P < 5e-8)
6. **Reproducibility**: Full command logging, parameter tracking, software versions

## Input Formats

| Format | Extension | Required Fields | Example |
|--------|-----------|-----------------|---------|
| PLINK binary | `.bed` + `.bim` + `.fam` | Standard PLINK format | `example.bed` |
| BGEN | `.bgen` | BGEN v1.2+ with sample info | `example.bgen` |
| Phenotype | `.txt` | FID, IID, trait column(s) | `phenotype_bin.txt` |
| Covariate | `.txt` | FID, IID, covariate columns | `covariates.txt` |

## Workflow

1. **Validate**: Check input files exist, detect format, verify binaries on PATH
2. **QC** (PLINK2): Variant missingness, sample missingness, MAF, HWE filtering; LD pruning for Step 1
3. **Step 1** (REGENIE): Whole-genome ridge regression on LD-pruned genotyped variants with LOCO
4. **Step 2** (REGENIE): Single-variant association with Firth correction (binary) or linear regression (quantitative)
5. **Post-GWAS**: Parse results, compute lambda GC, extract lead variants, generate plots
6. **Report**: Write report.md, result.json, summary statistics TSV, and reproducibility bundle

## CLI Reference

```bash
# Demo mode (REGENIE example data, binary trait Y1)
python skills/gwas-pipeline/gwas_pipeline.py --demo --output /tmp/gwas_demo

# Real data
python skills/gwas-pipeline/gwas_pipeline.py \
  --bed /path/to/data --pheno pheno.txt --covar covar.txt \
  --trait-type bt --trait Y1 --output results/

# Via ClawBio runner
python clawbio.py run gwas-pipe --demo
```

## Demo

```bash
python clawbio.py run gwas-pipe --demo
```

Expected output: A full GWAS report on REGENIE's official 500-sample, 1000-variant example dataset with binary trait Y1, including QC summary, REGENIE Step 1/2 output, Manhattan plot, QQ plot with lambda GC, and reproducibility bundle.

## Dependencies

**Required** (external binaries):
- `plink2` >= 2.0 — genotype QC and LD operations
- `regenie` >= 3.0 — two-step whole-genome regression

Install via conda: `CONDA_SUBDIR=osx-64 conda create -n clawbio-gwas -c conda-forge -c bioconda plink2 regenie`

**Python** (standard library + matplotlib):
- `matplotlib` >= 3.7 — Manhattan and QQ plots
- `numpy` >= 1.24 — QQ plot expected quantiles

## Safety

- **Local-first**: All computation runs locally via PLINK2/REGENIE subprocesses
- **Disclaimer**: Every report includes the ClawBio medical disclaimer
- **Audit trail**: Every PLINK2/REGENIE command logged to `reproducibility/commands.sh`
- **No hallucinated science**: All QC thresholds trace to Anderson et al. 2010 / REGENIE documentation

## Integration with Bio Orchestrator

**Trigger conditions** — the orchestrator routes here when:
- User mentions GWAS, association testing, Manhattan plot, or case-control study
- User provides genotype files (BED/BIM/FAM, BGEN, VCF) with a phenotype file

**Chaining partners**:
- `gwas-lookup`: Downstream — look up lead variants across federated databases
- `gwas-prs`: Downstream — compute polygenic risk scores from summary statistics
- `variant-annotation`: Downstream — annotate lead variants with VEP/ClinVar

## Citations

- [Mbatchou et al. (2021)](https://pubmed.ncbi.nlm.nih.gov/34017140/) — REGENIE: computationally efficient whole-genome regression. *Nature Genetics* 53:1097–1103
- [Chang et al. (2015)](https://pubmed.ncbi.nlm.nih.gov/25722852/) — Second-generation PLINK. *GigaScience* 4:7
- [Anderson et al. (2010)](https://pubmed.ncbi.nlm.nih.gov/21085122/) — Data quality control in genetic case-control association studies. *Nature Protocols* 5:1564–1573
