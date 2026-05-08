---
name: rare-disease-rnaseq
description: Blood RNA-seq expression-outlier detection for rare-disease diagnostics. Cases scored against a control reference panel; outliers ranked and filtered by a haploinsufficient disease-gene panel.
version: 0.1.0
tags: [rna-seq, rare-disease, outlier-detection, OUTRIDER, FRASER, diagnostic, blood, transcriptomics, haploinsufficiency]
trigger_keywords: [rare disease rnaseq, expression outlier, OUTRIDER, FRASER, blood rna-seq diagnostic, NGRL, undiagnosed, candidate diagnosis]
metadata:
  openclaw:
    requires:
      bins:
        - python3
      env: []
      config: []
    always: false
    emoji: "🩸"
    homepage: https://github.com/ClawBio/ClawBio
    os: [darwin, linux]
    install:
      - kind: uv
        package: pandas
        bins: []
      - kind: uv
        package: numpy
        bins: []
      - kind: uv
        package: matplotlib
        bins: []
---

# 🩸 Rare-Disease Blood RNA-seq Outlier Detection

Reproduces the diagnostic principle of the Genomics England NGRL paper (Blood-based RNA-Seq of 5,412 individuals, medRxiv 2026.03.19.26348811). For each case sample, scores per-gene expression against a control reference panel and flags candidates falling in a curated dosage-sensitive disease-gene panel.

## When To Use

- A WGS-negative or WGS-VUS rare-disease patient with a paired blood RNA-seq sample
- A clinical bioinformatician triaging candidate diagnoses before MDT review
- A population-biobank team building an ancestry-matched control reference for outlier calling (e.g. Qatar Biobank for Sidra paediatric cases)

## Method

Per-gene robust outlier scoring on log2(CPM+1):

1. Library-size normalise (CPM), log-transform
2. For each gene: compute median and MAD across the control panel
3. For each case-gene cell: modified z = 0.6745 (x − median) / MAD
4. Flag |z| ≥ threshold (default 3) and gene in disease panel
5. Rank by |z|, separate down-outliers (haploinsufficiency-consistent) from up-outliers

This implements the **diagnostic principle** of OUTRIDER (per-gene outlier vs control panel) without the autoencoder, so it runs in seconds with no R/Bioconductor stack. For clinical-grade calls swap to the full DROP pipeline (gagneurlab/drop) which adds OUTRIDER's denoising autoencoder, FRASER2 splicing outliers, and confounder correction. The skill's I/O contract is the same so the upgrade is drop-in.

## Input Contract

- Counts matrix (`.csv` or `.tsv`): rows = genes (HGNC symbol), columns = sample IDs
- Cases file (`.txt`): one case sample ID per line
- Controls file (`.txt`): one control sample ID per line (typically n ≥ 50)
- Disease panel (optional, `.csv` with `gene` and `mechanism` columns): defaults to a built-in 50-gene haploinsufficient panel

## Output Structure

```
rdoutlier_report/
├── report.md                     # per-case candidate diagnoses + clinical narrative
├── result.json                   # standard ClawBio envelope
├── figures/
│   └── case_outlier_heatmap.png  # z-scores across cases × top genes
├── tables/
│   ├── outlier_calls.csv         # all flagged outliers with z-score, direction, mechanism
│   └── per_gene_stats.csv        # control median + MAD per gene
└── reproducibility/
    ├── commands.sh
    ├── environment.yml
    └── checksums.sha256
```

## Demo

```bash
python clawbio.py run rdoutlier --demo
```

Generates 100 synthetic Gulf-ancestry control samples + 2 cases with injected outliers (FBN1 down, NF1 up) across a 200-gene panel. Demonstrates the diagnostic loop end-to-end in seconds.

## Production Path (Sidra / QBB Reference)

| Component             | Demo                          | Production                                    |
|-----------------------|-------------------------------|-----------------------------------------------|
| Aligner + quantifier  | none (synthetic counts)       | STAR + featureCounts (or Salmon)              |
| Outlier algorithm     | robust per-gene z-score       | OUTRIDER autoencoder + FRASER2 splicing       |
| Control panel         | 100 synthetic samples         | QBB n≈12K PAXgene blood RNA-seq               |
| Confounder correction | none                          | DROP pipeline (RIN, batch, hidden factors)    |
| Disease panel         | 50 haploinsufficient genes    | ClinGen haploinsufficient + PanelApp          |
| Return-of-result loop | report.md                     | Sidra MDT reflex from WGS-negative referrals  |

## Safety

- Local-only processing, no network calls in core pipeline
- Compatible with secure research environments (Genomics England RE pattern; Sidra clinical genomics environment)
- Disclaimer required on every report

## Disclaimer

ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions.
