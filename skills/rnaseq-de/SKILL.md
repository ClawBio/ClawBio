---
name: rnaseq-de
description: Differential expression analysis for bulk RNA-seq and pseudo-bulk count matrices with QC, PCA, and contrast testing.
version: 0.1.0
tags: [rna-seq, differential expression, bulk, pseudo-bulk, transcriptomics, DESeq2, PyDESeq2, QC, PCA]
trigger_keywords: [rna-seq, differential expression, bulk RNA, pseudo-bulk, volcano plot, ma plot, count matrix, DESeq2, pydeseq2]
metadata:
  openclaw:
    requires:
      bins:
        - python3
      env: []
      config: []
    always: false
    emoji: "🧬"
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
      - kind: uv
        package: scikit-learn
        bins: []
---

# 🧬 RNA-seq Differential Expression

This skill performs differential expression on bulk RNA-seq or pseudo-bulk count matrices.

## Core Capabilities

1. Input validation for count matrix and sample metadata
2. Pre-DE QC (library size, detected genes, low-count filtering)
3. PCA visualisation on normalized expression
4. Differential expression from formula + contrast
5. Volcano and MA plots
6. Markdown report with reproducibility files

## Input Contract

- Count matrix (`.csv` or `.tsv`): rows are genes, columns are samples, first column is gene identifier
- Metadata table (`.csv` or `.tsv`): one row per sample, must include `sample_id`
- Formula: e.g. `~ condition` or `~ batch + condition`
- Contrast: `factor,numerator,denominator` (e.g. `condition,treated,control`)

## Output Structure

```
rnaseq_de_report/
├── report.md
├── figures/
│   ├── pca.png
│   ├── volcano.png
│   └── ma_plot.png
├── tables/
│   ├── qc_summary.csv
│   ├── normalized_counts.csv
│   └── de_results.csv
└── reproducibility/
    ├── commands.sh
    ├── environment.yml
    └── checksums.sha256
```

## Usage

```bash
python rnaseq_de.py \
  --counts counts.csv \
  --metadata metadata.csv \
  --formula "~ batch + condition" \
  --contrast "condition,treated,control" \
  --output report_dir
```

## Safety

- Local-only processing
- Warn before overwriting existing output
- Report-level disclaimer required
