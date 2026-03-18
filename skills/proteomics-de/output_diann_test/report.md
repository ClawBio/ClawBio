# Proteomics Differential Expression Report

## Summary
- Samples: 6
- Proteins pre-filter: 10
- Proteins post-filter/imputation: 10
- Contrast: `treated vs control`
- Significant proteins: 0
  - Upregulated: 0
  - Downregulated: 0

## Preprocessing
- Imputation distribution: `figures/imputation_distribution.png`
- PCA (sample clustering): `figures/pca.png`

## Differential Expression Results
- Full results: `tables/de_results.csv`
- Imputed protein matrix: `tables/imputed_proteinGroups.csv`
- Volcano plot: `figures/volcano.png`

### Top 10 Upregulated Proteins (by log2 fold change)

| Protein ID | log2FoldChange | -log10(pvalue) |
|---|---:|---:|
| No upregulated proteins | | |

### Top 10 Downregulated Proteins (by log2 fold change)

| Protein ID | log2FoldChange | -log10(pvalue) |
|---|---:|---:|
| No downregulated proteins | | |

## Reproducibility
- Commands: `reproducibility/commands.sh`
- Environment: `reproducibility/environment.yml`
- Checksums: `reproducibility/checksums.sha256`

## Disclaimer
ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions.
