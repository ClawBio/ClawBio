# Proteomics Differential Expression Report

## Summary
- Samples: 6
- Proteins pre-filter: 2383
- Proteins post-filter/imputation: 2069
- Contrast: `treated vs control`
- Significant proteins: 138
  - Upregulated: 45
  - Downregulated: 93

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
| P16522 | 12.509 | 5.345 |
| Q12379 | 8.234 | 2.729 |
| P09798 | 8.031 | 3.217 |
| P42935 | 6.054 | 2.336 |
| P06633 | 5.861 | 1.452 |
| P53886 | 5.663 | 2.165 |
| Q3E747 | 5.629 | 1.775 |
| P14724 | 5.530 | 1.440 |
| P53912;P25608;P54007 | 5.492 | 2.144 |
| P31111 | 5.474 | 1.993 |

### Top 10 Downregulated Proteins (by log2 fold change)

| Protein ID | log2FoldChange | -log10(pvalue) |
|---|---:|---:|
| Q12341 | -9.621 | 4.294 |
| P38210 | -7.432 | 2.567 |
| P39004;P39003 | -7.085 | 1.839 |
| P38238 | -6.099 | 3.672 |
| P53978 | -5.936 | 3.098 |
| Q03774 | -5.816 | 2.230 |
| P38314 | -5.668 | 2.365 |
| P34163 | -5.450 | 1.919 |
| P15992 | -5.427 | 1.694 |
| Q08723 | -5.425 | 3.365 |

## Reproducibility
- Commands: `reproducibility/commands.sh`
- Environment: `reproducibility/environment.yml`
- Checksums: `reproducibility/checksums.sha256`

## Disclaimer
ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions.
