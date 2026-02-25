# HEIM Equity Report

**Date**: 2026-02-23 13:32 UTC
**Input**: `demo_populations.vcf`
**Checksum (SHA-256)**: `f847bb32a384508cd062bda6caac0979297dbef3fe14d7e645c2b21f1f9f7c68`
**Samples**: 50
**Populations**: 5
**Variants analysed**: 500

---

## HEIM Equity Score: 76.2/100 (Good)

![HEIM Gauge](figures/heim_gauge.png)

### Score Breakdown

| Component | Value | Weight | Description |
|-----------|-------|--------|-------------|
| Representation Index | 0.720 | 0.35 | Match to global population proportions |
| Heterozygosity Balance | 0.667 | 0.25 | Genetic diversity relative to theoretical max |
| FST Coverage | 1.000 | 0.2 | Fraction of pairwise comparisons computed |
| Geographic Spread | 0.714 | 0.2 | Continental groups represented (out of 7) |

### Key Findings

- **Most represented**: EUR (44.0%, 2.8x global proportion)
- **Least represented**: AMR (10.0%, 0.8x global proportion)
- **Mean observed heterozygosity**: 0.3335 (highest: AFR at 0.3543)

---

## Population Distribution

| Population | Count | Sample % | Global % | Ratio | Obs Het | Exp Het |
|------------|-------|-----------|-----------|-------|---------|---------|
| AFR | 8 | 16.0% | 17.0% | 0.94x | 0.3543 | 0.3338 |
| AMR | 5 | 10.0% | 13.0% | 0.77x | 0.3472 | 0.3064 |
| EAS | 7 | 14.0% | 22.0% | 0.64x | 0.3014 | 0.2788 |
| EUR | 22 | 44.0% | 16.0% | 2.75x | 0.3139 | 0.3075 |
| SAS | 8 | 16.0% | 26.0% | 0.62x | 0.3508 | 0.3220 |

![Ancestry Distribution](figures/ancestry_bar.png)

## Heterozygosity

![Heterozygosity](figures/heterozygosity.png)

## Pairwise FST

| Comparison | Hudson FST |
|------------|-----------|
| AFR vs AMR | 0.0773 |
| AFR vs EAS | 0.1011 |
| AFR vs EUR | 0.0590 |
| AFR vs SAS | 0.0619 |
| AMR vs EAS | 0.0898 |
| AMR vs EUR | 0.0425 |
| AMR vs SAS | 0.0692 |
| EAS vs EUR | 0.0439 |
| EAS vs SAS | 0.0782 |
| EUR vs SAS | 0.0467 |

![FST Heatmap](figures/fst_heatmap.png)


## Principal Component Analysis

![PCA](figures/pca_plot.png)

- PC1 explains 7.6% of variance
- PC2 explains 4.9% of variance
- Top 5 components explain 23.0% of total variance


---

## Methods

- **Tool**: ClawBio Equity Scorer v0.1.0
- **HEIM framework**: Health Equity Index for Minorities (Corpas, 2026)
- **Heterozygosity**: Observed = proportion of heterozygous genotypes per site, averaged across 500 variants. Expected = 2pq from population allele frequencies.
- **FST**: Nei's GST (HT-HS)/HT, ratio of averages across sites. Values floored at 0.
- **PCA**: scikit-learn PCA on mean-imputed genotype matrix (0/1/2 encoding).
- **Global reference**: Approximate continental proportions from the 1000 Genomes Project.

## Reproducibility

```bash
# Re-run this analysis
python equity_scorer.py --input demo_populations.vcf --output demo_report
```

**Input checksum**: `f847bb32a384508cd062bda6caac0979297dbef3fe14d7e645c2b21f1f9f7c68`

## References

- Corpas, M. (2026). ClawBio. https://github.com/manuelcorpas/ClawBio
- Hudson, R.R., Slatkin, M. & Maddison, W.P. (1992). Estimation of levels of gene flow from DNA sequence data. Genetics, 132(2), 583-589.
- The 1000 Genomes Project Consortium (2015). A global reference for human genetic variation. Nature, 526, 68-74.
