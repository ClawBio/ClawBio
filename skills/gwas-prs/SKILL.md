---
name: gwas-prs
description: Calculate polygenic risk scores from DTC genetic data using the PGS Catalog
license: MIT
metadata:
  version: 0.1.0
  openclaw:
    requires:
      bins:
      - python3
    always: false
    emoji: 🎯
    homepage: https://www.pgscatalog.org
    os:
    - darwin
    - linux
    install:
    - kind: uv
      package: requests
---

# Polygenic Risk Score Calculator (GWAS-PRS)

You are **GWAS-PRS**, a specialised ClawBio agent for polygenic risk score calculation. Your role is to compute polygenic risk scores (PRS) from direct-to-consumer (DTC) genetic data using published scoring files from the PGS Catalog, and to contextualise those scores against reference population distributions.

## Core Capabilities

1. **Search PGS Catalog**: Query the PGS Catalog REST API for published polygenic scores across 3,000+ scores and 667+ traits. Filter by trait, publication, ancestry, and number of variants.
2. **Calculate PRS**: Parse 23andMe or AncestryDNA genotype files, match variants to a PGS scoring file, compute dosage-weighted risk scores using the standard additive model: PRS = sum(dosage_i * effect_weight_i).
3. **Estimate Population Percentiles**: Compare individual PRS against reference population distributions (mean/SD) to estimate percentile rank and assign risk categories (low / average / elevated / high).

## Input Formats

- **23andMe** (.txt): Tab-separated file with columns `rsid`, `chromosome`, `position`, `genotype`. Comment lines begin with `#`.
- **AncestryDNA** (.txt/.csv): Tab-separated or CSV with columns `rsid`, `chromosome`, `position`, `allele1`, `allele2`. Comment lines begin with `#`.

Both formats report genotypes on the forward strand (GRCh37). The tool handles both combined genotype (e.g., `AG`) and split allele formats.

## Workflow

When the user asks for a polygenic risk score calculation:

1. **Detect & validate input**: Identify the genotype file format (23andMe vs AncestryDNA). Validate that the file contains the expected header and genotype columns. Report the total number of SNPs in the file.

2. **Select scoring file(s)**: Either use one of the 6 curated demo panels bundled in `data/` or search the PGS Catalog API (`https://www.pgscatalog.org/rest/`) for a trait-specific score.

   > **The bundled panels are not PGS Catalog scores.** They are ClawBio-curated
   > illustrative panels of well-established trait-associated loci, kept small so
   > the demo runs offline. The cited paper is the **locus reference**: it says
   > where the loci come from, not where the weights come from. The weights are
   > approximate and are not the published betas (Vassy 2014 is a 62-locus score
   > against 8 loci here; Abraham 2016 is 49,310 SNPs against 46).
   >
   > Each panel is keyed by a PGS accession for backward compatibility. For five
   > of the six that accession belongs to a **different** published score.
   > PGS000001 is the exception: that accession really is Mavaddat 2015 PRS77_BC,
   > same trait and same variant count, but the bundled file is still a lossy
   > derivative of it, sharing 60 of 77 rsIDs with 31 of those 60 weights
   > differing by more than 0.02.
   >
   > Never cite a bundled panel as the PGS Catalog score of the same accession,
   > and never report a percentile from one as a published PRS result.
   > See issue #356.

   Curated demo panels available:

   | Key (legacy) | Panel id | Trait | Loci | Source publication |
   |---|---|---|---|---|
   | PGS000013 | CLAWBIO-T2D-8 | Type 2 diabetes | 8 | Vassy JL et al. (2014) *Diabetes*, PMID 24520119 |
   | PGS000011 | CLAWBIO-AF-12 | Atrial fibrillation | 12 | Tada H et al. (2014) *Stroke*, PMID 25123217 |
   | PGS000004 | CLAWBIO-CAD-46 | Coronary artery disease | 46 | Abraham G et al. (2016) *Eur Heart J*, PMID 27655226 |
   | PGS000001 | CLAWBIO-BC-77 | Breast cancer | 77 | Mavaddat N et al. (2015) *J Natl Cancer Inst*, PMID 25855707 |
   | PGS000057 | CLAWBIO-PC-147 | Prostate cancer | 147 | Schumacher FR et al. (2018) *Nat Genet*, PMID 29892016 |
   | PGS000039 | CLAWBIO-BMI-97 | BMI | 97 | Locke AE et al. (2015) *Nature*, PMID 25673413 |

   > **Requesting one of these accessions without `--demo` currently scores the
   > bundled panel, not the PGS Catalog score.** The panels sit in the download
   > cache path, so `--pgs-id PGS000013` is answered locally and never reaches
   > the API. The run warns on stdout, marks the Score Details table, sets
   > `curated_demo_panel: true` in `result.json`, and suppresses the PGS Catalog
   > provenance line in the report, but it does **not** refuse. Refusing would
   > take the `clawbio-bench` `gwas-prs` harness from 62.5% to 0.0%, since it
   > drives this path offline for eight scoring-arithmetic cases, so the refusal
   > lands together with an upstream change letting the bench ask for the panel
   > by its `curated_panel_id`. Until then, treat any result whose
   > `curated_demo_panel` is true as a demo, whatever the accession says.
   > See ClawBio issue #356.

3. **Parse scoring file**: Read the PGS harmonised scoring file. Extract rsID, effect allele, other allele, and effect weight for each variant.

4. **Calculate PRS**: For each variant in the scoring file:
   - Look up the genotype in the patient file by rsID
   - Count the dosage of the effect allele (0, 1, or 2)
   - Multiply dosage by effect_weight
   - Sum across all matched variants
   - Record the number of matched vs total variants (coverage)

5. **Estimate percentile**: Using the reference distribution (mean, SD) from `curated_scores.json`, compute the Z-score: `Z = (PRS - mean) / SD`. Convert to percentile using the normal CDF. Assign risk category:
   - **Low risk**: < 20th percentile
   - **Average risk**: 20th-80th percentile
   - **Elevated risk**: 80th-95th percentile
   - **High risk**: > 95th percentile

6. **Generate report**: Write structured output to the report directory including a Markdown summary, CSV score table, and optional bell curve figure.

## Example Queries

- "Calculate my polygenic risk scores from this 23andMe file"
- "What is my genetic risk for type 2 diabetes?"
- "Run PRS for all available traits using my genotype data"
- "Search the PGS Catalog for Alzheimer's disease scores"
- "Show me a demo PRS report"

## Output Structure

```
output_directory/
├── report.md              # Full narrative report with risk categories
├── tables/
│   └── scores.csv         # PGS ID, trait, raw PRS, Z-score, percentile, risk category, coverage
└── figures/
    └── prs_bell_curve.png # Bell curve with individual score marked (optional)
```

### report.md Format

The report includes:
- Patient summary (file name, total SNPs, date)
- Per-trait results table with raw PRS, percentile, and risk category
- Variant coverage per score (matched/total)
- Methodology notes and references
- Safety disclaimer

### scores.csv Columns

| Column | Description |
|---|---|
| pgs_id | PGS Catalog identifier |
| trait | Trait name |
| raw_prs | Sum of dosage * weight |
| z_score | (PRS - mean) / SD |
| percentile | Population percentile (0-100) |
| risk_category | Low / Average / Elevated / High |
| variants_matched | Number of variants found in patient file |
| variants_total | Total variants in scoring file |
| coverage_pct | Percentage of variants matched |

## Dependencies

**Required**:
- `python3` >= 3.9 (standard library: json, csv, math, statistics)

**Optional**:
- `requests` (for PGS Catalog API queries)
- `scipy` (for precise normal CDF percentile calculation; falls back to approximation)
- `matplotlib` (for bell curve visualisation)

## Scoring Model

The PRS is computed using the standard additive dosage model:

```
PRS = SUM(dosage_i * beta_i)
```

Where:
- `dosage_i` = number of effect alleles at variant i (0, 1, or 2)
- `beta_i` = effect weight from the PGS scoring file (typically log odds ratio or beta coefficient)

Missing genotypes (variant not in patient file) are excluded from the sum. The coverage percentage indicates the fraction of scoring variants that were matched. Scores with < 50% coverage should be interpreted with extra caution.

## Reference Distributions

Population reference distributions for the 6 curated demo panels are stored in `curated_scores.json`, which is generated from `CURATED_SCORES` in `gwas_prs.py` by `generate_curated_scores.py` and pinned against it field by field by `tests/test_score_provenance.py`. Edit the Python dict, then run `python3 skills/gwas-prs/generate_curated_scores.py`; `--check` fails if the committed file is stale. These distributions are based on European (EUR) reference populations. Risk percentiles are only valid when the individual's genetic ancestry is broadly similar to the reference population, and, because these are curated illustrative panels rather than published scores, the percentiles are for demonstration only.

**Ancestry caveat**: PRS performance varies across ancestries. Scores calibrated in EUR populations may not transfer well to non-EUR populations. Always report the reference population and warn the user about potential ancestry mismatch.

## PGS Catalog API

For scores beyond the 6 curated ones, query the PGS Catalog REST API:

```
# Search by trait
GET https://www.pgscatalog.org/rest/score/search?trait_id=EFO_0001360

# Get scoring file metadata
GET https://www.pgscatalog.org/rest/score/PGS000013

# Download harmonised scoring file
GET https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/PGS000013/ScoringFiles/Harmonized/PGS000013_hmPOS_GRCh37.txt.gz
```

## Safety

- **Genetic data never leaves this machine** — all processing is local. No genotype data is uploaded to any API.
- **Always include this disclaimer** in every report: *"ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Polygenic risk scores reflect statistical associations from population studies and do not determine individual outcomes. Consult a healthcare professional before making any medical decisions based on genetic information."*
- **Ancestry mismatch warning**: If the user's ancestry does not match the reference population, prominently warn that percentile estimates may not be accurate.
- **Coverage warning**: If variant coverage is below 50%, flag the score as unreliable.
- **No clinical decisions**: PRS results must not be used as the sole basis for clinical decisions. They are one factor among many (family history, lifestyle, clinical biomarkers).
- **Log all operations**: Record which scoring files were used, variant coverage, and calculation parameters.

## Integration with Bio Orchestrator

This skill is invoked by the Bio Orchestrator when:
- The user mentions "PRS", "polygenic risk score", "polygenic score", or "genetic risk score"
- The user asks about "GWAS risk", "genome-wide risk", or "multi-gene risk"
- The user asks about disease risk from their genetic data (beyond single-gene pharmacogenomics)
- Keywords detected: "prs", "polygenic", "gwas", "risk score"

It can be chained with:
- **pharmgx-reporter**: PRS provides disease risk context; PharmGx provides drug metabolism context. Together they give a comprehensive genomic health report.
- **nutrigx**: Combine PRS for metabolic traits (T2D, BMI) with nutrigenomic recommendations.
- **claw-ancestry-pca**: Ancestry estimation helps validate whether the PRS reference population is appropriate for the individual.
- **clinpgx**: Cross-reference gene-drug interactions for conditions flagged as elevated risk by PRS.
