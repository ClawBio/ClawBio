---
name: fine-mapping
description: >-
  Statistical fine-mapping of GWAS loci using SuSiE (Sum of Single Effects)
  and Approximate Bayes Factors to identify credible sets and posterior inclusion
  probabilities (PIPs) for causal variant discovery.
version: 0.1.0
author: ClawBio
license: MIT
tags: [gwas, fine-mapping, susie, credible-sets, pip, causal-variants, statistics]
metadata:
  openclaw:
    requires:
      bins:
        - python3
      env: []
      config: []
    always: false
    emoji: "🎯"
    homepage: https://github.com/ClawBio/ClawBio
    os: [darwin, linux]
    install:
      - kind: pip
        package: numpy
        bins: []
      - kind: pip
        package: scipy
        bins: []
      - kind: pip
        package: pandas
        bins: []
      - kind: pip
        package: matplotlib
        bins: []
    trigger_keywords:
      - fine-mapping
      - finemapping
      - susie
      - credible set
      - posterior inclusion probability
      - PIP
      - causal variant
      - fine map
      - ABF
      - approximate bayes factor
      - FINEMAP
      - polyfun
      - fine map locus
      - causal SNP
---

# 🎯 SuSiE Fine-Mapper

You are **SuSiE Fine-Mapper**, a specialised ClawBio agent for statistical fine-mapping of GWAS loci. Your role is to identify credible sets of likely causal variants and compute per-variant posterior inclusion probabilities (PIPs) from GWAS summary statistics.

## Why This Exists

GWAS identifies associated loci, not causal variants. A single GWAS signal can contain dozens of correlated SNPs in high LD — fine-mapping colocalises the signal onto the minimal credible set of likely causal variants.

- **Without it**: Researchers must manually triage 10–200 correlated SNPs per locus with no principled prioritisation
- **With it**: A ranked credible set with PIPs and 95% credible set boundaries in seconds
- **Why ClawBio**: Runs locally without uploading individual-level data; implements ABF natively and wraps SuSiE (via polyfun) when available — no R dependency required

## Core Capabilities

1. **Approximate Bayes Factors (ABF)**: Single-causal-variant fine-mapping from z-scores alone; no LD matrix required
2. **SuSiE (Sum of Single Effects)**: Multi-signal fine-mapping with LD using the iterative Bayesian stepwise selection algorithm; wraps `polyfun` if installed, falls back to a pure-Python SuSiE implementation
3. **Credible sets**: 95% and 99% credible sets computed from PIPs; reports size, coverage, and lead variant
4. **Visualisation**: Locus PIP plot (colour-coded by LD r²), regional association plot overlaid with PIPs (optionally with a gene track fetched from Ensembl), credible set summary table
5. **LD computation**: Optionally computes LD from a PLINK .bed file in the locus window; or accepts a pre-computed LD matrix

## Input Formats

| Format | Extension | Required Fields | Example |
|--------|-----------|-----------------|---------|
| GWAS summary stats | `.tsv` / `.csv` / `.txt` | `rsid`, `chr`, `pos`, `beta`, `se` **or** `z` | `locus_sumstats.tsv` |
| Pre-computed LD matrix | `.npy` / `.tsv` | Square correlation matrix, row/col = variant order | `ld_matrix.npy` |
| PLINK genotype data | `.bed/.bim/.fam` | Standard PLINK binary format | `reference_panel` |
| Demo (built-in) | — | — | `--demo` |

Optional columns in sumstats: `p`, `maf`, `n`, `a1`, `a2`

## Workflow

When the user asks for fine-mapping:

1. **Parse**: Load sumstats TSV; detect z-score vs beta+se input; filter to locus window if `--chr`/`--start`/`--end` provided
2. **LD**: If `--ld` matrix supplied, load and validate dimensions match variants; if `--bfile` supplied, compute LD with `bed-reader`+numpy; if neither, run ABF (no LD needed)
3. **Fine-map**: Run ABF for single-signal or SuSiE for multi-signal; compute PIPs and credible sets
4. **Visualise**: Generate locus PIP plot; colour variants by LD r² to lead variant
5. **Report**: Write `report.md` with credible set tables, PIPs, methodology note, and reproducibility bundle

## CLI Reference

```bash
# ABF single-signal fine-mapping (no LD needed)
python skills/fine-mapping/fine_mapping.py \
  --sumstats locus.tsv --output /tmp/finemapping

# SuSiE multi-signal with pre-computed LD matrix
python skills/fine-mapping/fine_mapping.py \
  --sumstats locus.tsv --ld ld_matrix.npy --output /tmp/finemapping

# Filter to a specific locus window
python skills/fine-mapping/fine_mapping.py \
  --sumstats gwas_full.tsv --chr 1 --start 109000000 --end 110000000 \
  --ld ld_matrix.npy --output /tmp/finemapping

# Set maximum number of causal signals (SuSiE L parameter)
python skills/fine-mapping/fine_mapping.py \
  --sumstats locus.tsv --ld ld_matrix.npy --max-signals 5 --output /tmp/finemapping

# Add a gene track below the regional association plot (requires internet)
python skills/fine-mapping/fine_mapping.py \
  --sumstats locus.tsv --ld ld_matrix.npy --gene-track --output /tmp/finemapping

# Demo mode (synthetic 200-variant locus, two causal signals)
python skills/fine-mapping/fine_mapping.py --demo --output /tmp/finemapping_demo
```

## Demo

```bash
python skills/fine-mapping/fine_mapping.py --demo --output /tmp/finemapping_demo
```

Expected output: a report covering a synthetic 200-variant locus with two injected causal signals, SuSiE credible sets of ~3–8 variants each, per-variant PIP plot, and reproducibility bundle.

## Algorithm / Methodology

### Approximate Bayes Factors (ABF)

Used when no LD matrix is available (assumes variants are independent).

For each variant *i* with z-score *z_i* and prior variance *W*:

```
V_i  = 1 / n_eff    (if se available: V_i = se_i^2)
ABF_i = sqrt(V_i / (V_i + W)) * exp(z_i^2 * W / (2 * (V_i + W)))
PIP_i = ABF_i / sum(ABF_j)
```

Default prior: W = 0.04 (σ = 0.2 on log-OR scale; Wakefield 2009)

### SuSiE (Sum of Single Effects, Wang et al. 2020)

When an LD matrix **R** is provided:

1. Initialise L single-effect vectors **α_l** (L = number of expected causal signals, default 10)
2. Iterative Bayesian Stepwise Selection (IBSS):
   - For each effect l, compute residual z-scores removing all other effects
   - Update **α_l** via single-effect regression posterior: `α_l ∝ ABF(z_residual | R)`
   - Update posterior variance `μ_l²` and `σ_l²`
3. Converge when ELBO change < 1e-3 (max 100 iterations)
4. PIPs: `PIP_i = 1 - prod_l (1 - α_l_i)`
5. Credible sets: greedily add highest-PIP variants until cumulative PIP ≥ 0.95

**Key thresholds / parameters**:
- Prior W (ABF): 0.04 (source: Wakefield 2009, Am J Hum Genet)
- Credible set coverage: 95% (adjustable via `--coverage`)
- Max signals L: 10 (adjustable via `--max-signals`)
- Min purity (SuSiE CS filter): 0.5 average pairwise LD r² within set
- Convergence tolerance: ELBO Δ < 1e-3

## Example Queries

- "Fine-map the PCSK9 locus from my GWAS summary stats"
- "Run SuSiE on this locus with the LD matrix"
- "What's the credible set for rs562556?"
- "Compute PIPs for all variants in my GWAS locus file"
- "Run fine-mapping demo so I can see the output"
- "Which variants have PIP > 0.1 in this locus?"

## Output Structure

```
output_directory/
├── report.md                    # Primary markdown report
├── fine_mapping.json            # Machine-readable PIPs + credible sets
├── figures/
│   ├── pip_locus_plot.png       # Per-variant PIP coloured by LD r²
│   ├── regional_association.png # -log10(p) with lead variant highlighted (only if p-values present)
│   └── ld_heatmap.png           # LD r² heatmap with credible set annotations (only if LD matrix provided)
├── tables/
│   ├── pips.tsv                 # rsid, chr, pos, pip, cs_membership
│   └── credible_sets.tsv        # cs_id, size, coverage, lead_rsid, variants
└── reproducibility/
    ├── commands.sh              # Exact command to reproduce
    └── environment.yml          # Package versions
```

## Dependencies

**Required**:
- `numpy` >= 1.24 — array maths, LD matrix operations
- `scipy` >= 1.10 — statistical functions
- `pandas` >= 1.5 — sumstats parsing
- `matplotlib` >= 3.7 — locus plots

**Optional**:
- `polyfun` — SuSiE-inf (infinite-mixture prior); automatically used if importable
- `bed-reader` — LD computation from PLINK .bed files

## Safety

- **Local-first**: No data upload; all computation is on-machine
- **Disclaimer**: Every report includes the ClawBio medical disclaimer
- **Audit trail**: `reproducibility/commands.sh` logs exact inputs and parameters
- **No hallucinated science**: All parameters trace to cited papers; model outputs are probabilistic, not clinical diagnoses

## Integration with Bio Orchestrator

**Trigger conditions** — the orchestrator routes here when:
- Query contains "fine-map", "finemapping", "credible set", "PIP", "posterior inclusion"
- File has columns: `beta`/`z` + `se` (looks like GWAS summary stats)
- Query mentions SuSiE, FINEMAP, CAVIAR, ABF, polyfun

**Chaining partners** — this skill connects with:
- `gwas-lookup`: look up the lead variant before fine-mapping to confirm locus context
- `gwas-prs`: fine-mapped causal variants can be used as a more precise PRS variant set
- `vcf-annotator`: annotate the credible set variants with functional consequences

## Citations

- [Wang et al. (2020) JRSS-B](https://doi.org/10.1111/rssb.12388) — SuSiE algorithm
- [Wakefield (2009) Am J Hum Genet](https://doi.org/10.1016/j.ajhg.2008.12.010) — Approximate Bayes Factors for GWAS
- [Benner et al. (2016) Bioinformatics](https://doi.org/10.1093/bioinformatics/btw018) — FINEMAP comparison
- [Zou et al. (2022) Nature Genetics](https://doi.org/10.1038/s41588-022-01167-z) — SuSiE-inf / polyfun
