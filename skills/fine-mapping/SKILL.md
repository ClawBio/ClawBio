---
name: fine-mapping
description: Statistical fine-mapping of GWAS loci using SuSiE, SuSiE-inf, and Approximate Bayes Factors to identify credible
  sets and posterior inclusion probabilities (PIPs) for causal variant discovery. SuSiE-inf adds an infinitesimal polygenic
  component for improved calibration at well-powered loci.
license: MIT
metadata:
  version: 0.3.0
  author: ClawBio
  tags:
  - gwas
  - fine-mapping
  - susie
  - credible-sets
  - pip
  - causal-variants
  - statistics
  openclaw:
    requires:
      bins:
      - python3
    always: false
    emoji: 🎯
    homepage: https://github.com/ClawBio/ClawBio
    os:
    - darwin
    - linux
    install:
    - kind: pip
      package: numpy
    - kind: pip
      package: scipy
    - kind: pip
      package: pandas
    - kind: pip
      package: matplotlib
    - kind: pip
      package: sushie
    trigger_keywords:
    - fine-mapping
    - finemapping
    - susie
    - susie-inf
    - susieinf
    - infinitesimal fine-mapping
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
- **Why ClawBio**: Runs locally without uploading individual-level data; implements ABF natively and delegates SuSiE to the published [sushie](https://github.com/mancusolab/sushie) package — no R dependency required

## Core Capabilities

1. **Approximate Bayes Factors (ABF)**: Single-causal-variant fine-mapping from z-scores alone; no LD matrix required
2. **SuSiE (Sum of Single Effects)**: Multi-signal fine-mapping with LD, delegated to the sushie package (mancusolab/sushie, JAX-based) via its summary-statistics interface run with a single ancestry; requires the `fine-mapping` extra (`uv sync --extra fine-mapping`)
3. **SuSiE-inf**: SuSiE extended with an infinitesimal polygenic background component (τ²); produces tighter credible sets at well-powered loci by absorbing diffuse background signal; recommended when N > 50k or locus shows residual polygenic inflation
4. **Swappable benchmark**: `tests/benchmark/finemapping_benchmark.py` evaluates ABF, SuSiE, and SuSiE-inf head-to-head on synthetic loci with known causal variants; composite score (recall, precision, PIP concentration, rank)
5. **Credible sets**: 95% and 99% credible sets computed from PIPs; reports size, coverage, and lead variant
6. **Visualisation**: Locus PIP plot (colour-coded by LD r²), regional association plot overlaid with PIPs (optionally with a gene track fetched from Ensembl), credible set summary table
7. **LD computation**: Accepts a pre-computed LD matrix (`.npy` or `.tsv`)

## Input Formats

| Format | Extension | Required Fields | Example |
|--------|-----------|-----------------|---------|
| GWAS summary stats | `.tsv` / `.csv` / `.txt` | `rsid`, `chr`, `pos`, `beta`, `se` **or** `z` | `locus_sumstats.tsv` |
| Pre-computed LD matrix | `.npy` / `.tsv` | Square correlation matrix, row/col = variant order | `ld_matrix.npy` |
| Demo (built-in) | — | — | `--demo` |

Optional columns in sumstats: `p`, `maf`, `n`, `a1`, `a2`

## Workflow

When the user asks for fine-mapping:

1. **Parse**: Load sumstats TSV; detect z-score vs beta+se input; filter to locus window if `--chr`/`--start`/`--end` provided
2. **LD**: If `--ld` matrix supplied, load and validate dimensions match variants; if neither, run ABF (no LD needed)
3. **Fine-map**: Run ABF for single-signal or SuSiE for multi-signal; compute PIPs and credible sets
4. **Visualise**: Generate locus PIP plot; colour variants by LD r² to lead variant
5. **Report**: Write `report.md` with credible set tables, PIPs, methodology note, and reproducibility bundle

## CLI Reference

```bash
# ABF single-signal fine-mapping (no LD needed; no extra required)
python skills/fine-mapping/fine_mapping.py \
  --sumstats locus.tsv --output /tmp/finemapping

# SuSiE multi-signal with pre-computed LD matrix (sushie engine:
# install once with `uv sync --extra fine-mapping`)
uv run --extra fine-mapping python skills/fine-mapping/fine_mapping.py \
  --sumstats locus.tsv --ld ld_matrix.npy --output /tmp/finemapping

# Filter to a specific locus window
uv run --extra fine-mapping python skills/fine-mapping/fine_mapping.py \
  --sumstats gwas_full.tsv --chr 1 --start 109000000 --end 110000000 \
  --ld ld_matrix.npy --output /tmp/finemapping

# Set maximum number of causal signals (SuSiE L parameter)
uv run --extra fine-mapping python skills/fine-mapping/fine_mapping.py \
  --sumstats locus.tsv --ld ld_matrix.npy --max-signals 5 --output /tmp/finemapping

# Add a gene track below the regional association plot (requires internet)
uv run --extra fine-mapping python skills/fine-mapping/fine_mapping.py \
  --sumstats locus.tsv --ld ld_matrix.npy --gene-track --output /tmp/finemapping

# Demo mode (synthetic 200-variant locus, two causal signals)
uv run --extra fine-mapping python skills/fine-mapping/fine_mapping.py --demo --output /tmp/finemapping_demo
```

## Demo

```bash
uv run --extra fine-mapping python skills/fine-mapping/fine_mapping.py --demo --output /tmp/finemapping_demo
```

Expected output: a report covering a synthetic 200-variant locus with two injected causal signals, two single-variant SuSiE credible sets pinpointing the causal variants (indices 60 and 140), per-variant PIP plot, and reproducibility bundle.

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

When an LD matrix **R** is provided, the locus is fine-mapped by the sushie
package (`sushie.infer_ss.infer_sushie_ss` with a single ancestry), which
implements the SuSiE model with effect variances estimated by EM:

1. sushie fits L single effects (default 10) on the z-scores and LD matrix, run in float64 (jax x64) to avoid ELBO precision failures
2. sushie prunes effects to those forming valid credible sets (coverage threshold + purity); the adapter returns only these active signals, so null loci yield zero credible sets and no phantom PIPs
3. PIPs are sushie's `pip_cs` — computed over the kept signals: `PIP_i = 1 - prod_l (1 - α_l_i)`
4. Credible sets for the report: greedily add highest-α variants per signal until cumulative α ≥ 0.95, with the min-|r| purity flag (Wang 2020 §3.2)

### SuSiE-inf (Cui et al. 2024)

Extends SuSiE with an infinitesimal variance component τ² that captures diffuse polygenic signal. The residual precision matrix becomes:

```
Ω = (τ² · D² + σ² · I)⁻¹   in the LD eigenbasis
```

where D² are eigenvalues of X'X (n × LD eigenvalues). When τ²→0 the model reduces to standard SuSiE.

1. Eigendecompose LD once: `LD = V diag(d²/n) V'`
2. IBSS loop with Ω-weighted residuals instead of σ²-only residuals
3. Method-of-moments update for σ² and τ² each iteration
4. Credible sets via per-effect PIPs (p×L matrix) with purity filter

**When to prefer SuSiE-inf over SuSiE**:
- Large cohort (N > 50k): background polygenic signal is detectable
- Locus shows many nominally associated variants (diffuse signal)
- SuSiE returns very large credible sets (many variants absorbed as "sparse" effects)

**Key thresholds / parameters**:
- Prior W (ABF): 0.04 (source: Wakefield 2009, Am J Hum Genet)
- Credible set coverage: 95% (adjustable via `--coverage`)
- Max signals L: 10 (adjustable via `--max-signals`)
- Min purity (SuSiE/SuSiE-inf CS filter): 0.5 average pairwise LD r² within set
- Convergence tolerance (SuSiE engine): ELBO change < 1e-4 (sushie `min_tol`)

## Gotchas

1. **`--prior-variance` is a seed, not a fixed prior, under SuSiE.** The model will want to treat `w` as it does for ABF (a fixed Wakefield prior). Do not. sushie seeds its `effect_var` with `w` and then re-estimates it by EM every iteration, so two runs with different `w` usually converge to the same fit. Only ABF honours `w` exactly.
2. **`mu`/`mu2` from `run_susie` are not susieR z-unit moments.** The model will want to sanity-check `mu` against the single-effect shrinkage formula `r · z` with `r = w/(w + 1/n)`. Do not. sushie reports conditional posterior moments on its standardised effect-size scale; on a `z=[5,5,0], n=100` locus susieR-style `mu` is 4.0 while sushie's `post_mean` is ~0.2. Same quantity, different units — compare shapes and ordering, not magnitudes.
3. **Pruned signals are dropped, not zeroed.** `alpha`, `mu` and `mu2` contain only the signals sushie kept as credible sets at the requested `coverage` and `min_purity`. A null locus, or a locus whose only signal is spread over uncorrelated variants (purity 0), returns arrays with **zero rows** and all-zero PIPs. Do not index `alpha[0]` without checking `alpha.shape[0]` first.
4. **Non-convergence is a warning plus a flag, not an exception.** Hitting `max_iter` emits a `RuntimeWarning` and sets `converged: False`, mirroring `susieR::susie_rss`; finite PIPs are still returned. The model will want to report those PIPs as results. Do not — surface `converged` in the report and say the estimate is provisional.
5. **`coverage` and `min_purity` must lie strictly inside (0, 1).** sushie rejects the endpoints, so `--coverage 1.0` or `--min-purity 0` raise a `ValueError` before fitting. The old pure-Python engine accepted them; scripts that passed `1.0` need updating.

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

**SuSiE engine** (ABF works without it):
- `sushie` >= 0.20, < 0.21 — SuSiE inference (pulls jax, jaxlib, equinox, polars, glimix-core); install with `uv sync --extra fine-mapping` and run via `uv run --extra fine-mapping python ...`


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
- [Cui et al. (2024) Nature Genetics](https://doi.org/10.1038/s41588-023-01597-3) — SuSiE-inf: improving fine-mapping by modeling infinitesimal effects
