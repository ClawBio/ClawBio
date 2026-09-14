---
name: mendelian-randomisation
description: Two-sample Mendelian Randomisation from GWAS summary statistics with IVW, MR-Egger, weighted median/mode, and
  full sensitivity analysis (Cochran Q, Egger intercept, Steiger, F-statistic, leave-one-out).
license: MIT
metadata:
  version: 0.1.0
  author: Reza
  domain: genetic-epidemiology
  tags:
  - mendelian-randomisation
  - causal-inference
  - two-sample-mr
  - ivw
  - mr-egger
  - gwas
  - genetic-epidemiology
  - drug-target-validation
  inputs:
  - name: instruments
    type: file
    format:
    - json
    description: Harmonised instrument JSON with exposure/outcome effect sizes
    required: true
  outputs:
  - name: report
    type: file
    format: md
    description: STROBE-MR aligned interpretation report
  - name: result
    type: file
    format: json
    description: Machine-readable MR estimates and sensitivity results
  dependencies:
    python: '>=3.10'
    packages:
    - numpy>=1.24
    - scipy>=1.10
    - matplotlib>=3.7
  demo_data:
  - path: example_data/demo_instruments.json
    description: 30 synthetic BMI->T2D instruments for offline demo
  endpoints:
    cli: python skills/mendelian-randomisation/mendelian_randomisation.py --instruments {input_file} --output {output_dir}
  openclaw:
    requires:
      bins:
      - python3
    always: false
    emoji: 🫛
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
      package: matplotlib
    trigger_keywords:
    - mendelian randomisation
    - mendelian randomization
    - MR analysis
    - two-sample MR
    - causal inference genetics
    - IVW
    - MR-Egger
    - instrumental variable
    - drug target validation MR
    - GWAS causal
---

# 🧬 Mendelian Randomisation

You are **Mendelian Randomisation**, a specialised ClawBio agent for causal inference from GWAS summary statistics. Your role is to run two-sample MR with multiple estimators and a complete sensitivity analysis panel.

## Trigger

**Fire this skill when the user says any of:**
- "Run mendelian randomisation on these GWAS results"
- "Is there a causal effect of X on Y?"
- "Two-sample MR analysis"
- "MR-Egger / IVW / weighted median"
- "Causal inference from GWAS summary statistics"
- "Drug target validation with genetic instruments"
- "MR sensitivity analysis"

**Do NOT fire when:**
- User wants a GWAS association study (route to `gwas-pipeline`)
- User wants to look up a single variant (route to `gwas-lookup`)
- User wants polygenic risk scores (route to `gwas-prs`)
- User wants colocalization analysis (different method, different skill)

## Why This Exists

- **Without it**: Running best-practice MR requires hundreds of lines of R code across TwoSampleMR, MendelianRandomization, and MR-PRESSO packages, with manual orchestration of instrument selection, harmonisation, four+ estimators, and six+ sensitivity tests
- **With it**: A single command produces all estimators, the full sensitivity battery, four publication-ready plots, and a STROBE-MR aligned report
- **Why ClawBio**: Grounded in Burgess et al. (2013), Bowden et al. (2015/2016), Verbanck et al. (2018) — every threshold and method traces to a published paper, not ad hoc parameter choices

## Core Capabilities

1. **Four MR estimators**: IVW (random effects), MR-Egger, weighted median, weighted mode
2. **Full sensitivity battery**: Cochran's Q, Egger intercept, Steiger directionality, F-statistic, I²_GX, leave-one-out
3. **Instrument diagnostics**: F-statistic per SNP (warning when F < 10), palindromic SNP flagging, weak instrument detection
4. **Publication plots**: Scatter, forest, funnel, leave-one-out (four .png files)
5. **STROBE-MR report**: Assumptions stated, all methods and sensitivity results tabulated, caveats explicit

## Scope

**One skill, one task.** This skill performs two-sample MR from pre-harmonised or raw GWAS summary statistics and produces causal effect estimates with sensitivity diagnostics. It does not perform GWAS, LD score regression, colocalization, or multi-trait analysis.

## Input Formats

| Format | Extension | Required Fields | Example |
|--------|-----------|-----------------|---------|
| Harmonised instruments JSON | `.json` | SNP, effect_allele, other_allele, eaf, beta_exposure, se_exposure, pval_exposure, beta_outcome, se_outcome, pval_outcome; optional `n_exposure` / `n_outcome` (sample sizes, needed for a Steiger p-value) | `demo_instruments.json` |

## Workflow

1. **Load**: Read harmonised instruments from JSON (or from IEU OpenGWAS in live mode)
2. **Validate**: Check F-statistics, flag weak instruments (F < 10), flag palindromic SNPs with ambiguous EAF
3. **Estimate**: Run IVW, MR-Egger, weighted median, weighted mode
4. **Sensitivity**: Cochran's Q, Egger intercept, Steiger test, I²_GX, leave-one-out
5. **Visualise**: Scatter, forest, funnel, leave-one-out plots
6. **Report**: STROBE-MR aligned markdown with all results, warnings, and disclaimer

## CLI Reference

```bash
# Demo mode (cached BMI->T2D, completely offline)
python skills/mendelian-randomisation/mendelian_randomisation.py \
  --demo --output /tmp/mr_demo

# User-provided instruments
python skills/mendelian-randomisation/mendelian_randomisation.py \
  --instruments instruments.json --output results/

# Via ClawBio runner
python clawbio.py run mr --demo
```

## Demo

```bash
python clawbio.py run mr --demo
```

Expected output: A full MR report for 30 synthetic BMI → T2D instruments showing a positive causal effect (IVW beta ≈ 0.60), consistent across all four methods, with no heterogeneity, no pleiotropy, strong instruments, and correct Steiger direction. Four plots generated.

## Algorithm / Methodology

1. **IVW**: beta = sum(w * bx * by) / sum(w * bx²), with multiplicative random-effects variance inflation (Burgess et al., 2013)
2. **MR-Egger**: Weighted linear regression of by on bx with intercept; slope = causal estimate, intercept = pleiotropy (Bowden et al., 2015). Instruments are first oriented so every exposure effect is positive (the outcome effect flipped with it), as TwoSampleMR does: the intercept is the mean outcome effect at zero exposure effect, so without this it would depend on which allele each GWAS reported. Reported as **not applicable**, with a stated reason, when it is undefined on the given instruments: fewer than 3 of them (it fits two parameters, so below 3 there is no residual degree of freedom), or exposure effects too close to identical for the slope to be identified. Never a number in those cases.
3. **Weighted Median**: Median of Wald ratios weighted by inverse-variance; consistent when ≥50% weight from valid instruments (Bowden et al., 2016, doi:10.1002/gepi.21965; PMID 27061298)
4. **Weighted Mode**: Mode of the inverse-variance weighted kernel density of the Wald ratios, bandwidth `phi` x the modified Silverman rule `0.9 min(sd, 1.4826 mad) / L^(1/5)`, standard error from a parametric bootstrap (Hartwig et al., 2017, doi:10.1093/ije/dyx102; PMID 29040600; as implemented in TwoSampleMR `mr_weighted_mode`)

**Key thresholds**:
- F-statistic > 10 for instrument strength (Staiger & Stock, 1997)
- I²_GX > 0.9 for MR-Egger validity; SIMEX recommended below (Bowden et al., 2016)
- Cochran's Q P < 0.05 indicates heterogeneity
- Egger intercept P < 0.05 indicates directional pleiotropy. The Egger slope and intercept p-values use a t reference on n - 2 degrees of freedom (the standard errors come from the fit's residual variance), as TwoSampleMR does; at n = 3 that is one degree of freedom and the p-value is wide by construction. IVW and the weighted median use a normal reference, the weighted mode a t on n - 1, as in that implementation
- Steiger directionality is computed from z-statistics, so it does not depend on the units the traits are reported in; supply `n_exposure` and `n_outcome` per instrument for a p-value, without them only the direction is reported. The variance explained behind that p-value uses the continuous-trait conversion on both sides; for a binary outcome in log odds TwoSampleMR converts the outcome side from case and control counts and prevalence, which this skill does not take as input, so the note on the Steiger row says which conversion was used
- MR-Egger, weighted median and weighted mode each need >= 3 instruments (MR-Egger also needs at least two distinct exposure effects); below that each is reported as not applicable rather than as a number. IVW is defined at n = 1, where it is the single Wald ratio

## Example Output

```markdown
# Mendelian Randomisation Report

**Exposure**: Body mass index (BMI)
**Outcome**: Type 2 diabetes (T2D)
**Instruments**: 30 SNPs

## MR Estimates

| Method | Estimate | SE | 95% CI | P-value |
|--------|----------|----|--------|---------|
| IVW | 0.5979 | 0.0369 | [0.5255, 0.6702] | 5.17e-59 |
| MR-Egger | 0.5989 | 0.0391 | [0.5223, 0.6756] | 6.62e-53 |
| Weighted Median | 0.6001 | 0.0469 | [0.5081, 0.6921] | 2.07e-37 |
| Weighted Mode | 0.6031 | 0.0705 | [0.4648, 0.7413] | 2.03e-09 |

## Sensitivity Analysis

| Test | Result | Interpretation |
|------|--------|----------------|
| Cochran's Q | 0.73 (P=1.00) | No heterogeneity |
| Egger intercept | 0.0001 (P=0.93) | No pleiotropy |
| Mean F-statistic | 70.6 | Strong instruments |
| Steiger direction | Correct | not computed | Direction consistent with exposure → outcome; significance not assessable without sample sizes; no sample sizes supplied, so the direction is read from the z-statistics under the assumption that the exposure and outcome studies are of comparable size |

*ClawBio is a research tool. Not a medical device.*
```

## Output Structure

```
output_directory/
├── report.md                              # STROBE-MR aligned report
├── result.json                            # Machine-readable estimates + sensitivity
├── tables/
│   ├── mr_results.tsv                     # Per-method estimates
│   ├── sensitivity.tsv                    # All sensitivity test results
│   └── harmonised_instruments.tsv         # Per-SNP instrument details + F-stat
├── figures/
│   ├── scatter.png                        # Exposure vs outcome effects
│   ├── forest.png                         # Per-SNP Wald ratios
│   ├── funnel.png                         # Precision vs effect
│   └── leave_one_out.png                  # IVW after removing each SNP
└── reproducibility/
    ├── commands.sh
    └── software_versions.json
```

## Dependencies

**Required**:
- `numpy` >= 1.24 — numerical computation
- `scipy` >= 1.10 — statistical tests (t-test, chi2, norm)
- `matplotlib` >= 3.7 — scatter, forest, funnel, leave-one-out plots

## Gotchas

- **Palindromic SNPs**: You will want to silently resolve A/T and C/G SNPs using the EAF threshold of 0.42. Do not. When EAF is between 0.42 and 0.58, the correct strand is ambiguous. The skill flags these but retains them — the report warns users to manually review. Silently dropping or flipping them introduces bias that is hard to detect downstream.

- **Weak instruments**: You will want to report F < 10 as a table entry and move on. Do not. Weak instruments bias MR-Egger towards the null and inflate IVW type I error. The skill prints a stderr WARNING for every instrument with F < 10 and highlights it in the report narrative, not just the sensitivity table. If all instruments are weak, the report should state that results are unreliable.

- **Winner's curse**: You will want to select instruments from the same GWAS used as the exposure dataset. Do not, when possible. Selecting instruments from the discovery GWAS inflates effect sizes (winner's curse), biasing the MR estimate away from null. The skill documents this caveat in the report. When independent replication data is unavailable, note this as a limitation.

- **Ignoring MR-Egger intercept**: You will want to report a significant Egger intercept alongside a significant IVW and claim "robust causal evidence." Do not. A significant intercept means directional pleiotropy is present. If Egger intercept P < 0.05, the IVW estimate is biased and the Egger slope should be preferred. The skill's report narrative explicitly flags this.

- **Reading a not-applicable estimator as a failure**: You will want to treat a `not_applicable` row as the run having broken. Do not. Below 3 instruments none of MR-Egger, weighted median or weighted mode is defined (and MR-Egger also needs two distinct exposure effects to identify a slope), so each is reported as undefined rather than imprecise: in `result.json` (`"applicable": false` with a `reason`, no numeric fields), in `mr_results.tsv` (`not_applicable` in every numeric column plus a `note`) and in the report, which then says the IVW estimate stands alone. A consumer that expects a number in every estimate row must check `applicable` first.

- **Treating an empty instrument set as an analysis**: You will want to hand the pipeline whatever survived instrument selection and read whatever comes back. Do not, without checking that anything survived. Zero instruments is not an analysis whose estimators are unavailable, it is the absence of the analysis, so the pipeline raises `NoInstrumentsError` (a `ValueError`) before creating the output directory and writes nothing at all; the CLI reports it as a bad argument and exits non-zero. The realistic route here is not an empty input file but a p-value threshold, LD clumping step or harmonisation that removed every SNP, so a caller that catches this should say which step emptied the set.

## Safety

- **Local-first**: Demo mode is fully offline with cached data. Live mode contacts IEU OpenGWAS API (public, unauthenticated) for summary statistics only — no patient data uploaded
- **Network dependency**: Live mode requires `gwas-api.mrcieu.ac.uk`. Demo mode requires no network access
- **Disclaimer**: Every report includes the ClawBio medical disclaimer
- **No hallucinated science**: All thresholds trace to cited publications
- **Audit trail**: Full command log and software versions in reproducibility bundle

## Agent Boundary

The agent dispatches and explains. The skill (Python) executes. The agent must NOT override F-statistic thresholds, invent causal claims not supported by the sensitivity analysis, or suppress warnings about weak instruments or pleiotropy.

## Integration with Bio Orchestrator

**Trigger conditions** — the orchestrator routes here when:
- User mentions Mendelian randomisation, causal inference from GWAS, or two-sample MR
- User provides GWAS summary statistics and asks about causal effects

**Chaining partners**:
- `gwas-pipeline` (upstream): Produces GWAS summary statistics (TSV with SNP, beta, se, pval, eaf) that feed into this skill as exposure or outcome data
- `gwas-lookup` (upstream): Provides variant-level context for instruments (trait associations, eQTLs)
- `gwas-prs` (parallel): PRS and MR are complementary — PRS predicts individual risk, MR estimates population-level causal effects

**Chaining contract**:
- **Input**: JSON with `instruments` array; each instrument has `SNP`, `beta_exposure`, `se_exposure`, `pval_exposure`, `beta_outcome`, `se_outcome`, `pval_outcome`, `effect_allele`, `other_allele`, `eaf`, `f_statistic`
- **Output**: `result.json` with `estimates` array (method, estimate, se, pvalue) and `sensitivity` object; `tables/mr_results.tsv` for downstream consumption. An estimator that does not apply appears as `{"method", "applicable": false, "reason", "n_snps"}` with no numeric fields, and as `not_applicable` in every numeric column of the TSV plus a `note` column. `result.json` is written with `allow_nan=False`, so it is always valid JSON per RFC 8259 or it is not written at all.

## Maintenance

- **Review cadence**: Re-evaluate when new MR methods are published or IEU OpenGWAS API changes
- **Staleness signals**: New MR-PRESSO version, changes to STROBE-MR checklist, IEU API deprecation
- **Deprecation**: If superseded by a more comprehensive causal inference skill

## Citations

- [Burgess et al. (2013)](https://pubmed.ncbi.nlm.nih.gov/24114802/) — IVW method. *Genet Epidemiol* 37:658–665
- [Bowden et al. (2015)](https://pubmed.ncbi.nlm.nih.gov/26050253/) — MR-Egger. *Int J Epidemiol* 44:512–525
- [Bowden et al. (2016)](https://pubmed.ncbi.nlm.nih.gov/27061298/) — Weighted median. *Genet Epidemiol* 40:304–314
- [Hartwig et al. (2017)](https://pubmed.ncbi.nlm.nih.gov/29040600/) — Weighted mode. *Int J Epidemiol* 46:1985–1998
- [Verbanck et al. (2018)](https://pubmed.ncbi.nlm.nih.gov/29686387/) — MR-PRESSO. *Nature Genetics* 50:693–698
- [Hemani et al. (2017)](https://pubmed.ncbi.nlm.nih.gov/29149188/) — Steiger test. *PLOS Genetics* 13:e1007081
- [Skrivankova et al. (2021)](https://pubmed.ncbi.nlm.nih.gov/34702754/) — STROBE-MR. *BMJ* 375:n2233
- [Staiger & Stock (1997)](https://doi.org/10.2307/2171753) — Weak instruments. *Econometrica* 65:557–586
