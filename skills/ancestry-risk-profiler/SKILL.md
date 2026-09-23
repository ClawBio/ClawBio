---
name: ancestry-risk-profiler
description: >-
  Infers genetic super-population ancestry from a 23andMe/AncestryDNA file and
  computes ancestry-stratified odds ratios with an exploratory Ancestry Elevation
  Score (AES) showing where ancestry-specific GWAS effect sizes diverge from
  European reference estimates. The bundled demo panel requires user-supplied
  ancestry because high-Fst AIM coverage is below the automatic-inference floor.
license: MIT
metadata:
  version: "1.4.0"
  author: ClawBio
  domain: population-genetics
  tags:
    - ancestry
    - disease-risk
    - population-genetics
    - gwas
    - ancestry-stratified
  inputs:
    - name: genotype_file
      type: file
      format:
        - txt
      description: 23andMe or AncestryDNA raw data file
      required: false
  outputs:
    - name: ancestry_risk_report.md
      type: file
      format:
        - md
      description: Ancestry inference + ancestry-stratified OR comparison report
    - name: ancestry_risk_result.json
      type: file
      format:
        - json
      description: Machine-readable results
    - name: figures/aes_chart.png
      type: file
      format:
        - png
      description: Ancestry Elevation Score bar chart (exploratory)
  dependencies:
    python: ">=3.11"
    packages:
      - matplotlib>=3.6
  demo_data:
    - path: data/demo_patient_south_asian.txt
      description: Synthetic South Asian 23andMe profile with T2D and CAD risk alleles; run with --demo --ancestry SAS
  endpoints:
    cli: python skills/ancestry-risk-profiler/ancestry_risk_profiler.py --input {genotype_file} --output {output_dir}
  openclaw:
    requires:
      bins:
        - python3
    always: false
    emoji: "🧬"
    homepage: https://github.com/ClawBio/ClawBio
    os:
      - darwin
      - linux
    install:
      - kind: pip
        package: matplotlib
    trigger_keywords:
      - ancestry risk
      - population-stratified risk
      - South Asian diabetes risk
      - ancestry-aware variant risk
      - which diseases am I at risk for given my ancestry
      - ancestry elevation score
      - APOL1 African kidney
      - KCNQ1 East Asian diabetes
      - genetic super-population disease risk
---

# 🧬 Ancestry-Aware Disease Risk Profiler

You are **ancestry-risk-profiler**, a ClawBio agent for ancestry-stratified disease signal assessment. Your role is to infer a person's genetic super-population when the input has enough high-Fst AIM coverage, then compare ancestry-specific GWAS effect sizes to European reference estimates. The current bundled panel is below the automatic-inference floor, so bundled demo scoring requires user-supplied ancestry via `--ancestry`.

## Trigger

**Fire this skill when the user says any of:**
- "given my ancestry, what diseases am I at risk for?"
- "does my genetic background affect my disease risk?"
- "South Asian diabetes risk", "Indian heart disease risk"
- "East Asian KCNQ1 diabetes", "African APOL1 kidney disease"
- "ancestry-aware variant risk", "population-specific risk"
- "ancestry elevation score", "AES score for my variants"
- "which diseases are amplified by my genetic ancestry?"

**Do NOT fire when:**
- User asks for pharmacogenomics / drug interactions → use `pharmgx-reporter`
- User asks for standard PRS / polygenic risk scores → use `gwas-prs`
- User asks for general variant annotation → use `variant-annotation`
- User asks to look up a specific rsID → use `gwas-lookup`
- User asks about ethnicity, nationality, or cultural background (this skill infers genetic super-population, not those things)

## Why This Exists

- **Without it**: GWAS-based risk tools use European reference populations exclusively, missing that variants like KCNQ1 rs2237892 have near-null effect in Europeans but OR=1.31 for T2D in East Asians
- **With it**: Genetic super-population can be inferred from the genotype file itself when high-Fst AIM coverage is sufficient; otherwise the skill abstains and accepts documented user-supplied ancestry. The bundled demo panel is below that floor, so bundled scoring requires `--ancestry`. Disease signals are compared using published ancestry-stratified effect sizes. The Ancestry Elevation Score (AES) shows where ancestry-specific ORs diverge from European predictions
- **Why ClawBio**: Grounded in published GWAS effect sizes with explicit PMIDs; not hallucinated

## Core Capabilities

1. **Ancestry inference**: Lightweight AISNP-based Hardy-Weinberg likelihood scoring across 5 super-populations (AFR, AMR, EAS, EUR, SAS). Matched panel SNPs contribute to the likelihood, but automatic inference is allowed only when at least 30 matched SNPs have Wright Fst ≥ 0.3. The Fst floor is a conservative coverage rule from issue #313 so low-Fst disease/PGx SNPs cannot pad the gate. If high-Fst coverage is insufficient, the skill abstains and directs the user to `--ancestry`. Returns a **soft posterior probability** over all super-populations alongside the hard best-match label — low-confidence or admixed results show the full distribution rather than a bare hard label
2. **Ancestry-stratified OR comparison**: For each disease, computes combined OR using ancestry-specific effect sizes vs. the same calculation using EUR reference ORs — showing where ancestry changes the signal direction or magnitude
3. **Ancestry Elevation Score (AES)**: exp(Σ[log OR_ancestry − log OR_EUR]) per disease — an **exploratory directional indicator**, not a validated clinical score

## Scope

**One skill, one task.** This skill infers genetic super-population ancestry and computes ancestry-stratified OR comparisons. It does NOT:
- Compute absolute lifetime risk percentages (applying ORs on top of population baseline prevalence double-counts allele contributions already embedded in that baseline; use `gwas-prs` instead)
- Perform pharmacogenomics, full PRS, variant annotation, or clinical ACMG classification
- Report on self-reported ethnicity, cultural identity, or nationality

**Genetic ancestry vs. ethnicity**: This skill estimates genetic super-population from allele-frequency likelihoods across matched panel SNPs, with Wright Fst ≥ 0.3 used as the coverage gate for automatic inference. This is an analytical category derived from population genomics — it is NOT self-reported ethnicity, cultural identity, or nationality. Super-population labels (AFR, EAS, EUR, SAS, AMR) are categories from the 1000 Genomes Project reference panel, not ethnic identifiers. Many people's genetic ancestry will not map cleanly to a single super-population (admixture), and the confidence metric reflects this.

## Input Formats

| Format | Extension | Notes |
|--------|-----------|-------|
| 23andMe raw | `.txt` | Tab-separated, rsid/chr/pos/genotype columns |
| AncestryDNA raw | `.txt` | Comma-separated, RSID/CHROMOSOME/POSITION/ALLELE1/ALLELE2 |

## Workflow

1. **Parse** genotype file → extract {rsid: genotype} dict, skip `--` no-calls
2. **Count informative AISNP coverage** → compute Wright Fst for each matched panel SNP; count only Fst ≥ 0.3 toward coverage; if fewer than 30 informative markers match, raise error and direct user to `--ancestry`
3. **Infer ancestry** → compute Hardy-Weinberg log-likelihood from matched panel SNPs for each of 5 super-populations; assign best match + confidence (gap in log-likelihood units)
4. **Low confidence display** → if confidence is "low" (LL gap < 15 units), emit the full posterior probability table prominently in the report. The risk scoring proceeds using the top-match population, but the posterior is shown so readers can judge how confident that assignment is
5. **Load associations** → read curated `ancestry_risk_associations.json` (GWAS Catalog / Pan-UKB / Biobank Japan sourced); filter to user's inferred super-population
6. **Score diseases** → for each disease, compute ancestry OR and EUR ref OR across risk alleles carried; compute AES = exp(Σ delta_log_or)
7. **Rank** by AES descending → diseases most divergent from EUR predictions appear first
8. **Generate report** → markdown with OR comparison table, AES bar chart, variant detail, gwas-prs referral, disclaimer

## CLI Reference

```bash
# Standard run
python skills/ancestry-risk-profiler/ancestry_risk_profiler.py \
  --input <23andme_file.txt> --output <report_dir>

# Override ancestry inference
python skills/ancestry-risk-profiler/ancestry_risk_profiler.py \
  --input <23andme_file.txt> --ancestry SAS --output <report_dir>

# Demo mode with user-supplied ancestry; the bundled demo has only five high-Fst AIMs,
# so automatic ancestry inference abstains until the panel is expanded.
python skills/ancestry-risk-profiler/ancestry_risk_profiler.py \
  --demo --ancestry SAS --output /tmp/ancestry_risk_demo
```

## Example Output

```markdown
# Ancestry-Aware Disease Risk Profile

## 1. Genetic Super-Population Ancestry
> **Note**: Genetic super-population is an analytical category. It is **not**
> self-reported ethnicity, cultural identity, or nationality. Labels (AFR, EAS, EUR,
> SAS, AMR) are analytical categories from the 1000 Genomes Project — not ethnic identifiers.

| Field | Value |
|---|---|
| Genetic super-population | **SAS** — South Asian *(user-supplied)* |
| Confidence | user-supplied |
| Informative AISNPs matched (Fst >= 0.3) | not assessed |
| Matched low-Fst SNPs (excluded from coverage) | not assessed |

Ancestry was supplied with `--ancestry`. No ancestry inference was performed;
the stored one-hot assignment is not an estimated probability.

## 2. Ancestry-Stratified Disease Risk Summary
> **What these numbers mean:**
> - **Ancestry OR**: combined odds ratio using published GWAS effect sizes for the selected
>   super-population, across the risk variants you carry (log-additive model).
> - **EUR ref OR**: the same calculation using European reference effect sizes for those
>   same variants — allows direct comparison.
> - **AES** (Ancestry Elevation Score): exp(Σ[log OR_ancestry − log OR_EUR]) — an
>   **exploratory directional indicator** showing where ancestry-specific effect sizes
>   diverge from European estimates. AES > 1 = ancestry amplifies the signal; AES < 1 =
>   ancestry attenuates it. **AES has not been externally validated and is not a clinical
>   risk score.** Use it as a signal to explore further, not a probability.
>
> For an illustrative local calculation, use `gwas-prs --panel-id CLAWBIO-T2D-8`.
> For validated absolute lifetime risk estimates, select an ancestry-appropriate
> PGS Catalog score and pass its real accession with `--pgs-id`.

| Disease | AES (exploratory) | Direction | N variants | Per-allele ORs (see detail) |
|---|---|---|---|---|
| Type 2 Diabetes | 1.84 | 🔴 Elevated By Ancestry | 7 | rs7903146 1.40x, rs13266634 1.17x, rs2237892 1.14x, rs7756992 1.18x, rs1552224 1.34x, rs8042680 1.21x, rs5219 1.23x |
| Coronary Artery Disease | 1.06 | 🟡 Neutral | 1 | rs1333049 1.34x |

> **Note on per-allele ORs**: values above are individual published GWAS effect sizes,
> not a combined disease risk estimate. Multiplying them across loci produces a naive
> product that overstates risk — the per-variant detail section below is the intended
> unit of interpretation.
```

## Output Structure

```
output_directory/
├── ancestry_risk_report.md       # Primary report
├── ancestry_risk_result.json     # Machine-readable results
└── figures/
    └── aes_chart.png             # AES horizontal bar chart (optional)
```

## Scoring Methodology

**Ancestry-stratified OR** (log-additive model):
```
combined_or = exp( Σᵢ log(OR_ancestry_i) × dosage_i )
or_eur_combined = exp( Σᵢ log(OR_EUR_i) × dosage_i )
```

**Ancestry Elevation Score (AES)** — exploratory, not validated:
```
AES = exp( Σᵢ [ log(OR_ancestry_i) − log(OR_EUR_i) ] × dosage_i )
```
- AES > 1.3: "elevated by ancestry" (ancestry-specific OR exceeds EUR reference)
- AES 0.77–1.3: "neutral"
- AES < 0.77: "reduced by ancestry"

These thresholds are for display colouring only. AES has no published external validation and is an exploratory metric.

**Why no absolute lifetime risk %?** Applying these ORs to a population baseline prevalence (e.g., 26.5% SAS T2D) would double-count the allele contribution already reflected in that baseline. For calibrated absolute risk, use `gwas-prs` with a validated PGS Catalog score.

## Gotchas

- **The model will want to run this for any variant question.** Do not. Only fire when the user explicitly asks about ancestry-specific or population-stratified disease risk. For general PRS, use `gwas-prs`.
- **If fewer than 30 ancestry-informative markers (Wright Fst ≥ 0.3) match, the skill MUST abstain and direct the user to `--ancestry`.** Nassir et al. (2009) and Kosoy et al. (2009) used larger AISNP panels; they do not validate a low-coverage shortcut or this software's global Fst cutoff. Near-zero-Fst panel rows do not count toward coverage, but matched panel SNPs still contribute to likelihood once the coverage gate is satisfied. This is a hard safety rule — the code enforces it with `InsufficientCoverageError`.
- **Low confidence does not mean wrong ancestry — it means admixed or ambiguous signal.** If confidence is "low" (LL gap < 15 units), warn prominently and suggest the user specify `--ancestry`. Do not refuse to run, but make the limitation visible.
- **Dosage is additive per allele for most loci.** If a user is homozygous for a risk allele (dosage=2), the log-OR is doubled. This is the standard log-additive GWAS assumption.
- **APOL1 (rs73885319 G1, rs60910145 G2) is an exception — it is recessive.** A single heterozygous APOL1 allele does NOT confer the full OR. Risk requires two high-risk alleles (G1+G2, G1/G1, or G2/G2). The code uses `model: "recessive_compound"` and counts total alleles across both loci before applying the validated compound OR (~7x). Do not change APOL1 to additive.
- **Combined OR is a naive product.** `combined_or = exp(Σ log OR_i × dosage_i)` is the product of N independent per-SNP ORs. It is **not** a validated polygenic score. Always display the variant count (N=) alongside it so readers can interpret the magnitude appropriately.
- **The curated panel covers ~10 diseases.** Diseases not in the panel are simply not reported — do not extrapolate or add unsupported associations.
- **AES is exploratory.** Do not present it as a validated score, percentile, or clinical probability. Always use the word "exploratory" when explaining it to the user.
- **"Genetic super-population" ≠ "ethnicity".** Use precise language. Never say "your ethnicity" when you mean "your inferred genetic super-population."

## Safety

- **Local-first**: All computation runs on-device; no genotype data is uploaded
- **Disclaimer**: Every report includes the ClawBio medical disclaimer
- **Cited sources**: Every association entry carries a non-null PMID (enforced by a CI test). One entry (ALDH2 rs671 ESCC) uses PMID 22960999 (Wu et al. 2012 Nat Genet) with GWAS Catalog accession GCST001563; the entry note flags that the accession-to-PMID mapping was not independently confirmed via the Catalog. See `data/PROVENANCE.md` for full correction history.
- **LCT entries removed in v1.3.0**: All `rs4988235` Lactose Intolerance entries were removed because (a) PMID 14507249 cited as Enattah 2002 resolves to an unrelated bladder-cancer paper, and (b) the EUR `or=0.45` and non-EUR `or=3.2–6.8` encoded opposite outcome framings for the same allele, manufacturing spurious AES of 7–15x.
- **Soft posterior gating in v1.3.0**: When ancestry confidence is "low" and top posterior < 0.35, disease risk scoring is skipped entirely with an explanatory message. Between 0.35–0.6 a caveat note is shown alongside results. User-supplied `--ancestry` overrides the gate.
- **No hallucinated ORs**: All effect sizes trace to the bundled `ancestry_risk_associations.json`
- **AIM Fst gate (v1.4.0)**: Coverage counts only panel SNPs with Wright Fst ≥ 0.3, with a minimum of 30 informative markers. A 23andMe extract of 30 T2D SNPs must still abstain. The current bundled panel has five high-Fst markers, so automatic inference abstains unless the user supplies `--ancestry`. See issue #313.
- **APOL1 recessive model**: APOL1 G1/G2 use a compound-recessive model; per-allele log-additive OR is biologically wrong for this locus
- **Combined OR is naive**: `combined_or` is the product of independent per-SNP ORs (log-additive); it is NOT a validated aggregate risk score. `N=` in the report shows how many variants contribute so readers can judge the calculation
- **No absolute risk claims**: The Cornfield/baseline prevalence calculation has been removed to prevent double-counting

## Agent Boundary

The agent (LLM) dispatches and explains results. The skill (Python) executes the inference and scoring. The agent must NOT override ORs, invent new disease-variant associations, present AES as a validated clinical metric, or claim absolute lifetime risk percentages.

## Integration with Bio Orchestrator

**Trigger conditions**: routes here when:
- Query mentions genetic ancestry + disease risk together
- File is a 23andMe/AncestryDNA raw data file AND query is about disease risk stratified by population

**Chaining partners**:
- `gwas-prs`: for validated absolute risk scores with ancestry-appropriate PGS Catalog scores — always signpost this when users ask about lifetime risk
- `pharmgx-reporter`: after ancestry signal profiling, run pharmgx to add drug response context
- `profile-report`: ancestry-risk-profiler output feeds into the unified profile report

## Maintenance

- **Review cadence**: Update `ancestry_risk_associations.json` when major multi-ancestry GWAS meta-analyses are published (Pan-UKB updates, Global Biobank Meta-Analysis Initiative releases)
- **Staleness signals**: New population-specific GWAS not yet in panel; GWAS Catalog accumulates new non-EUR studies
- **Deprecation**: Archive if a more comprehensive multi-ancestry PRS tool with per-ancestry calibration supersedes this approach

## Citations

- Genovese et al. (2010) Science 329:841. PMID 20566908. APOL1 G1/G2 kidney disease (recessive compound model)
- Karczewski et al. (2020) Nature 581:434. PMID 32461654. gnomAD v3.1 allele frequencies for AISNP panel
- Nassir et al. (2009) BMC Genetics 10:39. [PMID 19630973](https://pubmed.ncbi.nlm.nih.gov/19630973/). DOI 10.1186/1471-2156-10-39. 93 AISNP panel for continental ancestry inference; does not validate a low-coverage implementation shortcut
- Kosoy et al. (2009) Human Mutation 30(1):69–78. [PMID 18683858](https://pubmed.ncbi.nlm.nih.gov/18683858/). DOI 10.1002/humu.20822. 128-marker AISNP set and 96/64/48/24-marker subsets; does not define this software's global Fst ≥ 0.3 gate
- Dubois et al. (2010) Nat Genet 42:295–302. PMID 20190752. Celiac disease GWAS (PTPN22 R620W)
- Wu et al. (2012) Nat Genet 44:1090–1093. PMID 22960999. GWAS Catalog GCST001563. ESCC GWAS in Chinese (ALDH2 rs671); two earlier wrong PMIDs (20686008, 22561518) were corrected, and the accession-to-PMID mapping is flagged for confirmation in the entry note
- Grant et al. (2006) Nat Genet 38:320–323. PMID 16415884. TCF7L2 T2D discovery
- Zeggini et al. (2007) Nat Genet 39:638–644. PMID 17463246. T2D replication
- Feder et al. (1996) Nat Genet 13:399–408. PMID 9068472. HFE hereditary haemochromatosis

**Removed citations**: Enattah et al. (2002) LCT lactase persistence — LCT rs4988235 entries removed in v1.3.0 (direction artifact + PMID 14507249 was wrong). PMID 22561518 (Wu 2012 ESCC) — resolves to Jin 2012 vitiligo GWAS.

See `data/PROVENANCE.md` for the full citation table with correction history.
