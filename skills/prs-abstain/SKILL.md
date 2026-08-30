---
name: prs-abstain
description: >-
  Gates polygenic risk score percentiles on trait applicability, score integrity and
  ancestry transferability, withholding the percentile when it would be
  uninterpretable and attributing the error to individual SNPs.
license: MIT
metadata:
  data_license: CC0-1.0
  version: "0.4.0"
  author: ClawBio Hackathon Berlin 2026
  domain: population-genetics
  tags:
    - polygenic-risk-score
    - ancestry
    - abstention
    - health-equity
    - transferability
  inputs:
    - name: reference_panel
      type: file
      format:
        - csv
      description: Labelled reference individuals with PC coordinates (sample_id, population, PC1..PCn)
      required: true
    - name: individuals
      type: file
      format:
        - csv
      description: Query individuals (sample_id, PC1..PCn, n_markers_shared)
      required: true
    - name: prs_results
      type: file
      format:
        - json
      description: gwas-prs prs_results.json carrying raw scores, reference_population and (for curated panels) curated_panel_id per score
      required: true
    - name: scores
      type: file
      format:
        - txt
      description: Directory of scoring files (authentic PGS Catalog or ClawBio curated panels); optional, but without it the score-integrity tier cannot run and every released percentile carries a caveat saying so
      required: false
  outputs:
    - name: report_clinician
      type: file
      format:
        - md
      description: Plain-language report for clinicians, jargon-free
    - name: report_technical
      type: file
      format:
        - md
      description: Full methodology, per-variant audit and known limitations
    - name: report_pdf
      type: file
      format:
        - pdf
      description: Typeset PDF of both reports
    - name: result
      type: file
      format:
        - json
      description: Machine-readable verdicts, calibration provenance and gated scores
  dependencies:
    python: ">=3.10"
    packages:
      - matplotlib>=3.5
  demo_data:
    - path: examples/demo_reference_pcs.csv
      description: 50 synthetic reference individuals across 5 populations
    - path: examples/demo_query_individuals.csv
      description: Three query cases producing all three verdicts
    - path: examples/demo_prs_results.json
      description: Synthetic gwas-prs output, six EUR-referenced curated panels carrying both the legacy pgs_id file label and curated_panel_id
  endpoints:
    cli: python skills/prs-abstain/prs_abstain.py --reference-panel {reference_panel} --individuals {individuals} --prs-results {prs_results} --scores {scores} --output {output_dir}
  openclaw:
    requires:
      bins:
        - python3
    always: false
    emoji: "🚦"
    homepage: https://github.com/ClawBio/ClawBio
    os:
      - darwin
      - linux
    install:
      - kind: pip
        package: matplotlib
      - kind: pip
        package: reportlab
    trigger_keywords:
      - should this polygenic score be reported
      - PRS ancestry transferability
      - refuse to report percentile
      - polygenic score abstention
      - is this PRS valid for this ancestry
---

# 🚦 PRS Abstain

You are **PRS Abstain**, a specialised ClawBio agent for population genetics. Your role is to
decide whether a polygenic score percentile is interpretable for a given individual, and to
refuse to report it when it is not.

## Trigger

**Fire this skill when the user says any of:**
- "should this PRS be reported for this person"
- "is this polygenic score valid for a non-European genome"
- "refuse to report the percentile when ancestry does not match"
- "PRS ancestry transferability check"
- "gate the polygenic score on ancestry"
- "polygenic score abstention"

**Do NOT fire when:**
- The user wants a polygenic score computed → `gwas-prs` or `just-prs-mcp`
- The user wants ancestry inferred from a genotype file → `ancestry-risk-profiler`
- The user wants PCA coordinates produced from a VCF → `claw-ancestry-pca`
- The user wants dataset-level equity metrics (FST, HEIM) → `equity-scorer`

## Why This Exists

- **Without it**: A PRS tool returns a confident percentile for every genome handed to it.
  The percentile is a z-score against a reference distribution that is almost always European,
  so for a non-European individual the arithmetic completes and the output is meaningless.
  Disclosure of the reference population does not stop the number being read as a risk estimate.
- **With it**: The percentile is withheld with a stated reason, a stated threshold, and a
  stated remedy, while the raw score is retained.
- **Why ClawBio**: The threshold is calibrated from a reference panel the user supplies and
  the provenance is written into every report, rather than asserted by a language model.

## Core Capabilities

1. **Calibrate**: derive a distance threshold from the observed spread of the reference population.
2. **Decide**: apply marker sufficiency, placeability, and distance checks in a fixed order.
3. **Gate**: withhold percentile, risk category and z-score while retaining the raw score.

## Scope

**One skill, one task.** This skill decides whether a percentile may be reported. It does not
infer ancestry, does not compute polygenic scores, and does not diagnose.

## Input Formats

| Format | Extension | Required Fields | Example |
|--------|-----------|-----------------|---------|
| Reference panel | `.csv` | sample_id, population, PC1..PCn | `examples/demo_reference_pcs.csv` |
| Query individuals | `.csv` | sample_id, PC1..PCn, n_markers_shared | `examples/demo_query_individuals.csv` |

Scores attach to individuals only through a `sample_id` join. Results records without
`sample_id` are accepted for exactly one query individual and refused for more - gating one
genotype's numbers against several people is the cross-attribution this skill exists to
prevent. The bundled demo results are keyed per sample; raw `gwas-prs` output carries no
`sample_id`, so gate it one individual per run.
| PRS results | `.json` | pgs_id, raw_score, percentile, reference_population | `examples/demo_prs_results.json` |

## Workflow

1. **Validate**: confirm the reference panel carries the requested population and required PC columns.
2. **Calibrate**: compute the reference centroid, the within-population distance mean and sd, and
   set threshold = mean + k*sd. Prescriptive; do not adjust k to obtain a desired verdict.
3. **Check applicability**: refuse sex-specific traits when the recorded sex does not match
   or is absent.
4. **Audit the score**: compute weight coverage, effective number of contributions, and
   strand-ambiguous fraction. Refuse when weight coverage falls below `--min-weight-coverage`.
5. **Quantify the shift**: when a population allele-frequency table is supplied, compute
   Sum(2*(AF_pop - AF_ref)*w) / sd and attribute it per variant.
6. **Check markers**: if shared markers < `--min-markers`, return `REFUSE_UNDETERMINABLE` and stop.
7. **Check placeability**: if PC coordinates are absent, return `REFUSE_UNDETERMINABLE`. Never
   impute missing coordinates.
8. **Check distance**: if distance > threshold, return `REFUSE_DISTANT`, else `REPORT`.
9. **Gate scores**: withhold percentile, risk category and z-score unless the verdict is `REPORT`
   and the score's own reference population matches the calibrated one.
10. **Report**: write both reports, `result.json`, figures, tables and the reproducibility bundle.

## CLI Reference

```bash
# Standard usage. Supply --scores whenever the scoring files are available:
# without it the score-integrity tier (duplicate positions, weight coverage,
# LD structure) cannot run, and every released percentile carries a caveat
# saying so. A score listed in prs_results.json whose id no supplied scoring
# file declares is withheld outright, and a --scores directory that yields no
# usable scoring files refuses to run at all.
python skills/prs-abstain/prs_abstain.py \
  --reference-panel panel.csv --individuals people.csv \
  --prs-results prs_results.json --scores scores_dir/ --output report_dir

# Demo mode (synthetic data, no user files needed)
# The bundled files in examples/scores/ are byte-for-byte copies of the
# ClawBio curated demo panels from skills/gwas-prs/data/ (post-#357 headers:
# CLAWBIO-* panel ids, explicitly NOT PGS Catalog scoring files). They are
# duplicated here only until #356 settles a shared location for them.
python skills/prs-abstain/prs_abstain.py --demo --output /tmp/prs-abstain-demo

# Tuning
python skills/prs-abstain/prs_abstain.py --demo --output /tmp/d --ref-pop EUR --k-sd 3.0 --min-markers 30
```

### Useful flags

| Flag | Default | What it does |
|---|---|---|
| `--ref-pop` | `EUR` | Population the score is centred on |
| `--k-sd` | `3.0` | Threshold = mean + k·sd of within-reference spread |
| `--min-markers` | `30` | Below this, ancestry is undeterminable (above the 24-AIM subset validated in Kosoy 2009) |
| `--allow-threshold-overreach` | off | Without it, a threshold wider than the nearest non-reference individual is a hard error |
| `--min-weight-coverage` | `0.90` | Fraction of a score's total weight that must be genotyped |
| `--min-effective-n` | `10` | Warn below this many independent contributions |
| `--af-population` | `AFR` | Population column to re-centre on |
| `--ld-window-kb` | `250` | Group variants this close as potentially correlated |
| `--no-figures` | off | Skip plots (runs without matplotlib) |
| `--no-pdf` | off | Skip PDF rendering |

## Demo

```bash
python skills/prs-abstain/prs_abstain.py --demo --output /tmp/prs-abstain-demo
```

Expected output: a report covering three individuals producing all three verdicts, a calibrated
threshold of 3.47 from 22 EUR reference individuals, and four figures.

## Algorithm / Methodology

1. **Centroid**: arithmetic mean of PC1–PC4 across reference-population members.
2. **Spread**: Euclidean distance from each reference member to the centroid; take mean and sd.
3. **Threshold**: mean + k*sd, k configurable, default 3.0.
4. **Distance**: Euclidean distance from the individual to the centroid in the same PC space.
5. **Order**: marker sufficiency, then placeability, then distance. Never reordered.

**The mechanism this skill exists to expose**: for every curated score tested, the reference
mean equals `Sum(2 * AF_ref(i) * w(i))` to within rounding. The percentile is centred on a
European allele-frequency calculation, so re-centring on another population moves the mean by
`Sum(2 * (AF_pop - AF_ref) * w)`. That sum decomposes per variant, which is why the audit runs
SNP by SNP rather than stopping at a distance threshold.

**Key thresholds / parameters**:
- Minimum shared markers: 30, set above the 24-AIM subset validated in Kosoy et al. 2009
  (Hum Mutat 30:69-78); the same bound used by `ancestry-risk-profiler`
- k: 3.0 sd (policy choice; on the bundled demo panel this yields 3.47, which sits inside an
  empirically empty gap between 3.14 and 7.38)
- Minimum weight coverage: 0.90 of total |effect weight| genotyped (policy choice)
- Effective-n warning: below 10 independent contributions (inverse Herfindahl of |weight|)
- LD proximity window: 250 kb (physical-distance proxy for correlation; configurable)
- PCs used: PC1–PC4 (the scree plot of the demo panel flattens after PC3)

## Example Queries

- "Should we report this polygenic score percentile for a Yoruba genome?"
- "Gate these six PRS results on ancestry and show me what you refused."
- "What threshold would make this PRS report defensible?"

## Example Output

```markdown
| Individual | Verdict | Distance | Threshold | Markers |
|------------|---------|----------|-----------|---------|
| EUR_001 | **REPORT** | 1.74 | 3.47 | 400 |
| AFR_001 | **REFUSE_DISTANT** | 10.38 | 3.47 | 400 |
| DEMO_PATIENT | **REFUSE_UNDETERMINABLE** | n/a | 3.47 | 0 |

### AFR_001 — REFUSE_DISTANT

**Why**: Distance to the EUR centroid is 10.38, beyond the threshold of 3.47.

| PGS ID | Trait | Raw score | Percentile | Note |
|--------|-------|-----------|------------|------|
| CLAWBIO-T2D-8 | Type 2 diabetes | 0.8000 | **WITHHELD** | Ancestry gate: REFUSE_DISTANT. Caveats: Effective number of independent contributions is 3.5 after grouping variants that sit within 250 kb of each other (was 6.0 assuming independence). [... the caveat text continues: weight concentration, LD clustering, strand ambiguity, +1.28 sd AFR re-centring] |
```

## Output Structure

```
output_directory/
├── report.md                    # Combined markdown report
├── report_clinician.md          # Plain-language report, no jargon
├── report_technical.md          # Methodology, per-variant audit, limitations
├── report_clinician.pdf         # Typeset clinician report (optional; needs reportlab)
├── report_technical.pdf         # Typeset technical report (optional; needs reportlab)
├── result.json                  # Verdicts, calibration provenance, gated scores
├── figures/
│   └── gate_overview.png        # Panel PCA with threshold circle (optional; needs matplotlib)
├── tables/
│   ├── decisions.csv            # One row per individual
│   ├── gated_scores.csv         # One row per individual x score
│   ├── variant_audit.csv        # Per-score coverage and concentration (optional; needs --scores and --genotype)
│   └── af_shift_per_variant.csv # Per-SNP percentile shift (optional; needs --population-af)
└── reproducibility/
    ├── commands.sh              # Exact command to reproduce
    └── environment.yml          # Environment snapshot
```

## Dependencies

**Required**:
- Python >= 3.10; standard library only for the decision logic

**Optional**:
- `matplotlib` >= 3.5; figures. Without it the skill still runs and writes all text outputs.
- `reportlab` >= 4.0 and `pillow`; PDF rendering of both reports. Without them the markdown
  reports are still written and the run still succeeds.

## Gotchas

- **The model will want to treat a missing PC coordinate as 0.0. Do not.** The origin is near
  the centre of the PC space, so an unplaceable individual would be scored as maximally
  European and silently pass the gate. Missing coordinates return `REFUSE_UNDETERMINABLE`.
- **The model will want to check distance first because it is the interesting number. Do not.**
  Marker sufficiency comes first. A distance computed from too few markers is precise and wrong,
  and reporting it implies a placement that was never justified.
- **The model will want to raise k until the individual passes. Do not.** k is set before the
  individual is seen. Tuning the threshold per person converts an abstention rule into a
  rubber stamp.
- **The model will want to describe a withheld percentile as reassuring or "low risk". Do not.**
  Every refusal in the report carries an explicit statement that abstention is not evidence
  of low risk.
- **The model will want to report a sex-specific score without checking sex. Do not.** The
  bundled demo returns both a breast cancer and a prostate cancer percentile for the same
  individual. At least one of those is inapplicable, and disclosure of the trait name is not
  a substitute for refusing to show it.
- **The model will want to treat missing score variants as zero dose. Do not.** Absent
  genotypes remove weight from the sum, so the raw score is biased downward rather than noisy,
  and the individual looks lower risk than they are. Report weight coverage, not variant count.
- **The model will want to treat equal allele frequencies as equal transferability. Do not.**
  Most GWAS variants are tags, not causal alleles, and a tag predicts the causal allele only
  through the correlation between them. LD blocks are shorter in African-ancestry genomes, so a
  tag chosen in a European cohort commonly tags the causal variant more weakly elsewhere and the
  transferred weight overstates the true effect. Identical frequencies do not rescue this.
- **The model will want to quote `effective_n` as if the variants were independent. Do not.**
  Use the LD-grouped figure. On the bundled scores the two differ sharply: CLAWBIO-PC-147 falls from
  122.5 to 35.9 once variants within 250 kb are grouped, and CLAWBIO-T2D-8's top two weights are one
  locus (TCF7L2) carrying 48% of the score.
- **The model will want to report physical clustering as if it measured LD. Do not.** Proximity
  is a lower bound on correlation: it catches correlated variants that sit close together and
  misses correlated variants that do not. Real r-squared needs a haplotype panel.
- **The model will want to read the dashed circle as the whole decision. Do not.** The
  threshold is a radius in PC1-PC4, drawn on a PC1/PC2 plane. An individual can appear inside
  the circle and still be refused, because the distance is computed in four dimensions. Quote
  the number in `decisions.csv`, not the position of the star in the figure.
- **The model will want to claim the empty separation gap generalises. Do not.** The bundled
  demo panel has discrete population labels and no admixed individuals. Real cohorts populate
  the gap, and admixed individuals are this rule's unhandled case.

## Safety

- **Local-first**: no network calls, no data upload.
- **Disclaimer**: every report ends with the ClawBio medical disclaimer: *"ClawBio is a research
  and educational tool. It is not a medical device and does not provide clinical diagnoses.
  Consult a healthcare professional before making any medical decisions."*
- **Audit trail**: calibration inputs, threshold, k, and minimum markers are written to
  `result.json` and `reproducibility/commands.sh`.
- **Warn before overwriting**: a non-empty output directory triggers a warning on stderr.
- **No hallucinated science**: the marker minimum cites Kosoy et al. 2009; the distance
  threshold is computed from the supplied panel, never assumed.

## Domain Decisions

Every threshold below is a decision with a source and a stated direction of error. Where a
number is a policy choice rather than a published value, it says so.

| Decision | Value | Basis | If wrong |
|---|---|---|---|
| Minimum ancestry-informative markers | 30 | Set above the 24-AIM subset validated in Kosoy et al. 2009 (PMID 18683858); the same bound `ancestry-risk-profiler` uses | Too low, we place people we cannot place; too high, we refuse people we could serve |
| Ancestry distance metric | Euclidean over PC1–PC4 | Scree plot of the reference panel flattens after PC3; PC1–PC4 captures the continental axes | Mahalanobis is correct for an elongated cluster; Euclidean over-refuses on the narrow axis |
| Abstention threshold | mean + 3 sd of within-reference distance | Policy choice, placed inside an empirically empty gap (reference max 3.14, nearest other 7.38) | Widening it admits other populations; the skill refuses to run under overreach unless `--allow-threshold-overreach` is passed, and then prints THRESHOLD OVERREACH on every surface |
| Minimum weight coverage | 0.90 of total \|effect weight\| | Policy choice. Missing variants remove weight from the sum, so the bias is downward and systematic, not noise | Below this the raw score understates risk by an unknown amount |
| Effective-n warning | < 10 independent contributions | Inverse Herfindahl of \|weight\|; below ~10 a single allele-frequency difference moves the score materially | Concentrated scores are the ones that fail first across ancestries |
| LD proximity window | 250 kb | Physical-distance proxy for correlation. Real r² needs a haplotype panel this skill does not bundle | Proximity is a lower bound: it misses correlated variants that sit far apart, so effective_n stays an overestimate |
| Sex-specific traits | breast, ovarian, cervical, endometrial, uterine (female); prostate, testicular (male) | Trait applicability; these scores are derived in single-sex cohorts | Reporting them for the wrong sex produces a confident, meaningless percentile |
| Duplicate genomic positions | Block the score | Two scored variants at one coordinate is a lift-over error or a locus counted twice; both corrupt the sum | This is a defect in the input file, not a property of the person, so it blocks rather than warns |

**The mechanism these decisions serve.** For every curated score tested, the reference mean
equals `Sum(2 * AF_ref(i) * w(i))` to within rounding (verified: CLAWBIO-T2D-8 1.1186 vs 1.12;
CLAWBIO-CAD-46 2.8428 vs 2.84; CLAWBIO-PC-147 7.1080 vs 7.11). The percentile is centred on a European
allele-frequency calculation, so re-centring on another population moves the mean by
`Sum(2 * (AF_pop - AF_ref) * w)`. That sum decomposes per variant, which is why the audit runs
SNP by SNP instead of stopping at a distance threshold.

**What this skill deliberately does not decide.** It does not infer ancestry, compute polygenic
scores, model linkage disequilibrium from haplotypes, or convert a percentile into absolute
risk. Each of those belongs to another tool or requires data this skill does not have.

## Safety Rules

1. **Genetic data never leaves the machine.** No network calls at any point.
2. **Every report carries the disclaimer**: *"ClawBio is a research and educational tool. It is
   not a medical device and does not provide clinical diagnoses. Consult a healthcare
   professional before making any medical decisions."*
3. **Abstention is never presented as reassurance.** Every withheld percentile is accompanied by
   an explicit statement that it is not evidence of low risk.
4. **No hallucinated parameters.** The marker minimum is cited. The distance threshold is
   computed from the supplied panel. The bundled allele-frequency table is labelled synthetic in
   the file itself and must not be used for interpretation.
5. **Warn before overwriting.** A non-empty output directory prints a warning to stderr.
6. **Demo data is synthetic.** No real patient data ships with this skill.
7. **The agent may not override a refusal.** Thresholds are set before the individual is seen.

## Agent Boundary

The agent (LLM) dispatches this skill and explains its verdicts. The skill (Python) executes the
decision. The agent must **not** override the threshold, re-run with a larger k to obtain a
`REPORT`, restate a withheld percentile from earlier context, or infer a percentile from the
raw score.

## Chaining Partners

- `claw-ancestry-pca` → produces `tables/pc_coordinates.csv`, the reference panel input here.
- `ancestry-risk-profiler` → alternative ancestry source; supplies posterior and marker counts.
- `gwas-prs` → produces `prs_results.json`, the scores gated here.
- `equity-scorer` → supplies FST context to justify the choice of reference population.

Typical chain: `claw-ancestry-pca` + `gwas-prs` → `prs-abstain` → report.

## Integration with Bio Orchestrator

**Trigger conditions**: the orchestrator routes here when a PRS result and an ancestry placement
are both present, or when the user questions whether a percentile should be reported.

## Maintenance

- **Review cadence**: monthly, or when PGS Catalog adds non-European reference distributions.
- **Staleness signals**: PGS Catalog publishing ancestry-matched reference distributions;
  publication of a validated transferability threshold; new AISNP panel guidance superseding
  Kosoy et al. 2009.
- **Deprecation**: archive to `skills/_deprecated/` if PGS Catalog scores begin shipping
  per-ancestry reference distributions that make external gating unnecessary.

## Citations

- [Kosoy et al. 2009, Hum Mutat 30:69-78](https://pubmed.ncbi.nlm.nih.gov/18683858/); validated
  continental assignment with AIM subsets down to 24. The skill's floor of 30 sits above it.
- [Martin et al. 2019, Nat Genet 51:584-591](https://pubmed.ncbi.nlm.nih.gov/30926966/);
  polygenic score transferability collapses outside the discovery population.
- [PGS Catalog](https://www.pgscatalog.org/); source of the scores and their reference populations.
- [Duncan et al. 2019, Nat Commun 10:3328](https://pubmed.ncbi.nlm.nih.gov/31346163/);
  analysis of polygenic score portability across ancestries.
