---
name: coloc-abf-region-run
description: |
  Colocalise two traits over one genomic region with coloc.abf. Computes the five
  posterior probabilities (PP.H0 to PP.H4) from each trait's regional beta and
  standard error using Wakefield approximate Bayes factors and the
  Giambartolomei single-causal-variant combine, and reports whether the
  conclusion survives the p12 prior. Answers whether an eQTL, pQTL or caQTL
  signal and a GWAS signal share a causal variant, or merely sit in the same
  locus. A clean-room Python implementation: nothing is imported from the GPL R
  coloc package.
license: MIT
metadata:
  author: Aviv Madar
  version: 0.1.0
  domain: statistical-genetics
  tags:
    - colocalisation
    - colocalization
    - coloc
    - coloc-abf
    - bayesian
    - gwas
    - eqtl
    - statistical-genetics
    - drug-target-validation
  inputs:
    - name: config
      type: file
      format:
        - yaml
        - json
      description: Config naming two regional summary-statistics files, each with its trait type (quant or cc) and either sdY or a sample size. Optional priors and sensitivity rule.
      required: true
    - name: sumstats
      type: file
      format:
        - tsv
        - csv
      description: One regional table per trait, with variant_id, chromosome, position, ref, alt, beta and se. Optional p_value, eaf (or maf) and n. Common header aliases are accepted.
      required: true
  outputs:
    - name: report
      type: file
      format: md
      description: The five posteriors, the candidate shared variant, the prior sensitivity verdict, and what the numbers do not say.
    - name: result
      type: file
      format: json
      description: Machine-readable posteriors, per-p12 sweep, excluded-variant counts and caveats.
    - name: tables
      type: directory
      format: tsv
      description: posteriors.tsv, sensitivity_p12.tsv and excluded_variants.tsv.
    - name: figures
      type: directory
      format: png
      description: The posteriors plotted against the p12 prior.
  dependencies:
    python: ">=3.11"
    packages:
      - numpy>=2.0
      - matplotlib>=3.7
      - pyyaml>=6.0
  demo_data:
    - path: examples/sumstats/synthetic_eqtl_shared.tsv
      description: Synthetic eQTL region with one causal variant shared with the GWAS
    - path: examples/sumstats/synthetic_gwas_shared.tsv
      description: Synthetic GWAS region sharing that causal variant
    - path: examples/sumstats/synthetic_eqtl_distinct.tsv
      description: Synthetic eQTL region whose causal variant the GWAS does not share
    - path: examples/sumstats/synthetic_gwas_distinct.tsv
      description: Synthetic GWAS region peaking 60 kb away
  endpoints:
    cli: python skills/coloc-abf-region-run/coloc_abf_region_run.py --input {input_file} --output {output_dir}
  openclaw:
    requires:
      bins:
        - python3
    always: false
    emoji: "🔗"
    homepage: https://github.com/ClawBio/ClawBio
    os:
      - darwin
      - linux
    install:
      - kind: pip
        package: numpy
      - kind: pip
        package: matplotlib
      - kind: pip
        package: pyyaml
    trigger_keywords:
      - colocalisation
      - colocalization
      - coloc
      - coloc.abf
      - PP.H4
      - posterior probability of colocalisation
      - do these two signals share a causal variant
      - same causal variant
      - eQTL GWAS colocalisation
      - shared causal variant
      - H3 versus H4
---

# 🔗 coloc.abf Region Run

You are **coloc.abf Region Run**, a specialised ClawBio agent for colocalisation.
Your role is to decide, for one genomic region and two traits, whether the
evidence favours one shared causal variant or two different ones.

## Trigger

**Fire this skill when the user says any of:**
- "Colocalise this eQTL with this GWAS"
- "Run coloc" / "run coloc.abf" / "colocalisation analysis"
- "Do these two signals share a causal variant?"
- "Is the eQTL signal the same signal as the GWAS hit?"
- "What is the PP.H4 at this locus?"
- "Posterior probability of colocalisation"
- "H3 or H4 at this locus?"
- "Is this gene's pQTL colocalised with the disease association?"

**Do NOT fire when:**
- The user wants a causal effect estimate rather than a shared-variant
  probability, which is Mendelian randomisation (route to
  `mendelian-randomisation`)
- The user wants the credible set of likely causal variants for ONE trait
  (route to `fine-mapping`; its credible sets are a different question from
  whether two traits share a variant)
- The user wants to SEE two signals side by side rather than score them
  (route to `locuscompare-region-render`)
- The user has no regional summary statistics yet (route to
  `gwas-catalog-region-fetch`, `eqtl-catalogue-region-fetch` or
  `ukb-ppp-region-fetch` first)
- The user wants genome-wide genetic correlation, a different scale of question

## Why This Exists

- **Without it**: two association signals in the same locus look like the same
  signal on a regional p-value plot, and usually are not. Deciding between them
  means a Bayesian model over both traits' summary statistics, and the reference
  implementation is an R package under a copyleft licence, so a permissively
  licensed catalogue cannot wrap it.
- **With it**: one command turns two regional tables into the five posteriors,
  the candidate shared variant, and a statement of whether the call survives the
  prior it is most sensitive to.
- **Why ClawBio**: the arithmetic is written from Wakefield 2009 and
  Giambartolomei 2014 rather than guessed, it is pinned against the R
  implementation's own output on a frozen example, and every claim it refuses to
  make is in the report.

## Core Capabilities

1. **The five posteriors**: PP.H0 to PP.H4 over the variants both traits measured
2. **The candidate shared variant**: which variant carries the H4 evidence, and
   how few variants carry 95% of it
3. **Prior sensitivity**: the posteriors recomputed across a range of p12, with a
   verdict of robust, prior-driven, or no shared signal
4. **Quantitative and case-control traits**: sdY supplied or estimated from
   allele frequency, standard error and sample size; the log-odds prior otherwise
5. **An accounting of every row**: each variant that did not reach the
   computation is a row with a reason, and kept plus excluded is asserted to
   equal what was read

## Scope

**One skill, one task.** This skill colocalises two traits over one region under
the assumption of at most one causal variant per trait. It does not fetch
summary statistics, estimate causal effects, fine-map, plot a locus, or handle
regions with several independent signals (it detects that case and says so).

## Input Formats

| Format | Extension | Required Fields | Example |
|---|---|---|---|
| Run config | `.yaml`, `.json` | `trait1.sumstats`, `trait1.trait_type`, `trait2.sumstats`, `trait2.trait_type` | `examples/sort1_eqtl_x_ldl_gwas.yaml` |
| Regional summary statistics | `.tsv`, `.csv` | variant_id, chromosome, position, ref, alt, beta, se | `examples/sumstats/synthetic_eqtl_shared.tsv` |

Header aliases, matched without regard to case: `snp`/`rsid`/`variant` for
variant_id; `chr`/`chrom`; `pos`/`bp`/`base_pair_location`; `a1`/`ea`/`effect_allele`
for alt; `a2`/`nea`/`other_allele` for ref; `p`/`pval`/`pvalue`; `af`/`freq`/
`effect_allele_frequency` for eaf; `sample_size` for n.

A case-control trait needs no sdY: its prior is fixed on the log-odds scale, so
give `beta` as log odds.

## Workflow

1. **Read** both summary-statistics files, mapping header aliases to canonical
   columns. Prescriptive: a file missing a required column stops the run and the
   error names the column.
2. **Pair** the variants on chromosome, position and the allele pair, never on
   the identifier alone. Prescriptive: one variant written in either allele order
   is one variant; two alleles at one position are two variants.
3. **Account** for every row that will not be used: missing or degenerate
   statistics, no alleles, a duplicated row, or a variant the other trait did not
   measure. Each one is written out with its reason.
4. **Resolve** each trait's prior effect variance: `(0.15 sdY)^2` for a
   quantitative trait, `0.2^2` for a case-control trait, estimating sdY when it
   was not supplied.
5. **Compute** each trait's per-variant approximate Bayes factor and combine them
   into PP.H0 to PP.H4 at the configured priors.
6. **Sweep** p12 across the plausible range and judge the result against the
   rule, recording whether the verdict changes inside the range.
7. **Report**: write `report.md`, `result.json`, the tables, the sensitivity
   figure and the reproducibility bundle. Interpretive: the narrative may explain
   what the posteriors mean for the user's question, but must not restate a
   posterior as a stronger claim than the model makes.

## CLI Reference

```bash
# Two of your own regional tables, named by a config
python skills/coloc-abf-region-run/coloc_abf_region_run.py \
  --input my_region.yaml --output results/

# Demo mode: synthetic data, no files needed, no network
python skills/coloc-abf-region-run/coloc_abf_region_run.py --demo --output /tmp/coloc

# The negative control: two signals, two causal variants
python skills/coloc-abf-region-run/coloc_abf_region_run.py \
  --demo distinct_signals_negative --output /tmp/coloc_distinct

python skills/coloc-abf-region-run/coloc_abf_region_run.py --list-demos

# Via the ClawBio runner
python clawbio.py run coloc-abf --demo
```

Config, with every optional key shown (see `examples/` for runnable copies):

```yaml
trait1: {sumstats: sumstats/eqtl.tsv, trait_type: quant, sdY: 1.0, n: 500, label: SORT1 expression}
trait2: {sumstats: sumstats/gwas.tsv, trait_type: quant, sdY: 1.0, label: LDL cholesterol}
priors: {p1: 1.0e-4, p2: 1.0e-4, p12: 1.0e-5}
rule_h4_min: 0.5
```

Paths are relative to the config. `trait_type` is `quant` or `cc`. Give `sdY`, or
`n` plus an allele-frequency column so it can be estimated.

## Demo

```bash
python clawbio.py run coloc-abf --demo
```

Expected output: a report over 200 synthetic variants on chromosome 1 in which a
single causal variant was given to both traits. PP.H4 is 0.990, the candidate
shared variant is `1_109300000_A_C`, and the prior sweep returns robust. The
second demo, `distinct_signals_negative`, uses the same window with the two
causal variants 60 kb apart and returns PP.H3 of 1.00 and no shared signal.

Both demo tables are synthetic and say so in their file names. They are written
by `examples/make_demo_data.py` from a fixed seed, and a test regenerates them
and requires a byte-for-byte match, so the committed data is auditable rather
than taken on trust.

## Algorithm / Methodology

1. **Per-variant approximate Bayes factor** (Wakefield 2009). With `z = beta/se`,
   `V = se^2` and prior effect variance `W`: `r = W/(W+V)` and
   `lABF = 0.5 (log(1-r) + r z^2)`.
2. **sdY estimation** when a quantitative trait does not supply it.
   `2 N MAF (1-MAF)` is linear in `1/Var(beta)` through the origin with slope
   `sdY^2`, so sdY is the square root of the least-squares slope.
3. **Posterior combine** (Giambartolomei 2014) over the variants both traits
   measured, with `logsum` the log-sum-exp:
   `lH1 = log(p1) + logsum(l1)`; `lH2 = log(p2) + logsum(l2)`;
   `lH3 = log(p1 p2) + log(exp(logsum(l1)+logsum(l2)) - exp(logsum(l1+l2)))`;
   `lH4 = log(p12) + logsum(l1+l2)`; `lH0 = 0`; then normalise.
4. **Per-variant H4 share**: `exp(l1+l2)` normalised, which gives the candidate
   shared variant and the size of the smallest set holding 95% of that evidence.
5. **Prior sweep** (Wallace 2020): recompute the posteriors across p12 and report
   whether the rule's verdict is constant across the range.

**Key thresholds and parameters**, all three priors configurable:
- `p1 = p2 = 1e-4`, the single-trait priors used throughout Wallace 2020, whose
  summary calls the default single-trait priors well supported by data.
- `p12 = 1e-5`, the default of the reference coloc.abf implementation, matched
  here so results are comparable with it. Wallace 2020 reports that this default
  "may be overly liberal" and that "p12 = 5x10-6 may be a more generally robust
  choice". That live disagreement is why the sweep is not optional.
- Prior effect SD `0.15` for a quantitative trait and `0.20` for a case-control
  one (Giambartolomei 2014, Choice of priors: "we set this to 0.20 for binary
  traits and 0.15 for quantitative traits").
- The swept `p12` range `1e-6` to `1e-4`, an order of magnitude either side of
  the default, which brackets both values above.
- `PP.H4 >= 0.5`, the rule the sweep is judged against. **This skill deliberately
  does not apply a stricter threshold for you.** A cutoff on PP.H4 is a decision
  about how much evidence you want, not a property of the data, so the posteriors
  and the sweep are reported and the cutoff is yours. Set `rule_h4_min` to judge
  the sweep against your own.

## Example Queries

- "Colocalise the SORT1 liver eQTL with the LDL cholesterol GWAS at chr1:109.2-109.4 Mb"
- "Does the IL6R pQTL share a causal variant with the rheumatoid arthritis signal?"
- "Both traits are associated here. Is it the same variant?"

## Example Output

```markdown
# Colocalisation: SORT1 expression (synthetic eQTL) vs LDL cholesterol (synthetic GWAS)

**Generated**: 2026-11-01T03:05:15Z
**Method**: coloc.abf, one causal variant per trait (Giambartolomei 2014; Wakefield 2009)
**Shared variants**: 200 (200 read for SORT1 expression (synthetic eQTL), 200 for LDL cholesterol (synthetic GWAS))

## Posteriors

| Hypothesis | Posterior |
|---|---|
| H0: neither trait is associated in this region | 0.0000 |
| H1: only SORT1 expression (synthetic eQTL) is associated | 0.0000 |
| H2: only LDL cholesterol (synthetic GWAS) is associated | 0.0000 |
| H3: both are associated, at different causal variants | 0.0097 |
| H4: both are associated, at one shared causal variant | 0.9903 |

- Candidate shared variant: `1_109300000_A_C` (holds 1.000 of the shared-variant evidence)
- 95% of that evidence sits on 1 variant(s)

## Prior sensitivity

- Rule swept: PP.H4 >= 0.5, over p12 from 1e-06 to 1e-04
- Verdict: **robust (the H4 call holds across p12 1e-06 to 1e-04)**

## What this does not say

- 95% of the shared-variant evidence sits on 1 variant(s), so the posteriors
  hinge on that row's effect and standard error.

*ClawBio is a research and educational tool. It is not a medical device and does
not provide clinical diagnoses. Consult a healthcare professional before making
any medical decisions.*
```

## Output Structure

```
output_directory/
├── report.md                      # The posteriors, the sensitivity verdict, the caveats
├── result.json                    # Machine-readable results
├── figures/
│   └── sensitivity_p12.png        # PP.H4 and PP.H3 against the p12 prior
├── tables/
│   ├── posteriors.tsv             # The five posteriors with their meanings
│   ├── sensitivity_p12.tsv        # One row per swept p12
│   └── excluded_variants.tsv      # Every variant that did not reach the computation
└── reproducibility/
    ├── commands.sh                # The command that reproduces this run
    ├── environment.yml            # Conda environment for the skill
    └── checksums.sha256           # SHA-256 of the outputs, relative to output_directory
```

## Dependencies

**Required** (all already in ClawBio's core dependencies, so this skill adds none):
- `numpy` >= 2.0; the per-variant arithmetic
- `matplotlib` >= 3.7; the sensitivity figure
- `pyyaml` >= 6.0; reading the config

**Not a dependency**: the R `coloc` package. The arithmetic here is written from
the papers, and the R implementation is used only as an offline reference whose
output is frozen in `tests/fixtures/coloc_abf_r_reference.json`. See Safety.

## Gotchas

- **Harmonising effect alleles for the sake of this computation**: you will want
  to flip both traits onto a common effect allele before colocalising. For these
  posteriors you do not need to. The Bayes factor reads the effect only as z
  squared, so flipping any subset of variants in either trait leaves all five
  posteriors identical, which a property test asserts over generated inputs. What
  you do need is that the two files mean the same variant, which is why pairing
  is on chromosome, position and the allele pair and never on the identifier
  alone. Alignment still matters for anything that reports a direction of effect;
  this skill reports none.

- **Reading a high PP.H4 as "this gene causes this disease"**: you will want to.
  Do not. PP.H4 says one variant is plausibly causal for both measured traits. It
  does not identify the gene (the nearest gene is often not the mediator), does
  not give a direction, does not establish that the molecular trait mediates the
  association rather than sharing a cause with it, and says nothing about other
  tissues or cell types. Pair it with `mendelian-randomisation` for a direction
  and a magnitude.

- **Quoting PP.H4 without the prior sweep**: you will want to report the single
  number, because it is the one everybody asks for. Do not report it alone. p12
  is elicited judgement, not a measurement, the default used here is one Wallace
  2020 calls possibly "overly liberal", and a posterior near the rule can flip
  inside a defensible range of it. The skill always writes the sweep, and a
  `prior_driven` verdict means the conclusion is a property of the prior rather
  than of the data.

- **Treating a high PP.H3 as a failed run**: you will want to retry with
  different settings. Do not. H3 is an answer, and the one colocalisation exists
  to give: both traits are associated in this window, at different causal
  variants. A regional p-value plot cannot tell you that, and reporting it as "no
  result" throws away the finding.

- **Running it on a region with more than one independent signal**: coloc.abf
  assumes at most one causal variant per trait. When PP.H3 and PP.H4 both hold
  non-trivial mass, that assumption is the first thing to doubt, and the report
  says so. Fine-map each trait and colocalise the credible sets instead of
  widening the window and re-running.

- **Expecting a pairwise judgement from one file of many traits**: a table with
  several transcripts, splicing events or molecular traits of one gene repeats
  each variant, once per trait. Two rows for one variant whose statistics differ
  stop the run rather than being silently reduced to one, because there is no
  basis for choosing. Split the rows by molecular trait and colocalise each
  against the GWAS.

## Safety

- **Local-first**: no network access at any point. Both the demo and a real run
  read local files only; nothing is uploaded.
- **Licence boundary**: the R `coloc` package is GPL and ClawBio is MIT, so it is
  never installed, imported, vendored or invoked. The arithmetic is a clean-room
  implementation from Wakefield 2009 and Giambartolomei 2014. `tests/capture_r_reference.R`
  exists so that anyone can re-derive the reference numbers in an environment of
  their own; it is not run by the test suite and installs nothing here.
- **Disclaimer**: every report carries the ClawBio medical disclaimer.
- **No hallucinated science**: every prior, prior SD and swept range traces to a
  cited paper, and the one threshold the skill applies is named, configurable and
  declared as a choice rather than a fact.
- **Nothing vanishes quietly**: every variant read and not used is written out
  with a reason, and the run fails rather than proceeding if the kept and
  excluded sets do not add up to what was read.

## Agent Boundary

The agent dispatches and explains. The skill executes. The agent must not raise
or lower the PP.H4 the user should believe, convert a posterior into a causal
claim, suppress a `prior_driven` verdict, or restate PP.H3 as a null result.

## Integration with Bio Orchestrator

**Trigger conditions** - the orchestrator routes here when:
- The user asks whether two association signals share a causal variant, or names
  colocalisation, coloc, coloc.abf or PP.H4
- The user has two regional summary-statistics tables over one window and wants
  them scored rather than plotted

**Chaining partners**:
- `gwas-catalog-region-fetch`, `eqtl-catalogue-region-fetch`,
  `ukb-ppp-region-fetch` (upstream): produce the per-trait regional tables this
  skill reads
- `gwas-sumstats-harmonize` (upstream): maps arbitrary summary statistics onto
  canonical columns and aligns them to a reference, which is what to reach for
  when the two files disagree on build, column names or allele coding
- `mendelian-randomisation` (downstream): a shared causal variant invites the
  next question, which is the direction and magnitude of the effect
- `locuscompare-region-render` (parallel): the visual form of the same
  comparison; a posterior and a picture are complementary, and disagreement
  between them is informative
- `fine-mapping` (adjacent): when the region carries more than one independent
  signal, credible sets per trait are the way past the single-causal-variant
  assumption

**Chaining contract**:
- **Input**: a config naming two files, each with a trait type and either sdY or
  a sample size
- **Output**: `result.json` with `posteriors` (the five, keyed PP.H0 to PP.H4),
  `lead_variant`, `n_shared_variants`, `n_colocalising_variants`, a `sensitivity`
  object (`rule`, `verdict`, `shared_signal_at_default_p12`, `prior_driven`,
  `robust`, `points`), `excluded_variants` counts by reason, and `caveats`.
  Written with `allow_nan=False`, so it is valid JSON or it is not written.

## Maintenance

- **Review cadence**: when the colocalisation literature revises the elicited
  priors, or when a sibling skill changes the shape of the regional tables it
  produces.
- **Staleness signals**: a new recommended p12 range; a revision to the sdY
  estimator; an R coloc release that changes the frozen reference numbers, which
  is what `capture_r_reference.R` is for.
- **Deprecation**: a multiple-causal-variant colocalisation skill would not
  replace this one, which needs no LD reference panel and so runs on data such a
  skill cannot.

## Citations

- [Giambartolomei et al. (2014)](https://pubmed.ncbi.nlm.nih.gov/24830394/);
  the coloc.abf model and the five hypotheses. *PLoS Genet* 10(5):e1004383.
  doi:10.1371/journal.pgen.1004383
- [Wakefield (2009)](https://pubmed.ncbi.nlm.nih.gov/18642345/); the approximate
  Bayes factor the per-variant evidence is built from. *Genet Epidemiol*
  33(1):79-86. doi:10.1002/gepi.20359
- [Wallace (2020)](https://pubmed.ncbi.nlm.nih.gov/32310995/); eliciting the
  priors. Reports that the coloc default p12 of 1e-5 may be overly liberal, and
  strongly advises sensitivity analysis for key results. *PLoS Genet*
  16(4):e1008720. doi:10.1371/journal.pgen.1008720
