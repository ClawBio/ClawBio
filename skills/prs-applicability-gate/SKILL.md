---
name: prs-applicability-gate
description: >-
  Deterministic gate for ONE already-computed polygenic score result and ONE person: returns SUPPORTED (raw score,
  standardised score and reference percentile may be shown), RAW_ONLY (raw score only) or ABSTAIN (nothing), with
  reason codes, a full rule trace and what would change the result. Never produces absolute risk.
license: MIT
metadata:
  version: "2.2.0"
  author: Rinat Rizvanov
  domain: genomics
  tags:
    - polygenic-score
    - prs
    - applicability
    - ancestry
    - pgs-catalog
    - safety-gate
  inputs:
    - name: gate_input
      type: file
      format:
        - json
      description: Normalised evidence for one candidate score and one person (schema prs-applicability-gate.input.v2)
      required: true
  outputs:
    - name: report
      type: file
      format:
        - md
      description: Status, allowed claims, full rule trace and what would change the result
    - name: result
      type: file
      format:
        - json
      description: ClawBio result envelope carrying the full deterministic decision (schema prs-applicability-gate.output.v2)
  dependencies:
    python: ">=3.11"
    packages:
      - pyyaml>=6.0
  demo_data:
    - path: examples/synthetic_A_supported.gate_input.json
      description: SYNTHETIC gate input (invented values) that is SUPPORTED
    - path: examples/synthetic_B_raw_only.gate_input.json
      description: SYNTHETIC gate input (invented values) that is RAW_ONLY, intermediate reference placement
    - path: examples/synthetic_C_abstain.gate_input.json
      description: SYNTHETIC gate input (invented values) that is ABSTAIN, sparse genotype and low scoreability
  endpoints:
    cli: python skills/prs-applicability-gate/prs_applicability_gate.py --input {gate_input} --output {output_dir}
  openclaw:
    requires:
      bins:
        - python3
    always: false
    emoji: "🚦"
    homepage: https://github.com/rinatrizvanov/prsguard-demo
    os:
      - darwin
      - linux
    install:
      - kind: pip
        package: pyyaml
    trigger_keywords:
      - can this polygenic score be interpreted
      - is this PRS percentile valid for me
      - PRS applicability
      - prs gate
      - should this PGS percentile be shown
---

# 🚦 PRS Applicability Gate

You are **PRS Applicability Gate**, a specialised ClawBio agent for polygenic-score interpretation safety. Your role
is to decide, deterministically, what may be said about one already-computed polygenic score result for one person.

## Trigger

**Fire this skill when the user says any of:**
- "can this polygenic score / PRS / PGS be interpreted for me (for this person)?"
- "is this PRS percentile valid / applicable / trustworthy for this sample?"
- "should the percentile be shown?", "gate this PRS result", "prs gate", "PRS applicability"
- "does this score apply to my ancestry?" (once a score and a reference placement exist)
- any time a PGS percentile is about to be reported to a person, even if the upstream steps looked clean

**Do NOT fire when:**
- the user wants a score *computed* (`gwas-prs` / `prs`, `just-prs-mcp`, `wgs-prs`)
- the user wants candidate scores *found, ranked or chosen* for a trait (PGS Catalog search is not this skill)
- the user wants ancestry *inferred* (`claw-ancestry-pca`, `ancestry-risk-profiler`)
- the user wants population equity metrics (`equity-scorer`)
- no gate input exists yet: the evidence must be assembled first (see Input Formats); the gate never fetches it

## Why This Exists

- **Without it**: a raw PGS is silently turned into a percentile against a reference population the person may
  not belong to, from a genotype file that may not contain the score's important variants, for a score never
  evaluated in anyone like them. Nothing looks broken; the number is confidently wrong.
- **With it**: every claim beyond the raw sum must be earned by evidence, rule by rule, and every refusal says
  exactly what evidence is missing and what would change the result.
- **Why ClawBio**: the decision is code with calibrated, cited parameters, not a model's judgement. The same input
  always gives the same answer, byte for byte.

## Core Capabilities

1. **Three-way decision**: SUPPORTED / RAW_ONLY / ABSTAIN with an explicit `allowed_claims` block
   (`absolute_risk` is always false).
2. **Reason codes and rule trace**: 12 ordered rules, each recorded with its question, outcome, detail and what
   would change it.
3. **Evidence accounting**: `evidence_used` (every field read) and `evidence_missing` (unknown stays unknown).
4. **Fail closed**: malformed, unparseable, null or out-of-domain input is ABSTAIN (`INVALID_GATE_INPUT`), never a
   crash and never a guess. An optional `input_digest` is re-verified, so an edited input is rejected.
5. **Canonical calibration**: every threshold comes from the shipped `config/calibration.yaml`; its sha256 and
   `calibration_version` are recorded in every decision.

## Scope

**One skill, one task.** Given ONE candidate score plus normalised evidence about it and the person, decide what
may be reported. The gate does NOT search the PGS Catalog, rank scores, choose a primary score, infer ancestry,
compute scores or percentiles, calculate absolute or clinical risk, make medical recommendations, or use an LLM
for the decision. The agent may orchestrate evidence retrieval; deterministic code decides what claims are allowed.

## Input Formats

| Format | Extension | Required Fields | Example |
|--------|-----------|-----------------|---------|
| Gate input v2 | `.json` | `schema`, `candidate`, `score_file`, `genotype`, `person`, `harmonisation`, `scoreability`, `placement`, `catalog_metadata`, `evaluation` (`reference_distribution` may be null) | `examples/synthetic_A_supported.gate_input.json` |

Fields the rules read (others are carried through but ignored). Evidence must be assembled upstream, for example
by PRSGuard's evidence builder; the gate never fetches or infers it.

| Block | Fields | Rule |
|---|---|---|
| `schema` | `"prs-applicability-gate.input.v2"` | G1 |
| `candidate` | `pgs_id` (string), `sex_specific` (`female` / `male` / null) | G1, G6 |
| `score_file` | `weight_type`, `ratio_weight_type`, `unsupported_features` (list), `variants_interactions`, `n_parse_problems` | G2 |
| `genotype` | `build` (null or `UNRESOLVED` fails), `build_evidence.method` | G3 |
| `person` | `sex` (`female` / `male` / null) | G6, G10 |
| `harmonisation` | `n_variants`, `n_matched`, `n_located` (integers >= 0); `allele_mismatch_fraction` in [0, 1], null exactly when `n_located` is 0; `weight_loss_by_status` (status to fraction, largest loss first) | G4, G5 |
| `scoreability` | `r`: correlation, in the placed reference group, of the score computable from this genotype with the complete published score ([-1, 1] or null); `method` | G5 |
| `placement` | `status` (`RESOLVED` / `INTERMEDIATE` / `UNSTABLE` / `UNRESOLVED`); `placement` (`AFR` / `AMR` / `EAS` / `EUR` / `SAS`, required when RESOLVED); `placement_stability` ([0, 1] or null); `detail` | G8 |
| `catalog_metadata` | `status` (`resolved`, `contradictory`, anything else = unavailable), `detail` | G7 |
| `evaluation` | `units`: list of `{code, pooled, n, percent_male, metrics: [{name, ci_lower, ci_upper, null}]}` | G9, G10 |
| `reference_distribution` | `available`, `reference_group`, `reference_n`, `n_intersection`, `reference_sensitive`, `reference_sensitive_pairs` | G11, G12 |
| `input_digest` (optional) | `sha256:` of the canonical JSON (sorted keys, compact) of every other field | G1 |

## Workflow

When the user asks whether a PGS result may be interpreted:

1. **Validate** (prescriptive): schema, types, numeric domains and `input_digest` (G1). Anything invalid is ABSTAIN
   with `INVALID_GATE_INPUT`; the remaining rules are not evaluated.
2. **Evaluate** (prescriptive): rules G2-G12 in order, every one recorded, so the trace of a valid input always
   lists G1-G12.
3. **Decide** (prescriptive): ABSTAIN if any ABSTAIN rule failed, else RAW_ONLY if any RAW_ONLY rule failed, else
   SUPPORTED. The primary reason is the first failing rule at the decisive severity.
4. **Report** (prescriptive): write `report.md`, `result.json` and `reproducibility/`.
5. **Explain** (agent, flexible wording, fixed content): state the status, `allowed_claims`, the primary reason
   and `what_would_change_result` in plain language. Never reinterpret, soften or override the status.

## CLI Reference

```bash
# One gate input
python skills/prs-applicability-gate/prs_applicability_gate.py \
  --input gate_input.json --output <report_dir>

# Demo mode (three SYNTHETIC inputs bundled with the skill)
python skills/prs-applicability-gate/prs_applicability_gate.py --demo --output /tmp/prs_gate_demo

# Via ClawBio runner (always the shipped calibration)
python clawbio.py run prs-gate --input gate_input.json --output <dir>
python clawbio.py run prs-gate --demo
```

Exit status is 0 whenever a decision was written, ABSTAIN included; usage errors (missing input file, unusable
`--config`) exit 2 without writing a decision.

`--config` (an alternative calibration file) exists only on the direct script, for expert calibration research.
Its results carry `config_canonical: false`, a NON-CANONICAL CALIBRATION banner in `report.md` and a warning on
stderr. It is deliberately **not** forwarded by `clawbio.py run prs-gate`, the interface an agent uses: the alias
has an empty `allowed_extra_flags`, so the runner drops the flag, and a test proves the shipped calibration is used.

## Demo

```bash
python clawbio.py run prs-gate --demo
```

Expected output: three decisions over SYNTHETIC gate inputs (invented values, no person, no genotype, no real PGS
Catalog record): `PGS_SYNTHETIC_A` SUPPORTED; `PGS_SYNTHETIC_B` RAW_ONLY (`TARGET_REFERENCE_UNRESOLVED`, an
intermediate placement between two reference clouds); `PGS_SYNTHETIC_C` ABSTAIN (`LOW_SCOREABILITY`,
`VARIANTS_MISSING`, `TARGET_REFERENCE_UNRESOLVED`, a 150-site genotype file).

## Algorithm / Methodology

| Rule | Question | Fails with | Effect |
|---|---|---|---|
| G1 INPUT_VALID | Well-formed, in range, digest intact? | INVALID_GATE_INPUT | ABSTAIN |
| G2 SCORE_FORMAT | Plain additive score on a log scale? | UNSUPPORTED_SCORE_FORMAT | ABSTAIN |
| G3 BUILD | Genotype build established, never assumed? | BUILD_UNRESOLVED | ABSTAIN |
| G4 ALLELES | Located variants' alleles reconcilable (<= `max_mismatch_fraction`)? | ALLELE_HARMONIZATION_FAILED | ABSTAIN |
| G5 SCOREABILITY | r(computable score, published score) >= `r_min`? | LOW_SCOREABILITY (+ PALINDROMIC_VARIANT_UNRESOLVED / VARIANTS_MISSING / DUPLICATE_OR_CONFLICTING_VARIANTS naming the largest loss), SCOREABILITY_UNVERIFIED | ABSTAIN |
| G6 SEX_SCORE | Sex-specific score used for that sex? | SEX_POPULATION_MISMATCH (ABSTAIN), SEX_NOT_PROVIDED (RAW_ONLY) | |
| G7 METADATA | PGS Catalog record resolved and consistent? | METADATA_CONTRADICTION, EVALUATION_METADATA_UNAVAILABLE | RAW_ONLY |
| G8 PLACEMENT | Person stably inside one reference group? | TARGET_REFERENCE_UNRESOLVED | RAW_ONLY |
| G9 EVALUATION | Single-ancestry evaluation in that group with a metric whose 95% CI lies entirely above the null? | NO_RELEVANT_EVALUATION, EVALUATION_NOT_INFORMATIVE | RAW_ONLY |
| G10 SEX_EVALUATION | Do those evaluations include the person's sex? | SEX_POPULATION_MISMATCH | RAW_ONLY |
| G11 REFERENCE_DISTRIBUTION | Reference distribution on the person's matched variants? | REFERENCE_DISTRIBUTION_UNAVAILABLE | RAW_ONLY |
| G12 REFERENCE_SENSITIVITY | Percentile robust across equally defensible reference populations? | REFERENCE_SENSITIVE | RAW_ONLY |

**Key parameters** (`config/calibration.yaml`, `calibration_version` `2026.09.26-3`; derivations and experiments in
[PRSGuard `docs/calibration.md`](https://github.com/rinatrizvanov/prsguard-demo/blob/v1.0.1/docs/calibration.md)):
- `r_min = 0.90`: policy anchored in the published metric (a reduced score keeps about r of the published per-SD
  association), with its consequences for percentile error measured by masking experiments in 1000 Genomes.
- `max_mismatch_fraction = 0.05`: detected allele mismatches estimate a similar rate of undetectable wrong calls,
  which attenuate the score roughly as r ~ 1 - e; every real, correctly built public file tested showed <= 0.4%.
- Evaluation rule: qualitative (at least one single-ancestry evaluation in the placed group whose metric CI lies
  entirely above its null); no sample-size or fraction cut-off. Pooled units and inverse associations never count.
- Placement parameters (`min_sites = 200`, `min_stability = 0.95`, cloud quantile 0.999) belong to the upstream
  placement step and are recorded in the config so that changing them changes `calibration_version`.

**What SUPPORTED means, and does not.** G9 establishes evidence of association in a relevant evaluation group. It
does not establish clinically useful discrimination (an AUROC of 0.55 with a CI above 0.5 passes) or calibration.
SUPPORTED is a research-prototype reportability state: the percentile may be shown with its intervals. It is not a
clinical recommendation and never an absolute risk. The rule assumes a higher score means a higher phenotype value
or risk; the PGS Catalog has no structured direction field, so below-null effects are never accepted as support.

## Example Queries

- "Can the percentile of my breast cancer polygenic score be trusted?"
- "Gate this PGS result for this sample"
- "Is it valid to show a percentile for this admixed genome?"

## Example Output

Excerpt of `report.md` from `python clawbio.py run prs-gate --demo`:

```markdown
# PRS Applicability Gate Report

**Mode**: demo, 3 SYNTHETIC gate inputs bundled with the skill (invented values; no person, no genotype, no real PGS Catalog record)
**Calibration**: `2026.09.26-3` (canonical; config sha256:727942ca…) · **Gate**: 2.2.0

| Input | Score | Status | Primary reason | Reason codes |
|---|---|---|---|---|
| synthetic_A_supported.gate_input.json | PGS_SYNTHETIC_A | **SUPPORTED** | - | - |
| synthetic_B_raw_only.gate_input.json | PGS_SYNTHETIC_B | **RAW_ONLY** | TARGET_REFERENCE_UNRESOLVED | TARGET_REFERENCE_UNRESOLVED |
| synthetic_C_abstain.gate_input.json | PGS_SYNTHETIC_C | **ABSTAIN** | LOW_SCOREABILITY | LOW_SCOREABILITY, VARIANTS_MISSING, TARGET_REFERENCE_UNRESOLVED |

## PRS applicability gate: PGS_SYNTHETIC_B

**Status: RAW_ONLY** (TARGET_REFERENCE_UNRESOLVED)

### What may be reported

- Raw score: yes
- Standardised score and percentile: no
- Absolute risk: never

### Rule trace

| Rule | Question | Outcome | Detail |
|---|---|---|---|
| G5 SCOREABILITY | Does the computable score represent the published score (r >= r_min)? | pass | r = 0.955 (reference_panel_correlation); r_min 0.9; 1152/1200 variants matched; largest loss: missing |
| G8 PLACEMENT | Is the person placed stably inside one reference group? | fail -> RAW_ONLY | INTERMEDIATE: between the AFR and EUR reference clouds (intermediate placement) |
| G9 EVALUATION | Is there evidence of association (95% CI above the null) in an evaluation of the person's group? … | not_applicable | no resolved reference group to match evaluations against |

### What would change the result

- TARGET_REFERENCE_UNRESOLVED: a reference panel that represents this person's genetic background (e.g. admixed references with local-ancestry-aware scoring); not something the gate can relax

*Research software, not a medical device. SUPPORTED is a research-prototype reportability state, not a clinical
recommendation … ClawBio is a research and educational tool. It is not a medical device and does not provide
clinical diagnoses. Consult a healthcare professional before making any medical decisions.*
```

## Output Structure

```
output_directory/
├── report.md              # summary table, then status, allowed claims, rule trace and what would change, per input
├── result.json            # ClawBio envelope: summary plus the full decision per input
└── reproducibility/
    ├── commands.sh        # exact command to reproduce the run
    ├── environment.yml    # pip/conda environment
    └── checksums.sha256   # sha256 of report.md, result.json and the bundle files
```

`result.json` keys: `skill`, `version`, `input_checksum`, `datasets` (sha256 per input file), `summary` (mode,
`synthetic_demo`, `status_counts`, `calibration_version`, `gate_version`, `config_sha256`, `config_canonical`, one
row per input), `data.decisions` (the full decision, schema `prs-applicability-gate.output.v2`: `status`,
`primary_reason`, `reason_codes`, `allowed_claims`, `evidence_used`, `evidence_missing`, `calibration_version`,
`rule_trace`, `what_would_change_result`, `provenance`, `disclaimer`), `chat_summary_lines`, `disclaimer`. No
completion timestamp and no local path is written, so identical runs give byte-identical `report.md` and `result.json`.

## Dependencies

**Required**: `pyyaml` >= 6.0 (calibration config; a ClawBio core dependency) and `clawbio.common` (disclaimer and
reproducibility helpers). Nothing else: the gate makes no network calls and uses no randomness.

## Gotchas

- **You will want to pass `--config`, or edit `config/calibration.yaml`, to turn RAW_ONLY into SUPPORTED.** Do not.
  Thresholds are calibrated and versioned; an alternative config is labelled NON-CANONICAL and is for calibration
  research, not for answering a user.
- **You will want to fill a missing field with a plausible default, or edit a failing value in the gate input.**
  Do not. Missing stays in `evidence_missing` and the affected rule fails; an edited input fails its digest.
- **You will want to re-run with a different reference population to get a percentile.** Do not. Reference choice
  is decided by placement upstream; G12 exists precisely to catch answers that depend on that choice.
- **You will want to call a RAW_ONLY score "roughly average".** Do not. RAW_ONLY means no relative statement is
  supported; the raw sum has no meaning on its own scale.
- **You will want to treat ABSTAIN as a failed run and retry.** Do not. ABSTAIN is a valid answer (exit status 0);
  report it with its reason and what would change it.
- **You will want to treat a pooled multi-ancestry evaluation as evidence for every group in it.** Do not; pooled
  units never count as group-specific evidence (G9).
- **You will want to read "percent of variants matched" as scoreability.** Do not; a few heavily weighted variants
  can matter more than hundreds of small ones. The gate uses r(full, reduced) measured with LD.
- **You will want to treat MAE/NR evaluations, an AUROC without a CI, or an inverse association (CI below the null,
  e.g. a case-only subtype comparison) as supporting evidence.** Do not; they are not.
- **You will want to describe a passed G9 as "the score performs well" or "is clinically useful".** Do not; it shows
  association only.
- **You will want to present a percentile as a risk.** Never. `absolute_risk` is always false.

## Safety

- **Local-first**: the gate reads one JSON file and writes files; no network, and its input holds no genotypes
  (counts, correlations, placement status and catalog metadata only).
- **Disclaimer**: every report and every decision carries the ClawBio medical disclaimer, *"ClawBio is a research
  and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare
  professional before making any medical decisions."*, after the statement that SUPPORTED is not a clinical
  recommendation.
- **Audit trail**: input digest, config sha256, calibration and gate versions in every decision; the
  reproducibility bundle records the exact command and checksums.
- **No hallucinated science**: every parameter traces to `config/calibration.yaml` and its documented derivation.
- **Fail closed**: anything the gate cannot read or validate is ABSTAIN, never a best guess.

## Agent Boundary

The agent (LLM) may gather evidence, call this skill, and explain its output; the skill (Python) decides. The agent
must NOT change thresholds, pass `--config`, edit the config or a gate input, re-run until SUPPORTED, choose a score
or a reference population by the personal result, fill in missing evidence, convert a raw score into clinical or
absolute risk, or override, soften or reinterpret the status.

## Integration with Bio Orchestrator

**Trigger conditions**: the orchestrator routes here when:
- a PGS result for a person exists and a standardised score or percentile is about to be shown
- the user asks whether a PRS percentile applies to them or to a sample

**Chaining partners**: this skill connects with:
- `gwas-prs`, `just-prs-mcp`, `wgs-prs`: upstream; they compute the score. Their outputs are not converted into a gate
  input automatically: an evidence builder must express harmonisation, scoreability and placement in schema v2.
- `claw-ancestry-pca`, `ancestry-risk-profiler`: ancestry context only. The gate needs a placement status against
  1000 Genomes super-populations (RESOLVED / INTERMEDIATE / UNSTABLE / UNRESOLVED); it never infers one.
- `profile-report`: downstream; may show a PRS only within the gate's `allowed_claims`.
- `equity-scorer`: context shown next to the decision, never an input to it.

## Maintenance

- **Review cadence**: re-run the calibration benchmarks when the reference panel, the placement method or the PGS
  Catalog schema changes; any change to `config/calibration.yaml` requires a new `calibration_version`.
- **Staleness signals**: PGS Catalog REST schema change; new evaluation deposits for the scores in use; a new
  reference panel with admixed or finer-grained groups.
- **Deprecation**: archive to `skills/_deprecated/` if a maintained PGS-applicability standard supersedes it.
- **Provenance**: originated in PRSGuard, one of the winning projects of the ClawBio Boston Hackathon (Challenge 3;
  team Rinat Rizvanov, Timur Rizvanov, Takato Honda, Bradley Sheppard). This v2 gate was rewritten, calibrated and
  tested afterwards by Rinat Rizvanov (PRSGuard gate 2.1.1). 2.2.0 adapts it to ClawBio (result envelope,
  reproducibility bundle, synthetic demo) and completes the trace (G10 is recorded when placement is unresolved);
  no threshold or decision changed.

## Citations

- [PGS Catalog](https://www.pgscatalog.org/) (Lambert et al., Nat Genet 2021); score metadata, evaluations,
  harmonised scoring files
- [1000 Genomes Project phase 3](https://www.internationalgenome.org/) (Nature 2015); reference super-populations
- [pgsc_calc](https://github.com/PGScatalog/pgsc_calc) (Lambert et al., Nat Genet 2024); default exclusion of
  strand-ambiguous variants that the gate's inputs follow
- Martin et al., Nat Genet 2019; Privé et al., AJHG 2022; reduced PGS accuracy across ancestries, the reason
  evaluations in the person's group are required
- Brown, Cai & DasGupta, Stat Sci 2001; Jeffreys interval that the upstream percentile step uses for finite-panel
  uncertainty (the source of `reference_sensitive`)
- [PRSGuard](https://github.com/rinatrizvanov/prsguard-demo); the full pipeline, calibration experiments and
  benchmarks this skill comes from
