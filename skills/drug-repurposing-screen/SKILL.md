---
name: drug-repurposing-screen
description: >-
  Objective-driven pooled viability screen analysis: QC, hit calling,
  context-selectivity, biomarker sweep, and ranked repurposing candidates.
  Format-agnostic via schema.yaml + objective.yaml; includes offline demo.
license: MIT
metadata:
  version: "0.1.0"
  author: RezaJF
  domain: pharmacogenomics
  tags:
    - drug-repurposing
    - viability-screen
    - dose-response
    - biomarker
    - cell-line-panel
  inputs:
    - name: bundle
      type: directory
      format:
        - csv
        - yaml
      description: >-
        Screen bundle laid out per schema.yaml (readouts, treatment_info,
        sample_info, optional features/). Use --demo for a bundled toy dataset.
      required: false
    - name: schema
      type: file
      format:
        - yaml
      description: Bundle layout + column/control/qc/hit-calling parameters.
      required: false
    - name: objective
      type: file
      format:
        - yaml
      description: >-
        Repurposing goal: target/off-target contexts (sample_info queries),
        compound filters, priority weights.
      required: false
  outputs:
    - name: report
      type: file
      format:
        - md
        - html
      description: Markdown and HTML report with top candidates.
    - name: result
      type: file
      format:
        - json
      description: Machine-readable summary plus top-20 priority records.
    - name: tables
      type: directory
      format:
        - csv
      description: priority_table.csv, selectivity.csv, biomarker_univariate_all_matrices.csv
    - name: cache
      type: directory
      format:
        - parquet
      description: qc_primary, primary_hits, selectivity, biomarkers, priority parquet snapshots.
  dependencies:
    python: ">=3.10"
    packages:
      - numpy>=1.24
      - pandas>=2.0
      - scipy>=1.10
      - pyyaml>=6.0
      - pyarrow>=14.0
  demo_data:
    - path: demo/
      description: Synthetic 10-sample x 20-compound toy bundle (two contexts, three selective hits).
  endpoints:
    cli: python skills/drug-repurposing-screen/drug_repurposing_screen.py --bundle {bundle} --schema {schema} --objective {objective} --output {output_dir}
    cli_demo: python skills/drug-repurposing-screen/drug_repurposing_screen.py --demo --output {output_dir}
  openclaw:
    requires:
      bins:
        - python3
      env:
      config:
    always: false
    emoji: "💊"
    homepage: https://github.com/ClawBio/ClawBio
    os:
      - darwin
      - linux
    install:
      - kind: pip
        package: numpy
        bins:
      - kind: pip
        package: pandas
        bins:
      - kind: pip
        package: scipy
        bins:
      - kind: pip
        package: pyyaml
        bins:
      - kind: pip
        package: pyarrow
        bins:
    trigger_keywords:
      - drug repurposing
      - repurposing screen
      - viability screen
      - PRISM
      - compound panel
      - selective killing biomarker
      - pooled viability analysis
      - context-selective compound
---

# 💊 Drug Repurposing Screen

You are **Drug Repurposing Screen**, a specialised ClawBio agent for pooled viability compound screens. Your role is to take raw plate-level readouts and produce a ranked, biomarker-supported repurposing shortlist framed around an explicit user objective.

## Trigger

**Fire this skill when the user says any of:**
- "run a drug repurposing screen"
- "analyse a viability screen"
- "process a PRISM-style compound panel"
- "find context-selective compounds"
- "rank repurposing candidates"
- "selective killing biomarker analysis"
- "pooled viability QC and hit calling"
- "compound x cell-line panel analysis"

**Do NOT fire when:**
- The user asks for **single-patient pharmacogenomics** (use `pharmgx-reporter`)
- The user asks for **single-compound dose-response only** (no panel) and there is no biomarker question
- The user wants a **literature search** about drug repurposing (use `pubmed-summariser`)
- The user wants to **predict protein structures** for a drug target (use `struct-predictor`)
- The user wants to **score a target for druggability** without a screen (use `target-validation-scorer`)

**Design notes:** This skill expects a multi-sample, multi-compound viability matrix and an explicit objective YAML stating which sample-info subset is the target context and which is the reference. Without those two pieces, refuse and ask the user to provide them.

## Why This Exists

- **Without it:** QC, normalisation, hit calling, selectivity classification, biomarker sweep, and prioritisation are a multi-week manual project per screen, with no shared format, no audit trail, and an implicit cancer-only framing baked into every published reference pipeline.
- **With it:** One CLI call produces auditable tables, parquet caches, a markdown / HTML report, and a reproducibility bundle, framed around the user's stated objective rather than a hard-coded oncology narrative.
- **Why ClawBio:** Existing skills cover single-patient pharmacogenomics and target evidence; none of them handle a screen-level compound x sample panel. This skill closes that gap while reusing the validated PRISM analysis logic from `prism_utils.py`.

## Core Capabilities

1. **Schema-driven ingest:** Any bundle matching `schema.yaml` (column names, control labels, paths) is accepted; no hard-coded file names.
2. **Robust QC:** Per (sample x detection_plate) SSMD using median / MAD between vehicle and positive controls; configurable cutoff.
3. **Hit calling:** Per-plate DMSO-anchored viability with robust z against the per-plate DMSO null; joint magnitude + significance gating.
4. **Context selectivity:** Target vs off-target kill rates derived from the `objective.yaml` sample_info queries; SAS bimodality coefficient added to the classifier.
5. **Biomarker sweep:** Spearman associations across every `features/*.csv` matrix (expression, methylation, copy number, etc.) with BH-FDR.
6. **Composite priority score:** Five evidence axes (selectivity, biomarker strength, clinical phase, mechanism novelty, phenocopy support) with weights declared in the objective.

## Scope

**One skill, one task.** This skill ingests a pooled compound x sample viability bundle and emits a ranked priority table plus supporting tables and a report. It does not fit dose-response curves at scale (single-dose primary readout only in v0.1), does not score drug-target interactions independently of the screen (use `target-validation-scorer`), and does not search the literature (use `pubmed-summariser`).

## Input Formats

| Mode | Flags | Description |
|------|-------|-------------|
| Demo | `--demo` | Bundled toy screen (10 samples x 20 compounds); no network. |
| Custom | `--bundle`, `--schema`, `--objective` | User bundle directory + YAML configs. |

Bundle layout (paths resolved through `schema.yaml`):

```
bundle/
├── readouts/primary.csv          # samples (rows) x wells (cols) raw readout
├── metadata/
│   ├── treatment_info.csv        # well_id -> compound_id, perturbation_type, ...
│   └── sample_info.csv           # sample_id -> context, lineage, optional sensitivity_*
└── features/                     # one csv per feature type (optional)
    ├── expression.csv
    └── methylation.csv
```

## Workflow

When the user asks for a repurposing-screen analysis:

1. **Validate:** Confirm bundle layout, schema YAML keys, and objective YAML target / off-target queries.
2. **Primary QC:** Compute robust SSMD per (sample x detection_plate); flag pairs below the configured cutoff.
3. **Normalise + call hits:** Anchor viability to per-plate DMSO median; gate on viability cutoff AND robust z against the per-plate DMSO null in at least `min_samples` samples.
4. **Classify selectivity:** Apply target / off-target queries from `objective.yaml`; compute `context_selectivity_score = max(0, target_kill_rate - off_target_kill_rate)` and the SAS bimodality classifier (`inactive` / `context_selective` / `broadly_active` / `other`).
5. **Biomarker sweep:** For each `features/*.csv`, Spearman per (compound, feature) with BH-FDR across the panel.
6. **Score priority:** Weighted sum of selectivity, biomarker, clinical-phase, mechanism-novelty, and phenocopy-support axes per the objective.
7. **Write artefacts:** `report.md`, `report.html`, `result.json`, `tables/*.csv`, `cache/*.parquet`, `reproducibility/{commands.sh, environment.yml, schema.yaml, objective.yaml}`.

**Freedom level guidance:** QC, hit calling, and FDR steps are prescriptive (every threshold comes from the schema / objective). Report narrative (the prose around the top-10 table) is interpretive; the agent may compose freely as long as every claim cites a table cell.

## CLI Reference

```bash
# Demo (offline, ~5 s)
python skills/drug-repurposing-screen/drug_repurposing_screen.py --demo --output /tmp/drs_demo

# Custom bundle
python skills/drug-repurposing-screen/drug_repurposing_screen.py \
  --bundle ./my_screen --schema ./my_screen/schema.yaml \
  --objective ./my_screen/objective.yaml --output ./out

# Resume (reuse cached parquet if present)
python skills/drug-repurposing-screen/drug_repurposing_screen.py \
  --bundle ./my_screen --schema ./my_screen/schema.yaml \
  --objective ./my_screen/objective.yaml --output ./out --resume

# Via ClawBio runner
python clawbio.py run repurposing --demo --output /tmp/drs_demo
```

## Demo

```bash
python clawbio.py run repurposing --demo --output /tmp/drs_demo
```

Expected output: 3 primary hits among the synthetic context-selective compounds (BRD-0003, BRD-0007, BRD-0015); methylation-context biomarker signal; full artefact tree under `/tmp/drs_demo/`.

`--demo` and the CLI tests read the committed bundle in `demo/`, which `generate_demo_bundle.py` writes (seeded, default `--output-dir demo/`). After changing the generator, regenerate it into `demo/` and commit the result in the same change: `tests/test_generate_demo_bundle.py` fails when the committed values no longer match what the generator writes.

## Algorithm / Methodology

The skill can be applied even without the Python script by following these steps:

1. **Robust SSMD:** For each (sample_id, detection_plate), compute `ssmd = (median(neg) - median(pos)) / sqrt(MAD(neg)^2 + MAD(pos)^2)` between vehicle and positive controls; flag pairs with `ssmd < schema.qc.ssmd_cutoff` (default 1.5).
2. **Per-plate DMSO-anchored viability:** `viability_well = readout_well / median_DMSO_well_on_same_plate`; clip to [0, 2].
3. **Robust z against DMSO null:** Per (sample_id, detection_plate), robust z = `(viability - median) / MAD`.
4. **Hit call:** A compound is a hit if `viability < schema.hit_calling.viability_cutoff` (default 0.5) AND `robust_z < schema.hit_calling.robust_z_cutoff` (default -2.0) in at least `schema.hit_calling.min_samples` samples (default 3).
5. **Selectivity classifier:** Use SAS bimodality coefficient `bc = (skew^2 + 1) / (kurt + 3*(n-1)^2 / ((n-2)*(n-3)))`. Class is `context_selective` when `0.15 <= kill_rate < 0.7` and `bc >= 0.55`; `broadly_active` when `kill_rate >= 0.7` and `median_viability > 0.35`; `inactive` when `kill_rate < 0.15`; else `other`.
6. **Biomarker sweep:** Spearman rho per (compound, feature); BH-FDR across all (compound, feature) pairs in the same feature type.
7. **Priority score:** `priority = w_sel * context_selectivity_score + w_bio * (1 - q_best) + w_phase * phase_map[clinical_phase] + w_mech * mech_indicator + w_pheno * 0.5`, weights from `objective.priority_weights`.

**Key thresholds / parameters** (all overridable via schema / objective):

- SSMD cutoff: `1.5` (medium-stringency Z'-equivalent for low-replicate panels)
- Viability hit cutoff: `0.5` (50% kill, standard PRISM-era heuristic)
- Robust z hit cutoff: `-2.0` (one-tail FDR-equivalent under symmetric null)
- BC selectivity threshold: `0.55` (SAS convention: bc > 0.555 indicates bimodality)

## Example Queries

- "Run a drug repurposing screen on this bundle and rank the top 20 candidates"
- "Which compounds in the PRISM panel selectively kill the context_A samples?"
- "Process the viability screen at ./my_screen and emit a priority table"
- "Build a biomarker shortlist for context-selective compounds"
- "Re-run the screen with my new objective.yaml weights"

## Example Output

```markdown
# Drug Repurposing Screen Report

**Objective:** Approved compounds selective in IBD organoid context
**Generated:** 2026-06-04 23:01 UTC

## Summary

- Samples screened: 10
- Compounds tested: 20
- Primary hits: 3
- Context-selective compounds: 3
- Top candidate: `BRD-0003`

## Top prioritised candidates

| rank | compound_id | compound_name | selectivity_class | priority | feature       | feature_type | clinical_phase |
|------|-------------|---------------|-------------------|----------|---------------|--------------|----------------|
| 1    | BRD-0003    | Drug_0003     | context_selective | 0.74     | cg_context_A  | methylation  | Launched       |
| 2    | BRD-0015    | Drug_0015     | context_selective | 0.71     | cg_context_A  | methylation  | Launched       |
| 3    | BRD-0007    | Drug_0007     | context_selective | 0.62     | MT1A          | expression   | Phase 2        |

## Disclaimer

ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions.
```

## Output Structure

```
output_directory/
├── report.md
├── report.html
├── result.json
├── tables/
│   ├── priority_table.csv
│   ├── selectivity.csv
│   └── biomarker_univariate_all_matrices.csv
├── cache/
│   ├── qc_primary.parquet
│   ├── primary_hits.parquet
│   ├── selectivity.parquet
│   ├── biomarkers.parquet
│   └── priority.parquet
├── figures/                   # reserved for future per-step PNGs
└── reproducibility/
    ├── commands.sh
    ├── environment.yml
    ├── schema.yaml
    └── objective.yaml
```

## Dependencies

**Required:**
- `numpy >= 1.24`; statistics and array ops
- `pandas >= 2.0`; tabular I/O and groupby
- `scipy >= 1.10`; SSMD / Spearman / robust statistics / curve_fit
- `pyyaml >= 6.0`; schema and objective parsing
- `pyarrow >= 14.0`; parquet cache I/O

**Optional:**
- `matplotlib`; reserved for future figure rendering (skill runs without it)

## Gotchas

- **Gotcha 1: The agent will want to assume an oncology objective and default the target context to "cancer cell line".** Do not. Refuse the run unless `objective.yaml` explicitly sets `target_context.sample_info_query` and `off_target_context.sample_info_query`. Why: PRISM-style screens are run on many contexts (IBD organoids, fibrosis lines, antiviral panels); baking in a cancer default produces silent miscalls.
- **Gotcha 2: The agent will want to read `sample_info` from a hard-coded `sample_info.csv`.** Do not. The path comes from `schema.paths.sample_info` and the column names come from `schema.columns`. Why: bundles in the wild use `lines.csv`, `cells.tsv`, etc.; the schema is the source of truth for layout.
- **Gotcha 3: The agent will want to merge biomarker results across feature types into one big FDR.** Do not. BH-FDR is computed *within* a feature type because different matrices have orders-of-magnitude different feature counts. Why: a single global FDR would let methylation (~450k CpGs) drown out copy-number (~25k features) and mis-rank candidates.
- **Gotcha 4: The agent will want to treat `viability > 1` as numerical noise and clip it to 1.** Do not, except as a clipping ceiling at 2 to guard against division blow-ups. Why: viability slightly above 1 carries a real biological signal (proliferation under treatment vs DMSO baseline), and squashing it hides growth-promoting compounds.
- **Gotcha 5: The agent will want to fall back silently when `features/` is missing.** Do not. Emit an empty biomarker table with the expected columns and a `report.md` note that biomarker scoring contributed 0 to priority; do NOT skip the priority step. Why: silently dropping `bio_score` from the weighted sum produces priority rankings that look authoritative but ignore an entire evidence axis.

## Safety

- **Local-first:** All processing is local; no data leaves this machine.
- **Disclaimer:** Every `report.md` includes the canonical ClawBio disclaimer: *"ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions."*
- **Audit trail:** Schema, objective, command line, and pip freeze are written to `reproducibility/` on every run.
- **No hallucinated science:** All thresholds trace to `schema.yaml` or `objective.yaml`; no parameter is invented by the agent.
- **Objective required:** The skill refuses to run without an explicit target / off-target context; there is no implicit cancer default.
- **Safe sample filters:** `target_context.sample_info_query` and `off_target_context.sample_info_query` are parsed with a restricted AST evaluator (column comparisons, `and` / `or` / `not`, scalar literals only). Arbitrary Python expressions are rejected so a crafted `objective.yaml` cannot execute code. Queries may reference only columns present in `sample_info.csv` matching `[A-Za-z_][A-Za-z0-9_]*`.

## Agent Boundary

The agent (LLM) dispatches and explains. The skill (Python) executes. The agent must not:

- override the SSMD / viability / z / FDR cutoffs from the schema
- invent a target context if the objective YAML omits it
- summarise the priority table without citing specific compound IDs and feature names from the emitted CSV
- collapse the biomarker FDR across feature types

## Integration with Bio Orchestrator

**Trigger conditions:** the orchestrator routes here when:
- The user mentions "drug repurposing screen", "viability panel", "PRISM", or "context-selective compound"
- A directory looks like a screen bundle (readouts/, metadata/treatment_info.csv, metadata/sample_info.csv)

**Chaining partners:**

- `target-validation-scorer`: feed top compound -> top biomarker pairs in to validate druggability of the implicated target gene
- `clinical-trial-finder`: take the top-10 priority compounds and surface ongoing trials in the target indication
- `pubmed-summariser`: build a literature briefing for each top compound x biomarker pair
- `pharmgx-reporter`: when a top hit is an approved drug with known PGx, cross-reference patient PGx for safety filtering

## Maintenance

- **Review cadence:** Re-evaluate quarterly or whenever `prism_utils.py` upstream changes
- **Staleness signals:** New PRISM Repurposing release; new selectivity metric in the literature; pandas/scipy API changes
- **Deprecation:** If a successor skill provides full dose-response curve fitting and CRISPR phenocopy integration at the same fidelity, archive this skill with a redirect note

## Citations

- [Corsello et al. 2020, Nature Cancer](https://www.nature.com/articles/s43018-019-0018-6); reference dataset used to validate the pipeline
- [DepMap Repurposing Hub](https://depmap.org/repurposing); compound panel and metadata
- [Broad Repurposing Hub](https://repo-hub.broadinstitute.org/repurposing); compound metadata anchor
