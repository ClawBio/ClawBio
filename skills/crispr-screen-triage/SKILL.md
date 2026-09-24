---
name: crispr-screen-triage
description: Deterministic CRISPR screen hit ranking from local guide-level count tables
license: MIT
metadata:
  version: "0.1.0"
  author: ClawBio
  domain: functional-genomics
  tags:
    - crispr
    - screen
    - triage
  inputs:
    - name: input_file
      type: file
      format:
        - csv
      description: Guide-level CRISPR count and annotation table
      required: true
  outputs:
    - name: report
      type: file
      format:
        - md
      description: Ranked hit report
    - name: result
      type: file
      format:
        - json
      description: Machine-readable triage results
  dependencies:
    python: ">=3.10"
    packages:
  demo_data:
    - path: demo_screen_counts.csv
      description: Synthetic twelve-guide, six-gene CRISPR screen table
  endpoints:
    cli: python skills/crispr-screen-triage/crispr_screen_triage.py --input {input_file} --output {output_dir}
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
    trigger_keywords:
      - CRISPR screen triage
      - guide count ranking
      - rank CRISPR hits
      - depleted guide screen
---

# CRISPR Screen Triage

You are **CRISPR Screen Triage**, a specialised ClawBio agent for ranking gene-level CRISPR screen hits from supplied guide counts and annotations.

## Trigger

**Fire this skill when the user says any of:**
- "triage CRISPR screen hits"
- "rank guide-level CRISPR counts"
- "rank depleted CRISPR genes"
- "score genes from a knockout screen"
- "which CRISPR hits should I follow up"

**Do NOT fire when:**
- The user asks for variant interpretation.
- The user asks for single-cell clustering.
- The user asks for clinical actionability.

## Why This Exists

- **Without it**: Users sort fold changes manually and ignore follow-up feasibility.
- **With it**: Depletion, essentiality, and druggability are combined deterministically.
- **Why ClawBio**: The score is transparent, local, and reproducible.

## Core Capabilities

1. **Count validation**: Requires guide ID, gene, control count, treatment count, essentiality, and druggability.
2. **Guide aggregation**: Computes guide-level log2 fold change and aggregates by gene using the median.
3. **Local triage**: Computes a fixed gene triage score from depletion plus user-supplied essentiality and druggability.
4. **Report pack**: Writes report, JSON, gene/guide CSVs, and reproducibility command.

## Scope

One skill, one task. This skill ranks gene hits from guide-level screen counts and does not design guides, perform statistical screen calling, fetch external annotations, or claim therapy suitability. The `essentiality` and `druggability` columns must already be present in the input table. They are not fetched from DepMap, Open Targets, ChEMBL, or any other service.

## Input Formats

| Format | Extension | Required Fields | Example |
|--------|-----------|-----------------|---------|
| CSV | `.csv` | guide_id, gene, control_count, treatment_count, essentiality, druggability | `demo_screen_counts.csv` |

`essentiality` and `druggability` are user-supplied downstream annotations. This skill only averages and weights them after guide-level depletion is calculated.

## Workflow

1. **Validate**: Confirm required columns and numeric counts/scores.
2. **Compute**: Calculate guide-level `log2((treatment + 1) / (control + 1))`.
3. **Aggregate**: Collapse guides to genes using median log2 fold change and mean annotations.
4. **Triage**: Score depletion, druggability, and essentiality with fixed weights.
5. **Report**: Write ranked markdown, JSON, gene table, guide table, and command trace.

## CLI Reference

```bash
python skills/crispr-screen-triage/crispr_screen_triage.py --input counts.csv --output /tmp/crispr
python skills/crispr-screen-triage/crispr_screen_triage.py --demo --output /tmp/crispr
python clawbio.py run crispr-triage --demo
```

## Demo

```bash
python clawbio.py run crispr-triage --demo
```

Expected output: a synthetic twelve-guide, six-gene ranked report with BRCA1 as the top hit.

## Algorithm / Methodology

1. **Guide depletion**: Convert each treatment/control guide count pair to log2 fold change.
2. **Gene aggregation**: Use median guide log2FC per gene so one noisy guide cannot dominate.
3. **Score**: `0.55 * max(0, -median_log2FC) + 0.25 * druggability + 0.20 * essentiality`.
4. **Priority**: High requires score >= 1.35 and median log2FC <= -1.0.
5. **Non-goal**: This is not a canonical statistical screen caller. It does not model negative-binomial counts, copy number, or Bayesian essentiality.

## Example Queries

- "Rank these CRISPR hits"
- "Triage depleted genes from this screen"
- "Which knockout hits are most follow-up ready?"

## Example Output

```markdown
# CRISPR Screen Triage Report

| Rank | Gene | Guides | Median log2FC | Priority |
|---:|---|---:|---:|---|
| 1 | BRCA1 | 2 | -2.66 | high |
```

## Output Structure

```
output_directory/
├── report.md
├── result.json
├── tables/
│   ├── triaged_genes.csv
│   └── guide_metrics.csv
└── reproducibility/
    ├── commands.sh
    ├── environment.yml
    └── checksums.sha256
```

## Dependencies

- Python 3.11+ and the standard library only.

## Gotchas

- **Do not treat the priority score as validation**: It is a triage score only.
- **Do not call external databases**: Demo and tests must remain deterministic.
- **Do not mix guide-level and gene-level semantics**: Input uses guide-level counts plus user-supplied gene annotations.

## Safety

- **Local-first**: No external APIs or uploads.
- **Disclaimer**: Every report includes the ClawBio medical disclaimer.
- **Audit trail**: Commands are written to `reproducibility/commands.sh`.

## Agent Boundary

The agent dispatches and explains. The Python skill scores and writes outputs.

## Integration with Bio Orchestrator

**Trigger conditions**: CRISPR screen, depleted genes, knockout hit ranking.

## Chaining Partners

- `target-validation-scorer`: downstream target evidence synthesis.
- `omics-target-evidence-mapper`: cross-omics support for top hits.

## Maintenance

- **Review cadence**: Re-evaluate weights quarterly.
- **Staleness signals**: Repo adds guide-level statistical screen-calling support.
- **Deprecation**: Archive if replaced by a full screen-analysis workflow.

## Author & Attribution

Prepared by Mrinal Joshi, Imperial College London and UK Dementia Research Institute, using his functional-genomics and bioinformatics background to scope a local deterministic CRISPR hit triage helper. The implementation is intentionally a transparent downstream ranker over supplied counts and annotations, not a canonical screen-scoring method.

## Citations

- ClawBio local deterministic triage rules in `crispr_screen_triage.py`; this skill does not claim a method-paper implementation.
