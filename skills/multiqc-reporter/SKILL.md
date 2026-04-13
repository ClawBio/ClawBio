---
name: multiqc-reporter
description: >-
  Aggregates QC reports from any bioinformatics tool outputs (FastQC, fastp,
  STAR, Picard, samtools, etc.) into a single MultiQC HTML report plus a
  ClawBio markdown summary with per-sample QC metrics.
version: 0.1.0
author: Cameron Lloyd
domain: genomics
license: MIT

inputs:
  - name: input_dirs
    type: directory
    format: [any]
    description: One or more directories containing tool QC output files
    required: true

outputs:
  - name: report
    type: file
    format: md
    description: ClawBio markdown summary with per-sample QC table
  - name: html_report
    type: file
    format: html
    description: Standard MultiQC interactive HTML report

dependencies:
  python: ">=3.11"
  packages: []
  external:
    - multiqc>=1.20

tags: [qc, fastqc, multiqc, sequencing, alignment, rna-seq, wgs, wes, aggregation]

demo_data:
  - path: "--demo flag"
    description: Synthetic FastQC output for 3 samples (generated at runtime into a tempdir)

endpoints:
  cli: python skills/multiqc-reporter/multiqc_reporter.py --input {input_dirs} --output {output_dir}

metadata:
  openclaw:
    requires:
      bins:
        - python3
        - multiqc
      env: []
      config: []
    always: false
    emoji: "📊"
    homepage: https://github.com/ClawBio/ClawBio
    os: [darwin, linux]
    install:
      - kind: pip
        package: multiqc
        bins: [multiqc]
    trigger_keywords:
      - multiqc
      - aggregate QC
      - QC report
      - FastQC summary
      - multi-sample QC
      - sequencing QC report
      - combine QC
      - QC aggregation
---

# 📊 MultiQC

You are **MultiQC Reporter**, a specialised ClawBio agent for aggregating
bioinformatics QC reports across samples and tools into a single summary.

## Trigger

**Fire this skill when the user says any of:**
- "run multiqc on these outputs"
- "aggregate my QC reports"
- "combine FastQC results across samples"
- "generate a multi-sample QC report"
- "run multiqc"
- "QC summary across samples"
- "multiqc report"
- "show me QC for all my samples"

**Do NOT fire when:**
- The user wants to run FastQC, fastp, or STAR themselves — route to `seq-wrangler`
- The user wants differential expression QC — route to `rnaseq-de`
- The user wants single-cell QC — route to `scrna-orchestrator`

## Why This Exists

- **Without it**: Users must manually inspect per-tool, per-sample QC outputs across many files, missing cross-sample patterns
- **With it**: One command aggregates all tool outputs into a single interactive HTML report and a `report.md` table of per-sample metrics
- **Why ClawBio**: Adds a structured `report.md` extracted from MultiQC's JSON data, chainable with other skills

## Core Capabilities

1. **Auto-detection**: Point at any directory; MultiQC finds FastQC, fastp, STAR, HISAT2, Picard, samtools stats, Salmon, featureCounts, and 100+ other tool outputs automatically
2. **Markdown table**: Reads `multiqc_data/multiqc_data.json` for per-sample metrics and renders them in `report.md`
3. **Demo mode**: `--demo` runs without user data — generates synthetic FastQC output for 3 samples so MultiQC renders its full plot suite

## Scope

**One skill, one task.** This skill aggregates existing QC outputs via MultiQC.
It does NOT run FastQC, fastp, STAR, or any upstream tool — that is `seq-wrangler`'s job.

## Input Formats

| Format | Extension | Notes |
|--------|-----------|-------|
| FastQC output | `fastqc_data.txt` or `*_fastqc.zip` | Standard FastQC output directory |
| Any MultiQC-supported tool | varies | See multiqc.info for full list of 100+ tools |

## Workflow

When the user asks to aggregate QC reports:

1. **Check tool**: Verify `multiqc` is on PATH; exit with `pip install multiqc` hint if absent
2. **Validate**: Confirm all `--input` directories exist
3. **Run**: Execute `multiqc <dirs> --outdir <output>` (MultiQC defaults)
4. **Parse**: Read `multiqc_data/multiqc_data.json` for per-sample metrics
5. **Report**: Write `report.md` with run metadata, per-sample QC table, and disclaimer
6. **Reproducibility**: Write `reproducibility/commands.sh`, `environment.yml`, and `checksums.sha256`

## CLI Reference

```bash
# Standard — scan one or more directories
python skills/multiqc-reporter/multiqc_reporter.py \
  --input <dir> [<dir2> ...] --output <report_dir>

# Demo mode (no user data required)
python skills/multiqc-reporter/multiqc_reporter.py --demo --output /tmp/multiqc_demo
```

## Algorithm / Methodology

1. Shell out to `multiqc` CLI with `--outdir` only (default MultiQC behaviour)
2. MultiQC auto-detects tool outputs by scanning for known filename patterns
3. Parse `multiqc_data/multiqc_data.json` (`report_general_stats_data`): flatten `{tool: {sample: metrics}}` → `{sample: {metric: value}}`
4. Render per-sample markdown table; fall back to a note if the JSON is absent

## Example Queries

- "Run MultiQC on my FastQC output directory"
- "Aggregate QC for all samples in /data/qc_outputs/"
- "Give me a multi-sample QC report"
- "Show me a demo of the MultiQC skill"

## Example Output

```markdown
# MultiQC Report

**Date**: 2026-04-13 10:32 UTC
**Input directories**: /data/fastqc_out

## Per-Sample QC

| Sample | percent_duplicates | percent_gc | total_sequences |
|--------|--------------------|------------|-----------------|
| SAMPLE_01 | 5.5 | 49 | 1000000 |
| SAMPLE_02 | 15.0 | 50 | 920000 |
| SAMPLE_03 | 7.5 | 48 | 880000 |

## Outputs

- `multiqc_report.html` — interactive HTML report
- `multiqc_data/` — raw data files

## Reproducibility

- `reproducibility/commands.sh` — replay this ClawBio MultiQC run
- `reproducibility/environment.yml` — suggested conda environment
- `reproducibility/checksums.sha256` — key outputs

---

*ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions.*
```

## Output Structure

```
output_dir/
├── report.md                        # ClawBio markdown summary
├── multiqc_report.html              # Standard MultiQC HTML
├── multiqc_data/
│   ├── multiqc_data.json            # Structured stats (default MultiQC output)
│   └── ...
├── reproducibility/
│   ├── commands.sh                  # Exact replay command
│   ├── environment.yml              # Suggested env (multiqc via pip)
│   └── checksums.sha256             # Output digests
```

## Dependencies

**External binary** (not a Python package import):
- `multiqc >= 1.20`; install with `pip install multiqc`

**Python** (repo-local `clawbio` package for reproducibility helpers):
- `subprocess`, `json`, `shutil`, `argparse`, `tempfile`, `math`
- `clawbio.common.reproducibility` — `commands.sh`, `environment.yml`, `checksums.sha256`

## Gotchas

- **You will want to parse tool-specific files directly.** Do not. MultiQC's auto-detection handles this; let it do its job. Parsing FastQC text yourself will miss 99 other supported tools.
- **`report_general_stats_data` metric keys are already short** (e.g. `percent_duplicates`, `percent_gc`) — no further processing needed. If the table looks empty, check that `multiqc_data/multiqc_data.json` exists and that `report_general_stats_data` is non-empty.
- **`--demo` creates files in a `tempfile.TemporaryDirectory` that is deleted after `run_multiqc` returns.** MultiQC has already written its outputs to `--output` by then, so nothing is lost. Don't move the `with` block boundary.
- **MultiQC exits 0 even if it found no recognised files** — it just produces an empty report. The skill does not treat this as an error; the user will see an empty table in `report.md` and an HTML report noting no modules were found.
- **Static PNG/SVG/PDF plots are not produced by this skill** — it never passes MultiQC `--export`. Interactive plots remain in `multiqc_report.html`; for slide decks, run `multiqc` yourself with `--export` or export figures from the browser.

## Safety

- **Local-first**: All processing is local; no data is uploaded
- **Disclaimer**: Every `report.md` includes the ClawBio medical disclaimer
- **No hallucinated metrics**: All values in the table come directly from `multiqc_data/multiqc_data.json`

## Agent Boundary

The agent (LLM) dispatches and explains results. The skill (Python + MultiQC CLI) executes.
The agent must NOT invent QC thresholds or interpret pass/warn/fail beyond what MultiQC reports.

## Integration with Bio Orchestrator

**Trigger conditions**: the orchestrator routes here when:
- User mentions "multiqc", "aggregate QC", "multi-sample QC report"
- Output directory from seq-wrangler, rnaseq-de, or scrna-orchestrator is provided alongside a request to summarise QC

**Chaining partners**:
- `seq-wrangler`: produces FastQC/fastp/BAM stats directories → feed into multiqc
- `rnaseq-de`: STAR/HISAT2 alignment logs → feed into multiqc for alignment QC
- `scrna-orchestrator`: STARsolo per-sample QC dirs → feed into multiqc
- `repro-enforcer`: folds the `reproducibility/` trio into pipeline-wide bundles

## Maintenance

- **Review cadence**: Re-evaluate when MultiQC releases a major version (check `multiqc --version`)
- **Staleness signals**: If per-sample tables are empty after a MultiQC upgrade, check whether `report_general_stats_data` still exists in `multiqc_data.json`
- **Deprecation**: Archive to `skills/_deprecated/` if MultiQC adds a native ClawBio integration

## Citations

- [MultiQC](https://multiqc.info); Ewels et al., Bioinformatics 2016 — multi-tool QC aggregation
