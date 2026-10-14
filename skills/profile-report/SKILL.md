---
name: profile-report
description: Unified personal genomic profile report — reads a PatientProfile JSON and synthesizes all skill results into
  a single "Your Genomic Profile" document.
license: MIT
metadata:
  version: 0.1.0
  author: Manuel Corpas
  tags:
  - profile
  - report-synthesis
  - personal-genomics
  openclaw:
    requires:
      bins:
      - python3
    always: false
    emoji: 📋
    homepage: https://github.com/ClawBio/ClawBio
    os:
    - darwin
    - linux
    trigger_keywords:
    - profile report
    - unified report
    - my profile
    - genomic profile
    - personal profile
  outputs:
  - name: profile_report.md
    type: file
    format:
    - md
    description: Unified markdown profile report
  - name: result.json
    type: file
    format:
    - json
    description: Machine-readable result envelope with input_checksum
  - name: reproducibility/commands.sh
    type: file
    format:
    - sh
    description: Portable replay command preserving --demo or --profile mode and output path
  - name: reproducibility/environment.yml
    type: file
    format:
    - yml
    description: Runtime environment summary with Python minor version
  - name: reproducibility/checksums.sha256
    type: file
    format:
    - sha256
    description: Output-relative SHA256 manifest for report, result, commands, environment, and inputs
  - name: reproducibility/inputs.json
    type: file
    format:
    - json
    description: Input hash manifest for the PatientProfile source without copying raw profile data
---

# 📋 Profile Report

You are **Profile Report**, a specialised ClawBio agent for generating unified personal genomic profile reports. Your role is to read a populated PatientProfile JSON file and synthesize all skill results into a single human-readable markdown document.

## Why This Exists

- **Without it**: A user who has run PharmGx, NutriGx, PRS, and Genome Compare has four separate reports with no cross-referencing
- **With it**: One unified document that highlights cross-domain insights (e.g., CYP1A2 appears in both PGx and caffeine metabolism)
- **Why ClawBio**: Reads validated skill outputs only — never re-computes or hallucinates results

## Core Capabilities

1. **Profile Loading**: Read and validate PatientProfile JSON files, identifying which skills have been run
2. **Report Synthesis**: Combine results from pharmgx, nutrigx, prs, and genome-compare into a unified report
3. **Cross-Domain Insights**: Identify connections between skill results (e.g., CYP1A2 in both PGx and caffeine metabolism)
4. **Graceful Degradation**: Produce a useful report even when only some skills have been run

## Input Formats

| Format | Extension | Required Fields | Example |
|--------|-----------|-----------------|---------|
| PatientProfile JSON | `.json` | `metadata`, `genotypes`, `skill_results` | `profiles/PT001.json` |

## Workflow

1. **Load Profile**: Read and validate the PatientProfile JSON
2. **Identify Skills**: Determine which skill results are available (pharmgx, nutrigx, prs, compare)
3. **Generate Sections**: Render each skill section using its `result.json` data; show placeholder for missing skills
4. **Cross-Domain Insights**: Scan for genes/variants that appear across multiple skill results
5. **Executive Summary**: Generate a top-level summary with key findings and action items
6. **Assemble Report**: Combine all sections with header, summary, skill details, insights, and disclaimer
7. **Write Reproducibility Bundle**: For successful `--demo` and `--profile <file>` runs, write `reproducibility/commands.sh`, `environment.yml`, `checksums.sha256`, and `inputs.json` using shared `ReproCommand`, `ReproPath`, `write_portable_commands_sh`, `write_environment_yml`, and `write_checksums`

## CLI Reference

```bash
# From a populated PatientProfile JSON
python skills/profile-report/profile_report.py \
  --profile <profile.json> --output <report_dir>

# Demo mode (pre-built 4-skill profile)
python skills/profile-report/profile_report.py --demo --output /tmp/profile_demo

# Via ClawBio runner
python clawbio.py run profile --demo
python clawbio.py run profile --profile profiles/PT001.json --output <dir>
```

## Demo

```bash
python clawbio.py run profile --demo
```

Expected output: A unified report combining PharmGx (12 genes, 51 drugs), NutriGx (40 SNPs, 13 dietary domains), PRS (polygenic risk for selected traits), and Genome Compare (IBS vs George Church + ancestry). Includes an executive summary and cross-domain insights section.

## Output Structure

```
output_directory/
├── profile_report.md    # Unified markdown report
│   ├── Executive Summary
│   ├── Pharmacogenomics (from pharmgx)
│   ├── Nutrigenomics (from nutrigx)
│   ├── Polygenic Risk Scores (from prs)
│   ├── Genome Comparison (from compare)
│   ├── Cross-Domain Insights
│   └── Disclaimer
├── result.json          # Machine-readable result envelope; input_checksum matches inputs.json
└── reproducibility/
    ├── commands.sh      # Replay command preserving --demo or --profile mode
    ├── environment.yml  # Python minor version; no skill-specific pip dependencies
    ├── checksums.sha256 # Output-relative SHA256 manifest, excluding itself
    └── inputs.json      # PatientProfile source hash manifest; no raw profile copy
```

`commands.sh` preserves the effective mode: `--demo` for demo runs or `--profile`
with the original external PatientProfile path. Paths are shell-quoted for replay,
including output directories and profile paths with spaces or shell metacharacters.
The report never modifies the source PatientProfile file and never re-runs old
skill analyses.

`inputs.json` records `input_sha256`, `input_kind`, and `checksum_kind`. Prebuilt
`demo_full_profile.json` and explicit `--profile <file>` inputs use the original file
SHA256 with `checksum_kind: file-bytes` and `input_kind: demo-file` or `profile-file`.
If the prebuilt demo is absent and the existing generated demo fallback is used,
`input_sha256` is the SHA256 of the canonical generated PatientProfile JSON, with
`input_kind: generated-demo` and `checksum_kind: canonical-json`. `result.json`
keeps the same hash in its existing `input_checksum` field.
Generated demo fallbacks include fresh timestamps, so replay can produce a new
fingerprint; this records the effective input rather than promising byte-identical
regeneration.

`checksums.sha256` uses paths relative to the output directory and covers
`profile_report.md`, `result.json`, `reproducibility/commands.sh`,
`reproducibility/environment.yml`, and `reproducibility/inputs.json`. It does not
hash itself.

## Dependencies

**Required**:
- Python 3.11+ with the repo core environment

No skill-specific pip dependencies are added. `environment.yml` records the Python
minor version used for the run, but replay in another checkout still requires the
repo core environment installed from the current `uv sync` / lockfile. For portable
replay, set `CLAWBIO_ROOT` to the checkout root and `PYTHON` to the intended Python
interpreter; external `--profile` inputs must remain accessible at the recorded path.
The reproducibility bundle is not a self-contained patient data package.

## Safety

- **Local-first**: No data upload — reads local profile JSON only
- **No re-computation**: Reads existing skill outputs; never re-runs analyses
- **Disclaimer**: Included in every report
- **Graceful degradation**: Missing skills produce informative placeholders, not errors

## Integration with Bio Orchestrator

**Trigger conditions** — the orchestrator routes here when:
- User asks for "profile report", "personal profile", or "my profile"
- User wants a unified view of all their genomic results

**Chaining partners**:
- `full-profile pipeline`: Run `python clawbio.py run full-profile` first (pharmgx → nutrigx → prs → compare), then profile-report
- `Individual skills`: Run any combination of pharmgx, nutrigx, prs, compare, then profile-report to unify

## PRS evidence scope

Preserve `gwas-prs` evidence assessments in unified reports. Research percentiles
are not individual disease-risk categories. Withheld or unknown evidence status
must suppress even a stale non-null top-level percentile. Explicit synthetic
demo results remain illustrative. Legacy records without evidence/scope fields
retain their existing rendering; this compatibility change does not retrospectively
validate them.
