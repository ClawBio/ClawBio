---
name: just-prs-mcp
description: >-
  Compute evidence-aware polygenic risk scores from a local VCF or WGS file
  through the validated just-prs engine and a pinned local just-prs MCP server.
license: MIT
metadata:
  version: "0.1.0"
  author: Anton Kulaga
  domain: genomics
  tags:
    - polygenic-risk-score
    - personal-genomics
    - model-quality
  inputs:
    - name: input_vcf
      type: file
      format:
        - vcf
        - vcf.gz
        - vcf.bgz
      description: Local single-sample VCF containing genotypes
      required: true
  outputs:
    - name: report
      type: file
      format:
        - md
      description: Evidence-aware PRS interpretation report
    - name: result
      type: file
      format:
        - json
      description: Machine-readable model results and uncertainty signals
    - name: scores
      type: file
      format:
        - csv
      description: Per-model scores and interpretation fields
  dependencies:
    python: ">=3.11"
    packages:
      - fastmcp>=3.4.4,<4
      - typer>=0.27.0,<1
  demo_data:
    - path: examples/demo_patient.vcf
      description: Synthetic three-variant GRCh38 VCF used with a cached response
  endpoints:
    cli: uv run --extra just-prs python skills/just-prs-mcp/just_prs_mcp_bridge.py --input {input_file} --trait {trait} --output {output_dir}
  openclaw:
    requires:
      bins:
        - uvx
    always: false
    emoji: "🧬"
    homepage: https://github.com/dna-seq/just-prs-mcp
    os:
      - darwin
      - linux
    install:
      - kind: uv
        command: uv sync --extra just-prs
    trigger_keywords:
      - compute PRS from VCF
      - WGS polygenic risk score
      - absolute genetic risk from VCF
      - just-prs
---

# just-prs MCP bridge

You are **just-prs MCP bridge**, a specialised ClawBio agent for evidence-aware
polygenic scoring of local VCF and WGS genotypes.

## Trigger

**Fire this skill when the user says any of:**
- "compute PRS from my VCF"
- "score this WGS genome for type 2 diabetes"
- "estimate absolute genetic risk from this VCF"
- "use just-prs on my genome"
- "compare the PRS models for this trait"

**Do NOT fire when:**
- The input is a 23andMe or AncestryDNA text export; use `gwas-prs`.
- The input is raw FASTQ/BAM requiring variant calling; use `wgs-prs` first.
- The user asks for one variant's disease association; use `gwas-lookup`.
- A bare VCF is supplied without PRS or absolute-risk intent.

## Why This Exists

- **Without it**: VCF users must manually identify PGS models, normalize data,
  run scores, inspect coverage, obtain ancestry-matched percentiles, and compare models.
- **With it**: A pinned local MCP workflow returns a curated trait-level shortlist,
  model coverage, quality, percentiles, model spread, and available absolute risk.
- **Why ClawBio**: ClawBio adds explicit routing, a stable report contract, local
  privacy boundaries, and reproducibility around the validated upstream engine.

## Core Capabilities

1. **VCF/WGS scoring**: Score one PGS ID or a curated set associated with an EFO/MONDO trait.
2. **Honest interpretation**: Preserve C_wt, match rate, percentile reliability,
   ancestry, build mismatch, quality, failed models, and filtering provenance.
3. **Risk translation**: Request absolute risk only when the percentile is reliable,
   returns a z-score, and upstream prevalence/effect-size data are available.
4. **Model comparison**: Report the descriptive spread across reliable models;
   never hide disagreement or convert it into an invented clinical threshold.

## Scope

**One skill, one task.** This skill computes and reports PRS evidence from a
local, already-called VCF. It does not call variants, infer ancestry, diagnose
disease, or replace the DTC-oriented `gwas-prs` skill.

## Input Formats

| Format | Extension | Required fields | Example |
|---|---|---|---|
| VCF 4.x | `.vcf` | `#CHROM`, `POS`, `REF`, `ALT`, sample `GT` | `examples/demo_patient.vcf` |
| Compressed VCF | `.vcf.gz`, `.vcf.bgz` | Same fields, bgzip-compatible | user-provided |

## Workflow

1. **Validate (prescriptive)**: Require one local VCF and exactly one selector:
   trait term, EFO/MONDO trait ID, or PGS ID.
2. **Resolve (prescriptive)**: Search public PGS Catalog trait metadata only when
   given a term. Stop on ambiguity and require `--trait-id`.
3. **Compute (prescriptive)**: Launch `just-prs-mcp==0.2.0` with local stdio in
   essentials mode. Pass the resolved local path, never VCF bytes.
4. **Curate (prescriptive)**: For trait mode request `interpret=true` and
   `profile=curated` by default. Preserve the upstream filter summary and failures.
5. **Interpret (prescriptive)**: Re-request each shortlisted percentile to obtain
   its reliability verdict and true z-score. Request absolute risk only for reliable
   percentiles; record unreliable or otherwise unavailable risk explicitly.
6. **Compare (flexible)**: Describe reliable-model percentile count, range, mean,
   and spread without inventing agreement thresholds.
7. **Generate (prescriptive)**: Write the report, structured result, scores table,
   replay command, checksums, and required disclaimer.

## CLI Reference

```bash
uv sync --extra just-prs

uv run --extra just-prs python skills/just-prs-mcp/just_prs_mcp_bridge.py \
  --input sample.vcf.gz \
  --trait "type 2 diabetes" \
  --superpopulation EUR \
  --output output/just-prs-t2d

uv run --extra just-prs python skills/just-prs-mcp/just_prs_mcp_bridge.py \
  --input sample.vcf.gz --pgs-id PGS000014 --output output/just-prs-single

uv run --extra just-prs python skills/just-prs-mcp/just_prs_mcp_bridge.py \
  --demo --output /tmp/just_prs_demo

uv run --extra just-prs clawbio.py run just-prs --demo
```

## Demo

Run:

```bash
uv run --extra just-prs clawbio.py run just-prs --demo
```

The demo is deterministic and offline. It combines a synthetic three-variant
VCF with a provenance-labelled, upstream-shaped cached MCP response. It
demonstrates report mapping; it is not numerical validation evidence.

## Algorithm / Methodology

1. Resolve a specific ontology trait or PGS Catalog score.
2. Compute `sum(effect_weight × dosage)` in upstream `just-prs`.
3. Retain matched/total variant coverage and weight-mass coverage, C_wt.
4. Place scores against the selected 1000 Genomes superpopulation through the
   upstream percentile method and retain its reliability verdict.
5. Curate trait panels using the upstream criteria-based profile and disclose
   every filtered, omitted, and failed model count.
6. Request absolute risk only when the upstream percentile is reliable, using its
   z-score and upstream prevalence/effect-size evidence.
7. Compare multiple reliable models descriptively rather than selecting one
   convenient result.

**Key parameters**:
- Superpopulations: AFR, AMR, EAS, EUR, SAS (1000 Genomes).
- Default profile: `curated` (criteria owned by `just-prs-mcp`, not ClawBio).
- Default returned models: 5, ranked by the upstream coverage-aware ordering.

## Example Queries

- "Compute my type 2 diabetes PRS from this GRCh38 VCF."
- "Use PGS000014 on this WGS callset and show absolute risk if available."
- "Do the reliable models agree on my coronary artery disease percentile?"

## Example Output

```markdown
# just-prs Polygenic Risk Report

## Model agreement
- Reliable models: **2**
- Verdict: **descriptive_spread_only**
- Reliable percentile range: **61.00–74.00**
- Descriptive percentile spread: **13.00**

## Score details
| PGS ID | Status | Percentile | Reliable | C_wt | Quality | Absolute risk |
|---|---|---:|---|---:|---|---|
| PGS000014 | scored | 74.00 | True | 94.0% | High | 18.0% |
| PGS000013 | scored | 61.00 | True | 91.0% | Normal | unavailable |
```

## Output Structure

```text
output_directory/
├── report.md
├── result.json
├── tables/
│   └── scores.csv
└── reproducibility/
    ├── commands.sh
    └── checksums.sha256
```

## Dependencies

**Required**:
- `uvx`; launches the isolated Python 3.13+ upstream server.
- `fastmcp` in the `just-prs` optional extra; Python 3.11-compatible client.
- `typer` in the `just-prs` optional extra; typed CLI.

**Optional**:
- A warm upstream cache; avoids repeat downloads but does not change results.

The upstream cache defaults to the platform-specific `just-prs` cache. Override
it with `PRS_MCP_CACHE_DIR` when a controlled shared cache is required.

## Validation Evidence

- Upstream `just-prs/tests/test_cross_engine.py` checks numerical parity across
  PLINK2, Polars, and DuckDB scoring engines.
- Upstream `test_scoring.py`, `test_vcf.py`, and `test_percentile.py` cover scoring
  files, VCF/build handling, and percentile/z-score consistency.
- This bridge tests its integration boundary: pinned real stdio MCP calls,
  FastMCP result decoding, local-path-only requests, report mapping, routing,
  packaging, unreliable-risk suppression, and offline demo reproducibility.
- ClawBio does not reimplement or claim independent validation of the upstream
  numerical engine; it reuses that evidence and tests its own adapter behavior.

## Gotchas

- **You will want to treat match rate as model coverage. Do not. Here is why.**
  Match rate counts variants equally; retain C_wt because effect-weight mass is
  the upstream scale-free honesty signal.
- **You will want to report the highest percentile as the answer. Do not. Here
  is why.** Multiple models can disagree because of coverage, ancestry, and
  model design; report the reliable-model spread and filtering provenance.
- **You will want to interpret a missing absolute-risk result as average risk.
  Do not. Here is why.** Missing prevalence or effect-size evidence means the
  estimate is unavailable, not normal.
- **You will want to treat an empty curated row list as a failed run. Do not.
  Here is why.** The engine may have scored models and then explicitly removed
  every one for weak evidence or coverage; inspect `n_filtered` and
  `filter_summary` so the exclusion is visible.
- Do not silently accept trait-search ambiguity; require a stable ontology ID.
- Do not interpret a build-mismatched score until the coordinate build is resolved.
- Do not use a VCF fixture with sample genotypes unless it declares the `GT`
  FORMAT header; real readers correctly omit an undeclared genotype field.

## Safety

- **Local-first**: The VCF remains on the machine. The bridge passes only its
  resolved path to a local stdio child process and has no upload path.
- **Credential isolation**: The child receives an allow-listed runtime, network,
  and cache environment; unrelated API keys and service credentials are not forwarded.
- **Network egress**: Live runs fetch only public PGS Catalog metadata, scoring
  files, and reference distributions. This is a documented dependency, not a
  claim of zero network access. The offline demo makes no network calls.
- **Consent boundary**: Never switch this bridge to hosted HTTP/SSE for personal
  genomic data. A remote server cannot access the local path and must not receive the VCF.
- **Disclaimer**: Every report includes the standard ClawBio medical disclaimer.
- **Audit trail**: Reports record versions, input checksum, replay command, and
  output checksums without copying genotype content into outputs.
- **No hallucinated science**: Scientific calculations and curation criteria
  remain upstream; ClawBio maps and explains the returned evidence.

## Agent Boundary

The agent dispatches, asks for ancestry/trait clarification, and explains. The
skill executes scoring and evidence retrieval. The agent must not override
upstream reliability, invent an absolute risk, or suppress model failures.

## Integration with Bio Orchestrator

**Trigger conditions**:
- A VCF/WGS input plus explicit PRS, polygenic-risk, or absolute-risk intent.
- An explicit request to use just-prs on a local genome.

**Chaining partners**:
- `wgs-prs`: Produces a called VCF from raw sequencing before this skill.
- `gwas-prs`: Handles 23andMe/AncestryDNA text inputs instead of VCF/WGS.
- `profile-report`: May consume `result.json` in a later compatibility PR.

## Maintenance

- **Review cadence**: Review monthly and on every `just-prs-mcp` release.
- **Staleness signals**: MCP schema drift, changed curation fields, a new
  `just-prs` coverage contract, or a PGS Catalog API change.
- **Deprecation**: Archive if upstream no longer supports local stdio or ClawBio
  adopts a generic MCP bridge with the same tested report contract.

## Citations

- [PGS Catalog](https://www.pgscatalog.org/); public score metadata and scoring files.
- [just-prs](https://github.com/dna-seq/just-prs); scoring, normalization,
  percentile, quality, and PLINK2 parity validation.
- [just-prs-mcp](https://github.com/dna-seq/just-prs-mcp); typed MCP contracts,
  curated trait workflow, and local stdio privacy boundary.
- [just-dna-lite FAQ](https://github.com/dna-seq/just-dna-lite/blob/main/docs/FAQ.md);
  local computation, user ownership, open access, and citizen-science context.

*ClawBio is a research and educational tool. It is not a medical device and
does not provide clinical diagnoses. Consult a healthcare professional before
making any medical decisions.*
