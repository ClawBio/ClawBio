---
name: nfcore-scrnaseq-wrapper
description: Wrapper skill for running nf-core/scrnaseq-style single-cell RNA-seq preprocessing from FASTQ with strict preflight, reproducibility outputs, and downstream handoff to ClawBio scRNA skills.
license: MIT
metadata:
  version: "0.1.0"
  author: ClawBio
  domain: transcriptomics
  tags:
    - scrna
    - single-cell
    - nextflow
    - nf-core
    - fastq
    - 10x
    - h5ad
    - preprocessing
  inputs:
    - name: samplesheet
      type: file
      format:
        - csv
      description: Valid scrnaseq samplesheet with sample, fastq_1, fastq_2 columns
      required: true
  outputs:
    - name: report
      type: file
      format:
        - md
      description: Wrapper run summary and downstream handoff recommendations
    - name: result
      type: file
      format:
        - json
      description: Structured result payload with detected outputs and provenance
  dependencies:
    python: ">=3.11"
    packages: []
  demo_data:
    - path: demo/README.md
      description: Demo mode uses the upstream scrnaseq test profile rather than bundled sample data
  endpoints:
    cli: python skills/nfcore-scrnaseq-wrapper/nfcore_scrnaseq_wrapper.py --input {samplesheet} --output {output_dir}
  openclaw:
    requires:
      bins:
        - python3
        - nextflow
        - java
      env: []
      config: []
    always: false
    emoji: "🧫"
    homepage: https://github.com/ClawBio/ClawBio
    os:
      - darwin
      - linux
    install: []
    trigger_keywords:
      - scrnaseq
      - nf-core scrnaseq
      - run scrnaseq from fastq
      - preprocess 10x fastqs
      - generate h5ad from single-cell fastq
      - single-cell preprocessing
      - nextflow scrna pipeline
      - 10x chromium fastq pipeline
      - starsolo upstream processing
      - alevin-fry fastq to counts
      - run nextflow scrnaseq
      - upstream single-cell pipeline
      - fastq to h5ad single cell
      - 10x genomics fastq pipeline
---

# 🧫 nfcore-scrnaseq-wrapper

You are **nfcore-scrnaseq-wrapper**, a specialised ClawBio agent for upstream single-cell RNA-seq preprocessing from FASTQ using the `scrnaseq` Nextflow pipeline.

## Trigger

**Fire when:**
- User wants to run `scrnaseq` from raw FASTQ files
- User asks to preprocess 10x Chromium single-cell data
- User wants to execute `nf-core/scrnaseq`
- User wants to generate `.h5ad` from raw single-cell FASTQs
- User asks for primary scRNA preprocessing (FASTQ → h5ad)
- User mentions `simpleaf`, `STARsolo`, or `kb-python` as the aligner for upstream processing

**Do NOT fire when:**
- User already has an `.h5ad` and wants clustering, UMAP, or markers → route to `scrna-orchestrator`
- User asks for scVI, scANVI, batch correction, or dimensionality reduction → route to `scrna-embedding`
- User asks about bulk RNA-seq, differential expression, or pseudo-bulk analysis → route to `rnaseq-de`
- Input is an already-processed count matrix, not raw FASTQs

## Scope

One skill, one task: run upstream scRNA preprocessing from FASTQ using the `scrnaseq` Nextflow pipeline and produce canonical outputs for downstream ClawBio skills.

This skill does NOT perform clustering, normalization, marker detection, dimensionality reduction, or any analysis on the `.h5ad` it produces.

## Why This Exists

Running upstream scRNA pipelines manually is error-prone.

- **Without it**: Users hand-build samplesheets, guess reference combinations, miss backend issues, and struggle to find the right `.h5ad` for downstream analysis.
- **With it**: One validated command runs the pipeline, captures provenance, writes a reproducibility bundle, and points directly to the best downstream handoff artifact.
- **Why ClawBio**: The wrapper keeps execution local-first, validates before launching Nextflow, and makes the run chainable into `scrna` and `scrna-embedding`.

## Core Capabilities

1. **Strict Preflight**: Validate Java, Nextflow, backend, samplesheet, FASTQs, and references before execution.
2. **Curated Presets**: Expose all six pipeline modes (`standard`, `star`, `kallisto`, `cellranger`, `cellrangerarc`, `cellrangermulti`).
3. **Controlled Execution**: Always run with `-params-file`, fixed pipeline source, and explicit reproducibility artifacts.
4. **Output Resolution**: Detect MultiQC, pipeline_info, `.h5ad`, `.rds`, and choose a canonical `preferred_h5ad` when possible.
5. **Downstream Handoff**: Recommend the next command for `scrna-orchestrator` (auto via `--run-downstream`); `scrna-embedding` can follow as a second step.

## Input Formats

| Format | Extension | Required Fields | Example |
|--------|-----------|-----------------|---------|
| Samplesheet | `.csv` | `sample`, `fastq_1`, `fastq_2` | `samplesheet.csv` |
| Samplesheet CSV | `.csv` | `sample`, `fastq_1`, `fastq_2`; optional `expected_cells` and preset-specific columns | `--input samplesheet.csv` |
| Demo mode | n/a | none | `python clawbio.py run scrnaseq-pipeline --demo` |

## Workflow

When the user asks for scRNA preprocessing from FASTQ:

1. **Validate**: Check the selected preset, samplesheet structure, FASTQ accessibility, references, Java, Nextflow, and backend.
2. **Normalize**: Write a validated samplesheet copy into the reproducibility bundle.
3. **Configure**: Build one effective `params.yaml` and a fixed Nextflow command.
4. **Execute**: Run the upstream `scrnaseq` pipeline using the local checkout when available.
5. **Parse**: Detect MultiQC, pipeline_info, `.h5ad`, `.rds`, and CellBender-derived outputs.
6. **Generate**: Write `report.md`, `result.json`, provenance JSON files, and reproducibility files.
7. **Hand off**: Recommend the next ClawBio command using the preferred `.h5ad` when available.

## Algorithm / Methodology

The wrapper executes a strictly ordered 7-step pipeline. Each step is gated: a failure at any step raises a structured `SkillError` with an `error_code` and a `fix` hint, and no subsequent step runs.

1. **Pipeline source resolution** (`pipeline_source.py`): Determine whether to use a local sibling `scrnaseq/` checkout (pinned commit, audit-safe) or fall back to the remote pipeline tag. Local checkouts always take priority.

2. **Samplesheet validation** (`samplesheet_builder.py`): Parse the user-supplied CSV, resolve all FASTQ paths against the CSV's parent directory (not CWD), normalize whitespace in sample names to underscores, verify readability and FASTQ extensions, reject FASTQ basenames with whitespace, preserve upstream/preset-specific columns, reject exact duplicate FASTQ rows, enforce consistent `expected_cells` and `seq_center` for repeated sample rows, reject `expected_cells < 1`, and write a normalized copy with absolute POSIX FASTQ paths to `reproducibility/samplesheet.valid.csv`.

3. **Preflight** (`preflight.py`): Verify Java (≥17), Nextflow (≥25.04.0), and the selected backend (`docker`, `conda`, `singularity`, or `apptainer`). For Docker, run `docker info` and gate on its exit code. Compare version tuples after zero-padding to 3 elements (avoids `(25, 4) < (25, 4, 0)` false negatives). All subprocess calls have a 30-second timeout.

4. **Params construction** (`params_builder.py`): Translate the preset + CLI flags into a `params.yaml` that Nextflow consumes via `-params-file`. All paths use `.as_posix()` for forward-slash consistency across platforms.

5. **Command build + execution** (`command_builder.py`, `executor.py`): Construct the `nextflow run` command with `-params-file` and `-work-dir <output>/upstream/work`, then launch it via `subprocess.Popen` with stdout and stderr piped directly to log files on disk — never buffered in RAM. On `TimeoutExpired`, the process is killed and `EXECUTION_FAILED` is raised.

6. **Output parsing** (`outputs_parser.py`): Scan the output tree for MultiQC HTML, `pipeline_info/`, `.h5ad` (combined matrix preferred over per-sample), `.rds`, and CellBender-derived files. Set `handoff_available = True` only when a `preferred_h5ad` is confirmed on disk.

7. **Provenance + reporting** (`provenance.py`, `reporting.py`): Write JSON provenance bundles, a SHA-256 checksum manifest, `environment.yml`, a reproducibility shell script, `report.md`, and `result.json`. SHA-256 is computed only on `is_file()` paths — never on directories such as `star_index/`.

## Presets

| Preset | Aligner | Use case |
|--------|---------|---------|
| `standard` | simpleaf (alevin-fry) | Default for 10x GEX; fast, memory-efficient |
| `star` | STARsolo | Best FASTQ QC metrics; supports RNA velocity (`--star-feature "Gene Velocyto"`) |
| `kallisto` | kb-python / BUStools | Pseudo-alignment; fastest; lamanno/nac RNA velocity via `--kb-workflow` |
| `cellranger` | CellRanger | CellRanger v2/v3 compatibility; requires CellRanger binary in PATH |
| `cellrangerarc` | CellRanger ARC | Multiome (GEX + ATAC); accepts prebuilt `--cellranger-index` or reference-build inputs |
| `cellrangermulti` | CellRanger Multi | GEX + VDJ + feature barcoding in one run; `--cellranger-multi-barcodes` is required for CMO/multiplexed assignments and FFPE probe-set demultiplexing |

Each preset requires at least one reference option (in priority order): `--genome <iGenomes_name>` OR a pre-built index (`--star-index`, `--simpleaf-index`, etc.) OR `--fasta` + `--gtf`.

## CLI Reference

```bash
# Standard usage
python skills/nfcore-scrnaseq-wrapper/nfcore_scrnaseq_wrapper.py \
  --input samplesheet.csv --output ./scrnaseq_wrapper_run

# Preflight check only (no Nextflow execution)
python skills/nfcore-scrnaseq-wrapper/nfcore_scrnaseq_wrapper.py \
  --input samplesheet.csv --output ./scrnaseq_wrapper_run --check

# Demo mode (uses nf-core test profile; forces star preset; Docker must be running)
python skills/nfcore-scrnaseq-wrapper/nfcore_scrnaseq_wrapper.py \
  --demo --output ./scrnaseq_wrapper_demo

# Via ClawBio runner
python clawbio.py run scrnaseq-pipeline --input samplesheet.csv --output ./scrnaseq_wrapper_run
python clawbio.py run scrnaseq-pipeline --demo --output ./scrnaseq_wrapper_demo

# STARsolo with iGenomes reference
python skills/nfcore-scrnaseq-wrapper/nfcore_scrnaseq_wrapper.py \
  --input samplesheet.csv --output ./run --preset star --genome GRCh38

# STARsolo RNA velocity
python skills/nfcore-scrnaseq-wrapper/nfcore_scrnaseq_wrapper.py \
  --input samplesheet.csv --output ./run --preset star \
  --star-feature "Gene Velocyto" --star-ignore-sjdbgtf \
  --fasta /refs/hg38.fa --gtf /refs/hg38.gtf

# Kallisto RNA velocity (NAC workflow)
python skills/nfcore-scrnaseq-wrapper/nfcore_scrnaseq_wrapper.py \
  --input samplesheet.csv --output ./run --preset kallisto \
  --kb-workflow nac --fasta /refs/hg38.fa --gtf /refs/hg38.gtf

# CellRanger ARC (Multiome)
python skills/nfcore-scrnaseq-wrapper/nfcore_scrnaseq_wrapper.py \
  --input samplesheet.csv --output ./run --preset cellrangerarc \
  --genome GRCh38 \
  --cellrangerarc-reference GRCh38-2020-A-2.0.0 \
  --cellrangerarc-config /refs/arc_config.json \
  --motifs /refs/motifs.pfm

# CellRanger Multi (multiplexed)
python skills/nfcore-scrnaseq-wrapper/nfcore_scrnaseq_wrapper.py \
  --input samplesheet.csv --output ./run --preset cellrangermulti \
  --cellranger-index /refs/refdata-gex-GRCh38 \
  --gex-cmo-set /refs/cmo_set.csv \
  --fb-reference /refs/feature_barcodes.csv
```

### Key flags

| Flag | Type | Description |
|------|------|-------------|
| `--preset` | string | Aligner preset (default: `standard`) |
| `--profile` | string | Execution backend: `docker`, `conda`, `singularity`, `apptainer` (default: `docker`) |
| `--genome` | string | iGenomes shortcut (`GRCh38`, `mm10`, etc.) — mutually exclusive with `--fasta`/`--gtf` |
| `--fasta` | path | Genome FASTA (resolved to absolute POSIX path) |
| `--gtf` | path | GTF annotation |
| `--star-feature` | enum | STARsolo feature type: `Gene`, `GeneFull`, `Gene Velocyto` |
| `--star-ignore-sjdbgtf` | flag | Do not use GTF for SJDB (use with RNA velocity) |
| `--kb-workflow` | enum | Kallisto workflow: `standard`, `lamanno`, `nac` |
| `--kb-t1c` | path | cDNA t2c file for RNA velocity (lamanno/nac) |
| `--kb-t2c` | path | Intron t2c file for RNA velocity (lamanno/nac) |
| `--simpleaf-umi-resolution` | enum | UMI strategy: `cr-like`, `cr-like-em`, `parsimony`, `parsimony-em`, `parsimony-gene`, `parsimony-gene-em` |
| `--skip-fastqc` | flag | Skip FastQC (when QC was done externally) |
| `--skip-multiqc` | flag | Skip MultiQC report generation |
| `--skip-cellbender` | flag | Disable CellBender ambient RNA removal. Upstream never runs CellBender for `cellrangerarc`, regardless of this flag |
| `--skip-emptydrops` | flag | Deprecated compatibility alias for `--skip-cellbender`; replay commands use `--skip-cellbender` |
| `--skip-cellranger-renaming` | flag | Skip automatic sample renaming in CellRanger |
| `--skip-cellrangermulti-vdjref` | flag | Skip mkvdjref in cellrangermulti only when there is no VDJ data or a prebuilt `--cellranger-vdj-index` is supplied |
| `--save-reference` | flag | Save the built reference index for future reuse |
| `--save-align-intermeds` | flag | Save alignment intermediate BAM files |
| `--email` | string | Email for pipeline completion notification |
| `--multiqc-title` | string | Custom title for the MultiQC report |
| `--seq-center` | string | Sequencing center name for BAM read group tag |
| `--expected-cells` | int | Override expected cells count (≥1) |
| `--protocol` | string | Chemistry/protocol. Required and non-`auto` for `standard`, `star`, and `kallisto`; `standard`/Simpleaf rejects `smartseq`; `cellranger` accepts only `auto` or `10XV1`-`10XV4`; ARC accepts only `auto` |
| `--resume` | flag | Nextflow resume (checksum-verified against prior run) |
| `--run-downstream` | flag | Opt in to `scrna_orchestrator` handoff after pipeline completion |
| `--skip-downstream` | flag | Compatibility flag; downstream handoff is skipped unless `--run-downstream` is set |

## Example Queries

The following are concrete user phrases that should route to this skill:

- *"I have 10x Chromium FASTQs from three donors. How do I get an h5ad I can cluster?"*
- *"Run nf-core/scrnaseq on my samplesheet and give me a reproducibility bundle."*
- *"I want to preprocess single-cell data from FASTQ using STARsolo. Can you set it up?"*
- *"Generate a combined_filtered_matrix.h5ad from these FASTQ files using alevin-fry."*
- *"We have paired-end FASTQ files from a 10x Genomics Chromium v3 experiment. Run the upstream scRNA pipeline."*
- *"Demo mode — show me what the scrnaseq wrapper produces without any input files."*

## Demo

```bash
python clawbio.py run scrnaseq-pipeline --demo --output ./outputs/scrnaseq_demo
```

Expected output:
- wrapper-level `report.md`
- wrapper-level `result.json`
- `logs/`
- `provenance/`
- `reproducibility/`
- upstream `results/` tree under `upstream/`

## Output Structure

```text
output_directory/
├── report.md
├── result.json
├── logs/
│   ├── stdout.txt
│   └── stderr.txt
├── upstream/
│   └── results/
├── reproducibility/
│   ├── samplesheet.valid.csv
│   ├── commands.sh
│   ├── params.yaml
│   ├── environment.yml
│   ├── checksums.sha256
│   └── manifest.json
└── provenance/
    ├── upstream.json
    ├── skill.json
    ├── invocation.json
    ├── inputs.json
    ├── outputs.json
    ├── runtime.json
    └── preflight.json
```

## Dependencies

**Required**:
- `nextflow`
- `java` 17+
- one supported backend (`docker`, `conda`, `singularity`, or `apptainer`)

**Optional**:
- local sibling checkout of `scrnaseq` for pinned local execution

## Safety

- **Local-first**: User FASTQs and outputs remain local.
- **Strict preflight**: Do not launch Nextflow if validation fails.
- **Controlled surface**: No arbitrary passthrough flags in v1.
- **Disclaimer**: Every report includes the ClawBio disclaimer.
- **No hallucinated outputs**: Report only artifacts that were actually detected on disk.

## Example Output

After a successful run, `result.json` contains:

```json
{
  "skill": "scrnaseq-pipeline",
  "version": "0.1.0",
  "summary": {
    "preset": "star",
    "aligner_effective": "star",
    "pipeline_source_kind": "local_checkout",
    "pipeline_version_or_commit": "abc1234abcde",
    "profile": "docker",
    "resume_used": false,
    "preferred_h5ad": "output/upstream/results/star/mtx_conversions/combined_filtered_matrix.h5ad",
    "handoff_available": true,
    "samples_detected": 2,
    "cellbender_used": false,
    "output_artifacts": {
      "preferred_h5ad": "output/upstream/results/star/mtx_conversions/combined_filtered_matrix.h5ad",
      "multiqc_report": "output/upstream/results/multiqc/star/multiqc_report.html",
      "pipeline_info_dir": "output/upstream/results/pipeline_info",
      "h5ad_candidates": [
        "output/upstream/results/star/mtx_conversions/combined_filtered_matrix.h5ad"
      ],
      "rds_candidates": []
    }
  }
}
```

And `report.md` closes with:
```
## Next Steps

- `python clawbio.py run scrna --input output/upstream/.../combined_filtered_matrix.h5ad --output <dir>`
- `python clawbio.py run scrna-embedding --input output/upstream/.../combined_filtered_matrix.h5ad --output <dir>`
```

## Gotchas

- **Preflight runs before any Nextflow call.** If Java, Nextflow, or Docker are missing, the pipeline never starts and you get a structured JSON error with `error_code` and a `fix` hint. Do not try to skip preflight.
- **`--resume` enforces strict compatibility.** The wrapper checks that the stored manifest matches the current preset, profile, and pipeline source. If any differ, it raises `INVALID_RESUME_STATE` rather than silently continuing with incompatible state.
- **`--demo` runs the pipeline's real `test` profile, not bundled sample data.** It exercises the full Nextflow code path with the upstream test fixtures. Docker must be running. Output reflects what `nf-core/scrnaseq -profile test,docker` produces.
- **`preferred_h5ad` may be empty.** If no combined matrix is found and there are multiple per-sample files, `handoff_available` is `false`. Always check `result.json` before chaining to `scrna-orchestrator` or manually invoking `scrna-embedding`.
- **No arbitrary Nextflow passthrough.** You cannot add `-c`, `--outdir`, or custom Nextflow flags. All pipeline configuration flows through the preset system and `params.yaml`. Attempting to inject flags directly is blocked by the wrapper's CLI contract.
- **`--genome` conflicts with any explicit reference flag.** Providing `--genome` alongside any of `--fasta`, `--gtf`, `--star-index`, `--simpleaf-index`, `--kallisto-index`, `--cellranger-index`, `--transcript-fasta`, or `--txp2gene` raises `CONFLICTING_REFERENCES` in preflight. Use either `--genome <shortcut>` or explicit reference/index flags — never both.
- **Protocol compatibility is enforced before Nextflow starts.** `standard`/Simpleaf, `star`, and `kallisto` cannot rely on upstream `protocol=auto`, so provide an explicit protocol such as `10XV3`, `dropseq`, or a supported custom chemistry string. `standard`/Simpleaf rejects `smartseq`; use `star` or `kallisto` for Smart-seq data. `cellranger` accepts only `auto` or `10XV1`-`10XV4`; `cellrangerarc` accepts only `auto`.
- **`--skip-emptydrops` is deprecated upstream.** The wrapper still accepts it for old commands, but normalizes it to the canonical `skip_cellbender` parameter in `params.yaml` and replay scripts.
- **CellBender follows upstream behavior.** CellBender runs by default when supported and is disabled by `--skip-cellbender` (or deprecated `--skip-emptydrops`). The upstream workflow never runs CellBender for `cellrangerarc`, even when `skip_cellbender=false`.
- **Schema-incompatible paths fail early.** The wrapper writes `params.input` as a relative bundle path and rejects FASTA paths that do not match the nf-core/scrnaseq 4.1.0 schema (`.fa`, `.fna`, `.fasta`, and gzip variants, with no whitespace).
- **RNA velocity requires two coordinated flags.** For STARsolo, pass `--star-feature "Gene Velocyto"` AND `--star-ignore-sjdbgtf` together — the former selects the feature type, the latter suppresses GTF-based SJDB construction (required for velocity matrices). For Kallisto, use `--kb-workflow lamanno` or `nac` with `--kb-t1c` and `--kb-t2c`.
- **`cellrangerarc` follows the upstream ARC rules.** It accepts a prebuilt `--cellranger-index`, or builds from `--genome` / `--fasta` + `--gtf`. `--motifs` is optional. If you provide a custom `--cellrangerarc-config`, you must also provide `--cellrangerarc-reference` (and vice versa), otherwise preflight fails with `INVALID_PRESET_CONFIGURATION`.
- **`cellrangermulti` validates requirements from the samplesheet.** The main samplesheet must include `feature_type`. If `feature_type=ab`, preflight requires `--fb-reference`; if `feature_type=cmo`, preflight requires `--cellranger-multi-barcodes`; if `feature_type=vdj` and `--skip-cellrangermulti-vdjref` is set, preflight requires `--cellranger-vdj-index`. CMO, FFPE probe-set demultiplexing, and OCM assignment are treated as mutually exclusive multiplexing strategies.
- **`--save-align-intermeds` defaults to false in `params.yaml`.** This overrides heavy upstream defaults unless the user explicitly requests alignment intermediates.
- **`--check` writes `check_result.json` to the output directory.** This file is ignored by subsequent full runs. Running `--check` and a full run to the same `--output` directory is safe — `check_result.json` will be overwritten by `result.json` on success.
- **`--demo` always forces `--preset star` and `--skip-cellbender`.** The nf-core test profile ships STAR-compatible data only. If you request another preset with `--demo`, the wrapper warns and overrides it.
- **Local checkout must be a sibling directory.** The wrapper looks for `../scrnaseq` relative to the ClawBio repo root — both repos must share the same parent folder. If `scrnaseq/` is anywhere else, the wrapper silently falls back to downloading the remote pipeline from GitHub (slower, not pinned to a commit). To use the local path, ensure `ClawBio/` and `scrnaseq/` are siblings.

## Agent Boundary

The agent dispatches and explains; this skill executes.

**The agent is responsible for:**
- Interpreting the user's preprocessing intent and choosing the right preset
- Verifying that `handoff_available` is true in `result.json` before routing to `scrna-orchestrator` (automatic) or `scrna-embedding` (manual follow-up)

**This skill is responsible for:**
- Validating the environment and inputs before any execution
- Running the Nextflow pipeline with controlled, reproducible parameters
- Writing all provenance, checksums, and reproducibility artifacts
- Reporting the detected `preferred_h5ad` for downstream handoff

## Integration with Bio Orchestrator

**Trigger conditions**:
- user wants to run `scrnaseq` from FASTQs
- user asks for upstream scRNA preprocessing
- user wants `.h5ad` generation from raw single-cell FASTQs

**Do NOT route here when**:
- the user already has `.h5ad` and wants clustering, markers, or UMAP
- the user wants scVI/scANVI or batch correction downstream

**Chaining partners**:
- `scrna-orchestrator`: downstream clustering and markers from `preferred_h5ad`
- `scrna-embedding`: downstream scVI/scANVI from `preferred_h5ad`

## Maintenance

**Review cadence**: Review this skill whenever `nf-core/scrnaseq` releases a new major version (currently tracked at `nf-co.re/scrnaseq/changelog`). Check that `NEXTFLOW_MIN_VERSION` in `schemas.py` and the `SUPPORTED_PRESETS` map in `schemas.py` remain correct.

**Staleness signals**:
- `preflight.py` raises `NEXTFLOW_VERSION_TOO_OLD` on a version that the current nf-core release explicitly supports — update `NEXTFLOW_MIN_VERSION` in `schemas.py`.
- New aligners appear in the upstream `scrnaseq` pipeline but are not in `SUPPORTED_PRESETS` — add them to `schemas.py` and update this SKILL.md and the test suite.
- The VirtioFS macOS workaround (`stageInMode = "copy"`) is only necessary while Apple Silicon runs Docker via QEMU. If Apple releases a native arm64 Docker runtime that eliminates VirtioFS deadlocks, remove `_write_macos_docker_config` and its tests.
- The local sibling checkout path (`../scrnaseq`) is a convention for this thesis environment. Downstream deployments should configure `CLAWBIO_SCRNASEQ_PATH` or equivalent.

**Deprecation criteria**: Deprecate this skill if nf-core/scrnaseq releases a Python SDK that exposes equivalent preflight, params, and provenance APIs — at that point the wrapper adds less value than it costs in maintenance.

## Citations

- [nf-core/scrnaseq](https://nf-co.re/scrnaseq)
- [Nextflow](https://www.nextflow.io/)
- [Alevin-fry / Simpleaf](https://simpleaf.readthedocs.io/)
