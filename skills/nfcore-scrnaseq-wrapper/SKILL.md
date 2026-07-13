---
name: nfcore-scrnaseq-wrapper
description: Wrapper skill for running nf-core/scrnaseq 4.1.0 upstream single-cell RNA-seq preprocessing from FASTQ with strict preflight, reproducibility outputs, and downstream handoff to ClawBio scRNA
  skills.
license: MIT
metadata:
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
  dependencies:
    python: '>=3.11'
    packages:
    - pyyaml
  endpoints:
    cli: python clawbio.py run scrnaseq-pipeline --input {samplesheet} --output {output_dir}
  openclaw:
    requires:
      bins:
      - python3
      - nextflow
      - java
    always: false
    emoji: 🧫
    homepage: https://github.com/ClawBio/ClawBio
    os:
    - darwin
    - linux
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
  author: ClawBio
  demo_data:
  - path: demo/README.md
    description: Demo mode uses the upstream nf-core/scrnaseq -profile test dataset rather than bundled FASTQs
  inputs:
  - name: samplesheet
    type: file
    format:
    - csv
    description: 'nf-core/scrnaseq samplesheet CSV with required columns: sample, fastq_1, fastq_2'
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
  version: 0.1.0
---

# nfcore-scrnaseq-wrapper

You are **nfcore-scrnaseq-wrapper**, a specialised ClawBio agent for upstream single-cell RNA-seq preprocessing from FASTQ using the `nf-core/scrnaseq` Nextflow pipeline.

## Trigger

**Fire when:**
- User wants to run `scrnaseq` from raw FASTQ files
- User asks to preprocess 10x Chromium single-cell data
- User wants to execute `nf-core/scrnaseq`
- User wants to generate `.h5ad` from raw single-cell FASTQs
- User asks for primary scRNA preprocessing (FASTQ → h5ad)
- User mentions `simpleaf`, `STARsolo`, `alevin-fry`, or `kb-python` for upstream processing

**Do NOT fire when:**
- User already has an `.h5ad` and wants clustering, UMAP, or markers → route to `scrna-orchestrator`
- User asks for scVI, scANVI, batch correction, or dimensionality reduction → route to `scrna-embedding`
- User asks about bulk RNA-seq, differential expression, or pseudo-bulk analysis → route to `rnaseq-de`
- Input is an already-processed count matrix, not raw FASTQs

## Scope

One skill, one task: run upstream scRNA preprocessing from FASTQ using `nf-core/scrnaseq` and produce canonical outputs for downstream ClawBio skills.

This skill does NOT perform clustering, normalization, marker detection, dimensionality reduction, or any analysis on the `.h5ad` it produces.

## Why This Exists

- **Without it**: Users hand-build samplesheets, guess reference combinations, miss backend issues, and struggle to locate the correct `.h5ad` for downstream analysis.
- **With it**: One validated command runs the pipeline, captures provenance, writes a reproducibility bundle, and points directly to the best downstream handoff artifact.
- **Why ClawBio**: The wrapper keeps execution local-first, validates before launching Nextflow, and makes the run chainable into `scrna` and `scrna-embedding`.

## Core Capabilities

1. **Strict Preflight**: Validate Java, Nextflow, backend, samplesheet, FASTQs, and references before execution.
2. **Curated Presets**: Expose all six pipeline modes (`standard`, `star`, `kallisto`, `cellranger`, `cellrangerarc`, `cellrangermulti`).
3. **Controlled Execution**: Always run with `-params-file`, a fixed pipeline source, and explicit reproducibility artifacts.
4. **Output Resolution**: Detect MultiQC, pipeline_info, `.h5ad`, `.rds`, and select a canonical `preferred_h5ad` when possible.
5. **Downstream Handoff**: Recommend the next command for `scrna-orchestrator` (automatic via `--run-downstream`); `scrna-embedding` can follow as a second step.

## Input Formats

| Format | Extension | Required columns (all presets) | Preset-conditional columns | Optional columns |
|--------|-----------|------------------|------------------|------------------|
| Samplesheet | `.csv` | `sample`, `fastq_1`, `fastq_2` | `sample_type` + `fastq_barcode` (required for `cellrangerarc`); `feature_type` (required for `cellrangermulti`) | `expected_cells`, `seq_center` |
| Demo mode | n/a | none — test profile provides its own data | — | — |

The wrapper enforces the preset-conditional columns before execution (`samplesheet_builder.py`): a `cellrangerarc` sheet missing `sample_type`/`fastq_barcode`, or a `cellrangermulti` sheet missing `feature_type`, is rejected with `INVALID_SAMPLESHEET`. Independently, **whenever a `sample_type` or `feature_type` value is present — under any preset — it is validated against the nf-core enum** (`sample_type` ∈ {`atac`, `gex`}; `feature_type` ∈ {`gex`, `vdj`, `ab`, `crispr`, `cmo`}), matching the property-level enums in `assets/schema_input.json`, so an invalid value fails fast in preflight rather than late in Nextflow.

## Workflow

1. **Validate**: Check the selected preset, samplesheet structure, FASTQ accessibility, references, Java, Nextflow, and backend.
2. **Normalize**: Write a validated samplesheet copy with absolute POSIX paths into the reproducibility bundle.
3. **Configure**: Build one effective `params.yaml` and a fixed Nextflow command.
4. **Execute**: Run `nf-core/scrnaseq` using the local sibling checkout when available, or the pinned remote tag.
5. **Parse**: Detect MultiQC, pipeline_info, `.h5ad`, `.rds`, and CellBender-derived outputs.
6. **Generate**: Write `report.md`, `result.json`, provenance JSON files, and reproducibility artifacts.
7. **Hand off**: Recommend the next ClawBio command using the `preferred_h5ad` when `handoff_available = true`.

## Algorithm / Methodology

The wrapper executes a strictly ordered 7-step pipeline. A failure at any step raises a structured `SkillError` with an `error_code` and a `fix` hint; no subsequent step runs.

1. **Pipeline source resolution** (`pipeline_source.py`): Prefer a local sibling `scrnaseq/` checkout (pinned commit, audit-safe). Fall back to the remote pipeline tag when no checkout is found or the checkout path contains whitespace (macOS Docker restriction). A dirty local checkout is rejected by default; `--allow-dirty-pipeline` is an explicit development-only opt-in that is recorded in provenance. Use `--require-local-pipeline` when fallback to the remote pipeline would be unacceptable.

2. **Samplesheet validation** (`samplesheet_builder.py`): Parse the CSV, resolve FASTQ paths relative to the CSV parent directory, normalize sample-name whitespace to underscores, verify readability and FASTQ extensions, reject FASTQ basenames with whitespace, enforce consistent `expected_cells` (≥1) and `seq_center` for repeated sample rows, reject exact duplicate FASTQ rows, and write a normalized copy to `reproducibility/samplesheet.valid.csv` with local FASTQ paths resolved to absolute POSIX (remote URIs, accepted only with `--allow-remote-inputs`, are passed through unchanged). The normalized replay sheet contains only headers accepted by the pinned 4.1.0 `assets/schema_input.json`; unsupported metadata columns are reported and omitted. In particular, `protocol` is supplied through the global `--protocol` parameter and is never forwarded as a samplesheet column.

3. **Preflight** (`preflight.py`): Verify Java (≥17) and Nextflow (≥25.04.0). Compare version tuples after zero-padding to 3 elements (avoids false negatives such as `(24, 4) < (24, 4, 0)`). For `docker`, run `docker info` and gate on exit code. For `conda`/`mamba`, locate the binary. Cell Ranger presets are rejected with conda/mamba unless `--allow-conda-cellranger` is supplied with a trusted site config. For `singularity`/`apptainer`, accept either binary interchangeably. For `wave` and `gpu`, no binary check is needed (Nextflow-native features). Safe institutional profile components are accepted for HPC/site profiles, every `-c/--config` file must exist before execution, and configs are treated as trusted Groovy code. All preflight subprocess calls have a 60-second timeout (`_SUBPROCESS_TIMEOUT` in `preflight.py`); the git probes in `pipeline_source.py` use a 10-second timeout.

4. **Params construction** (`params_builder.py`): Translate the preset + CLI flags into a `params.yaml` consumed by Nextflow via `-params-file`. All file paths use `.as_posix()` for forward-slash consistency across platforms. `igenomes_ignore` is automatically set to `true` whenever an explicit **genome** reference (`fasta`, `gtf`, `transcript_fasta`, `txp2gene`, or any prebuilt index) is provided — auxiliary files (barcode whitelist, CMO/probe/feature sets, primers, multi-barcode samplesheets) never trigger it, so they remain compatible with `--genome` (suppresses nf-schema DNS validation of the default iGenomes S3 URL). Skip flags are only written when `true`, keeping `params.yaml` minimal. In `--demo` mode no reference/protocol params are written at all (the test profile owns them).

5. **Command build + execution** (`command_builder.py`, `executor.py`): Construct the `nextflow run` command with `-params-file`, validated `-c/--config` files, and a work directory that defaults to `<output>/upstream/work` but may be overridden by `--work-dir` (including object-store URIs for cloud executors), then launch via `subprocess.Popen` with stdout and stderr piped to log files on disk — never buffered in RAM. On `TimeoutExpired`, the process is killed and `EXECUTION_FAILED` is raised. On `KeyboardInterrupt`, the child process tree is terminated before the interrupt is re-raised.

6. **Output parsing** (`outputs_parser.py`): Scan the upstream results tree for MultiQC HTML, `pipeline_info/`, aligner output directories, `.h5ad` (CellBender/filtered preferred over generic combined/raw), `.rds`, CellBender-derived files, and an `official_outputs` manifest for documented nf-core output families. Required outputs are validated before success artifacts are written; `handoff_available` is set to `true` only when a `preferred_h5ad` is confirmed on disk.

7. **Provenance + reporting** (`provenance.py`, `reporting.py`): Write JSON provenance bundles, a SHA-256 checksum manifest (files only — never directories), `environment.yml`, a portable `commands.sh`, `report.md`, and `result.json`.

## Presets

| Preset | Aligner | Use case |
|--------|---------|---------|
| `standard` | simpleaf (alevin-fry) | Default for 10x GEX; fast, memory-efficient |
| `star` | STARsolo | Best FASTQ QC metrics; supports RNA velocity (`--star-feature "Gene Velocyto"`) |
| `kallisto` | kb-python / BUStools | Pseudo-alignment; fastest; lamanno/nac RNA velocity via `--kb-workflow` |
| `cellranger` | CellRanger | CellRanger v2/v3 compatibility; CellRanger is provided by the nf-core container under `docker`/`singularity` (no host binary needed). Not available under `-profile conda` (10x licensing keeps it off bioconda) |
| `cellrangerarc` | CellRanger ARC | Multiome (GEX + ATAC); accepts prebuilt `--cellranger-index` or reference-build inputs |
| `cellrangermulti` | CellRanger Multi | GEX + VDJ + feature barcoding; `--cellranger-multi-barcodes` required for CMO/FFPE multiplexing |

Each preset requires at least one reference option: `--genome <iGenomes_shortcut>` OR a pre-built index (`--star-index`, `--simpleaf-index`, etc.) OR `--fasta` + `--gtf`. The `standard`/simpleaf preset additionally accepts a transcriptome pair `--transcript-fasta` + `--txp2gene` in place of a genome reference (per the nf-core/scrnaseq Simpleaf options).

## nf-core/scrnaseq 4.1.0 Compatibility Policy

This wrapper targets nf-core/scrnaseq `4.1.0`. It is not a free-form passthrough. Parameters are grouped as:

- **Supported upstream parameters:** input/output, aligner/preset, reference/index, skip, CellRanger, CellRanger ARC, CellRanger Multi, selected MultiQC/reporting options.
- **Wrapper policy parameters:** `--preset`, `--check`, `--run-downstream`, `--skip-downstream`, `--expected-cells`, `--timeout-hours`, `--work-dir`, `--allow-remote-inputs`, `--allow-dirty-pipeline`, `--require-local-pipeline`, `--allow-pipeline-version-override`, `--trust-config-params`, `--allow-conda-cellranger`, and `-c`/`--config`/`--nextflow-config`; these are ClawBio conveniences and are not nf-core parameters.
- **Deprecated compatibility aliases:** `skip_emptydrops` is accepted only as `--skip-emptydrops` and translated to `skip_cellbender: true`; the deprecated upstream parameter is never written.
- **Intentionally unsupported upstream parameters:** `custom_config_version`, `custom_config_base`, `config_profile_name`, `config_profile_description`, `config_profile_contact`, `config_profile_url`, `version`, `plaintext_email`, `max_multiqc_email_size`, `hook_url`, `validate_params`, `pipelines_testdata_base_path`, `help`, `help_full`, `show_hidden`.

Unsupported parameters are either hidden/institutional metadata, interactive help/version flags, or options that would weaken the wrapper's fixed validation/reproducibility policy.

## Input & Reference Path Policy

**Local-first by default.** Samplesheet FASTQs and reference/index inputs must be local paths unless you explicitly opt in. Remote URIs (`s3://`, `gs://`, `https://`, `ftp://`, …) are **rejected** at preflight with `REMOTE_INPUT_NOT_ALLOWED`, so genetic data and references stay on the local machine and no accidental cloud fetch happens. This guarantee is enforced by the code, not just advertised (`preflight._check_remote_inputs`).

**Opt-in for remote inputs.** Pass `--allow-remote-inputs` to permit remote samplesheet inputs and reference paths (parity with `nfcore-sarek-wrapper` / `nfcore-rnaseq-wrapper`, which share the same flag). When enabled, remote URIs are passed through **verbatim** (Nextflow resolves and stages them; only the FASTQ/FASTA basename is validated) and preflight emits a **runtime WARNING** listing every path that will be fetched over the network, so cloud access is always visible. The object-store `--work-dir` is a separate setting and is not gated.

**Local** paths are still validated eagerly at preflight so they fail fast with a clear error instead of a late Nextflow error:
- A supplied local reference/index path (`--fasta`, `--gtf`, `--star-index`, …) that does not exist raises `MISSING_REFERENCE` (`preflight.py`).
- A local FASTQ that does not exist (or is not a regular file) raises `MISSING_FASTQ` (`samplesheet_builder.py`).

Readability is never pre-checked: Nextflow reads inputs in the true execution context (often a root container under the default Docker profile), so a launcher-side `os.access(R_OK)` probe would false-block valid runs (`errors.py`).

## CLI Reference

```bash
# Standard real-data usage (explicit protocol and reference are required)
python skills/nfcore-scrnaseq-wrapper/nfcore_scrnaseq_wrapper.py \
  --input samplesheet.csv --output ./scrnaseq_run \
  --preset star --protocol 10XV3 --genome GRCh38

# Preflight check only (no Nextflow execution)
python skills/nfcore-scrnaseq-wrapper/nfcore_scrnaseq_wrapper.py \
  --input samplesheet.csv --output ./scrnaseq_run --check \
  --preset star --protocol 10XV3 --genome GRCh38

# Demo mode (runs the upstream nf-core test profile; forces star preset; uses the
# selected backend — default --profile docker, which must be running)
python skills/nfcore-scrnaseq-wrapper/nfcore_scrnaseq_wrapper.py \
  --demo --output ./scrnaseq_demo

# Via ClawBio runner
python clawbio.py run scrnaseq-pipeline --input samplesheet.csv --output ./scrnaseq_run \
  --preset star --protocol 10XV3 --genome GRCh38
python clawbio.py run scrnaseq-pipeline --demo --output ./scrnaseq_demo

# STARsolo with local FASTA+GTF (STAR index built by the pipeline)
python skills/nfcore-scrnaseq-wrapper/nfcore_scrnaseq_wrapper.py \
  --input samplesheet.csv --output ./run --preset star --protocol 10XV3 \
  --fasta /refs/hg38.fa --gtf /refs/hg38.gtf

# STARsolo with prebuilt STAR index
python skills/nfcore-scrnaseq-wrapper/nfcore_scrnaseq_wrapper.py \
  --input samplesheet.csv --output ./run --preset star --protocol 10XV3 \
  --star-index /refs/star_index

# STARsolo RNA velocity (star requires an explicit --protocol like every star/standard/kallisto run)
python skills/nfcore-scrnaseq-wrapper/nfcore_scrnaseq_wrapper.py \
  --input samplesheet.csv --output ./run --preset star --protocol 10XV3 \
  --star-feature "Gene Velocyto" --star-ignore-sjdbgtf \
  --fasta /refs/hg38.fa --gtf /refs/hg38.gtf

# Simpleaf (standard) with UMI resolution override
python skills/nfcore-scrnaseq-wrapper/nfcore_scrnaseq_wrapper.py \
  --input samplesheet.csv --output ./run --preset standard --protocol 10XV3 \
  --simpleaf-umi-resolution cr-like-em --genome GRCh38

# Kallisto RNA velocity (NAC workflow)
python skills/nfcore-scrnaseq-wrapper/nfcore_scrnaseq_wrapper.py \
  --input samplesheet.csv --output ./run --preset kallisto --protocol 10XV3 \
  --kb-workflow nac --fasta /refs/hg38.fa --gtf /refs/hg38.gtf

# Air-gapped cluster: local iGenomes mirror
python skills/nfcore-scrnaseq-wrapper/nfcore_scrnaseq_wrapper.py \
  --input samplesheet.csv --output ./run --preset star --protocol 10XV3 \
  --genome GRCh38 --igenomes-base /mnt/local_igenomes

# CellRanger Multi (CMO multiplexing)
python skills/nfcore-scrnaseq-wrapper/nfcore_scrnaseq_wrapper.py \
  --input samplesheet.csv --output ./run --preset cellrangermulti \
  --cellranger-index /refs/refdata-gex-GRCh38 \
  --gex-cmo-set /refs/cmo_set.csv \
  --cellranger-multi-barcodes /refs/multi_barcodes.csv
```

### Key flags

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `--input` | path | — | Samplesheet CSV (`sample,fastq_1,fastq_2`). Required unless `--demo` |
| `--output` | path | — | Output directory for results and the reproducibility bundle (required) |
| `--demo` | flag | — | Run the upstream nf-core `test` profile; forces `--preset star` and `skip_cellbender` |
| `--check` | flag | — | Run preflight only and exit (no Nextflow execution) |
| `--preset` | string | `standard` | Aligner preset |
| `--aligner` | string | — | nf-core/scrnaseq aligner alias for `--preset`: `simpleaf` maps to `standard`; the other values map to same-named presets |
| `--profile` | string | `docker` | Execution backend or comma-separated nf-core profile list: `docker`, `conda`, `mamba`, `singularity`, `apptainer`, `podman`, `shifter`, `charliecloud`, `wave`, `gpu`, `debug`, `arm64`, `emulate_amd64`, `test`, `test_full`, `test_cellrangermulti`, `test_multiome` |
| `--pipeline-version` | string | `4.1.0` | Remote `nf-core/scrnaseq` tag or commit (used when no local sibling checkout is found) |
| `--allow-dirty-pipeline` | flag | — | Development-only opt-in to run a modified local sibling `scrnaseq/` checkout; rejected by default for production reproducibility |
| `--require-local-pipeline` | flag | — | Require a verifiable local sibling `scrnaseq/` checkout; fail instead of falling back to the remote pipeline |
| `--allow-pipeline-version-override` | flag | — | Allow a `--pipeline-version` other than the pinned `4.1.0` contract (warned, recorded in provenance; validations stay 4.1.0) |
| `--trust-config-params` | flag | — | Allow `-c/--config` files that set `params.*` (otherwise blocked); detected overrides are recorded in provenance |
| `-c` / `--config` | path | — | Additional Nextflow config file. May be repeated; files are validated, applied to the live run, copied into the reproducibility bundle, and replayed by `commands.sh` |
| `--protocol` | string | `None` | Chemistry/protocol forwarded to the aligner. Required for `standard`, `star`, and `kallisto`; omitted only preserves CellRanger auto-detection for `cellranger`, `cellrangerarc`, and `cellrangermulti`. Explicit `auto` is invalid for `standard`, `star`, and `kallisto`. `smartseq` is valid for `star` and `kallisto` only. Other protocol strings are mapped when known or passed through by nf-core |
| `--genome` | string | — | iGenomes shortcut (`GRCh38`, `mm10`, etc.) — mutually exclusive with `--fasta`/`--gtf` and all index flags |
| `--igenomes-base` | string | — | Base URL or local path for iGenomes (default `s3://ngi-igenomes/igenomes/`). Use for local mirrors or air-gapped clusters. A **local** base path is existence-checked in preflight when `--genome` is set (remote `s3://`/`https://` bases are deferred to Nextflow) |
| `--igenomes-ignore` | flag | — | Do not load the iGenomes reference config. Set automatically whenever an explicit genome reference is supplied; only needed manually in unusual setups |
| `--fasta` | path | — | Genome FASTA (`.fa`, `.fna`, `.fasta`, `.gz` variants; no whitespace in path) |
| `--gtf` | path | — | Gene annotation GTF |
| `--star-index` | path | — | Prebuilt STAR genome index directory |
| `--simpleaf-index` | path | — | Prebuilt simpleaf/alevin-fry index |
| `--kallisto-index` | path | — | Prebuilt kallisto index |
| `--cellranger-index` | path | — | Prebuilt CellRanger or CellRanger ARC reference |
| `--transcript-fasta` | path | — | Transcriptome FASTA for simpleaf |
| `--txp2gene` | path | — | Transcript-to-gene mapping for simpleaf |
| `--barcode-whitelist` | path | — | Custom barcode whitelist (per-aligner format) |
| `--star-feature` | enum | — | STARsolo feature type: `Gene`, `GeneFull`, `Gene Velocyto` |
| `--star-ignore-sjdbgtf` | flag | — | Do not use GTF for SJDB construction (required for `Gene Velocyto`) |
| `--seq-center` | string | — | Sequencing center name for BAM read group tag |
| `--simpleaf-umi-resolution` | enum | — | UMI resolution strategy for alevin-fry: `cr-like`, `cr-like-em`, `parsimony`, `parsimony-em`, `parsimony-gene`, `parsimony-gene-em` |
| `--kb-workflow` | enum | — | Kallisto workflow: `standard`, `lamanno`, `nac` |
| `--kb-t1c` | path | — | cDNA transcripts-to-capture file for RNA velocity (lamanno/nac). Required only with a prebuilt `--kallisto-index`; auto-generated from `--fasta`/`--gtf` |
| `--kb-t2c` | path | — | Intron transcripts-to-capture file for RNA velocity (lamanno/nac). Required only with a prebuilt `--kallisto-index`; auto-generated from `--fasta`/`--gtf` |
| `--skip-cellbender` | flag | — | Disable the CellBender ambient RNA removal subworkflow |
| `--skip-emptydrops` | flag | — | Deprecated compatibility alias for `--skip-cellbender`; the wrapper writes `skip_cellbender: true` and never writes deprecated upstream `skip_emptydrops` |
| `--skip-fastqc` | flag | — | Skip FastQC quality control |
| `--skip-multiqc` | flag | — | Skip MultiQC report generation |
| `--skip-cellranger-renaming` | flag | — | Skip automatic sample renaming in CellRanger modules |
| `--skip-cellrangermulti-vdjref` | flag | — | Skip mkvdjref in cellrangermulti (when VDJ data is absent or a prebuilt `--cellranger-vdj-index` is supplied) |
| `--save-reference` | flag | — | Save the built reference index for future reuse |
| `--save-align-intermeds` / `--no-save-align-intermeds` | flag | — | Forward `save_align_intermeds: true`/`false`; when neither is given the upstream nf-core/scrnaseq 4.1.0 default (`true`) is preserved, so intermediate BAMs are published by default — pass `--no-save-align-intermeds` on large runs to save disk |
| `--expected-cells` | int | — | Override expected cell count for a single-sample samplesheet; multi-sample runs must set `expected_cells` per row. The wrapper enforces `≥1` (stricter than the upstream `integer` schema, which has no minimum) since a non-positive count is meaningless |
| `--timeout-hours` | float | `12` | Wall-clock cap for the Nextflow run, in hours. Use `0` to disable the cap for long HPC/cloud runs whose walltime is enforced by the scheduler. Via the ClawBio runner the runner's own timeout also applies |
| `--work-dir` | string | `<output>/upstream/work` | Nextflow work directory. Local paths are resolved before execution; object-store URIs such as `s3://...` or `gs://...` are passed through for cloud executors |
| `--allow-conda-cellranger` | flag | — | Allow Cell Ranger presets with conda/mamba only when a trusted site config supplies Cell Ranger |
| `--email` | string | — | Email address for pipeline completion notification |
| `--email-on-fail` | string | — | Email address for pipeline failure notification |
| `--multiqc-title` | string | — | Custom title for the MultiQC report |
| `--multiqc-config` | path | — | Custom MultiQC config YAML |
| `--multiqc-logo` | path | — | Custom MultiQC logo image |
| `--multiqc-methods-description` | path | — | Custom MultiQC methods-description YAML |
| `--publish-dir-mode` | enum | — | Forwarded nf-core publish mode: `symlink`, `rellink`, `link`, `copy`, `copyNoFollow`, or `move` |
| `--trace-report-suffix` | string | — | Suffix for Nextflow trace/report/timeline filenames |
| `--monochrome-logs` | flag | — | Disable ANSI colors in nf-core logs |
| `--resume` | flag | — | Nextflow resume (checksum-verified against prior manifest; preset/profile/source/work-dir must match) |
| `--run-downstream` | flag | — | Opt in to `scrna_orchestrator` handoff after pipeline completion |
| `--skip-downstream` | flag | — | Force-skip the downstream handoff even if `--run-downstream` is given (handoff is off by default) |
| `--cellrangerarc-config` | path | — | Config JSON for CellRanger ARC index construction |
| `--cellrangerarc-reference` | string | — | Reference genome name used inside the CellRanger ARC config |
| `--motifs` | path | — | Motif file (e.g. JASPAR) for CellRanger ARC |
| `--cellranger-vdj-index` | path | — | Prebuilt CellRanger VDJ reference |
| `--gex-frna-probe-set` | path | — | Probe set CSV for FFPE fixed RNA profiling (`cellrangermulti`) |
| `--gex-target-panel` | path | — | Target panel CSV for targeted GEX (`cellrangermulti`) |
| `--gex-cmo-set` | path | — | CMO reference CSV for multiplexed samples (`cellrangermulti`) |
| `--gex-barcode-sample-assignment` | path | — | Barcode-to-sample assignment override CSV (`cellrangermulti`). **Not** an OCM selector — OCM mode is encoded via the `ocm_ids` column of `--cellranger-multi-barcodes` |
| `--fb-reference` | path | — | Feature-barcode reference CSV for antibody capture (`cellrangermulti`) |
| `--vdj-inner-enrichment-primers` | path | — | V(D)J cDNA enrichment primer sequences (`cellrangermulti`) |
| `--cellranger-multi-barcodes` | path | — | Multiplexed sample samplesheet for CMO/FFPE demultiplexing (`cellrangermulti`) |

## Output Structure

```text
output_directory/
├── report.md                         # Wrapper run summary
├── result.json                       # Structured result payload
├── check_result.json                 # Written only with --check (preflight-only mode); no upstream/ is produced
├── logs/
│   ├── stdout.txt                    # Nextflow stdout
│   └── stderr.txt                    # Nextflow stderr
├── upstream/
│   └── results/                      # nf-core/scrnaseq output tree
│       ├── fastqc/                   # Per-read FastQC reports
│       ├── multiqc/                  # MultiQC HTML and data
│       │   └── multiqc_report.html
│       ├── pipeline_info/            # Execution report, timeline, trace, DAG
│       └── <aligner>/                # Aligner-specific outputs
│           └── mtx_conversions/      # AnnData (.h5ad), SCE (.rds), Seurat (.rds)
│               │  # Per nf-core/scrnaseq 4.1.0 conf/modules.config, the concatenated
│               │  # matrices sit at the TOP of mtx_conversions/ while each per-sample
│               │  # matrix is nested one level deeper under <sample>/ (MTX_TO_H5AD/
│               │  # CONCAT_H5AD/ANNDATA_BARCODES rewrite non-combined files to
│               │  # `${meta.id}/${filename}`). The wrapper scans BOTH depths and does
│               │  # not hard-code filenames: it ranks whatever .h5ad it finds to pick
│               │  # one preferred_h5ad (combined > per-sample; within a group
│               │  # cellbender_filter > filtered > plain > raw, matched on the filename
│               │  # suffix). combined_* filter-suffixed variants are version-dependent.
│               ├── combined_matrix.h5ad                       # documented concatenated matrix (top level)
│               ├── combined_filtered_matrix.h5ad              # version-dependent: filtering ran
│               ├── combined_cellbender_filter_matrix.h5ad     # version-dependent: CellBender ran → top preference
│               └── <sample>/                                  # per-sample matrices nested one level deeper
│                   ├── <sample>_raw_matrix.h5ad
│                   └── <sample>_filtered_matrix.h5ad          # conditional: filtering ran
├── reproducibility/                  # Reproducibility + provenance bundle (single directory)
│   ├── samplesheet.valid.csv         # Normalized samplesheet (absolute POSIX paths); named samplesheet.demo.csv with --demo
│   ├── params.yaml                   # Effective Nextflow parameters
│   ├── nextflow_configs/             # Written only when -c/--config is used: copies of the supplied config files, replayed by commands.sh
│   ├── commands.sh                   # Portable replay script
│   ├── environment.yml               # Conda environment spec (for reference)
│   ├── checksums.sha256              # SHA-256 for in-bundle artifacts only; `sha256sum -c` passes from output dir
│   ├── manifest.json                 # Run metadata: preset, profile, versions, checksums
│   ├── macos_docker.config           # macOS+Docker workarounds (VirtioFS, ARM64, STAR FIFOs)
│   ├── remap_paths.py                # Helper for replaying on a different machine
│   ├── compatibility_policy.json     # Copied policy snapshot (resume/update rules)
│   ├── pinned_versions.json          # Copied pinned-versions snapshot
│   ├── inputs.json                   # Samplesheet + fastq/reference paths and digests (reference_checksums)
│   ├── invocation.json               # Timestamp, preset, profile, pipeline version
│   ├── preflight.json                # Java/Nextflow/backend info
│   ├── upstream.json                 # Pipeline source resolution details
│   ├── outputs.json                  # Detected artifacts
│   ├── runtime.json                  # Execution timing
│   └── skill.json                    # Skill name and version
├── provenance/                       # Written only with --run-downstream
│   └── handoff.json                  # Downstream orchestrator path + checksum + outcome
└── scrna_analysis/                   # Written only with --run-downstream: scrna_orchestrator output
```

## Example Output

`result.json` (abbreviated):
```json
{
  "skill": "scrnaseq-pipeline",
  "version": "0.1.0",
  "summary": {
    "preset": "star",
    "aligner_effective": "star",
    "pipeline_source_kind": "remote_repo",
    "pipeline_version_or_commit": "4.1.0",
    "profile": "docker",
    "preferred_h5ad": "<output>/upstream/results/star/mtx_conversions/combined_filtered_matrix.h5ad",
    "handoff_available": true,
    "samples_detected": 2,
    "cellbender_used": false
  }
}
```

`report.md` closes with:
```
## Next Steps
- python clawbio.py run scrna --input <preferred_h5ad> --output <dir>
- python clawbio.py run scrna-embedding --input <preferred_h5ad> --output <dir>
```

## Gotchas

- **Preflight runs before any Nextflow call.** If Java, Nextflow, or the backend are missing or too old, the pipeline never starts and you get a structured JSON error with `error_code` and a `fix` hint. Nextflow ≥25.04.0 is required.
- **Conda/Mamba profiles need network access at runtime.** Preflight only verifies the `conda`/`mamba` binary exists — it deliberately does not probe network connectivity (a flaky check that would also reject valid offline caches). With `-profile conda`, Nextflow resolves environments from bioconda/conda-forge at runtime, so an offline or proxied host can make a task fail mid-run with a confusing resolver error. Cell Ranger presets are rejected with conda/mamba by default because Cell Ranger is not distributed via bioconda; use `--allow-conda-cellranger` only when a trusted site config provides Cell Ranger. For fully offline/reproducible execution prefer `docker`/`singularity` (pinned containers), or pre-warm the conda package cache before running.
- **`--genome` conflicts with any explicit genome-reference flag.** Providing `--genome` alongside `--fasta`, `--gtf`, `--transcript-fasta`, `--txp2gene`, or any prebuilt index raises `CONFLICTING_REFERENCES` in preflight. Use either `--genome <shortcut>` or explicit genome flags — never both. Auxiliary files (barcode whitelist, CMO/probe/feature sets, primers, multi-barcode samplesheets) are **not** genome references and are compatible with `--genome`.
- **`igenomes_ignore` is set automatically.** Whenever an explicit **genome** reference (`fasta`, `gtf`, `transcript_fasta`, `txp2gene`, or any prebuilt index) is provided, the wrapper writes `igenomes_ignore: true` to suppress nf-schema DNS validation of the default iGenomes S3 URL. Auxiliary files never trigger it. You do not need to set this manually. Use `--igenomes-base` only for local iGenomes mirrors.
- **Protocol compatibility is enforced before Nextflow starts.** `standard`, `star`, and `kallisto` require an explicit `--protocol`, because nf-core/scrnaseq 4.1.0 documents `auto` as CellRanger-only. Explicit `auto` is rejected for those presets. `smartseq` denotes **Smart-seq3** and is accepted for `star` and `kallisto` only; **`smartseq2` is rejected for every preset** because nf-core/scrnaseq 4.1.0 does not support Smart-seq2. `cellranger` accepts only `auto` and `10XV1`-`10XV4`; `cellrangerarc` accepts only `auto`; `cellrangermulti` is samplesheet-driven. Unknown custom protocol strings are passed through only for `standard`, `star`, and `kallisto`, where nf-core documents custom values.
- **`save_align_intermeds` defaults to ON in 4.1.0.** Upstream nf-core/scrnaseq 4.1.0 sets `save_align_intermeds: true` by default, so aligner intermediate BAMs are published into the results tree unless you pass `--no-save-align-intermeds`. On real (non-`--demo`) datasets these intermediates can be tens to hundreds of GB; disable them when you only need the count matrices.
- **Deprecated upstream parameters:** `skip_emptydrops` is deprecated in nf-core/scrnaseq 4.1.0. The wrapper accepts `--skip-emptydrops` only as a compatibility alias and emits `skip_cellbender: true`.
- **`--demo` forces preset=star and skip_cellbender=true.** The nf-core upstream `test` profile ships STAR-compatible data and explicitly disables CellBender (which does not work on small test datasets). If a different preset is requested with `--demo`, the wrapper warns and overrides it. **`--demo` also requires network access**: the test profile fetches its FASTQs and references from remote GitHub URLs (nf-core's design — no local test data is bundled), so on an offline/sandboxed host set `NXF_OFFLINE` and the wrapper fails fast at preflight with `DEMO_REQUIRES_NETWORK` instead of a cryptic Nextflow `does not exist` abort. This does not violate local-first (which governs your *genetic data* — never uploaded); `--demo` only *downloads* nf-core's public test data. For a fully offline run, use your own local `--input` and references.
- **`--demo` is fully hermetic — every pipeline flag is ignored.** The `test` profile owns the entire pipeline configuration: input, references, protocol, and all QC/skip/tuning/save/reporting knobs. Because a `-params-file` value overrides profile config in Nextflow, the wrapper writes **only** the four forced essentials into `params.yaml` in demo mode — `outdir`, `aligner` (`star`), `igenomes_ignore`, and `skip_cellbender` — and nothing else. Any other flag you pass with `--demo` (e.g. `--genome`, `--fasta`, an index, `--igenomes-base`, `--protocol`, `--star-feature`, `--skip-fastqc`, `--save-reference`, `--seq-center`, `--publish-dir-mode`, …) is ignored and listed in a WARNING. Demo output validation likewise requires FastQC and MultiQC regardless of any `--skip-*` you pass, since those flags are not written. Drop `--demo` to run on your own inputs.
- **The Nextflow run has a 12 h wall-clock cap by default.** Large multi-sample STARsolo/CellRanger runs on full genomes can exceed this; raise it with `--timeout-hours <n>` or pass `--timeout-hours 0` to disable the cap entirely (recommended on HPC/cloud where the scheduler enforces walltime). When the run targets an object-store work directory (`--work-dir s3://…`/`gs://…`) or an institutional/site profile while the cap is still active, preflight prints a WARNING reminding you to pass `--timeout-hours 0`, since the scheduler — not the wrapper — should bound such runs. When neither the cap fires nor the job is killed, behaviour is unchanged.
- **Required outputs are checked after Nextflow exits.** The wrapper expects `pipeline_info/`, the effective aligner directory, MultiQC unless `--skip-multiqc`, FastQC HTML/ZIP reports unless `--skip-fastqc`, and at least one `.h5ad` matrix. FastQC is a hard requirement for **every** aligner — including the Cell Ranger family (`cellranger`/`cellrangerarc`/`cellrangermulti`) — because nf-core/scrnaseq 4.1.0 runs FASTQC on the shared input-read channel (`ch_fastq`) *before* any aligner-specific branching and publishes `fastqc/` for all of them (the output docs state "FastQC is applied to all aligners' input reads"). A missing `fastqc/` tree on a non-skipped run is therefore a genuine failure, not a tolerated gap. All six presets run the pipeline's `MTX_TO_H5AD` conversion (CellRanger ARC/Multi reuse the CellRanger template), so a completed run with **zero `.h5ad` is treated as a failure** (`EXPECTED_OUTPUTS_NOT_FOUND`) rather than a silent success. A present-but-ambiguous selection (e.g. several per-sample matrices, no combined) does not fail the run; it is signalled by `handoff_available = false`, and `--run-downstream` prints a warning instead of silently returning.
- **`preferred_h5ad` may be absent.** If no combined matrix is produced and there are multiple per-sample files, `handoff_available` is `false`. Always check `result.json` before chaining to `scrna-orchestrator` or `scrna-embedding`.
- **No arbitrary Nextflow parameter passthrough.** Pipeline *parameters* flow only through the preset system and `params.yaml`; no custom `--param` flags can be injected. Validated `-c/--config` files are allowed for infrastructure/HPC configuration and are copied into the reproducibility bundle. Note this is a **trust boundary, not a sandbox**: a Nextflow config is executable Groovy (it can set `process.beforeScript`, `process.shell`, etc.), so a `-c` file is as trusted as code you run yourself. The wrapper validates that each config exists and lints it for `params.*` overrides (see the next gotcha), but otherwise does not sandbox its contents — only pass configs you authored or trust.
- **The work directory defaults to `<output>/upstream/work`.** This preserves local resume and bundle portability for workstation/HPC runs. Managed cloud-batch executors that need an object-store work directory can pass `--work-dir s3://...` or `--work-dir gs://...`; the value is recorded in provenance and replayed by `commands.sh`. Local custom paths are resolved before execution and are less portable than the default.
- **The published results directory (`outdir`) is intentionally the local `<output>/upstream/results`.** This is a wrapper policy, not an nf-core limitation: the wrapper parses the results tree to detect the `preferred_h5ad`, MultiQC, and `pipeline_info/`, and to hand off downstream, so results must land on the local filesystem. This is still compatible with cloud-batch executors — the work directory may be remote (`--work-dir s3://...`) while Nextflow copies the published results back to the local `outdir`. There is no `--outdir` override to an object store, because the wrapper could not then parse outputs or chain downstream; for a remote-published run, invoke the upstream pipeline directly. The `--output` directory must also be **outside** the ClawBio source tree — an output inside the repository is rejected at preflight with `OUTPUT_DIR_INSIDE_REPO` (parity with nfcore-rnaseq/nfcore-sarek) so multi-gigabyte pipeline artifacts never pollute the checkout.
- **`--pipeline-version` is pinned to the 4.1.0 contract.** The wrapper's parameter set, protocol matrix and output validation are written for nf-core/scrnaseq 4.1.0, so a different `--pipeline-version` is rejected unless you pass `--allow-pipeline-version-override` (a warned, provenance-recorded opt-in; validations stay 4.1.0).
- **`-c/--config` files may not set `params.*`.** nf-core advises against setting parameters via `-c`. The wrapper blocks any config that assigns parameters outside the audited `params.yaml` — dotted (`params.aligner = …`), bracket (`params['aligner'] = …`), whole-map (`params = […]`), or a `params { … }` block (including the `{` on the next line). Pass `--trust-config-params` to allow it (detected overrides are recorded in provenance). Note the config is still trusted Groovy and is not sandboxed; this lint only guarantees `params.yaml` stays the single parameter source.
- **`--resume` enforces strict compatibility.** The wrapper checks that the stored manifest matches the current preset, profile, pipeline source, effective `params.yaml` checksum, and Nextflow work directory. Mismatches raise `INVALID_RESUME_STATE`.
- **RNA velocity requires coordinated flags and is validated in preflight.** For STARsolo: `--star-feature "Gene Velocyto"` requires `--star-ignore-sjdbgtf`. For Kallisto `--kb-workflow lamanno`/`nac`: when the index is built from `--fasta`/`--gtf`, `kb ref` generates the capture files, so `--kb-t1c`/`--kb-t2c` are **not** required (matching the documented fasta+gtf example); they are required only alongside a prebuilt `--kallisto-index`, and supplying exactly one of the two is always rejected. Invalid combinations raise `INVALID_PRESET_CONFIGURATION` before Nextflow starts.
- **`cellrangerarc` config and reference are paired.** Providing `--cellrangerarc-config` without `--cellrangerarc-reference` (or vice versa) raises `INVALID_PRESET_CONFIGURATION`. `--motifs` is optional and independent.
- **`cellrangermulti` validates documented preflight constraints.** `feature_type=ab` **or** `feature_type=crispr` requires `--fb-reference` (in nf-core/scrnaseq 4.1.0 both Antibody Capture and CRISPR Guide Capture share the single `fb_reference` feature reference — there is no separate CRISPR reference param, so a crispr run without it has no feature reference and fails inside Cell Ranger); `feature_type=vdj` with `--skip-cellrangermulti-vdjref` requires `--cellranger-vdj-index`; `--gex-frna-probe-set` (FFPE) requires `--cellranger-multi-barcodes`; `feature_type=cmo` (CMO multiplexing) also requires `--cellranger-multi-barcodes`. The three multiplexing modes nf-core documents as mutually exclusive are encoded **per multiplexed sample** in the `--cellranger-multi-barcodes` samplesheet via the `probe_barcode_ids` (**FFPE**), `cmo_ids` (**CMO**) and `ocm_ids` (**OCM**) columns. The wrapper parses that samplesheet and rejects any physical `sample` that populates more than one of those columns (`INVALID_PRESET_CONFIGURATION`, with `conflicting_modes`). It additionally rejects the run-level FFPE↔CMO flag combination (`--gex-frna-probe-set` with `feature_type=cmo`/`--gex-cmo-set`). Note `--gex-barcode-sample-assignment` is a cell/tag-calling override, **not** an OCM selector, so it is never treated as a mode. `cellrangermulti` FASTQ filenames are **not** validated by the wrapper (unlike `cellranger`/`cellrangerarc`): Cell Ranger Multi maps libraries through its own multi samplesheet and `[libraries]` config, so the per-file 10x naming is left to Cell Ranger to validate.
- **FASTA schema validation.** The FASTA path must match `^\S+\.fn?a(sta)?(\.gz)?$` (the nf-core/scrnaseq 4.1.0 schema). Paths with whitespace or non-standard extensions are rejected in preflight.
- **Local checkout must be a sibling directory.** The wrapper looks for `../scrnaseq` relative to the ClawBio repo root. If the checkout path contains whitespace (common on macOS), the wrapper warns and falls back to the remote pipeline. Ensure `ClawBio/` and `scrnaseq/` share the same parent folder.
- **Local checkout version is verified.** A sibling checkout is used only when its exact git tag matches `--pipeline-version` or its commit equals `--pipeline-version`. Git-less or unverifiable local directories fall back to the remote pinned pipeline.
- **The `validation.*` / `validationSchemaIgnoreParams` warnings on a run or a bundle replay come from the pipeline, not from the wrapper. Do NOT try to silence them by adding params.** nf-core/scrnaseq 4.1.0's own `nextflow.config` declares `plugins { id 'nf-schema@2.5.1' }` and a matching `validation { defaultIgnoreParams = ["genomes"]; monochromeLogs = params.monochrome_logs }` scope, and it carries a `monochromeLogs = null` param that upstream itself annotates `// TODO temporary workaround a warning` (with a matching ignored lint check in its `.nf-core.yml`). nf-schema surfaces that as `WARN: The following invalid input values have been detected: * --validationSchemaIgnoreParams: genomes`. The wrapper writes **none** of these keys — its `params.yaml` carries only the audited CLI surface — and the bundle's `commands.sh` runs the pinned pipeline unmodified, so the same warning appears on a plain upstream `nextflow run nf-core/scrnaseq -r 4.1.0` too. It is harmless upstream noise: the plugin loads correctly (`Plugins declared=[nf-schema@2.5.1]`) and validation still runs. Writing `validation`/`monochromeLogs` into `params.yaml` to quiet it would override the pipeline's own config and is exactly the wrong fix — the wrapper's "omit = trust upstream default" rule applies. A regression test pins that these keys never enter `params.yaml`.
- **macOS Docker workaround is applied automatically.** On macOS with Docker, `macos_docker.config` is written to the reproducibility bundle and passed to Nextflow. It sets `stageInMode = "copy"` (avoids VirtioFS EDEADLK), `--platform linux/amd64` (Rosetta emulation), and routes STAR `_STARtmp` to the container's `/tmp` (avoids VirtioFS FIFO limitation). These three workarounds apply to **every** macOS+Docker run. A `resourceLimits` block (`cpus 4, memory 15.GB`) is written **only for `--demo` runs**, with `time` raised to `4.h` (the upstream `conf/test.config` caps tasks at `1.h`, too short for STAR genome generation under Apple-Silicon emulation). These ceilings are test-machine sized and are never applied to real datasets (a human STAR index needs far more than 15 GB). One consequence: a `--demo` run is capped at `4.h` on macOS+Docker but at the upstream `1.h` on Linux. **On non-macOS docker, a real (non-demo) run instead writes a host-scaled `process.resourceLimits` config** (physical RAM minus headroom, bounded by the Docker `MemTotal` when known) so the local executor does not abort with `Process requirement exceeds available memory` when a process's production request — e.g. STARsolo's index build — exceeds this machine. This follows nf-core's resourceLimits guidance (limits should match the machine maximum) and carries no 15 GB ceiling. It is written to `reproducibility/resource_limits.config` and shipped in the bundle; `commands.sh` re-applies it on non-macOS replay hosts so a from-scratch reproduction on the generating machine does not re-abort on memory. Override with your own `-c process { resourceLimits = [...] }` if needed. Output directories under `/tmp` emit a WARNING — use a path under HOME.

## Safety

- **Local-first by default**: user FASTQs/references and outputs stay on the local filesystem. Remote input/reference URIs are rejected (`REMOTE_INPUT_NOT_ALLOWED`) unless `--allow-remote-inputs` is explicitly passed, which also logs a runtime warning naming every path fetched over the network. `--allow-remote-inputs` relaxes only the wrapper's own preflight check: remote FASTQ/reference URIs are then written into the normalized samplesheet/`params.yaml` **verbatim** and staged natively by Nextflow at run time. The wrapper does not download them itself, so remote inputs require outbound network access and are incompatible with `NXF_OFFLINE` — under offline mode Nextflow's own file-existence validation (nf-schema) still runs and will fail on the remote paths.
- **Strict preflight**: Nextflow is never invoked if validation fails.
- **No hallucinated outputs**: Only artifacts confirmed on disk are reported.
- **Disclaimer**: Every report includes the ClawBio medical disclaimer.

## Reproducibility Scope

The bundle guarantees **configuration- and provenance-level** reproducibility, not bit-for-bit identical scientific results. Be explicit about the boundary:

**Guaranteed (captured in the bundle):**
- Pipeline pinned by `-r <version>` (or a git-verified local checkout commit).
- Nextflow engine pinned via `export NXF_VER=<version>` in `commands.sh`.
- Effective `params.yaml` + SHA-256; `--resume` rejects a mismatched params checksum.
- `checksums.sha256` (relative labels — `sha256sum -c` passes after the bundle is copied to any folder), reference-file digests in `inputs.json`, and a portable self-anchoring `commands.sh` + `remap_paths.py`. Because `commands.sh` self-anchors and `params.yaml` stores output-relative paths, **relocating the output directory needs no rewrite** — unlike the rnaseq/sarek bundles you do not pass `remap_paths.py --output-dir`; it is accepted for CLI parity and simply confirms this. Use `remap_paths.py --old/--new` and `--refs-old/--refs-new` only when the FASTQ/reference paths themselves move.

**NOT guaranteed (outside the wrapper's control):**
- **Bit-for-bit identical outputs.** Aligner thread-count can reorder records, and **CellBender is stochastic** (random-seed driven), so the same inputs+versions may not yield byte-identical `.h5ad`.
- **Container immutability.** Images are pinned by the tag nf-core ships, not by an immutable digest; a re-published tag could change the binaries.
- **`conda`/`mamba` is not offline-reproducible** — environments resolve from channels at run time (prefer `docker`/`singularity`).
- **External reference bytes are not bundled** — `--genome` pulls from iGenomes S3; local refs are checksummed for provenance but not copied into the bundle.

To approach bit-level reproduction: use `docker`/`singularity` (not conda), pin container digests, fix tool seeds where the pipeline allows, and archive the exact reference files alongside the bundle.

## Agent Boundary

The agent dispatches and explains; this skill executes.

**Agent**: Interpret the user's preprocessing intent, choose the preset, and verify that `handoff_available` is `true` in `result.json` before routing to downstream skills.

**Skill**: Validate environment and inputs, run the pipeline with controlled parameters, write all provenance and reproducibility artifacts, and report the detected `preferred_h5ad`.

## Chaining Partners

| Skill | When to chain |
|-------|--------------|
| `scrna-orchestrator` | After a successful run, pass `preferred_h5ad` for clustering, QC, and markers |
| `scrna-embedding` | Pass `preferred_h5ad` for scVI/scANVI batch integration and latent embeddings |
| `multiqc-reporter` | Re-aggregate QC across multiple wrapper runs |

## Maintenance

**Review cadence**: After each nf-core/scrnaseq major release. Check `NEXTFLOW_MIN_VERSION` (`schemas.py`), `SUPPORTED_PRESETS`, `SUPPORTED_PROFILES`, and this SKILL.md for accuracy.

**Staleness signals**:
- Preflight rejects a Nextflow version that the current pipeline supports → update `NEXTFLOW_MIN_VERSION` in `schemas.py` and `reproducibility/pinned_versions.json`.
- New aligners appear upstream but are absent from `PRESET_ALIGNERS` → add to `schemas.py` and update tests.
- New reference/index inputs appear upstream → add them to `GENOME_REFERENCE_FIELDS` or `AUXILIARY_PATH_FIELDS` in `schemas.py` (the single source both `params_builder` and `preflight` consume). Never re-introduce per-module copies.
- The pinned `STAR_ALIGN_BASE_EXT_ARGS` in `nfcore_4_1_0_contract.py` drifts from `conf/modules.config` → `test_pinned_star_args_match_sibling_checkout_if_present` fails when a sibling checkout is present; update the constant to match.
- The VirtioFS macOS workaround (`stageInMode = "copy"`) is only necessary while Apple Silicon runs Docker via QEMU. Remove `build_macos_docker_config` / `write_macos_docker_config` when a native arm64 Docker runtime eliminates VirtioFS deadlocks.

**Deprecation criteria**: Deprecate if nf-core/scrnaseq releases a Python SDK with equivalent preflight, params, and provenance APIs.

## Citations

- [nf-core/scrnaseq 4.1.0](https://nf-co.re/scrnaseq/4.1.0)
- [nf-core/scrnaseq usage](https://nf-co.re/scrnaseq/4.1.0/docs/usage/)
- [nf-core/scrnaseq parameters](https://nf-co.re/scrnaseq/4.1.0/parameters/)
- [nf-core/scrnaseq output](https://nf-co.re/scrnaseq/4.1.0/docs/output/)
- [Nextflow](https://www.nextflow.io/)
- [Alevin-fry / Simpleaf](https://simpleaf.readthedocs.io/)
- [STARsolo](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md)
- [kb-python / BUStools](https://www.kallistobus.tools/)
