# nfcore-scrnaseq-wrapper

ClawBio wrapper for running the upstream `scrnaseq` Nextflow pipeline from FASTQ inputs with strict preflight, reproducibility artifacts, provenance, and explicit handoff to downstream ClawBio scRNA skills.

## Scope

- Upstream scRNA preprocessing from FASTQ via Nextflow.
- Curated presets for `standard`/simpleaf, `star`, `kallisto`, `cellranger`, `cellrangerarc`, and `cellrangermulti`.
- Validated samplesheet input with absolute normalized FASTQ paths.
- Reproducibility bundle, provenance bundle, `report.md`, and `result.json`.
- Canonical `.h5ad` detection for downstream `scrna` and `scrna-embedding`.
- Minimum nf-core output validation before reporting a successful run.

## Out Of Scope

- Clustering, marker detection, normalization, scVI/scANVI, or other downstream analysis.
- Downstream chaining is **off by default**; it runs only as an explicit opt-in via `--run-downstream` (and `--skip-downstream` force-skips it). The wrapper never chains automatically without that flag.
- Free-form Nextflow passthrough flags beyond validated `-c`/`--config`/`--nextflow-config` files.
- Running without preflight.

## Quick Start

```bash
python clawbio.py run scrnaseq-pipeline --demo --output "$HOME/clawbio-runs/scrnaseq_demo"
```

For real data:

```bash
python clawbio.py run scrnaseq-pipeline \
  --input samplesheet.csv \
  --output "$HOME/clawbio-runs/scrnaseq_run" \
  --preset star \
  --protocol 10XV3 \
  --genome GRCh38
```

Use `--check` to run preflight without launching Nextflow. Checkout-mode runs without `--output` default to `~/clawbio-runs/`; explicit destinations must be outside the repository.

## nf-core/scrnaseq 4.1.0 Compatibility Policy

This wrapper targets nf-core/scrnaseq `4.1.0`. It is not a free-form passthrough. Parameters are grouped as:

- **Supported upstream parameters:** input/output, aligner/preset, reference/index, skip, CellRanger, CellRanger ARC, CellRanger Multi, selected MultiQC/reporting options.
- **Wrapper policy parameters:** `--preset`, `--check`, `--run-downstream`, `--skip-downstream`, `--expected-cells`, `--timeout-hours`, `--work-dir`, `--allow-remote-inputs`, `--allow-dirty-pipeline`, `--require-local-pipeline`, `--allow-pipeline-version-override`, `--trust-config-params`, `--allow-conda-cellranger`, and `-c`/`--config`/`--nextflow-config`; these are ClawBio conveniences and are not nf-core parameters.
- **Deprecated compatibility aliases:** `skip_emptydrops` is accepted only as `--skip-emptydrops` and translated to `skip_cellbender: true`; the deprecated upstream parameter is never written.
- **Intentionally unsupported upstream parameters:** `custom_config_version`, `custom_config_base`, `config_profile_name`, `config_profile_description`, `config_profile_contact`, `config_profile_url`, `version`, `plaintext_email`, `max_multiqc_email_size`, `hook_url`, `validate_params`, `pipelines_testdata_base_path`, `help`, `help_full`, `show_hidden`.

Unsupported parameters are either hidden/institutional metadata, interactive help/version flags, or options that would weaken the wrapper's fixed validation/reproducibility policy.

## Input & Reference Path Policy

**Local-first by default.** Samplesheet FASTQs and reference/index inputs must be local paths; remote URIs (`s3://`, `gs://`, `https://`, `ftp://`, …) are rejected at preflight (`REMOTE_INPUT_NOT_ALLOWED`) so genetic data and references stay on the local machine. Pass `--allow-remote-inputs` to opt in — remote URIs are then passed through verbatim (Nextflow stages them; only the basename is validated) and a runtime warning lists every path fetched over the network. The same flag is shared by all three nf-core wrappers. Local references and FASTQs are always existence-checked at preflight so missing local paths fail fast.

## Protocol Policy

For `standard`, `star`, and `kallisto`, `--protocol` is required because nf-core/scrnaseq 4.1.0 documents `auto` as CellRanger-only. Explicit `auto` is rejected for those presets. `smartseq` is accepted for `star` and `kallisto` only. `cellranger` accepts only `auto` and `10XV1`-`10XV4`; `cellrangerarc` accepts only `auto`; `cellrangermulti` is samplesheet-driven. Unknown custom protocol strings are passed through only for `standard`, `star`, and `kallisto`, where nf-core documents custom values.

## Output Policy

After Nextflow exits, the wrapper validates `pipeline_info/`, the effective aligner directory, MultiQC unless `--skip-multiqc` was used, FastQC unless `--skip-fastqc` was used, and at least one `.h5ad`. FastQC HTML/ZIP reports are a hard requirement for **every** aligner — including the Cell Ranger family — because nf-core/scrnaseq 4.1.0 runs FASTQC on the shared input-read channel before aligner branching and publishes `fastqc/` for all of them. Missing required outputs raise `EXPECTED_OUTPUTS_NOT_FOUND`. Present-but-ambiguous `.h5ad` files do not fail the run, but `handoff_available=false` is reported and `--run-downstream` prints a warning. `result.json` exposes an `official_outputs` manifest for documented nf-core output families, including FastQC, reference genome, aligner outputs, CellBender, and matrix conversions.

## Production Safeguards

- Cell Ranger presets fail fast with `conda`/`mamba` unless `--allow-conda-cellranger` is supplied for a trusted site configuration.
- `--work-dir` supports local paths and object-store URIs for cloud executors; the default remains `<output>/upstream/work`.
- `--resume` requires the same preset, profile, pipeline source, effective params checksum, and work directory recorded in the prior manifest.
- `--require-local-pipeline` fails instead of falling back to the remote nf-core pipeline when a local checkout is required.
- RNA velocity pairings are validated before execution: STARsolo `Gene Velocyto` requires `--star-ignore-sjdbgtf`, and Kallisto `lamanno`/`nac` requires both capture files.
- `--check` reports the official nf-core parameter surface, supported upstream parameters, and intentionally unsupported parameters.

## Validation and resume

nf-schema checks every parameter against the pinned `nextflow_schema.json` when Nextflow starts; the wrapper refuses only NaN and infinity, which nf-schema cannot represent, and a failed run names each value nf-schema refused. Malformed or ambiguous samplesheets fail before normalized snapshots are written. An atomic launch manifest supports resuming failed launches; resumed `--check` also validates params and samplesheet checksums.

Local Docker runs, including demos, cap both per-process requests and the local executor budget to detected host/Docker CPU and memory capacity. Site/cloud profiles and object-store work directories keep their own resource policy. Explicit `--nextflow-config` files have final precedence. Resource caps cannot make a dataset fit insufficient RAM.

Replay uses the smaller of recorded limits and current Docker capacity, on macOS and Linux, through a temporary config; the original checksummed config is preserved. User config snapshots retain final precedence, except for container images: a replay without a Nextflow cache applies `containers.config` last when a complete, locally verified container record is available. Registry-only observations do not establish the original image and withhold the pin config.

See the [compatibility audit](../../docs/nfcore-wrapper-audit.md) for upstream contracts, fixes, test results, and execution limits.
