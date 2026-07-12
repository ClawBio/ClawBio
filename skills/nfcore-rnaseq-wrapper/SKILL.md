---
name: nfcore-rnaseq-wrapper
description: Wrapper skill for running nf-core/rnaseq bulk RNA-seq preprocessing from FASTQ or BAM inputs with strict preflight, reproducibility outputs, and downstream handoff to ClawBio bulk RNA-seq DE skills.
license: MIT
metadata:
  version: "0.1.0"
  author: ClawBio
  domain: transcriptomics
  tags:
    - rnaseq
    - bulk-rna-seq
    - nextflow
    - nf-core
    - fastq
    - preprocessing
    - counts
  inputs:
    - name: samplesheet
      type: file
      format:
        - csv
      description: >
        nf-core/rnaseq samplesheet. Required columns: sample, fastq_1, strandedness.
        FASTQ mode may add fastq_2. BAM reprocessing mode preserves the original FASTQ
        columns and adds genome_bam and/or transcriptome_bam plus percent_mapped; use
        it only with --skip-alignment. Optional metadata columns: seq_platform, seq_center.
      required: false  # required for real runs; not for --demo or self-contained nf-core test profiles (the only universally required CLI arg is --output)
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
      description: Structured result payload with detected count matrices and provenance
  dependencies:
    python: ">=3.10"
    packages:
  demo_data:
    - path: demo/README.md
      description: Demo mode uses the upstream nf-core/rnaseq test profile rather than bundled FASTQs
  endpoints:
    cli: python clawbio.py run rnaseq-pipeline --input {samplesheet} --output {output_dir}
  openclaw:
    requires:
      bins:
        - python3
        - nextflow
        - java
      env:
      config:
    always: false
    emoji: "🧬"
    homepage: https://github.com/ClawBio/ClawBio
    os:
      - darwin
      - linux
    install:
    trigger_keywords:
      - bulk RNA-seq preprocessing
      - nf-core rnaseq
      - run rnaseq from fastq
      - preprocess RNA-seq FASTQs
      - FASTQ to count matrix
      - STAR Salmon RNA-seq pipeline
      - RSEM RNA-seq pipeline
      - HISAT2 RNA-seq alignment
      - bowtie2 salmon prokaryotic rnaseq
---

# 🧬 nfcore-rnaseq-wrapper

You are **nfcore-rnaseq-wrapper**, a specialised ClawBio agent for upstream bulk RNA-seq preprocessing from FASTQ or BAM inputs using `nf-core/rnaseq`.

## Trigger

**Fire when:**
- User wants to run `nf-core/rnaseq`
- User asks for bulk RNA-seq preprocessing from raw FASTQ files
- User wants FASTQ to gene-count matrix, Salmon counts, RSEM counts, or MultiQC outputs
- User mentions STAR/Salmon, STAR/RSEM, HISAT2, or Bowtie2/Salmon as upstream bulk RNA-seq routes
- User asks for a reproducible Nextflow wrapper before downstream differential expression

**Do NOT fire when:**
- User already has a count matrix and wants differential expression -> route to `rnaseq-de`
- User has single-cell FASTQs or wants `.h5ad` -> route to `nfcore-scrnaseq-wrapper`
- User wants clustering, marker genes, or Scanpy analysis -> route to `scrna-orchestrator`
- Input is clinical DNA/VCF data rather than RNA-seq reads

## Scope

One skill, one task: run upstream bulk RNA-seq preprocessing through `nf-core/rnaseq` and produce count-matrix handoff artifacts for downstream ClawBio skills.

This skill does not perform differential expression. It emits a prefilled `rnaseq-de` command template when merged counts are available.

## Why This Exists

- **Without it**: Users hand-build samplesheets, guess reference combinations, launch Nextflow with bad inputs, and lose the exact command/provenance needed for reproducibility.
- **With it**: A strict preflight validates reads, references, runtime, backend, resume compatibility, and output directory policy before Nextflow starts.
- **Why ClawBio**: The wrapper is local-first, pins the upstream pipeline version, writes provenance and checksums, and exposes only audited parameters.

## Core Capabilities

1. **Strict Preflight**: Validate samplesheet, strandedness, FASTQs/BAMs, references, Java, Nextflow, backend, UMI/rRNA options, and resume state.
2. **Audited Execution**: Run `nf-core/rnaseq` v3.26.0 through `-params-file` with deterministic work/result directories.
3. **Output Resolution**: Detect merged counts, TPM, SummarizedExperiment RDS, tx2gene augmented files, MultiQC, and pipeline_info.
4. **Reproducibility Bundle**: Write `commands.sh`, `params.yaml`, `manifest.json`, checksums, `environment.yml`, and seven provenance JSON files.
5. **Downstream Handoff**: Emit a template for `python clawbio.py run rnaseq --counts ...` when a merged count matrix is available.

## Aligners

| `--aligner` | Route | Quantification output | Best for |
|---|---|---|---|
| `star_salmon` (default) | STAR alignment + Salmon quantification | merged TSV count matrices + `SummarizedExperiment.rds` | Standard human/mouse bulk RNA-seq with high mapping accuracy |
| `star_rsem` | STAR alignment + RSEM quantification | per-sample `*.genes.results` + merged matrix + RDS | Encode-style isoform-level analyses |
| `hisat2` | HISAT2 alignment only (no quantification) | BAM only — `handoff_available=false` unless `--pseudo-aligner` is also set | Alignment-only workflows; add `--pseudo-aligner salmon` to re-enable downstream DE handoff |
| `bowtie2_salmon` | Bowtie2 alignment + Salmon quantification | merged TSV count matrices + RDS | Prokaryotic transcriptomes (combine with `--prokaryotic`) |

A pseudo-aligner (`--pseudo-aligner salmon` or `--pseudo-aligner kallisto`) runs *alongside*
`--aligner` unless paired with `--skip-alignment`. Each route may use either `--genome <iGenomes>`
(optionally with additive annotation/transcriptome overrides such as `--gtf` *or* `--gff`,
`--additional-fasta`, `--transcript-fasta`, `--gene-bed`, `--splicesites`, `--salmon-index`, or
`--kallisto-index`) *or* a fully explicit `--fasta`/`--gtf`(/`--gff`) reference plus optional
pre-built `--*-index` paths. You may not provide both `--genome` and your own genome `--fasta`
or a genome-level index (`--star-index`/`--rsem-index`/`--hisat2-index`/`--bowtie2-index`).
If both `--gtf` and `--gff` are supplied, the wrapper keeps `--gtf` and drops `--gff` with a
warning — matching nf-core/rnaseq, which uses the GTF and ignores the GFF when both are given.
For new analyses nf-core/rnaseq recommends supplying explicit `--fasta`/`--gtf` directly; the
iGenomes `--genome` catalogue is supported here for legacy compatibility and convenience.

## Input Formats

| Format | Extension | Required Fields | Example |
|--------|-----------|-----------------|---------|
| Samplesheet | `.csv` | `sample`, `fastq_1`, `strandedness`; optional `fastq_2` | `samplesheet.csv` |
| BAM reprocessing samplesheet | `.csv` | `sample`, `fastq_1`, `strandedness`, plus `genome_bam` and/or `transcriptome_bam`; use with `--skip-alignment` | `samplesheet_with_bams.csv` |
| Demo mode | n/a | none | `python clawbio.py run rnaseq-pipeline --demo` |

## Workflow

1. **Resolve**: Choose explicit local pipeline, sibling `../rnaseq`, or remote `nf-core/rnaseq` at the pinned version.
2. **Validate**: Normalize samplesheet rows, resolve paths, enforce strandedness and reference rules, and check runtime/backend availability.
3. **Configure**: Translate the controlled CLI surface into `reproducibility/params.yaml`.
4. **Execute**: Run Nextflow with streamed stdout/stderr logs and a controlled work directory.
5. **Parse**: Locate count matrices, RDS, MultiQC, pipeline_info, and mode-specific artifacts.
6. **Report**: Write `report.md`, `result.json`, provenance JSON, checksums, and replay commands.
7. **Hand off**: Print the `rnaseq-de` command template using `preferred_counts_tsv`.

## CLI Reference

```bash
# Preflight only; no Nextflow execution
python clawbio.py run rnaseq-pipeline \
  --input samplesheet.csv --output ./rnaseq_check --check \
  --genome GRCh38

# Demo mode using upstream test profile
python clawbio.py run rnaseq-pipeline --demo --output ./rnaseq_demo

# STAR + Salmon default route
python clawbio.py run rnaseq-pipeline \
  --input samplesheet.csv --output ./rnaseq_run \
  --aligner star_salmon --genome GRCh38

# Explicit FASTA/GTF reference
python clawbio.py run rnaseq-pipeline \
  --input samplesheet.csv --output ./rnaseq_run \
  --fasta /refs/genome.fa --gtf /refs/genes.gtf

# RSEM route
python clawbio.py run rnaseq-pipeline \
  --input samplesheet.csv --output ./rsem_run \
  --aligner star_rsem --genome GRCh38

# Contaminant screening with Kraken2 + Bracken
python clawbio.py run rnaseq-pipeline \
  --input samplesheet.csv --output ./rnaseq_run \
  --genome GRCh38 \
  --contaminant-screening kraken2_bracken \
  --kraken-db /refs/kraken2_db --bracken-precision G

# Auto-handoff to rnaseq-de when all flags are provided
python clawbio.py run rnaseq-pipeline \
  --input samplesheet.csv --output ./rnaseq_run \
  --genome GRCh38 --run-downstream \
  --metadata metadata.csv --formula "~ batch + condition" \
  --contrast "condition,treated,control"

# Prokaryotic transcriptomes via Bowtie2+Salmon
python clawbio.py run rnaseq-pipeline \
  --input samplesheet.csv --output ./prok_run \
  --aligner bowtie2_salmon --fasta /refs/genome.fa --gtf /refs/genes.gtf \
  --profile docker --prokaryotic

# ARM architecture (Apple M-series, AWS Graviton) — composes -profile docker,arm64
python clawbio.py run rnaseq-pipeline \
  --input samplesheet.csv --output ./rnaseq_arm \
  --genome GRCh38 --profile docker --arm

# BAM reprocessing from nf-core samplesheet_with_bams.csv output
python clawbio.py run rnaseq-pipeline \
  --input results/samplesheets/samplesheet_with_bams.csv \
  --output ./rnaseq_reprocess \
  --skip-alignment

# Wrapper runtime controls (parity with scrnaseq/sarek):
#   --timeout-hours N   wall-clock cap (default 12h; 0 disables for HPC/cloud)
#   --work-dir PATH     Nextflow work dir (local path or object-store URI; default <output>/upstream/work)
#   --nextflow-config / -c / --config   extra Nextflow config file(s), repeatable
#   --allow-pipeline-version-override    run a non-3.26.0 --pipeline-version at your own risk
#   --allow-remote-inputs               opt in to remote inputs/refs (default local-first)
python clawbio.py run rnaseq-pipeline \
  --input samplesheet.csv --output ./rnaseq_run \
  --genome GRCh38 --aligner star_salmon \
  --timeout-hours 0 --work-dir s3://my-bucket/rnaseq/work
```

## Demo

```bash
python clawbio.py run rnaseq-pipeline --demo --output /tmp/rnaseq_demo
```

Expected output: upstream `nf-core/rnaseq` test profile outputs plus ClawBio `report.md`, `result.json`, `provenance/`, and `reproducibility/`.

## Algorithm / Methodology

The wrapper uses a gated 7-step flow. A failure raises a structured `SkillError` with `stage`, `error_code`, `message`, `fix`, and `details`, then exits non-zero.

Key methods:
- Local samplesheet paths are resolved against the samplesheet directory and written as absolute POSIX paths; remote URIs (`s3://`, `https://`, ... — only accepted with `--allow-remote-inputs`) are passed through unchanged.
- `params.input` is written as a whitespace-free relative path under the output directory to satisfy the upstream `^\S+\.csv$` schema.
- References must use either `--genome`, `--fasta --gtf`, or `--fasta --gff`.
- `--genome` accepts additive annotation/transcriptome overrides (`--gtf` *or* `--gff`, `--gene-bed`, `--transcript-fasta`, `--additional-fasta`, `--splicesites`, `--salmon-index`, `--kallisto-index`) — matching nf-core/rnaseq — but is mutually exclusive with a genome `--fasta` or a genome-level index (`--star-index`/`--rsem-index`/`--hisat2-index`/`--bowtie2-index`).
- `--gtf` and `--gff` together are not rejected: nf-core/rnaseq uses the GTF and ignores the GFF when both are given, so the wrapper drops `--gff` (with a warning) and proceeds with `--gtf`, matching upstream in every reference mode (`--genome`, explicit `--fasta`, prebuilt indices).
- HISAT2 alignment-only mode sets `handoff_available=false`.
- Per-sample quantification mode does not auto-chain to `rnaseq-de`.

## Example Queries

- "Run nf-core/rnaseq on these FASTQs"
- "Preprocess bulk RNA-seq FASTQ files into a count matrix"
- "Run STAR Salmon and prepare counts for DESeq2"
- "Check my RNA-seq samplesheet before running Nextflow"

## Example Output

```markdown
# nf-core/rnaseq Wrapper Report

## Summary
- Aligner: `star_salmon`
- Samples: `5`

## Outputs
- Preferred counts TSV: `/run/upstream/results/star_salmon/salmon.merged.gene_counts_length_scaled.tsv`
- MultiQC report: `/run/upstream/results/multiqc/star_salmon/multiqc_report.html`

## Next Steps
python clawbio.py run rnaseq --counts <preferred_counts_tsv> --metadata <your_metadata.csv> ...
```

## Output Structure

```
output/
├── report.md
├── result.json
├── logs/
├── upstream/
│   ├── results/
│   │   ├── samplesheets/
│   │   │   └── samplesheet_with_bams.csv   # only when --save-align-intermeds; use with --skip-alignment for BAM reprocessing
│   │   ├── star_salmon/                    # star_salmon aligner outputs
│   │   │   ├── *.markdup.sorted.bam        # sorted, deduplicated BAMs (one per sample)
│   │   │   ├── log/                        # STAR alignment logs (*.Log.final.out, *.SJ.out.tab)
│   │   │   ├── salmon.merged.*.tsv         # merged gene/transcript count matrices
│   │   │   └── salmon.merged.*.rds         # SummarizedExperiment objects
│   │   └── ...
│   └── work/
├── provenance/
└── reproducibility/
    ├── samplesheet.valid.csv   # demo run → samplesheet.demo.csv; test profile → samplesheet.noinput.csv
    ├── params.yaml
    ├── commands.sh
    ├── remap_paths.py
    ├── manifest.json
    ├── environment.yml
    └── checksums.sha256
```

## Dependencies

**Required**
- Python >=3.10
- Java >=17
- Nextflow >=25.04.3
- One execution backend: Docker, Singularity, Apptainer, Podman, Conda/Mamba, Shifter, or Charliecloud

## Gotchas

- `strandedness` is required per row and must be `auto`, `forward`, `reverse`, or `unstranded`.
- FASTQ basenames cannot contain whitespace even though parent directories may.
- FASTQ basenames must end in `.fq`, `.fastq`, `.fq.gz`, or `.fastq.gz` (all four are accepted by the nf-core/rnaseq schema). Only the basename must be whitespace-free; parent directory paths may contain spaces.
- FASTQ and BAM samplesheet entries may be local paths or remote URIs such as `s3://...`/`https://...`. Local paths are normalized and existence-checked; remote URIs are preserved unchanged and left for Nextflow to stage.
- `--genome` may be combined with additive annotation/transcriptome overrides (`--gtf` *or* `--gff`, `--gene-bed`, `--transcript-fasta`, `--additional-fasta`, `--splicesites`, `--salmon-index`, `--kallisto-index`) — this matches nf-core/rnaseq and supports common cases such as ERCC spike-ins (`--genome GRCh38 --additional-fasta ercc.fa`) or overriding the dated iGenomes annotation (`--genome GRCh38 --gtf custom.gtf`). It is rejected only with a second genome **sequence** source (`--fasta`) or a genome-level index (`--star-index`/`--rsem-index`/`--hisat2-index`/`--bowtie2-index`), which would be ambiguous. If both `--gtf` and `--gff` are supplied, `--gff` is dropped with a warning and `--gtf` is used (matching nf-core/rnaseq). Names not in the built-in iGenomes catalogue emit a preflight warning but do not block execution — this is expected when using a user-defined genome catalogue (pass it via `--nextflow-config my_genomes.config`). If you intended an iGenomes entry, check the exact spelling and case (e.g. `GRCh38`, `GRCm38`).
- GENCODE autodetection (setting `gencode: true` from `gene_type`/`havana_gene` markers in the GTF) only inspects **local** `--gtf` files; for remote (`s3://`/`https://`) GTFs it is skipped silently — pass `--gencode` explicitly in that case. Autodetection scans only the **first 10 feature records** of the GTF (gzip is detected case-insensitively, e.g. `.gtf.gz` and `.gtf.GZ`); if your GENCODE markers appear later in the file, pass `--gencode` explicitly.
- `--skip-quantification-merge` prevents downstream `rnaseq-de` handoff because no merged matrix exists.
- `--aligner hisat2` is alignment-only for this handoff contract.
- `--with-umi` requires a barcode pattern unless `--skip-umi-extract` is set. Conversely, UMI options (`--umitools-bc-pattern`, `--umi-dedup-tool`, etc.) set *without* `--with-umi` are inert — preflight warns so a run is not mistaken for UMI-deduplicated when it is not.
- `--output` must be **outside** the ClawBio source tree. An output directory inside the repository is rejected at preflight with `OUTPUT_DIR_INSIDE_REPO`, so multi-gigabyte pipeline artifacts never pollute (or get committed to) the checkout — choose a path under your analysis workspace. This matches the nfcore-sarek and nfcore-scrnaseq wrappers.
- On macOS Docker, use an output directory under the home directory rather than `/tmp`. The wrapper writes a macOS Docker compatibility config whose per-process memory ceiling is derived from host RAM (75% share, floored at 8 GB, capped at 15 GB) and then capped to 90% of the actual Docker VM memory (`docker info`, when available) so a container process is never OOM-killed by requesting more than the VM has. Its per-process `time` ceiling tracks `--timeout-hours` (default 12, floored at 1 h) so raising the wrapper timeout does not leave processes capped at 12 h.
- The local Nextflow run is killed after `--timeout-hours` (default 12). Raise it for large cohorts (e.g. `--timeout-hours 48`) so a long but healthy run is not terminated, or pass `--timeout-hours 0` to disable the cap entirely for long HPC/cloud runs whose walltime is enforced by the scheduler (negative values are rejected). On a timeout the wrapper terminates Nextflow's process group, but containers started by the Docker/Singularity daemon are not in that group and may keep running — the timeout error reminds you to check for and remove leftover containers (e.g. `docker ps`).
- Reference paths (`--fasta`/`--gtf`/`--gff`/`--transcript-fasta`/`--additional-fasta`/`--gene-bed`) must resolve to a path **without whitespace** — the nf-core schema pattern `^\S+` rejects spaces. Preflight catches a whitespace-containing resolved path early with a precise `REFERENCE_PATH_HAS_WHITESPACE` error (mirroring the samplesheet input guard) instead of letting Nextflow abort late. Move or symlink the reference into a space-free directory.
- `--check` validates that Nextflow is present but defers the `>=25.04.3` version gate to the real run; it emits a warning so a passing check is not mistaken for confirmation of a compatible Nextflow version.
- Results are written under a **relative** `upstream/results` because the wrapper launches Nextflow with `cwd=<output>`; the relative path keeps the nf-core `^\S+$` `outdir` schema valid even when `--output` contains spaces (common on macOS). This is a deliberate **local-first** design. Running against cloud executors that require an absolute publish path (e.g. `outdir` on `s3://`/`gs://`) is outside the wrapper's audited surface.
- The wrapper exposes the **audited scientific parameter surface** of nf-core/rnaseq 3.26.0. A few cosmetic/notification options (`--plaintext_email`, `--max_multiqc_email_size`, `--monochrome_logs`, `--trace_report_suffix`, `--custom_config_*`) are intentionally not exposed. Non-parametric runtime settings (executor, resource limits, institutional config) are supplied through `--nextflow-config`.
- A sibling `../rnaseq` checkout is auto-detected and used, but its `manifest.version` must be `3.26.0` (the version this wrapper's validations are pinned to). A different version is rejected unless `--allow-pipeline-version-override` is passed; an unparseable manifest version is warned, not blocked.
- `--rseqc-modules` is validated against the eight nf-core/rnaseq 3.26.0 module names; a typo is rejected at preflight instead of failing later inside Nextflow.
- `--contaminant-screening kraken2`/`kraken2_bracken` requires `--kraken-db`, and `--contaminant-screening sylph` requires `--sylph-db`; local database paths are existence-checked before Nextflow starts, while URI schemes such as `s3://` and `https://` are passed through for Nextflow to stage. `--bracken-precision` only applies to `kraken2_bracken` and is warned (no effect) otherwise.
- Transcriptome-only pseudo-quantification (`--skip-alignment` + `--pseudo-aligner salmon`/`kallisto` + `--transcript-fasta` or a prebuilt `--salmon-index`/`--kallisto-index` + `--gtf`/`--gff`) is accepted without a genome `--fasta`. A pseudo-aligner running *alongside* a genome aligner still requires the genome reference.
- Fully prebuilt references need no `--fasta`: a genome index matching the aligner (`--star-index`/`--hisat2-index`/`--bowtie2-index`, or `--rsem-index` for `star_rsem`) plus `--gtf`/`--gff` and, for the Salmon routes, a transcript source (`--transcript-fasta` or `--salmon-index`) is accepted. A bare genome index without a transcript source (Salmon routes) or without `--rsem-index`/`--fasta` (RSEM) is rejected because quantification cannot run.
- `--pseudo-aligner-kmer-size` must be an **odd integer in 1..31** (Salmon and Kallisto both encode the index k-mer in a 64-bit word, so 31 is their shared hard cap; pipeline default 31). Preflight rejects an even or out-of-range value with `INVALID_PRESET_CONFIGURATION` instead of letting the pseudo-aligner indexing step crash. Lower it for short reads (<50 bp).
- Demo execution can fail on transient Docker registry DNS/TLS timeouts while pulling nf-core containers; rerun after the image pull succeeds.
- `--prokaryotic`, `--rapid-quant`, and `--arm` are profile-modifier flags. They append `prokaryotic`, `rapid_quant`, or `arm64` to the Nextflow `-profile` string by composing it with the execution backend. Use `--profile docker --prokaryotic` (composes `-profile docker,prokaryotic`). `--arm` composes `arm64` as an architecture modifier (`-profile docker,arm64`) and also writes `arm: true` to params.yaml — `arm` is a real hidden boolean parameter in the nf-core/rnaseq 3.26.0 schema (`"Use ARM architecture containers."`).
- BAM reprocessing samplesheets must preserve the official FASTQ columns: `sample`, `fastq_1`, `strandedness`, plus at least one of `genome_bam` or `transcriptome_bam`. Use the nf-core-generated `samplesheet_with_bams.csv` with `--skip-alignment`. Rows with BAMs and an empty `fastq_1` are rejected because they no longer match the audited nf-core/rnaseq 3.26.0 samplesheet contract. **Reprocess with the same `--aligner` used to generate the BAMs**: nf-core/rnaseq cannot mix quantifier types between BAM generation and reprocessing (BAMs from `star_salmon` must be reprocessed with `star_salmon`, `star_rsem` with `star_rsem`). The wrapper defaults to `star_salmon`, so pass `--aligner star_rsem` explicitly when reprocessing RSEM BAMs; preflight emits a reminder warning whenever BAM reprocessing is detected. The `samplesheet_with_bams.csv` you reprocess from is **only produced when the original alignment run used `--save-align-intermeds`** — nf-core/rnaseq creates it solely in that case, so add `--save-align-intermeds` to the run whose BAMs you intend to reprocess later.
- `--ribo-database-manifest` is preflight-checked when it is a local path; missing files or directories are rejected before Nextflow starts. URI schemes are preserved unchanged in `params.yaml`.
- `--use-parabricks-star` requires `--aligner star_salmon`; `--use-sentieon-star` requires a STAR-based aligner (`star_salmon` or `star_rsem`); `--use-gpu-ribodetector` requires `--remove-ribo-rna --ribo-removal-tool ribodetector`.
- Downstream `rnaseq-de` handoff is opt-in via `--run-downstream`. It **launches** `rnaseq-de` only when `--run-downstream` is set *and* `--metadata`, `--formula`, and `--contrast` are all provided. With `--run-downstream` but any of those three missing, only a copy-paste template `reproducibility/rnaseq_de_handoff.sh` is written. **Without `--run-downstream` (the default, including `--demo`), no handoff is launched and no template file is written** — the `report.md` "Next Steps" section still shows the suggested `rnaseq-de` command. `--skip-downstream` suppresses the template even when `--run-downstream` is set.
- `--rseqc-modules` runs a default set of 7 modules. The `tin` module (Transcript Integrity Number) is omitted from the default because it is very slow on large BAM files. Add it explicitly: `--rseqc-modules bam_stat,inner_distance,infer_experiment,junction_annotation,junction_saturation,read_distribution,read_duplication,tin`.
- `--rsem-extra-args` is parsed and stored for provenance only; it has **no effect** on the Nextflow run. nf-core/rnaseq ≥3.14 removed `extra_rsem_quant_args` from the schema. Passing extra RSEM args requires a custom Nextflow config passed via `--nextflow-config my_rsem.config`.
- `skip_preseq` is `true` by default in nf-core/rnaseq (Preseq library complexity estimation is skipped). Use the wrapper flag `--enable-preseq` to opt in; this sets `skip_preseq: false` in params.yaml. Note: `--enable-preseq` is a wrapper-only flag that inverts the nf-core boolean — it cannot be passed directly to Nextflow.
- `--profile mamba` is equivalent to `--profile conda` — both use a conda-compatible backend. The wrapper accepts either spelling.
- `--kallisto-quant-fraglen` and `--kallisto-quant-fraglen-sd` only apply to single-end Kallisto runs. Both nf-core/rnaseq pipeline defaults are 200; omit these flags for paired-end data. Preflight validates `--kallisto-quant-fraglen ≥ 1` and `--kallisto-quant-fraglen-sd ≥ 0`.
- `--min-trimmed-reads` must be ≥ 0 (pipeline default: 10000). Preflight rejects negative values. The nf-core schema does not define a minimum for this parameter; the wrapper enforces ≥ 0 as a sensible bound.
- **Omit = trust upstream default.** Several string parameters are intentionally absent from `params.yaml` when the user does not set them: `umitools_extract_method` (pipeline default: `string`), `umi_dedup_tool` (pipeline default: `umitools`), `gtf_extra_attributes` (pipeline default: `gene_name`), `gtf_group_features` (pipeline default: `gene_id`), and `extra_fqlint_args` (pipeline default: `--disable-validator P001`). Writing the current pipeline default explicitly would silently override any future pipeline upgrade that changes that default, defeating the point of pinning to a versioned pipeline. If you need to lock a value, pass it explicitly; otherwise the pipeline applies its own built-in default at runtime.
- Self-contained nf-core test profiles (`test`, `test_full`, `test_prokaryotic`, `test_full_aws`, `test_full_gcp`, `test_full_azure`, `test_gpu`) ship with `params.input` in their profile config and do not require `--input`. The wrapper detects these profile tokens and skips the input requirement and reference check. `test_full*` profiles use `genome='GRCh37'` via iGenomes — the wrapper does **not** set `igenomes_ignore: true` (nor `aligner`, unless you pass `--aligner` explicitly) for these, letting the profile config own them. `--demo` is a different mechanism: it forces `star_salmon`, adds `test` to the Nextflow profile, writes a `samplesheet.demo.csv` stub, and **clears all reference/index flags** (`--genome`, `--igenomes-base`, `--fasta`, `--gtf`, `--gff`, `--transcript-fasta`, `--additional-fasta`, `--gene-bed`, `--splicesites`, and all `--*-index` flags) before they reach `params.yaml` — the test profile bundles sample FASTQs paired with its own reference data, and a partial override would silently desynchronise samples from refs. Self-contained test profile runs produce `samplesheet.noinput.csv` instead so provenance audits can distinguish them. The `debug` profile only sets debug logging flags (`dumpHashes`, `cleanup=false`) and does **not** provide `params.input` — it still requires `--input`.
- **`--demo` requires network access.** It runs the upstream nf-core `-profile test`, whose sample FASTQs and reference FASTA/GTF are fetched from remote GitHub URLs (nf-core's design — the wrapper does not bundle local test data). On an offline/sandboxed host set `NXF_OFFLINE`, and the wrapper fails fast at preflight with `DEMO_REQUIRES_NETWORK` and a clear message, instead of a cryptic Nextflow `does not exist` abort during schema validation. This does **not** violate the local-first guarantee, which governs your *genetic data* (never uploaded); `--demo` only *downloads* nf-core's public test data. For a fully offline run, use a real analysis with your own local `--input` samplesheet and references.
- **nf-core-native (snake_case) flag spellings are accepted via the launcher.** You can paste an upstream nf-core command's parameters verbatim (`--gene_bed`, `--transcript_fasta`, …): `clawbio.py run rnaseq-pipeline` treats `_` and `-` as equivalent when matching the flag allowlist and forwards the wrapper's hyphenated spelling. No manual underscore-to-hyphen conversion is needed.
- **Host-limited memory is auto-capped on docker runs; IPv6-only networks are an environment issue — read the failure hint.** On a docker backend the wrapper writes a `process.resourceLimits` config scaled to this host (on macOS: host-RAM share capped to the Docker VM; on Linux/other: physical RAM minus headroom) so a real run does not abort with `Process requirement exceeds available memory` when an nf-core default request — e.g. `MAKE_TRANSCRIPTS_FASTA` — is larger than your machine (`--demo` is exempt: `-profile test` carries its own limits). If it still aborts (a non-docker backend, or one process that genuinely needs more RAM than the host has), override with your own `-c` config, e.g. `process { resourceLimits = [ memory: '12.GB', cpus: 4 ] }`; do **not** delete resource labels to force it through. On an IPv6-only / NAT64 host the JVM prefers IPv4 and downloads fail with `Network is unreachable`; export `NXF_OPTS='-Djava.net.preferIPv6Addresses=true'` and re-run. The wrapper inherits your environment and never overrides `NXF_OPTS`.
- **Replaying a bundle in place is idempotent — including a `--demo` bundle. Do not "fix" a replay by deleting the output directory.** Unlike the sarek/scrnaseq bundles (which replay Nextflow directly, and Nextflow tolerates a populated output dir), the rnaseq `commands.sh` re-invokes the *wrapper*, whose preflight rejects a non-empty `--output` with `OUTPUT_DIR_NOT_EMPTY`. So `commands.sh` carries a guard that adds `--resume` when the target output dir already holds a completed run of this bundle (`reproducibility/manifest.json` present); a fresh or `remap_paths.py --output-dir`-relocated directory has no manifest and runs clean. `--demo` bundles get the same guard: Nextflow's `-resume` is orthogonal to `-profile test` (nf-core documents no incompatibility), the demo samplesheet stub is content-stable so its checksum matches on replay, and the run's work tree (`upstream/work`) and Nextflow session cache (`.nextflow/`) both live under the output dir. Resuming across the demo/real boundary is still blocked — `demo` is compared against the manifest like `aligner`/`profile`/`arm`.
- **Relocating the bundle to a new output directory: use `remap_paths.py --output-dir <new-path>`.** The rnaseq bundle bakes the `--output` directory into `commands.sh` (its replay re-invokes the wrapper), so after moving the output tree run `python3 reproducibility/remap_paths.py --output-dir <new-path>` to rewrite it (it keeps the replay guard's manifest path in sync). Use `--old/--new` for relocated FASTQs and `--refs-old/--refs-new` for relocated references in `params.yaml`. The sarek bundle exposes the same `--output-dir`; the scrnaseq bundle self-relocates (its `commands.sh` self-anchors) and accepts `--output-dir` only for parity, as a no-op that confirms no rewrite is needed.

## Safety

- No patient data is bundled.
- Demo mode uses upstream test profile data.
- The wrapper does not upload data.
- **Local-first by default**: remote samplesheet inputs and reference paths are rejected (`REMOTE_INPUT_NOT_ALLOWED`) unless `--allow-remote-inputs` is explicitly passed, which also logs a runtime warning naming every path fetched over the network. The object-store `--work-dir` is not gated. `--allow-remote-inputs` relaxes only the wrapper's own preflight check: remote FASTQ/reference URIs are then written into the normalized samplesheet/`params.yaml` **verbatim** and staged natively by Nextflow at run time. The wrapper does not download them itself, so remote inputs require outbound network access and are incompatible with `NXF_OFFLINE` — under offline mode Nextflow's own file-existence validation (nf-schema) still runs and will fail on the remote paths.
- The wrapper does not pass arbitrary unvalidated Nextflow *parameters* via `--params-file`: only the audited CLI surface is translated to `params.yaml`. `--nextflow-config` forwards user-supplied `-c` config file(s) for trusted runtime settings such as process, executor, profiles, labels, institutional module tuning, and `params.genomes` custom genome catalogues. Configs that define `params` in any form — block (`params { … }`), property (`params.x`), assignment (`params = …`), subscript (`params['x']`), or map-merge (`params << …`) — are rejected so they cannot bypass the audited parameter surface (the documented `params.genomes` catalogue is the sole exception). Every locally-resolvable `includeConfig` target is audited recursively under the same rule; includes the wrapper cannot read (remote URIs, `${…}`-interpolated paths, or missing files) are surfaced as preflight **warnings** rather than silently trusted, so unaudited surface is always visible.
- `--resume` is rejected when the pipeline source/version, profile, aligner, pseudo-aligner, `--demo`/`--prokaryotic`/`--arm` modifiers, params checksum, or samplesheet checksum drift. `--demo` is part of that contract because it composes the upstream `test` profile, which supplies both its own samplesheet and its own bundled references — resuming across the demo/real boundary would swap both underneath the run.

> ClawBio is a research and educational tool. It is not a medical device and does not provide
> clinical diagnoses. Consult a healthcare professional before making any medical decisions.

## Agent Boundary

Use this skill to produce upstream bulk RNA-seq preprocessing outputs. Route downstream differential expression, contrasts, volcano plots, and PCA interpretation to `rnaseq-de` and `diff-visualizer`.

## Chaining Partners

- `rnaseq-de`: bulk/pseudo-bulk differential expression from `preferred_counts_tsv`
- `diff-visualizer`: plots from downstream DE results
- `multiqc-reporter`: optional QC aggregation/reporting follow-up
- `bio-orchestrator`: routes inbound bulk RNA-seq preprocessing requests to this wrapper

## Maintenance

Pinned upstream: `nf-core/rnaseq` v3.26.0. Before changing the default version, audit `nextflow.config`, `assets/schema_input.json`, `nextflow_schema.json`, `docs/output.md`, and changed module configs, then update tests and `reproducibility/pinned_versions.json`.

## Citations

- [nf-core/rnaseq 3.26.0](https://nf-co.re/rnaseq/3.26.0)
- [nf-core/rnaseq usage](https://nf-co.re/rnaseq/3.26.0/docs/usage/)
- [nf-core/rnaseq parameters](https://nf-co.re/rnaseq/3.26.0/parameters/)
- [nf-core/rnaseq output](https://nf-co.re/rnaseq/3.26.0/docs/output/)
- [Nextflow](https://www.nextflow.io/)
- [STAR](https://github.com/alexdobin/STAR)
- [Salmon](https://salmon.readthedocs.io/)
- [RSEM](https://github.com/deweylab/RSEM)
- [HISAT2](https://daehwankimlab.github.io/hisat2/)
- [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
