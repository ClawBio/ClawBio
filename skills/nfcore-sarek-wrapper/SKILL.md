---
name: nfcore-sarek-wrapper
description: ClawBio wrapper around nf-core/sarek 3.8.1 covering mapping through annotation for germline, tumor-only, and somatic paired analyses.
license: MIT
metadata:
  version: "0.1.0"
  author: ClawBio
  domain: genomics
  tags:
    - sarek
    - variant-calling
    - germline
    - somatic
    - tumor-only
    - wes
    - wgs
    - nextflow
    - nf-core
  inputs:
    - name: samplesheet
      type: file
      format:
        - csv
        - tsv
        - yaml
        - yml
        - json
      description: >
        nf-core/sarek samplesheet. Accepted formats: .csv, .tsv, .yaml, .yml,
        .json. Always-required columns: patient, sample. Step-dependent columns:
        for mapping, lane plus one input mode — fastq_1 + fastq_2,
        OR spring_1 (+ optional spring_2), OR bam (uBAM); for markduplicates/
        prepare_recalibration/variant_calling, bam+bai or cram+crai; for
        recalibrate, additionally table; for annotate, vcf (+ optional
        variantcaller). Optional metadata columns: sex (XX/XY/NA), status
        (0=normal, 1=tumor), contamination (varlociraptor). A row may not mix
        FASTQ, BAM/CRAM, and VCF input modes.
      required: true
  outputs:
    - name: report
      type: file
      format:
        - md
      description: Wrapper run summary, per-tool VCF inventory, and downstream handoff recommendations
    - name: result
      type: file
      format:
        - json
      description: Structured result payload with detected CRAMs, per-tool VCFs, MultiQC HTML, and provenance
  dependencies:
    python: ">=3.11"
    packages:
  demo_data:
    - path: demo/README.md
      description: Demo mode uses the upstream nf-core/sarek -profile test dataset rather than bundled FASTQs
  endpoints:
    cli: python clawbio.py run sarek-pipeline --input {samplesheet} --output {output_dir}
  openclaw:
    requires:
      bins:
        - python3
        - nextflow
        - java
      env:
      config:
    always: false
    emoji: 🧬
    homepage: https://github.com/ClawBio/ClawBio
    os:
      - darwin
      - linux
    install:
    trigger_keywords:
      - sarek
      - germline variant calling
      - somatic variant calling
      - tumor-normal pair
      - mutect2
      - strelka
      - haplotypecaller
      - ascat
      - WES variant calling
      - WGS variant calling
      - VEP annotation
      - SnpEff
      - GATK pipeline
      - nf-core
      - Nextflow
---

# nfcore-sarek-wrapper

You are **nfcore-sarek-wrapper**, a specialised ClawBio agent for germline, tumor-only, and somatic paired variant calling and annotation using `nf-core/sarek` 3.8.1.

## Trigger

**Fire when:**
- User wants to run `nf-core/sarek`
- User asks for germline variant calling from FASTQ, BAM, or CRAM
- User asks for somatic / tumor-normal paired variant calling
- User asks for tumor-only variant calling
- User mentions GATK HaplotypeCaller, Mutect2, Strelka, ASCAT, ControlFREEC, Manta, TIDDIT, MSIsensor2, MSIsensor-pro, FreeBayes, DeepVariant, or Sentieon (TNscope, Haplotyper, DNAscope)
- User wants WES or WGS variant calling with strict preflight, reproducibility outputs, and downstream handoff
- User asks to annotate VCFs with VEP or SnpEff
- User mentions UMI consensus calling with fgbio for germline/somatic variants

**Do NOT fire when:**
- User has FASTQ for bulk RNA-seq → route to `nfcore-rnaseq-wrapper`
- User has FASTQ for single-cell RNA-seq → route to `nfcore-scrnaseq-wrapper`
- User already has an annotated VCF and wants ACMG/AMP interpretation → route to `clinical-variant-reporter`
- User wants a clinical PDF report from a WES markdown summary → route to `wes-clinical-report-en` or `wes-clinical-report-es`
- User asks about PharmGx, PRS, methylation, or pharmacogenomics

## Scope

One skill, one task: orchestrate `nf-core/sarek` 3.8.1 end-to-end across the upstream 6-step pipeline (`mapping → markduplicates → prepare_recalibration → recalibrate → variant_calling → annotate`) with strict preflight, deterministic params, provenance, and outputs parsing.

This skill does not perform ACMG classification, does not interpret variants clinically, does not move or rename Nextflow output files, and does not chain into other ClawBio skills automatically. Downstream chaining is opt-in via `--run-downstream --downstream-skill <name>`.

## Why This Exists

- **Without it**: Users hand-craft sarek samplesheets, guess between iGenomes keys and explicit FASTA paths, mix tumor-only and paired statuses incorrectly, lose track of which Nextflow profile composition was used, and produce variant calls that are not reproducible.
- **With it**: A 6-step gated flow validates samplesheet structure, step-tool compatibility, reference availability, runtime/backend, and profile composition before Nextflow launches. Every run emits `params.yaml`, `commands.sh`, `manifest.json`, and a checksums bundle.
- **Why ClawBio**: Local-first, pinned to nf-core/sarek 3.8.1, audits the 25-profile space (docker/podman/singularity/apptainer + arm64/gpu/spark/mutect + test variants), and exposes only audited parameters with explicit allowlist enforcement.

## Core Capabilities

1. **Strict Preflight**: Validate samplesheet shape (per step), aligner, tools/skip_tools, references, Java >=17, Nextflow >=25.10.2, backend, UMI options, and resume-state drift.
2. **Profile Composition**: Compose docker/singularity/etc. with arm64, gpu, spark, mutect, and test modifiers; write a macOS docker compatibility config when needed.
3. **Audited Execution**: Run `nf-core/sarek` 3.8.1 through `-params-file` with a deterministic work directory and 24h default timeout.
4. **Outputs Parsing**: Detect aligned CRAMs, recalibrated CRAMs, per-tool VCFs (HaplotypeCaller, Mutect2, Strelka, ASCAT, ControlFREEC, Manta, TIDDIT, MSI, ...), annotated VCFs (SnpEff, VEP, merge, bcftools, SnpSift), and MultiQC.
5. **Reproducibility Bundle**: Write `commands.sh`, `params.yaml`, `manifest.json`, checksums, `environment.yml`, and provenance JSON under `reproducibility/`.
6. **Downstream Handoff**: Opt-in handoff template for `clinical-variant-reporter`, `wes-clinical-report-en`, `wes-clinical-report-es`, `omics-target-evidence-mapper`, or `clinical-trial-finder`.

## Steps

| `--step` | Inputs required | Best for |
|---|---|---|
| `mapping` (default) | `lane` plus one of: `fastq_1`+`fastq_2`, `spring_1`(+ optional `spring_2`), or `bam` (uBAM) | Standard FASTQ/Spring/uBAM-to-VCF runs |
| `markduplicates` | Aligned `bam`+`bai` or `cram`+`crai` | Restart from alignment |
| `prepare_recalibration` | Deduplicated BAM/CRAM | Pre-BQSR restart |
| `recalibrate` | BAM/CRAM + `table` | Restart at BQSR apply |
| `variant_calling` | Recalibrated BAM/CRAM | Tool re-run without realignment |
| `annotate` | `vcf` (+ optional `variantcaller`) | Annotate existing variant calls |

## Input Formats

The samplesheet may be `.csv`, `.tsv`, `.yaml`, `.yml`, or `.json` (CSV/TSV are
delimited; YAML/JSON are a top-level list of row records). File-column values
must be local paths by default (local-first): remote URLs (`https://`, `s3://`,
`gs://`, `ftp://`, …) — and remote reference paths — are rejected at preflight
(`REMOTE_INPUT_NOT_ALLOWED`) unless you pass `--allow-remote-inputs`, which also
logs a runtime warning naming every path fetched over the network. (The public
iGenomes mirror base and the object-store `--work-dir` are not gated.)

| Mode | Required Fields | Example |
|--------|-----------------|---------|
| Mapping (FASTQ) | `patient`, `sample`, `lane`, `fastq_1`, `fastq_2` | `samplesheet.csv` |
| Mapping (Spring) | `patient`, `sample`, `lane`, `spring_1` (+ optional `spring_2`) | `samplesheet_spring.csv` |
| Mapping (uBAM) | `patient`, `sample`, `lane`, `bam` | `samplesheet_ubam.csv` |
| BAM/CRAM restart | `patient`, `sample`, plus `bam`+`bai` or `cram`+`crai` | `samplesheet_bam.csv` |
| Recalibrate restart | above plus `table` | `samplesheet_recal.csv` |
| Annotate | `patient`, `sample`, `vcf` (+ optional `variantcaller`) | `samplesheet_vcf.csv` |
| Demo mode | none | `python clawbio.py run sarek-pipeline --demo` |

Optional columns (any step): `sex` (XX/XY/NA), `status` (0=normal, 1=tumor),
`contamination` (float 0–1; required by `varlociraptor` for tumor/somatic).

Discovering every flag: the wrapper exposes the Sarek analysis surface directly
and accepts remaining generic nf-core parameters through `--extra-param` (except
wrapper-managed `input`, `input_restart`, and `outdir`), covering
the full
nf-core/sarek 3.8.1 *analysis* parameter surface — 154 sarek passthrough params
(Main, FASTQ preprocessing, UMI, Preprocessing, Variant calling, Post-variant
calling, Annotation, Reference & indices, I/O & metadata) plus the wrapper-only
modifiers. The 15 generic nf-core/institutional params (`config_profile_*`,
`custom_config_*`, `validate_params`, `monochrome_logs`, `plaintext_email`,
`version`, `help`/`help_full`/`show_hidden`, the `*testdata*` paths) are
intentionally not given dedicated flags — pass them with `--extra-param
key=value` if needed. `python clawbio.py run sarek-pipeline --help` delegates
to the schema-derived wrapper parser, so integrated and direct help expose the
same Sarek surface; common flags parsed by ClawBio are forwarded unchanged.

## Workflow

1. **Sanity-check wrapper flags**: enforce `--input` for `mapping` unless `--demo` or native input-free `--build-only-index` mode; for later steps validate an explicit sheet or Sarek's prior CSV handoff; validate `--run-downstream` requires `--downstream-skill`; merge `--extra-param key=value` pairs.
2. **Compose profile**: merge user backend (docker/singularity/...) with `--arm`, `--gpu`, `--spark-profile`, `--mutect-profile`, and `--demo` (`test`) tokens.
3. **Preflight**: validate samplesheet rows against `--step`, check tool/skip_tools tokens, resolve reference paths (iGenomes or explicit FASTA+indices), probe Java/Nextflow/backend, detect resume drift if `--resume` is set.
4. **Build params**: assemble the effective `params.yaml` from CLI flags + extras + step-dependent defaults; clear all reference flags when `--demo` is set.
5. **Execute Nextflow**: launch with composed profile, `-params-file params.yaml`, deterministic `-work-dir`, streamed stdout/stderr.
6. **Parse outputs**: detect aligned/recalibrated CRAMs, per-tool VCFs (§1–§6 layout), annotated VCFs, MultiQC, and pipeline_info.
7. **Write provenance + report**: emit `report.md` and `result.json` at the output root (matching nfcore-rnaseq/scrnaseq), and under `reproducibility/` emit `commands.sh`, `params.yaml`, the normalized samplesheet, `environment.yml`, `checksums.sha256`, and seven JSON files (`manifest.json`, `parameters.json`, `samplesheet.json`, `pipeline_source.json`, `tool_versions.json`, `outputs.json`, `compatibility_policy.json`).

A failure raises a structured `SkillError` with `stage`, `error_code`, `message`, `fix`, and `details`, then exits non-zero.

## CLI Reference

```bash
# Preflight only; no Nextflow execution
python clawbio.py run sarek-pipeline \
  --input samplesheet.csv --output ./sarek_check --check \
  --genome GATK.GRCh38 --tools haplotypecaller

# Demo mode using upstream -profile test
python clawbio.py run sarek-pipeline --demo --output /tmp/sarek_demo

# Germline WES with custom targeted reference resources
python clawbio.py run sarek-pipeline \
  --input samplesheet.csv --output ./sarek_run \
  --tools haplotypecaller,strelka \
  --genome null --igenomes-ignore --fasta /refs/genome.fa \
  --known-indels /refs/known_indels.vcf.gz \
  --wes --intervals exome_targets.bed

# Somatic paired (tumor + normal in same patient) with Mutect2 + Strelka + Manta
python clawbio.py run sarek-pipeline \
  --input samplesheet_paired.csv --output ./sarek_somatic \
  --tools mutect2,strelka,manta,vep \
  --genome GATK.GRCh38

# Tumor-only with Mutect2 + PON
python clawbio.py run sarek-pipeline \
  --input samplesheet_tumor_only.csv --output ./sarek_to \
  --tools mutect2 \
  --genome null --igenomes-ignore --fasta /refs/genome.fa \
  --known-indels /refs/known_indels.vcf.gz \
  --pon /refs/pon.vcf.gz --pon-tbi /refs/pon.vcf.gz.tbi \
  --germline-resource /refs/af-only.vcf.gz --germline-resource-tbi /refs/af-only.vcf.gz.tbi

# Explicit FASTA reference (non-default genome build)
python clawbio.py run sarek-pipeline \
  --input samplesheet.csv --output ./sarek_run \
  --genome null --igenomes-ignore \
  --fasta /refs/genome.fa --fasta-fai /refs/genome.fa.fai --dict /refs/genome.dict \
  --bwa /refs/bwa/

# ARM (Apple M-series, AWS Graviton) — composes -profile docker,arm64
python clawbio.py run sarek-pipeline \
  --input samplesheet.csv --output ./sarek_arm \
  --profile docker --arm --genome GATK.GRCh38

# Opt-in downstream handoff to clinical-variant-reporter
python clawbio.py run sarek-pipeline \
  --input samplesheet.csv --output ./sarek_run \
  --tools haplotypecaller,vep \
  --genome GATK.GRCh38 \
  --run-downstream --downstream-skill clinical-variant-reporter

# Wrapper runtime controls (parity with scrnaseq/rnaseq):
#   --timeout-hours N   wall-clock cap (default 24h; 0 disables for HPC/cloud)
#   --work-dir PATH     Nextflow work dir (local path or object-store URI; default <output>/upstream/work)
#   --nextflow-config / -c / --config   extra Nextflow config file(s), repeatable
#   --allow-pipeline-version-override    run a non-3.8.1 --pipeline-version at your own risk
#   --allow-remote-inputs               opt in to remote inputs/refs (default local-first)
python clawbio.py run sarek-pipeline \
  --input samplesheet.csv --output ./sarek_run \
  --genome GATK.GRCh38 --tools haplotypecaller \
  --timeout-hours 0 --work-dir s3://my-bucket/sarek/work
```

## Demo

```bash
python clawbio.py run sarek-pipeline --demo --output /tmp/sarek_demo
```

Expected output: upstream `nf-core/sarek -profile test` outputs (synthetic small dataset) under `upstream/results/`, `report.md` and `result.json` at the output root, and the ClawBio `reproducibility/` bundle (params/commands/samplesheet snapshots, provenance JSON, `environment.yml`, `checksums.sha256`).

## Algorithm / Methodology

Key methods:
- Local data paths inside a samplesheet are resolved against its directory and
  written as absolute POSIX paths; remote data URLs are passed through unchanged.
  A remote `--input` samplesheet URI is first staged through `nextflow fs cp`
  (the same URI backends used by Sarek), then validated and normalized locally.
- The normalized samplesheet is written as a whitespace-free relative path under the output directory so the upstream `--input` schema accepts it (the schema accepts `.csv`, `.tsv`, `.yaml`, `.yml`, `.json`).
- Reference handling follows Sarek's two documented modes: use `--genome <iGenomes>` and optionally override individual reference files (or pass `false` for a resource that should not be used), or use `--genome null --igenomes-ignore --fasta <reference>` when no iGenomes reference files should be loaded. Optional FASTA indices and tool resources may be supplied in either mode when appropriate.
- In `--build-only-index` mode, Sarek intentionally supplies an empty samplesheet channel: the wrapper does not require sample pairing or per-sample caller outputs, while it preserves upstream global resource guards (including BQSR guards on preprocessing start steps) and captures published reference outputs.
- Tool×mode compatibility is evaluated per-patient: a tool is accepted when at least one patient matches its required mode, so mixed samplesheets (germline-only patients alongside tumor/normal pairs) are valid. Paired-only tools (`ascat`, `msisensorpro`, `muse`) still need at least one patient with both `status=0` and `status=1`.
- Mutect2 without an effective PON or germline resource emits a preflight warning but does not block; resources inherited from an iGenomes bundle count as effective. The bundled `GATK.GRCh38` PON still emits a recommendation to use a project-specific PON.
- Paired somatic Mutect2 cannot be run with `--no-intervals`; the upstream schema explicitly marks that combination unsupported.
- `--snv-consensus-calling` requires `--normalize-vcfs`, as enforced by the upstream workflow before post-variant processing.
- `--use-gatk-spark markduplicates` is incompatible with header/positional UMI dedup (`--umi-in-read-header` or `--umi-location`); it is fine with `--umi-read-structure` (fgbio consensus runs upstream).
- ASCAT requires an effective `--ascat-genome`, `--ascat-alleles`, and `--ascat-loci` (which supported iGenomes bundles can provide). With `--wes`, custom `--ascat-alleles`, `--ascat-loci`, `--ascat-loci-gc`, and `--ascat-loci-rt` resources are recommended; the wrapper warns because Sarek documents its iGenomes ASCAT resources as unsuitable for WES.

## Example Queries

- "Run nf-core/sarek for germline variant calling on these WES FASTQs"
- "Call somatic variants from this tumor-normal pair with Mutect2 and Strelka"
- "Annotate this VCF with VEP using sarek"
- "Tumor-only Mutect2 with PON for our WES cohort"
- "Restart sarek at the recalibration step"

## Example Output

```
output/                                       # the --output directory
├── .nextflow/                                # Nextflow cache/history (framework-created; excluded from checksums)
├── .nextflow.log                            # Nextflow launch log (framework-created; excluded from checksums)
├── upstream/
│   ├── results/                              # Nextflow --outdir
│   │   ├── csv/                              # handoff CSVs (mapped/markduplicates/recalibrated/variantcalled); legacy fallback: preprocessing/csv/
│   │   ├── preprocessing/
│   │   │   ├── mapped/                       # §1 aligned CRAMs (one per sample/lane)
│   │   │   ├── markduplicates/               # §2 deduplicated CRAMs
│   │   │   └── recalibrated/                 # §3 BQSR-recalibrated CRAMs
│   │   ├── variant_calling/
│   │   │   ├── haplotypecaller/              # §4 germline VCFs
│   │   │   ├── mutect2/                      # somatic / tumor-only VCFs
│   │   │   ├── strelka/                      # somatic + germline VCFs
│   │   │   ├── manta/                        # SV VCFs
│   │   │   ├── ascat/                        # CNV / purity / ploidy
│   │   │   ├── controlfreec/                 # CNV
│   │   │   ├── tiddit/                       # SV
│   │   │   ├── bcftools/                     # mpileup caller output
│   │   │   ├── msisensor2/                   # MSI
│   │   │   └── msisensorpro/                 # MSI (paired MSIsensorPro)
│   │   ├── annotation/<variantcaller>/<sample_or_pair>/   # §5 SnpEff/VEP/merge/bcftools/SnpSift annotated VCFs
│   │   ├── multiqc/                          # §6 MultiQC HTML + data
│   │   ├── pipeline_info/
│   │   └── reports/
│   └── work/                                 # Nextflow work directory
├── report.md                                # human-readable run summary (output root)
├── result.json                              # machine-readable run summary (output root)
├── check_result.json                        # written only with --check (preflight-only mode); parallel to result.json
├── logs/                                    # Nextflow stdout.txt / stderr.txt (real runs only; excluded from checksums)
└── reproducibility/                          # replay + provenance bundle
    ├── samplesheet.valid.csv                 # or samplesheet.demo.csv in --demo mode
    ├── params.yaml
    ├── commands.sh
    ├── remap_paths.py
    ├── environment.yml
    ├── checksums.sha256
    ├── compatibility_policy.json
    ├── parameters.json
    ├── samplesheet.json
    ├── pipeline_source.json
    ├── tool_versions.json
    ├── outputs.json                          # omitted if outputs parsing was skipped
    ├── manifest.json
    ├── macos_docker.config                   # written only on macOS + docker backend
    └── sarek_downstream_handoff.{sh,json}    # written only when --run-downstream is set
```

`report.md`, `result.json`, and `logs/` sit at the output root — the same layout
as the nfcore-rnaseq and nfcore-scrnaseq wrappers — so a consumer finds
`<output>/result.json` for any of the three pipelines. The `reproducibility/`
directory holds the portable replay + provenance bundle.

## Output Structure

Under the `output/` root the wrapper writes two child directories — `upstream/` (the Nextflow `results/` tree plus its `work/` directory) and `reproducibility/` (the portable replay + provenance bundle: the params/commands/samplesheet snapshots, the seven JSON provenance files, `environment.yml`, `checksums.sha256`, and — macOS + docker only — `macos_docker.config`) — alongside the run-summary files `report.md` and `result.json`. A real run also writes a root-level `logs/` directory (Nextflow `stdout.txt`/`stderr.txt`), and `--check` writes `check_result.json` at the root (parallel to `result.json`); both placements match the nfcore-rnaseq and nfcore-scrnaseq wrappers. Nextflow itself additionally writes its own hidden bookkeeping in the launch directory — `.nextflow/` (cache/history) and `.nextflow.log` — because the wrapper runs Nextflow with `cwd = output_dir` so the relative `input`/`outdir` paths resolve; both are excluded from `checksums.sha256` (`.nextflow` directory and any `.log` file are skipped). There is no separate top-level `provenance/` directory; all provenance JSON is co-located in `reproducibility/`. The `reproducibility/` tree, the root `logs/` directory, and the root summaries `report.md`/`result.json`/`check_result.json` are all excluded from `checksums.sha256`, so execution logs and wrapper summaries never enter the manifest. The §1–§6 layout in `outputs_parser.py` corresponds to: §1 mapped, §2 markduplicates, §3 recalibrated, §4 per-tool variant calls, §5 annotation, §6 MultiQC.

**Cross-machine / cross-OS portability.** The bundle stores absolute data/reference paths (required by Nextflow) but ships a stdlib-only `remap_paths.py` to rebase them on any host: `--old/--new` rewrites samplesheet data paths, `--refs-old/--refs-new` rewrites reference/index paths in `params.yaml` (and `commands.sh` if any were added there), `--output-dir <new-path>` rewrites the baked `--output` in `commands.sh` when you relocate the run, and `--verify` confirms every path resolves before replay. (The scrnaseq bundle self-relocates and needs no `--output-dir`; it accepts the flag only for parity.) URIs (`s3://`, `https://`, …) and the `false` disable sentinel are preserved. All bundle files use POSIX paths and `utf-8`/`\n`, so macOS↔Linux replay is byte-stable. The recommended replay path is a self-contained `bash commands.sh` (no environment variable required): it self-anchors via `BASH_SOURCE`, pins the Nextflow engine with `NXF_VER`, and applies the macOS-only Docker config through a `uname`-gated `-c reproducibility/macos_docker.config`, so the same bundle replays identically on Linux and macOS.

## Dependencies

**Required**
- Python >=3.11
- Java >=17
- Nextflow >=25.10.2
- One execution backend: Docker, Singularity, Apptainer, Podman, Conda/Mamba, Shifter, or Charliecloud

## Gotchas

- **ASCAT for WES should use custom resources.** Sarek allows `--wes --tools ascat` with iGenomes but warns that its default ASCAT resources are not suited for WES. The wrapper preserves that warning; for production WES runs, use `--genome null --igenomes-ignore` and provide target-appropriate `--ascat-genome`, `--ascat-alleles`, `--ascat-loci`, `--ascat-loci-gc`, and `--ascat-loci-rt`.
- **SnpEff needs a database identifier.** `snpeff` and `merge` require an effective `snpeff_db`: either pass `--snpeff-db` or use an iGenomes genome that supplies it. This is still required when downloading/building annotation caches because Sarek uses the database identifier to select/download the cache.
- **Mutect2 without PON or germline resource emits unreliable calls.** The model will want to skip the PON/germline resource for "simplicity". Do not. Mutect2 without a panel-of-normals returns the union of true somatic calls and recurrent technical artifacts, and without a germline resource no germline-based filtering is applied. The wrapper warns when neither explicit nor inherited resources are effective, and warns separately when the generic `GATK.GRCh38` PON is used; pass a technically matched `--pon`/`--pon-tbi` and appropriate `--germline-resource` for production runs.
- **MarkDuplicatesSpark cannot do UMI-based deduplication.** The model will want to set `--use-gatk-spark markduplicates` alongside header/positional UMIs (`--umi-in-read-header` or `--umi-location`) to "speed up dedup". Do not — MarkDuplicatesSpark cannot perform UMI-aware deduplication, so the wrapper rejects that combination (matching sarek). Note: `--use-gatk-spark markduplicates` IS compatible with `--umi-read-structure`, because fgbio collapses UMIs into consensus reads *before* deduplication.
- **iGenomes defaults to GATK.GRCh38 — a custom `--fasta` disables it automatically.** nf-core/sarek 3.8.1 defaults `genome = 'GATK.GRCh38'` in `nextflow.config`. When you pass `--fasta <reference>` **without** `--genome`, the wrapper sets `igenomes_ignore=true` for you (emitting a WARNING) so Sarek does not load the default GATK.GRCh38 iGenomes bundle and fail while validating its ~20 remote `s3://ngi-igenomes/...` reference paths as local files. You therefore no longer need to add `--igenomes-ignore` by hand for a custom build (passing it explicitly is still honoured, and combining `--fasta` with an explicit `--genome` is treated as a documented partial override and left untouched). Sarek documents FASTA as the only required custom-reference file and can build missing indices: the `.fai` index and `.dict` GATK needs are built automatically from your `--fasta` — you do not have to pre-generate or pass them (do so only if you already have them, and then they must be co-located and consistent with the FASTA). For a partial override, keep `--genome <iGenomes>` and pass only the replacement resource files, as described in the official usage guide.
- **BQSR requires known sites unless skipped.** Base recalibration (`baserecalibrator`) runs by default on the preprocessing start steps (`mapping`, `markduplicates`, `prepare_recalibration`, `recalibrate`), independently of the variant-calling `--tools` you request, and needs `--dbsnp` or `--known-indels` — supplied explicitly or inherited from an iGenomes `--genome`. The model will want to silently skip BQSR to "make it run" on a small or non-model genome. Do not. The wrapper does not inject any hidden default: a start step without those resources is rejected at preflight (`MISSING_REFERENCE`), and the error now states that baserecalibrator runs by default so the requirement is not surprising when you only asked for, say, `--tools haplotypecaller`. To run without known sites (e.g. a non-model organism or a synthetic test genome), make the choice explicit by adding `baserecalibrator` to `--skip-tools` (the nf-core-native `--skip_tools` spelling is accepted too). This mirrors nf-core/sarek 3.8.1 exactly (the upstream `-profile test` ships its own known-sites resources, so `--demo` is exempt).
- **Tool×mode pairing is inferred from the `status` column — there is no `--tumor-only` flag.** The model will want to write an all-`status=1` samplesheet and still ask for paired tools. Do not. The wrapper's preflight enforces each tool's supported modes against the samplesheet:
  - **Paired-only** (need both a `status=0` normal and a `status=1` tumor under the same `patient`): `ascat`, `msisensorpro`, `muse`.
  - **Tumor-only** (need a tumor; no normal required): `lofreq`, `msisensor2`.
  - **Tumor-only OR paired** (auto-routed): `mutect2`, `sentieon_tnscope`, `controlfreec`.
  - **Any mode** (germline + tumor-only + somatic): `manta`, `tiddit`, `cnvkit`, `freebayes`, `varlociraptor`.
  - **Germline, tumor-only, or paired-normal routing**: `mpileup`. **Germline or paired**: `strelka`. **WGS germline or paired only**: `indexcov` (`--wes` does not execute it upstream).
  - **Needs a normal/germline sample**: `haplotypecaller`, `deepvariant`, `sentieon_haplotyper`, `sentieon_dnascope`. These run on the normal sample of a tumor/normal pair as well as on standalone normals.
  - The check is **per-patient**, like sarek: a tool is accepted when *at least one* patient supplies its required input. Mixed samplesheets and paired cohorts are valid — for a paired patient, `haplotypecaller` runs on the normal while `mutect2` runs on the tumor/normal pair.
- **`--demo` clears all reference flags.** The model will want to combine `--demo` with `--genome GATK.GRCh38` or `--fasta`. Do not. The `--demo` flag forces the upstream `-profile test` dataset, which ships its own tiny reference; the wrapper clears every reference-path flag (genome, igenomes_base, fasta, intervals, dbsnp, known_indels, known_snps, germline_resource, pon, ...) before they reach `params.yaml` and ignores `--input`.
- **`--demo` requires network access.** The upstream `-profile test` fetches its sample FASTQs and reference files from remote GitHub URLs (nf-core's design — the wrapper does not bundle local test data). On an offline/sandboxed host set `NXF_OFFLINE`, and the wrapper fails fast at preflight with `DEMO_REQUIRES_NETWORK` and a clear message, instead of a cryptic Nextflow `does not exist` abort during schema validation. This does **not** violate the local-first guarantee, which governs your *genetic data* (never uploaded); `--demo` only *downloads* nf-core's public test data. For a fully offline run, use a real analysis with your own local `--input` and references.
- **`--resume` is rejected when the manifest drifts.** Changes to pipeline source, profile composition, step, aligner, tools, skip_tools, analysis_mode, wes, joint_germline, joint_mutect2, the complete emitted parameter checksum, reference fingerprints, or samplesheet checksum invalidate the Nextflow work directory. The wrapper refuses to resume in those cases — re-run with a fresh `--output` directory.
- **`--output` must be outside the ClawBio source tree.** An output directory inside the repository is rejected at preflight with `OUTPUT_DIR_INSIDE_REPO` so multi-gigabyte pipeline artifacts never pollute (or get committed to) the checkout — choose a path under your analysis workspace. This matches the nfcore-rnaseq and nfcore-scrnaseq wrappers.
- **`arm64` profile implicitly enables Wave.** The model will want to use `--arm` on offline / air-gapped machines. Do not without preparation. Sarek's `arm64` profile sets `wave.enabled=true` + `wave.freeze=true` + `wave.strategy='conda,container'` so containers are rebuilt for arm via the upstream Wave service. If outbound traffic to `wave.seqera.io` is blocked, the run stalls. On offline Apple Silicon, instead build images once via Wave on a connected host, mirror to a local registry, and add a custom `-c` config pointing `process.container` to the mirror.
- **Sentieon callers and aligner require a license.** The model will want to use `--aligner sentieon-bwamem` or `--tools sentieon_haplotyper,sentieon_dnascope,sentieon_tnscope` without setting `SENTIEON_LICENSE` (or `SENTIEON_LICENSE_BASE64`) first. Do not. Sarek does not vendor a Sentieon license; the upstream task fails late with a cryptic Sentieon binary error. Export `SENTIEON_LICENSE=host:port` (server) or `SENTIEON_LICENSE_BASE64=...` (offline) before invoking `clawbio.py run sarek-pipeline`, and document the license source in your run notes.
- **nf-core-native (snake_case) flag spellings are accepted.** You can paste an upstream nf-core command's parameters verbatim: `--skip_tools`, `--fasta_fai`, `--known_indels`, … all work as aliases of the hyphenated wrapper flags, both via `clawbio.py run sarek-pipeline` and when invoking the wrapper directly. There is no need to convert underscores to hyphens by hand.
- **`split_fastq` uses the nf-core default (50 M reads/chunk) in normal runs, but `--demo` (`-profile test`) disables it — results can differ on tiny datasets.** nf-core/sarek 3.8.1 defaults `split_fastq = 50000000`; the upstream `test` profile overrides it to `0` (splitting off). The wrapper faithfully passes the pipeline default — it does **not** silently inject `0` — so a normal run and a `--demo` run of the *same* tiny dataset can legitimately produce different variant sets: FASTQ splitting adds a `fastp/` split-merge stage whose read-group handling can perturb calls on degenerate, very-small inputs (a known artifact of running production splitting on test-scale data, not a wrapper bug). This never matters at real WGS/WES depth. If you are reproducing a `-profile test` result, benchmarking against demo output, or running a deliberately small dataset, pass `--split-fastq 0` to match the test profile exactly. Do **not** change the default for production cohorts — splitting parallelises mapping and is beneficial at scale. The divergence is a *profile* difference (`split_fastq`, and the test profile's own `--tools strelka` / reference choices), **not** a CPU/thread effect: nf-core/sarek pins `bwa mem -K 100000000` on every mapping process, which fixes the per-batch base count and makes alignment bit-identical regardless of the thread count — so the host-scaled `resourceLimits` (e.g. 28 vs the demo's 4 CPUs) does **not** change the calls. Do not attribute a demo-vs-normal difference to core count.
- **Host-limited memory is auto-capped on docker runs; IPv6-only networks are an environment issue — read the failure hint.** On a docker backend the wrapper writes a `process.resourceLimits` config scaled to this host (physical RAM minus headroom, and the Docker VM's `MemTotal` on macOS) so a real run does not abort with `Process requirement exceeds available memory` when an nf-core default request is larger than your machine (`--demo` is exempt: `-profile test` carries its own limits). If it still aborts (e.g. a non-docker backend, or a single process that genuinely needs more RAM than the host has), override with your own `-c` config, e.g. `process { resourceLimits = [ memory: '12.GB', cpus: 4 ] }` — do **not** delete resource labels to force it through. On an IPv6-only / NAT64 host the JVM prefers IPv4 and downloads fail with `Network is unreachable`; export `NXF_OPTS='-Djava.net.preferIPv6Addresses=true'` and re-run. The wrapper inherits your environment and never overrides `NXF_OPTS`, so your setting takes effect.
- **A demo replay can fail transiently on the first attempt — just re-run `bash commands.sh`, it auto-resumes.** The upstream `-profile test` stages its references from remote GitHub URLs, and nf-core/sarek's own `conf/test.config` sets `igenomes_base`/`modules_testdata_base_path` with a trailing slash (`…/test-datasets/modules/data/`), so joined reference URLs can contain a `//` that some CDNs momentarily 404, and concurrent remote staging occasionally races. This is upstream/network behaviour, not a wrapper or bundle defect — the pinned command is byte-identical every time. The generated `commands.sh` handles it: on a first run it starts clean, and because it adds `-resume` whenever a prior Nextflow session (`.nextflow/`) exists, simply invoking `bash commands.sh` **again** resumes from the cached tasks and completes the run without recomputing. Only investigate further if a second `-resume` attempt reproduces the *same* failing task deterministically.
- **`prov.*` / `validation.*` config keys belong to the pipeline's own plugins, not to the wrapper. Do NOT try to silence them by adding params or configs.** nf-core/sarek 3.8.1's `nextflow.config` declares `plugins { … nf-prov@1.2.2, nf-schema@2.6.1 }` and, right alongside them, the matching `prov { enabled = true; formats { bco { file = "${params.outdir}/pipeline_info/manifest_${params.trace_report_suffix}.bco.json" } } }` and `validation { defaultIgnoreParams = […] }` scopes. Those scopes are *owned by the plugins*, so Nextflow only reports them as unrecognised when the plugin is not active — an engine/environment condition, not a bundle defect. The wrapper writes **none** of these keys: `params.yaml` carries only the audited CLI surface, and `commands.sh` runs the pinned pipeline unmodified, so a plain upstream `nextflow run nf-core/sarek -r 3.8.1` behaves identically. Injecting `prov`/`validation` into `params.yaml` or a `-c` config to quiet a warning would override the pipeline's own provenance configuration and is exactly the wrong fix. A regression test pins that these keys never enter `params.yaml`.
- **Each run writes a fresh, timestamped set of `pipeline_info/` reports — a replay adds files rather than overwriting them.** Upstream sarek names its execution report, timeline, trace, DAG and nf-prov BCO manifest with `${params.trace_report_suffix}`, which defaults to the run's start timestamp. So re-running `commands.sh` in place leaves the previous run's reports next to the new ones. This is upstream nf-core behaviour, and it is deliberate — those files are the audit trail of *each* execution. It does not invalidate `checksums.sha256`, which pins only the bundle's own artefacts (samplesheet, `params.yaml`, logs, and the analysis outputs), never the timestamped report files. Do not "clean" it by pinning `trace_report_suffix`: that would make a replay silently overwrite the original run's execution evidence.
- **The nf-schema "positional argument 'nextflow' has been detected" warning is a benign upstream false positive.** nf-core/sarek 3.8.1's schema plugin inspects `workflow.commandLine`, whose first token is literally `nextflow` (argv[0]), and mislabels it as a positional argument. The wrapper's command is well-formed (`nextflow run nf-core/sarek -r 3.8.1 -profile … -params-file …`, identical in shape to the rnaseq/scrnaseq wrappers, which use newer plugin versions that do not warn). It does not affect the run — do **not** try to "fix" it by altering the command.

## Safety

- No patient data is bundled.
- Demo mode uses upstream `-profile test` data.
- The wrapper does not upload data.
- **Local-first by default**: remote samplesheet inputs and reference paths are rejected (`REMOTE_INPUT_NOT_ALLOWED`) unless `--allow-remote-inputs` is explicitly passed, which also logs a runtime warning naming every path fetched over the network. The iGenomes mirror base and the object-store `--work-dir` are not gated. `--allow-remote-inputs` relaxes only the wrapper's own preflight check: remote FASTQ/reference URIs are then written into the normalized samplesheet/`params.yaml` **verbatim** and staged natively by Nextflow at run time. The wrapper does not download them itself, so remote inputs require outbound network access and are incompatible with `NXF_OFFLINE` — under offline mode Nextflow's own file-existence validation (nf-schema) still runs and will fail on the remote paths.
- The wrapper does not pass arbitrary unvalidated Nextflow parameters; unknown native keys via `--extra-param` are passed through but tracked in provenance, while `input`, `input_restart`, and `outdir` remain wrapper-managed so normalized input and output provenance cannot be bypassed.
- `--resume` is rejected on manifest drift (pipeline source, profile, step, aligner, tools, skip_tools, analysis_mode, wes, joint_germline, joint_mutect2, effective params checksum, reference fingerprints, samplesheet checksum).

> ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions.

## Agent Boundary

Agent dispatches and explains; skill executes Nextflow. The agent must not invent variant-calling thresholds, fabricate reference paths, or override `--resume` drift errors. If a flag combination is rejected by preflight, the agent should explain the failure and suggest the corrected invocation — never bypass the check.

## Chaining Partners

This skill emits a downstream handoff template (under `reproducibility/sarek_downstream_handoff.sh`) when invoked with `--run-downstream --downstream-skill <name>`. Supported partners:

- `clinical-variant-reporter`: ACMG/AMP classification from VEP-annotated VCFs.
- `wes-clinical-report-en`: render an English-language WES clinical PDF from this run's markdown report.
- `wes-clinical-report-es`: render a Spanish-language WES clinical PDF from this run's markdown report.
- `omics-target-evidence-mapper`: aggregate target-level evidence across multi-omic sources using the gene list from annotated VCFs.
- `clinical-trial-finder`: match variants/genes to ClinicalTrials.gov and EUCTR trials.

The handoff is opt-in. Without `--run-downstream`, no template is written and no follow-on skill is launched.

## Maintenance

- **Review cadence**: quarterly, or whenever upstream `nf-core/sarek` ships a new release.
- **Staleness signals**: new sarek release (3.8.x → 3.9.x), new caller added upstream, schema_input.json changes, container manifest changes, deprecated profile tokens, ASCAT or VEP cache version bumps.
- **Deprecation**: a sarek major version bump that breaks the schema (3.x → 4.x) requires a new wrapper branch with regenerated `_SAREK_PASSTHROUGH_PARAMS`, regenerated preflight rules, and full test re-baselining. Before bumping the pinned version, audit `nextflow.config`, `assets/schema_input.json`, `nextflow_schema.json`, `docs/output.md`, and changed module configs.

## Citations

- nf-core/sarek 3.8.1: https://github.com/nf-core/sarek
- GATK Best Practices: https://gatk.broadinstitute.org/hc/en-us/sections/360007226651-Best-Practices-Workflows
- VEP: https://www.ensembl.org/info/docs/tools/vep/index.html
