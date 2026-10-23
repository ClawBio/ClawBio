# nfcore-sarek-wrapper

ClawBio wrapper around `nf-core/sarek` 3.8.1 covering mapping through annotation for germline, tumor-only, and somatic paired variant-calling analyses, with strict preflight, reproducibility artifacts, provenance, and explicit opt-in handoff to downstream ClawBio interpretation skills.

## Scope

- Upstream variant calling via Nextflow across the 6-step pipeline (`mapping → markduplicates → prepare_recalibration → recalibrate → variant_calling → annotate`).
- Germline, tumor-only, and somatic tumor-normal paired analyses.
- Supported callers: HaplotypeCaller, Mutect2, Strelka, ASCAT, ControlFREEC, Manta, TIDDIT, MSIsensor-pro, FreeBayes, DeepVariant, and Sentieon (TNscope/Haplotyper/DNAscope).
- Annotation with VEP and/or SnpEff.
- Validated samplesheet input (CSV/TSV/YAML/JSON) with normalized local paths (local-first by default; remote URLs require `--allow-remote-inputs`).
- Reproducibility bundle, provenance bundle, `report.md`, and `result.json`.

## Out Of Scope

- ACMG/AMP classification or clinical variant interpretation.
- Moving, renaming, or post-processing Nextflow output files.
- Automatic downstream chaining; handoff is opt-in via `--run-downstream --downstream-skill <name>`.
- Arbitrary free-form Nextflow passthrough flags (any non-promoted sarek parameter goes through `--extra-param key=value`; custom config via `-c/--config`).

## Quick Start

```bash
# See every supported flag (177-flag surface: 154 Sarek passthrough + 23 wrapper controls)
python skills/nfcore-sarek-wrapper/nfcore_sarek_wrapper.py --help

# Synthetic demo using the upstream -profile test dataset
python clawbio.py run sarek-pipeline --demo --output "$HOME/clawbio-runs/sarek_demo"
```

For real data:

```bash
python clawbio.py run sarek-pipeline \
  --input samplesheet.csv \
  --output "$HOME/clawbio-runs/sarek_run" \
  --genome GATK.GRCh38 \
  --tools haplotypecaller,vep
```

Use `--check` to run preflight without launching Nextflow. Checkout-mode runs without `--output` default to `~/clawbio-runs/`; explicit destinations must be outside the repository.

`clawbio.py run sarek-pipeline` promotes the common flags to first-class options; any other sarek parameter is passed through with `--extra-param key=value`.

See `SKILL.md` in this directory for trigger phrases, samplesheet schema per step, the full CLI reference, gotchas (ASCAT-WES, Mutect2 PON, GATK Spark + UMI, iGenomes default, somatic status pairing, `--demo` reference clearing), profile composition rules, downstream chaining, and the §1–§6 output layout.

## Validation and resume

nf-schema checks every parameter against the pinned `nextflow_schema.json` when Nextflow starts; the wrapper refuses only NaN and infinity, which nf-schema cannot represent, and a failed run names each value nf-schema refused. Malformed or ambiguous samplesheets fail before normalized snapshots are written. An atomic launch manifest supports resuming failed launches; resumed `--check` also validates params and samplesheet checksums.

Local Docker runs, including demos, cap both per-process requests and the local executor budget to detected host/Docker CPU and memory capacity. Site/cloud profiles and object-store work directories keep their own resource policy. Explicit `--nextflow-config` files have final precedence. Resource caps cannot make a dataset fit insufficient RAM.

Replay uses the smaller of recorded limits and current Docker capacity, on macOS and Linux, through a temporary config; the original checksummed config is preserved. User config snapshots retain final precedence, except for container images: a replay without a Nextflow cache applies `containers.config` last, so every process runs the digest it originally ran.

See the [compatibility audit](../../docs/nfcore-wrapper-audit.md) for upstream contracts, fixes, test results, and execution limits.
