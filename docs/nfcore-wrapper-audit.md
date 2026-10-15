# nf-core wrapper compatibility audit

This audit covers the three implemented Nextflow wrappers: `nfcore-rnaseq-wrapper`,
`nfcore-scrnaseq-wrapper`, and `nfcore-sarek-wrapper`. It compares wrapper behavior
with the exact upstream releases below, rather than with a moving default branch.
The starting ClawBio revision was `d67bd59af4c8c878e01a5f8deb8cd9a22e2e80b9`.
The completed audit was subsequently integrated with ClawBio `main` at
`9405cfd7dcbe2c86046091d58edacb808755c3f7` on 2026-10-14. The launcher audit
logging changes, updated locked dependencies, and newly registered skills were
retained; the catalog was regenerated from the combined tree. The focused matrix
below was repeated after this integration. The rebuilt wheel also passed installed
launcher success/failure audit logging and default-output routing checks.

## Alignment with PR #623

The audit now incorporates [PR #623](https://github.com/ClawBio/ClawBio/pull/623)
at `eb17be3604ec1cb5b4376d5ef9b9f9a3126232c2`. Its shared package is the
foundation; the competing `clawbio/common/nfcore.py` file has been removed.
The shared execution, command, path, failure-reporting, import-isolation, bundle
runtime, and bundle emitter modules are unchanged from that revision. The new
runner timeout, parser-environment, portable-remapping, checksum, demo-reporting,
and warning-relay fixes are retained.

Complementary audit helpers live in the relevant package modules: `parameters`
for strict mappings and scalar contracts, `source` for literal manifest versions,
`runtime` for Docker CPU probes, `cli` for manifest work-directory spelling,
`errors` for JSON-safe details, and `artifacts` for strict CSV input and atomic
prepared manifests. `resources` contains local Docker eligibility and the portable
replay guard. Wrappers still own scientific defaults, reference/input policies,
manifest payloads, checksums, and resource calculations. See
[nf-core wrapper architecture](nfcore-wrapper-architecture.md).

Sarek uses #623's `_run_with_samplesheet` boundary; fresh launches also stage
snapshots until complete emitted-parameter validation succeeds. A refused resume
preserves existing bundle files; a rejected fresh plan does not commit normalized
snapshots. The bulk
wrapper's test setup no longer replaces the public package with a stub, so shared,
wrapper, and public API tests can run together in either collection order. Tests
that probe ignored demo flags use schema-valid sample values and retain their
exhaustive warning assertions.

## Upstream contracts

| Wrapper | Pipeline release | Minimum Nextflow | Parameters in scalar contract |
| --- | --- | --- | --- |
| Bulk RNA-seq | [rnaseq 3.26.0](https://github.com/nf-core/rnaseq/tree/3.26.0) | 25.04.3 | 133 |
| Single-cell RNA-seq | [scrnaseq 4.1.0](https://github.com/nf-core/scrnaseq/tree/4.1.0) | 25.04.0 | 66 |
| Variant calling | [sarek 3.8.1](https://github.com/nf-core/sarek/tree/3.8.1) | 25.10.2 | 172 |

Each wrapper now includes `reproducibility/upstream_scalar_contract.json`, recording
the official schema URL, release, SHA-256, and scalar constraints. The generator
reads locally downloaded release sources and can check drift without network access:

```bash
python scripts/generate_nfcore_scalar_contracts.py --upstream-dir /path/to/releases --check
```

The release directory must contain `rnaseq-3.26.0`, `scrnaseq-4.1.0`, and `sarek-3.8.1`.
All three generated contracts matched the downloaded official schemas byte for byte.
The contracts validate emitted native types, finite numbers, numeric bounds,
enumerations, patterns, string lengths, and the numeric alternatives used by Sarek's
`split_fastq`. Existing pipeline-specific relationship and reference checks remain
in preflight. This is a projection of the parameter schemas, not a replacement for
nf-schema or its runtime file checks.

## Corrected behavior

### Input and parameter validation

- Reject duplicate or empty CSV headers, ragged records, malformed quoting, and
  invalid encodings before creating normalized samplesheet snapshots. Accept a
  UTF-8 BOM without treating it as part of the first column name.
- Reject duplicate YAML/JSON mapping keys and malformed Sarek parameter files with
  structured errors. Quoted booleans cannot silently enable scientific flags.
- Validate all emitted parameters against the pinned scalar constraints, including
  reference filename patterns and unsupported enum values.
- Preserve remote single-cell MultiQC configuration URIs instead of resolving them
  as local filesystem paths.
- Convert recognized single-cell protocol aliases to upstream's case-sensitive
  keys while preserving supported custom protocol strings.
- Apply Sarek's custom-FASTA iGenomes policy to parameter-file inputs as well as CLI
  flags. Demo mode also clears input/reference overrides from parameter files so
  the test profile owns its samples and reference data.
- Parse local pipeline manifest versions without confusing quoted braces, URLs, or
  comments. Sarek now checks local versions against the supported release policy.
- Route the three launcher aliases' checkout defaults to `~/clawbio-runs/` instead
  of a destination the wrappers reject inside the source tree. Preserve explicit
  destinations and other skills' output policy; correct the documented examples.
- Align the bulk wrapper's Python requirement with the repository's Python 3.11 floor.
- Forward `--work-dir` through the root launcher for all three pipelines, including
  resumed checks. Preserve Sarek's `--trace-report-suffix` instead of silently
  dropping it when the launcher parses the option.

### Resume and execution safety

- Write an atomic launch manifest before starting Nextflow. An interrupted or failed
  first launch therefore has a compatibility record for a subsequent resume.
- Reject missing, malformed, non-object, or incompatible resume manifests with
  structured errors. Track the work directory in new bulk/Sarek manifests.
- Validate both emitted parameters and normalized samplesheets in resumed check
  mode. Single-cell resumes now reject changed samplesheet contents.
- Stage Sarek samplesheets until preflight and resume checks pass, preserving the
  original snapshot when a proposed resume is rejected.
- Validate Sarek output-directory failures before creating artifacts.
- Reject nonfinite and negative timeout values. Read only a bounded tail of large
  execution logs and terminate process groups on keyboard interruption.
- Keep nonfinite values in error details representable in standards-compliant JSON.

### Local Docker capacity and replay

- Bound CPU and memory requests to the active Docker daemon as well as the host;
  a minimum-memory floor can no longer exceed the detected machine capacity.
- Configure the local executor's total CPU/memory budget alongside per-process
  `resourceLimits`, preventing concurrent tasks from oversubscribing a smaller VM.
  These safeguards include demo runs.
- Apply launcher capacity only to recognized local Docker profiles. Site profiles,
  cloud test profiles, and object-store work directories retain their configured
  resource policy. Explicit user config files retain final precedence.
- Replay single-cell/Sarek bundles on macOS and Linux with the smaller of recorded
  resource limits and current Docker capacity. Generate an ephemeral replay config
  and remove it on exit; preserve the checksummed original.
- Anchor Sarek replay to the bundle before reading resource files, and distinguish
  automatic resource configs from user configs with similar filenames.
- Write Sarek's environment through the shared reproducibility writer and reuse
  the parsed versions captured in preflight. Raw `nextflow -v` banners must not
  become conda pins or contradict replay/manifest versions.

The resource changes follow Nextflow's [local executor configuration](https://docs.seqera.io/nextflow/executor/local)
and [config precedence](https://docs.seqera.io/nextflow/config). A resource ceiling
limits requests; it cannot make an analysis fit a machine with insufficient RAM.

## Validation

The existing wrapper suite initially passed with **1,349 passed, 2 skipped**.
The package-integration suite also covers the shared modules, launcher routing,
runner/public package APIs, catalog generation, and config/flag forwarding.
Regression tests were added before fixes and exercised through failing and passing
runs. The final focused suite includes Hypothesis parsing/protocol/timeout checks,
mocked launch/failure/resume cases, shell replay tests, official schema comparisons,
and three real Nextflow scheduler tests.

| Python | Final focused result |
| --- | --- |
| 3.11.17 | 1,742 passed, 2 skipped |
| 3.12.15 | 1,742 passed, 2 skipped |
| 3.13.13 | 1,742 passed, 2 skipped |

The two existing skips require optional sibling pipeline checkouts. Both comparisons
(STAR arguments and the protocol matrix) also passed separately against the
downloaded scrnaseq 4.1.0 release on all three Python versions. Set
`CLAWBIO_NFCORE_UPSTREAM_DIR` to enable the source-contract comparison:

```bash
CLAWBIO_NFCORE_UPSTREAM_DIR=/path/to/releases python -m pytest \
  tests/test_nfcore_contract_generation.py tests/test_nfcore_resource_replay.py \
  tests/test_nfcore_wrapper_audit.py tests/test_nfcore_packaging.py \
  clawbio/common/tests/test_nfcore_*.py skills/nfcore-rnaseq-wrapper/tests \
  skills/nfcore-scrnaseq-wrapper/tests skills/nfcore-sarek-wrapper/tests \
  clawbio/tests/test_nfcore_default_output.py \
  clawbio/tests/test_pipeline_config_forwarding.py \
  clawbio/tests/test_pipeline_flag_spelling.py \
  clawbio/tests/test_extra_args_allowlist.py clawbio/tests/test_runner.py \
  clawbio/tests/test_package_api.py clawbio/tests/test_packaging.py \
  tests/test_generate_catalog.py -q
```

Using Nextflow **25.10.4** and Java **21.0.7**, all three generated scheduler configs
serialized a real three-task workflow constrained to one task's memory budget;
execution traces confirmed that the tasks did not overlap. During the original
pre-package audit, all three pinned pipeline sources also started successfully with wrapper-generated parameters and
`-profile test,docker --help`, then completed `-preview` with those same params.
Preview exercised upstream validation and workflow construction while skipping all
process execution. The single-cell preview reported
`validationSchemaIgnoreParams: genomes`, which originates in the pinned upstream
`conf/test.config`; that key is absent from wrapper-generated params.

After package adoption, the bulk RNA and Sarek Nextflow 25.10.4 previews also
completed through the shared command/execution adapters with all processes
skipped. The initial single-cell preview rejected one public FASTQ URL during
upstream file validation; a separate HTTP check of that exact URL returned 200.
This remote-input check is separate from the passing local source-contract tests.

After package adoption, a Nextflow 26.04.3 preview through the shared command and
execution adapters parsed the unmodified Sarek release with the forced v1 parser,
then timed out while downloading the official `nf-amazon@3.9.1` plugin. The registry
returned HTTP 200; a thread dump located the wait in the plugin downloader's socket
read, and a separate transfer also timed out. This is not recorded as a successful
preview or biological workflow run. Release-environment override, timeout, SIGINT,
SIGTERM and process-group cleanup are covered by the passing shared tests.

The post-package repository regression suite uses the locked dependencies with
all optional extras and excludes tests marked `integration`, `network`, or `slow`:
**8,104 passed, 63 skipped, 20 deselected, 2 expected failures, and 1 failure**.
Its only failure is
`skills/proteomics-clock/tests/test_proteomics_clock.py::test_mortality_to_years`,
which also failed on the original revision with the same dependency environment
(expected 5–15, observed 3.81255). The affected file is unchanged. An earlier run
hit DNS errors in unrelated live-API tests and one IGV timing failure; the final
run completed all of those tests successfully. The IGV case also passed alone.
The pre-package baseline at `ab0fd04` had 7,934 passed and the same proteomics
failure. No new repository failures were found after the package integration.

The integration wheel was built from the source distribution. It contains the
shared package and all three scalar contracts, with no competing monolithic module.
All shared module bytes match the tested checkout. All three installed package
remappers run with `python -I -S`, without site packages or `PYTHONPATH`.
Installed outside the checkout with the locked dependency environment available,
all three console/direct help commands and valid parameter emission/invalid-enum
rejection checks passed. With Docker stopped, all three installed preflights
returned structured `DOCKER_NOT_RUNNING` errors.

All three skill specifications passed `agentskills validate`. Contract-generation
checks, targeted unused-name lint, Python compilation, and `git diff --check` passed.

## Complete container validation (2026-10-15)

Full public-test workflows were subsequently executed on an isolated Linux x86-64
validation checkout with Nextflow **26.04.3**, Java **17.0.20.1**, Python **3.12.3**,
and Docker. The server exposes 28 CPUs and 62 GiB RAM. An explicit final-precedence
config requests an 8-CPU/20-GiB local executor budget and per-task resource limits;
at most two heavy validation workflows run at once. Other users' workloads,
containers, images, and output directories were left untouched.

| Pipeline | Original full workflow | Fresh bundle-only workflow |
| --- | --- | --- |
| rnaseq 3.26.0 | 234 tasks succeeded | 234 tasks succeeded, 25m 20s |
| sarek 3.8.1 | 9 tasks cached + 14 succeeded on recovery | 23 tasks succeeded |
| scrnaseq 4.1.0 | 19 tasks succeeded, 17m 25s | 19 tasks succeeded, 19m |

The original launches used `d11bbe1`; the launcher fixes were exercised through
`645df2f`, and the latest shared-package integration was validated at `d722647`.
These source checkouts remained frozen while their runs were active. Official
pipeline assets were used without editing their cached source. Only public nf-core
test data were downloaded. The initial Sarek run encountered a Docker parent-snapshot
error. Pulling its exact pinned GATK image again restored execution; the wrapper's
prepared manifest allowed recovery with `--resume` and reuse of nine completed tasks.

For every wrapper, validation included a real resumed workflow, a refused changed-work-directory
resume/check with the completed result and bundle preserved byte for byte,
in-place `commands.sh` execution from an unrelated working directory, a fresh
output containing only a copied `reproducibility/` directory, and execution after
renaming the complete original output so its former path no longer existed.
No work directory or Nextflow task database was copied into the fresh bundle-only
output. The official Nextflow asset/plugin and Docker image caches remained available.
Bulk replay requires its documented ClawBio installation/checkout dependency.
Standalone remappers were executed with `python -I -S`.

All three moved runs subsequently passed real resumes through the final integrated
checkout, with unchanged checked scientific outputs and valid current artifact
paths in `result.json`. Changed work-directory requests were refused without
changing stored snapshots. Sarek also refused a changed `--trace-report-suffix`,
confirming that the corrected root launcher forwards it. Every final checksum
and standalone path-verification check passed.

Scientific checks compare feature/sample identities as well as values:

- Sarek: both Strelka VCF record sets, including sample headers and eight reported
  variants, plus both mosdepth coverage-summary tables matched after every successful replay.
  Non-scientific VCF `##` execution metadata were excluded from comparison.
- Single-cell: four raw/filtered MatrixMarket outputs matched after canonicalizing
  gene and barcode identities. All six H5AD files were read with the upstream Scanpy
  container and matched by cell identity, gene identity, and counts, including both
  combined matrices. The combined filtered object contains 8,672 cells, 1,100 genes,
  13,847 nonzero entries, and 15,919 counts.
- Bulk RNA: 22 quantification/count/TPM tables were compared, including the preferred
  length-scaled gene counts. Feature IDs, sample IDs, and nonnumeric annotations
  matched. All four gene/transcript SummarizedExperiment RDS files loaded successfully
  and contained finite, nonnegative counts for five samples and 125 features.
  Fresh and relocated runs show native Salmon estimation variation; bitwise
  scientific equality is therefore not claimed for recomputed RNA quantification.

For bulk RNA, the fresh replay's preferred STAR/Salmon length-scaled counts differ
by at most **0.062291 counts**, with maximum per-sample L1-relative difference
**3.470e-6**. The relocated run differs by at most **0.045636 counts** and
**1.934e-6**. Raw STAR/Salmon gene counts differ by at most **0.007** and **0.002**
respectively. Pseudoalignment TPMs differ by at most **132.806438** in the fresh run.

This variation was investigated with **eight direct native controls**: two repeats
of each mapping/alignment Salmon task for `RAP1_IAA_30M_REP1` and `WT_REP2`, using
unchanged upstream task scripts, identical staged inputs, Salmon **1.10.3**, and the
same six threads. These controls also vary. The worst-case WT_REP2 mapping control
reproduces a **137.955965 TPM** difference, exceeding the fresh replay's observed
132.806438. Alignment controls also reproduce fractional-read variation. No arbitrary
numeric tolerance was used to silently pass differing scientific outputs. Cached
RNA resumes and in-place replay retain identical checked quantification tables.

### Snapshot integrity and upstream path constraints

`checksums.sha256` records a particular output snapshot. Direct Nextflow replay can
regenerate MultiQC reports, plot metadata, and execution reports. `--verify` checks
required files and input/reference paths; it is not a content-digest check. The
validator recorded checksum differences before repair, preserved each previous
manifest outside the bundle, deliberately made the checksum manifest missing, and
then exercised `--repair-bundle`. A second repair was a no-op, and both the resulting
`sha256sum -c` check and standalone path verification passed. Existing complete
snapshots are deliberately never silently rehashed by `--repair-bundle`.

Shell/YAML quoting tests do not guarantee that every upstream program accepts every
filesystem name. Nextflow 26.04.3 expands an unescaped literal `$` in its generated
work-directory task script; a deliberate stress test failed with `literal: unbound
variable`. Real relocated RNA and single-cell runs succeeded with spaces and an
apostrophe after removing the dollar sign. Sarek's GATK BaseRecalibrator repeatedly
failed in that quoted work-directory spelling with `AnnotatedIntervalCodec` errors;
ordinary-path relocation succeeded. A separate full run through the corrected
root launcher also succeeded with a quoted output directory and a plain external
`--work-dir`; its VCFs and coverage summaries matched the baseline, its manifest
recorded the requested work directory, and its checksums passed. Use a plain scratch
directory for Sarek, for example `--work-dir /analysis/sarek-work`, and retain the
failure logs when diagnosing an upstream tool. No cached Nextflow or nf-core source was patched to bypass these
constraints.

The earlier macOS VM had insufficient memory/storage for complete workflows; those
attempts remain resource-limit observations. The successful server runs above
supersede the earlier absence of complete biological pipeline validation. Public
test workflows validate these pinned releases and tested profiles, not every
production reference, optional caller, site scheduler, or future upstream release.
