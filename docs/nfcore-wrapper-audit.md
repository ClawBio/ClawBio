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
at `dc7951ef71fc8cd4dd64fb1ec3158f1180ce4914`. Its shared package is the
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
snapshots until complete emitted-parameter validation succeeds. Both refused
resumes and rejected fresh plans preserve unaccepted input snapshots. The bulk
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
| 3.11.17 | 1,730 passed, 2 skipped |
| 3.12.15 | 1,730 passed, 2 skipped |
| 3.13.13 | 1,730 passed, 2 skipped |

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

## End-to-end limits

The available Docker VM exposed **2 CPUs and approximately 1.91 GiB RAM**. A full
rnaseq demo attempted during the audit reached upstream tasks but repeatedly
returned exit 137 in BBMap/BBSplit, consistent with memory pressure, and was stopped.
The scheduler fixes were subsequently validated with the real-engine tests above.
A subsequent run used 4 CPUs, approximately 7.74 GiB Docker RAM, and Rosetta.
Deliberate interruption preserved its launch manifest; resumed check mode passed,
and the resumed run executed BBSplit successfully. It later failed with
`No space left on device` on the host. Temporary work files were removed and the
original Docker CPU/memory/emulation settings and stopped state were restored.
No complete biological pipeline run is claimed for any of the three wrappers.
Full container-based demos and scientific-output parity still require a host with
adequate RAM, storage, network access, and the applicable execution backend.
