# nf-core wrapper compatibility audit

This audit covers the three implemented Nextflow wrappers: `nfcore-rnaseq-wrapper`,
`nfcore-scrnaseq-wrapper`, and `nfcore-sarek-wrapper`. It compares wrapper behavior
with the exact upstream releases below, rather than with a moving default branch.
The starting ClawBio revision was `d67bd59af4c8c878e01a5f8deb8cd9a22e2e80b9`.

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
- Write Sarek's environment through the shared reproducibility writer.

The resource changes follow Nextflow's [local executor configuration](https://docs.seqera.io/nextflow/executor/local)
and [config precedence](https://docs.seqera.io/nextflow/config). A resource ceiling
limits requests; it cannot make an analysis fit a machine with insufficient RAM.

## Validation

The existing wrapper suite initially passed with **1,349 passed, 2 skipped**.
Regression tests were added before fixes and exercised through failing and passing
runs. The final focused suite includes Hypothesis parsing/protocol/timeout checks,
mocked launch/failure/resume cases, shell replay tests, official schema comparisons,
and three real Nextflow scheduler tests.

| Python | Final focused result |
| --- | --- |
| 3.11.17 | 1,468 passed, 2 skipped |
| 3.12.15 | 1,468 passed, 2 skipped |
| 3.13.13 | 1,468 passed, 2 skipped |

The two existing skips require optional sibling pipeline checkouts. Both comparisons
(STAR arguments and the protocol matrix) also passed separately against the
downloaded scrnaseq 4.1.0 release on all three Python versions. Set
`CLAWBIO_NFCORE_UPSTREAM_DIR` to enable the source-contract comparison:

```bash
CLAWBIO_NFCORE_UPSTREAM_DIR=/path/to/releases python -m pytest \
  tests/test_nfcore_contract_generation.py tests/test_nfcore_resource_replay.py \
  tests/test_nfcore_wrapper_audit.py skills/nfcore-rnaseq-wrapper/tests \
  skills/nfcore-scrnaseq-wrapper/tests skills/nfcore-sarek-wrapper/tests -q
```

Using Nextflow **25.10.4** and Java **21.0.7**, all three generated scheduler configs
serialized a real three-task workflow constrained to one task's memory budget;
execution traces confirmed that the tasks did not overlap. All three pinned
pipeline sources also started successfully with wrapper-generated parameters and
`-profile test,docker --help`, then completed `-preview` with those same params.
Preview exercised upstream validation and workflow construction while skipping all
process execution. The single-cell preview reported
`validationSchemaIgnoreParams: genomes`, which originates in the pinned upstream
`conf/test.config`; that key is absent from wrapper-generated params.

The repository regression suite uses the locked dependencies with all optional
extras and excludes tests marked `integration`, `network`, or `slow`: **7,318 passed,
76 skipped, 15 deselected, 2 expected failures, and 1 failure**. Its one
remaining failure,
`skills/proteomics-clock/tests/test_proteomics_clock.py::test_mortality_to_years`,
was reproduced on the original revision in the same environment (expected 5–15,
observed 3.81255). The affected file is unchanged by this audit. Initial broad runs
exhausted disk space; temporary audit outputs were cleared and dependency setup
was adjusted to reuse the installed environment before the clean run.

All three skill specifications passed `agentskills validate`. Contract-generation
checks, targeted unused-name lint, Python compilation, and `git diff --check` passed.

## End-to-end limits

The available Docker VM exposed **2 CPUs and approximately 1.91 GiB RAM**. A full
rnaseq demo attempted during the audit reached upstream tasks but repeatedly
returned exit 137 in BBMap/BBSplit, consistent with memory pressure, and was stopped.
The scheduler fixes were subsequently validated with the real-engine tests above.
No complete biological pipeline run is claimed for any of the three wrappers.
Full container-based demos and scientific-output parity still require a host with
adequate RAM, storage, network access, and the applicable execution backend.
