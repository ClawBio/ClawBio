# nfcore-rnaseq-wrapper — Changelog

All notable changes to the `nfcore-rnaseq-wrapper` ClawBio skill are documented
here. The format roughly follows [Keep a Changelog](https://keepachangelog.com)
and the wrapper version is tracked in `SKILL.md` YAML frontmatter.

## [Unreleased] — 0.1.0

### Fixed

- **`remap_paths.py`: combining actions no longer discards one of them silently.**
  `main()` dispatched through a chain of early returns, so only the first matching flag
  ever ran. `--output-dir <new> --verify` therefore rewrote nothing and then reported the
  bundle "ready to replay" while `commands.sh` still pointed at the **old** output
  directory — a confident green light on a stale bundle, which is worse than an error.
  (The scrnaseq bundle had the mirror image: `--output-dir` short-circuited first and
  `--verify` was dropped, so the requested readiness check never ran yet still exited 0.)

  Requested actions now compose and all of them run, in the order the module docstring
  already documented: `--repair-bundle`, then the remaps (`--old/--new`,
  `--refs-old/--refs-new`), then `--output-dir`, and `--verify` **last** — so `--verify`
  always reports on the bundle as it now stands, never on a pre-rewrite snapshot. The
  first failing step sets the exit code, but later steps still run so a single invocation
  surfaces every problem.

  Half a prefix pair (`--old` without `--new`, or `--refs-old` without `--refs-new`) is
  now a hard error instead of falling through to the help text and exiting non-zero with
  no explanation — the same silent-drop class of bug.

  `main()` also takes `argv`/`bundle_dir` so this dispatch is directly testable; regression
  tests pin that `--output-dir … --verify` both rewrites and verifies.
- **`clawbio.py run <pipeline> --help` now shows the wrapper's own help for every nf-core
  pipeline, not just sarek.** The launcher declares only a subset of each wrapper's flags as
  argparse options and forwards the rest through a per-skill allowlist, so rendering the
  launcher's own help hides flags that genuinely work. `--allow-remote-inputs` was the case
  in point: supported by all three wrappers, but discoverable via `--help` only for
  `sarek-pipeline`, because help delegation was hardcoded to that one skill. It now applies
  to every entry in `_NFCORE_PIPELINE_SKILLS`, which is what the delegation was introduced
  to do -- keep the launcher's help from drifting away from the wrapper parser.

### Fixed

- **A `--demo` run's reproducibility bundle is now replayable.** After a successful
  `--demo` run, `bash reproducibility/commands.sh` re-invoked the wrapper against the
  now-populated output directory and died with `OUTPUT_DIR_NOT_EMPTY`. The error's own
  suggested fix (`--resume`) could not work either, because `_apply_demo_overrides`
  force-cleared `--resume` for every demo run — leaving the bundle in a dead end.

  Root cause was that force-clear, which had no basis in the upstream documentation:
  Nextflow's `-resume` is orthogonal to `-profile test` and nf-core documents no
  incompatibility between them. Everything a resume needs already persisted inside the
  demo run's own output directory — the work tree (`upstream/work`) and the Nextflow
  session cache (`.nextflow/`) — and the demo samplesheet stub is content-stable, so its
  checksum matches on replay.

  `--demo` no longer touches `--resume`, and the idempotent-replay guard already emitted
  into `commands.sh` for real runs is now emitted for demo bundles too: it adds
  `--resume` when the target output directory already holds a completed run of this
  bundle (`reproducibility/manifest.json` present), and stays out of the way for a fresh
  or `remap_paths.py --output-dir`-relocated directory.

### Added

- **`demo` is now part of the resume compatibility contract.** The reproducibility
  manifest records `demo`, and preflight compares it like `aligner`/`profile`/`arm`.
  `--demo` composes the upstream `test` profile, which brings *both* its own samplesheet
  and its own bundled references, so resuming across the demo/real boundary would swap
  both underneath the run. This previously surfaced (if at all) as an opaque
  `params_checksum` mismatch; it is now rejected precisely, naming the `demo` field.
  Manifests written before this key existed are treated as `demo: false`, so older
  output directories still resume.

### Documentation

- **Bundle replay documented as a Gotcha.** SKILL.md now states that an in-place replay
  is idempotent (including for `--demo`), explains why the rnaseq bundle needs the guard
  when the sarek/scrnaseq bundles do not (rnaseq's `commands.sh` re-invokes the wrapper;
  theirs invoke Nextflow directly, which tolerates a populated output dir), and warns
  against "fixing" a replay by deleting the output directory.
- **`remap_paths.py --output-dir` documented as a Gotcha.** SKILL.md now explains that
  relocating the bundle's output directory uses `python3 reproducibility/remap_paths.py
  --output-dir <new-path>` (the rnaseq bundle bakes `--output` into `commands.sh`),
  alongside `--old/--new` and `--refs-old/--refs-new`, and notes the scrnaseq bundle
  self-relocates and accepts `--output-dir` only for parity. Resolves the cross-wrapper
  CLI inconsistency being undocumented.
- **`--allow-remote-inputs` semantics clarified.** SKILL.md now states explicitly
  that the flag relaxes only the wrapper's own local-first preflight check: remote
  FASTQ/reference URIs are written into the normalized samplesheet/`params.yaml`
  verbatim and staged natively by Nextflow at run time. The wrapper does not
  download them, so remote inputs require outbound network access and are
  incompatible with `NXF_OFFLINE` (under which Nextflow's own nf-schema file-existence
  validation still runs and fails on the remote paths). Shared wording across the
  three wrappers.

### Fixed

- **Output-dir "not empty" check now ignores the whole `reproducibility/` bundle.** The
  check used a per-file allowlist (`_ALLOWED_REPRO_FILES` = samplesheet/params/manifest
  only), so an output directory containing a *complete* reproducibility bundle — with the
  wrapper's own `commands.sh`, `checksums.sha256`, `environment.yml`, `remap_paths.py` and
  provenance JSON — was rejected with `OUTPUT_DIR_NOT_EMPTY` on a non-resume run. The
  `reproducibility/` directory is entirely wrapper-generated (never user data), so it is
  now ignored wholesale (added to `_IGNORED_ROOT_NAMES`), matching nfcore-sarek and
  nfcore-scrnaseq. This removes a fragile allowlist that had to grow with every new bundle
  file, and makes the three wrappers accept the same pre-existing output-dir contents.
  Genuine result artifacts at the output root (`report.md`, `result.json`, `upstream/`,
  `logs/`) still block a non-resume re-run, and the incomplete-prior-run guidance is
  unchanged.
- **Host memory auto-cap now also applies on Linux (not only macOS).** The
  `process.resourceLimits` config that stops Nextflow's local executor from aborting a
  real run with `Process requirement exceeds available memory` (when an nf-core default
  process request is larger than the host) was gated behind `platform.system() ==
  "Darwin"`, so a normal Linux docker run got no cap and failed on any host smaller than
  the pipeline's production requirements. A docker run on a non-macOS host now writes a
  portable `resourceLimits` config scaled to the machine — the smaller of physical RAM
  and the Docker runtime `MemTotal`, minus a few GB of headroom, with no 15 GB macOS-VM
  ceiling — per nf-core's resourceLimits guidance (values match the machine maximum).
  `--demo` is untouched (it relies on `-profile test`'s own limits) and the macOS path
  is byte-for-byte unchanged. The config carries only the resourceLimits block, none of
  the macOS-only workarounds (`--platform`, `stageInMode = 'copy'`).
- **IPv6/NAT64 hint now actually shows (scan `.nextflow.log`).** The
  `EXECUTION_FAILED` environment-hint scanner only read `logs/stdout.txt` and
  `logs/stderr.txt`. When Nextflow fails while parsing the config (before the pipeline
  starts), the top-level stderr often shows only "Unable to parse config file" while the
  real cause ("Network is unreachable") is recorded only in `.nextflow.log` — so the
  already-implemented IPv6/NAT64 hint (`NXF_OPTS='-Djava.net.preferIPv6Addresses=true'`)
  never appeared on the hosts that need it. `.nextflow.log` (in the Nextflow launch cwd)
  is now scanned too. Applied identically to all three wrappers (the hint function is
  shared verbatim). Covered by a new test.
- **`commands.sh` is regenerated atomically (no self-corruption on in-place replay).**
  An in-place `--resume` replay re-invokes the wrapper, which regenerates the very
  `commands.sh` that bash is still executing. The previous truncate-and-rewrite corrupted
  bash's mid-run read (a trailing shell error such as `…: No such file or directory`).
  `commands.sh` is now composed in memory and written once via an atomic
  `os.replace` (new `write_text_lf_atomic`), so an open reader keeps the original inode
  intact. Verified with an open-file-descriptor test.
- **`commands.sh` locates the ClawBio checkout automatically.** The script walked up
  from its own directory looking for `skills/`, but the bundle always lives OUTSIDE the
  repo (the wrapper forbids `--output` inside it), so it never found the checkout and
  required a manual `CLAWBIO_REPO` — contradicting the header. It now bakes the
  generating checkout as the `CLAWBIO_REPO` default (`REPO_ROOT="${CLAWBIO_REPO:-…}"`),
  so a same-machine replay needs no manual setup while remaining overridable, and the
  misleading "from anywhere inside the repository clone" header wording is corrected.
  Sarek/scRNA-seq don't need this (they replay Nextflow directly). Covered by new tests.
- **Config re-bundling on replay is idempotent (no `config_01_config_01_…` growth).**
  A `--nextflow-config` file is staged into `reproducibility/nextflow_configs/` as
  `config_NN_<name>`. On an in-place `--resume` replay, `commands.sh` re-invokes the
  wrapper with `--nextflow-config` already pointing at that staged copy
  (`${SCRIPT_DIR}/nextflow_configs/config_01_<name>`), and the wrapper copied it again
  under a fresh prefix (`config_01_config_01_<name>`), accumulating a prefix on every
  replay. Staging now detects that the source already lives in this bundle's
  `nextflow_configs/` directory and references it in place instead of re-copying.
  Sarek/scRNA-seq are unaffected (they replay Nextflow directly and do not re-stage
  configs). Covered by a new test.
- **`commands.sh` replay uses a portable interpreter at the source.** The shared
  `clawbio/common/portable_commands` template that builds `commands.sh` emitted a bare
  `python "$SKILL_SCRIPT"`; this wrapper previously rewrote it to `${PYTHON:-python3}`
  with a post-generation patch. The interpreter is now `"${PYTHON:-python3}"` directly in
  the template, so every skill that builds a bundle inherits the fix, and this wrapper's
  now-redundant post-generation patch has been removed. Behaviour is unchanged; the
  existing `commands.sh` interpreter test still passes.
- **`commands.sh` in-place replay is now idempotent.** Because this wrapper's
  `commands.sh` re-invokes the wrapper (which re-runs preflight, unlike the Sarek/scRNA-seq
  bundles that replay Nextflow directly), a plain re-run in the original `--output`
  directory failed with `OUTPUT_DIR_NOT_EMPTY`. `commands.sh` now guards the replay:
  when the target output directory already holds a completed run of this bundle
  (`reproducibility/manifest.json` present) it adds `--resume`; a fresh or
  `remap_paths.py --output-dir`-relocated output directory has no manifest and runs
  normally. The guard is omitted for `--demo` and for runs that already baked `--resume`.
  `remap_paths.py --output-dir` rewrites the guard's manifest path in lockstep with the
  `--output` flag. Covered by new tests.
- **Reproducibility helper instructions now invoke `python3`, matching the sibling
  wrappers.** `remap_paths.py` ships a `#!/usr/bin/env python3` shebang, but its own
  usage/help text and error hints — and the "Portability notice" header emitted into
  `commands.sh` by `reporting.py` — told users to run it as a bare `python
  remap_paths.py`. On modern macOS and many Linux distributions only `python3` exists
  in `PATH` (PEP 394), so the suggested command failed with `python: command not
  found`. All such instructions now say `python3`, consistent with the
  nfcore-sarek-wrapper and nfcore-scrnaseq-wrapper reproduction guides. A static guard
  test (`tests/test_portability_remap_interpreter.py`, mirroring the sibling wrappers)
  prevents regressions.
- **Generated `rnaseq_de_handoff.sh` now uses a portable interpreter.** The downstream
  handoff script — executed on a possibly-fresh machine via `bash rnaseq_de_handoff.sh`
  — invoked a bare `python "${CLAWBIO_REPO}/clawbio.py"`, which fails on python3-only
  systems. It now uses `"${PYTHON:-python3}"`, mirroring the `commands.sh` replay patch
  (defaults to `python3`, honours a `PYTHON` override). Covered by a new test.
- **`--allow-remote-inputs` is now replayed by `commands.sh`.** A run launched with
  `--allow-remote-inputs` (remote FASTQ/reference URIs) produced a bundle whose
  `commands.sh` omitted the flag, so re-running it failed preflight with
  `REMOTE_INPUT_NOT_ALLOWED` even though the bundle's samplesheet/`params.yaml` already
  carried the remote URIs. Because rnaseq's `commands.sh` re-invokes the wrapper (and
  thus re-runs preflight, unlike sarek/scrnaseq which replay Nextflow directly), the
  flag is now recorded in `_BOOLEAN_FLAGS` and emitted whenever the user opted in, so
  the replay faithfully reproduces the original run.
- **User `-c` configs are now copied into the bundle and replayed portably.** A
  `--nextflow-config`/`-c` file living outside the output directory was baked into
  `commands.sh` as a host-specific absolute path, so the reproduction script failed
  out-of-the-box on another machine (or after the config moved). Each config is now
  copied into `reproducibility/nextflow_configs/config_NN_<name>` and referenced via
  `${SCRIPT_DIR}`, matching the nfcore-scrnaseq-wrapper bundle. Remote config URIs are
  passed through verbatim.
- **Reproducibility bundle is now byte-stable (LF-only) on every OS.** Every bundle
  artifact the wrapper writes itself — `report.md`, `params.yaml`, `manifest.json`
  and the provenance JSONs, `check_result.json`, the error-path `result.json`,
  `commands.sh` (the repo-root/`python3` patches and the appended portability
  notice), the demo/no-input samplesheet stubs, the macOS Docker config, the
  `rnaseq_de_handoff.sh` template, and the standalone `remap_paths.py` rewrites and
  repair-bundle stubs — now routes through the shared `write_text_lf` choke-point
  (via a self-contained `_write_text_lf` in the bundle-shipped `remap_paths.py`,
  which cannot import ClawBio at replay). Previously these used raw
  `Path.write_text`, which emits CRLF on Windows and would corrupt `commands.sh`
  under bash and change every `checksums.sha256` digest. This brings the wrapper
  into line with the nfcore-scrnaseq and nfcore-sarek wrappers — whose bundle
  writers already routed through the same choke-point — and with the single-writer
  contract documented in `clawbio/common/textio.py`. Added
  `tests/test_bundle_portability.py` as the cross-OS regression guard (parity with
  the scrnaseq wrapper).
- **Detected Java/Nextflow versions preserve zero-padded components in reports.**
  The version string shown in `report.md` and `result.json` was reconstructed from
  the integer comparison tuple, so `int("04") = 4` turned `26.04.3` into `26.4.3` —
  not a real Nextflow release, and inconsistent with `manifest.json`, which kept the
  correct string. `_check_java`/`_check_nextflow` now record the version exactly as
  reported by the tool (via `_detected_version_string`) for display, using the
  integer tuple only for the minimum-version comparison. Parity with
  nfcore-scrnaseq-wrapper (which already did this) and nfcore-sarek-wrapper.
- **Config-parse failures now point at `NXF_OFFLINE`.** On Nextflow 26.x the
  nf-core `nextflow.config`'s `includeConfig ... ? <url> : '/dev/null'` line fails
  to parse when the remote `nfcore_custom.config` cannot be fetched. The executor's
  `EXECUTION_FAILED` fix now detects `Unable to parse config file` /
  `ConfigParseException` and suggests `NXF_OFFLINE=true` for a fully local run (or
  confirming outbound HTTPS/DNS). Shared verbatim across the three wrappers.

- **`EXPECTED_OUTPUTS_NOT_FOUND` now explains an all-skipped run.** nf-core/rnaseq
  drops samples below `min_trimmed_reads` (schema default 10000) before
  quantification, so a run with only tiny inputs completes "successfully" yet
  produces no merged count matrix — which the wrapper correctly flags. When the
  captured log shows nf-core "completed with skipped sample(s)", the error fix now
  points at `--min-trimmed-reads` (and the MultiQC/pipeline_info skipped-sample
  report) instead of leaving the cause unexplained. Diagnosed from the actual log
  signal, so it only appears when relevant.
- **`result.json` carries the shared `ok`/`status` contract.** A successful run's
  envelope now includes `status: "ok"` and `ok: true`, and a failed run includes
  `status: "error"` alongside the existing `ok: false` — a minimal discriminator
  shared with nfcore-sarek/scrnaseq. Implemented via opt-in `status`/`ok`
  parameters on the shared `clawbio.common.report.write_result_json` helper, so
  the ~40 other skills that use it are unaffected (the keys appear only when a
  caller passes them).

- **`clawbio run rnaseq-pipeline` now forwards `-c`/`--config` Nextflow config
  files.** The launcher (`clawbio/cli.py`) accepted `-c`/`--config` but forwarded
  only `--nextflow-config`, silently dropping configs supplied with the short or
  `--config` spelling. The three spellings are now a single repeatable option and
  every entry is forwarded as `--nextflow-config` (which the wrapper accepts as an
  alias), so config files reach Nextflow regardless of spelling.
- **`--timeout-hours` and `--allow-pipeline-version-override` are now forwardable.**
  Both are real wrapper flags (and were already exposed for the sibling pipelines)
  but were absent from the launcher allowlist, so the extra-args filter dropped
  them before they reached the wrapper. `--allow-pipeline-version-override` is also
  recorded as value-free so the filter cannot consume the following token.
- **Output directory inside the repository now raises the dedicated
  `OUTPUT_DIR_INSIDE_REPO` code** instead of the misleading
  `OUTPUT_DIR_NOT_WRITABLE` (which implied a permissions problem), matching
  nfcore-sarek and nfcore-scrnaseq.
- **`--demo` under `NXF_OFFLINE` now fails fast with a clear message.** Demo mode
  runs nf-core's upstream `-profile test`, whose FASTQs and references are remote
  GitHub URLs; on an offline/sandboxed host the nf-schema plugin previously aborted
  with a cryptic `does not exist`. Preflight now detects `NXF_OFFLINE` + demo and
  raises `DEMO_REQUIRES_NETWORK` with an actionable fix. Docs (SKILL.md, AGENTS.md)
  clarify that `--demo` downloads only nf-core public test data — no user/genetic
  data is uploaded — so it is compatible with the local-first guarantee.
- **macOS + Docker `/tmp` guard is now accurate and robust (demo-mode parity).**
  The preflight warning previously said output under `/tmp` "may be slow or
  unreliable due to VirtioFS behavior", which mis-describes the actual failure:
  Colima does not share `/tmp` into its VM, so a `/tmp` work-dir hard-fails with
  `.command.run: No such file or directory`. The warning now states the real
  cause and fix (move `--output` under HOME), uses a resolve-based
  `is_under_tmp` check (shared in `schemas.py`) instead of a brittle string
  prefix, and the executor appends the same actionable hint to `EXECUTION_FAILED`
  — matching nfcore-sarek and nfcore-scrnaseq.
- **nf-core-native snake_case flag spellings are now accepted via the launcher.**
  nf-core parameters are snake_case (`--gene_bed`, `--transcript_fasta`, …) while
  the wrapper exposes them as hyphenated flags. A user copying an upstream nf-core
  command previously had the token silently dropped by the launcher's INT-001
  allowlist filter, which matches exact (hyphenated) tokens. `clawbio run
  rnaseq-pipeline` now canonicalises `_`↔`-` when matching the allowlist and
  forwards the wrapper's canonical hyphen spelling (which its parser already
  accepts). Scoped to the three nf-core pipeline skills; other skills keep
  exact-match filtering.

### Added

- **Environment post-failure hints on `EXECUTION_FAILED`.** When a run fails, the
  executor scans the captured Nextflow logs and appends an actionable hint if it
  finds a known environment signature — diagnosed from the actual error text, so no
  resource thresholds are invented. `Process requirement exceeds available memory`
  (nf-core's default request larger than the host, e.g. `MAKE_TRANSCRIPTS_FASTA`)
  yields a hint to cap resources via a `-c` config using `process.resourceLimits`;
  `Network is unreachable` / a Java connection exception (common on IPv6-only /
  NAT64 hosts, where the JVM prefers IPv4) yields a hint to verify outbound
  DNS/HTTPS and to set `NXF_OPTS='-Djava.net.preferIPv6Addresses=true'`. Shared
  verbatim across the three wrappers.
- **`--allow-remote-inputs` opt-in (local-first by default).** Remote samplesheet
  inputs and reference paths (`s3://`, `gs://`, `https://`, `ftp://`, …) are now
  rejected at preflight (`REMOTE_INPUT_NOT_ALLOWED`) unless the flag is passed, in
  which case a runtime warning names every path fetched over the network. The
  object-store `--work-dir` is not gated. Shared verbatim with
  `nfcore-scrnaseq/sarek`.
- **Control-flag parity with the sibling wrappers.** Added `--work-dir` (Nextflow
  work directory override; accepts a local path or an object-store URI for cloud
  executors; was hardcoded to `<output>/upstream/work`) and `-c`/`--config` as
  aliases of `--nextflow-config`. The executor/command-builder now accept a
  `Path | str` work dir so remote URIs pass through verbatim.

### Changed

- **`--timeout-hours 0` now disables the wall-clock cap** (was rejected at
  preflight). The default still applies a 12h cap; `0` is an explicit opt-out for
  long HPC/cloud runs — parity with `nfcore-scrnaseq/sarek`. `_resolve_timeout_seconds`
  returns `int | None`, the executor accepts `None`, negative values are still
  rejected, and the macOS compatibility config falls back to the default ceiling
  when the cap is disabled. A `[provenance]` stage log was added so all three
  wrappers share the identical stage-prefix set.

### Added

- **Order-independent test suite (shared mechanism).** `tests/conftest.py` now
  also carries the canonical-object bare-module isolation block shared verbatim
  with the sibling wrappers (in addition to this skill's existing per-file
  guards), so the cross-skill isolation mechanism is identical across all three.
- **Converged CLI presentation with the sibling wrappers.** Added a startup
  banner (`--no-banner`), `-v/--verbose`, and sarek-style staged progress logs
  (`[preflight]`/`[execute]`/`[outputs]`/`[report]`/`[done]`) plus a
  human-readable boxed error on stdout — while keeping the machine-readable JSON
  error on stderr and `result.json` on disk. `--verbose`/`--no-banner` are
  registered in the `clawbio.py` runner allowlist and the replay-drift guards.
- **Robust `main()` entrypoint.** `main()` now returns exit code `130` on
  `KeyboardInterrupt` (the SIGINT convention) instead of letting the interrupt
  propagate as a traceback. Matches `nfcore-sarek-wrapper` and
  `nfcore-scrnaseq-wrapper`.

### Notes

- This wrapper is the reference implementation for the cross-wrapper
  **input-readability policy**: the launcher never pre-checks FASTQ/BAM
  readability (`os.access(R_OK)`) because Nextflow reads the data in the true
  execution context (often a root container under the default Docker profile),
  where a launcher-side probe would false-block valid runs. Existence is
  validated (`MISSING_FASTQ`); readability is deferred to Nextflow staging. The
  rationale is documented inline in `errors.py`.
