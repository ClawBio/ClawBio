# nfcore-sarek-wrapper — Changelog

All notable changes to the `nfcore-sarek-wrapper` ClawBio skill are documented
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

### Documentation

- **The `prov.*` / `validation.*` config keys are documented as pipeline-owned, with a
  regression guard.** A review reported the bundle's `commands.sh` emitting
  unrecognised-config warnings for `prov.enabled` and `prov.formats.bco.file`, and
  attributed them to the wrapper. They are not the wrapper's: nf-core/sarek 3.8.1's own
  `nextflow.config` declares `plugins { … nf-prov@1.2.2, nf-schema@2.6.1 }` and, alongside
  them, the matching `prov { enabled = true; formats { bco { file = "…/pipeline_info/manifest_${params.trace_report_suffix}.bco.json" } } }`
  and `validation { defaultIgnoreParams = [...] }` scopes. Those scopes are owned by the
  plugins, so Nextflow reports them as unrecognised only when a plugin is not active —
  an engine/environment condition, not a bundle defect.

  Verified: the wrapper writes none of these keys (`params.yaml` carries only the audited
  CLI surface) and `commands.sh` runs the pinned pipeline unmodified, so a plain upstream
  `nextflow run nf-core/sarek -r 3.8.1` behaves identically. No wrapper change is warranted.

  Documented as a Gotcha in SKILL.md, and pinned by a new regression test
  (`test_params_never_carry_upstream_plugin_config_scopes`) asserting the wrapper never
  writes `prov` / `validation` / `validationSchemaIgnoreParams` / `monochromeLogs` into
  `params.yaml`. Injecting those to silence a warning would override the pipeline's own
  provenance configuration — the wrong fix.
- **Timestamped `pipeline_info/` reports documented.** Upstream sarek names its execution
  report, timeline, trace, DAG and nf-prov BCO manifest with `${params.trace_report_suffix}`
  (the run's start timestamp), so re-running `commands.sh` in place adds a new set of report
  files next to the previous run's rather than overwriting them. That is upstream nf-core
  behaviour and is deliberate — those files are each execution's audit trail — and it does
  not invalidate `checksums.sha256`, which pins only the bundle's own artefacts. Recorded as
  a Gotcha with an explicit warning not to "clean" it by pinning `trace_report_suffix`, which
  would make a replay silently overwrite the original run's execution evidence.

### Fixed

- **`--demo` no longer discards `--resume`, so an interrupted demo can be resumed.**
  `_apply_demo_overrides` force-cleared `--resume` on every demo run. Because the upstream
  `-profile test` streams its reads and references from remote URLs, a demo is quite likely
  to be interrupted mid-run (a flaky download, a cancelled session) — and that left the user
  in a dead end: the half-populated output directory made preflight raise
  `OUTPUT_DIR_NOT_EMPTY`, whose fix text says to "re-run with `--resume` against a compatible
  manifest", which the demo override then silently threw away. The only escape was to delete
  the directory and recompute from scratch.

  The force-clear had no basis in the upstream documentation: Nextflow's `-resume` is
  orthogonal to `-profile test`, and nf-core documents no incompatibility between them. It is
  now left alone, matching `nfcore-scrnaseq-wrapper`, which never cleared it.

  Demo/real drift stays blocked, with no new field needed: `params["profile"]` records the
  *composed* profile (which carries the `test` token for a demo) and is already compared via
  `_RESUME_TRACKED_FIELDS`.

### Documentation

- **`split_fastq` gotcha: demo-vs-normal difference is a profile effect, not a thread/CPU
  effect.** The existing gotcha now records that a demo-vs-normal variant difference cannot
  come from the host-scaled `resourceLimits` CPU count: nf-core/sarek 3.8.1 pins
  `bwa mem -K 100000000` on every mapping process (`conf/modules/aligner.config`), which
  fixes the per-batch base count and makes alignment bit-identical regardless of thread
  count. The divergence is driven by the `-profile test` overrides already documented
  (`split_fastq = 0`, and the test profile's own `--tools strelka` / reference choices),
  so a run at 28 CPUs reproduces the same calls as the demo's 4. Documentation only.
- **Transient first-replay failures documented (upstream/network, auto-recovered).** A new
  gotcha explains that a demo replay can fail on the *first* attempt because nf-core/sarek's
  own `conf/test.config` uses a trailing-slash `igenomes_base`/`modules_testdata_base_path`
  (`…/test-datasets/modules/data/`) that can join into `//` URLs some CDNs momentarily 404,
  and remote staging occasionally races — upstream behaviour, not a wrapper/bundle defect.
  The generated `commands.sh` already recovers it: it adds `-resume` whenever a prior
  `.nextflow/` session exists, so re-invoking `bash commands.sh` resumes from cached tasks
  and completes without recomputing. Documentation only (the auto-`-resume` mechanism was
  added in an earlier change).
- **`remap_paths.py --output-dir` documented.** The portability section now lists
  `--output-dir <new-path>` (rewrites the baked `--output` in `commands.sh` when a run is
  relocated) alongside `--old/--new` and `--refs-old/--refs-new`, and notes that the
  scrnaseq bundle self-relocates and only accepts `--output-dir` for parity.
- **`split_fastq` gotcha added.** SKILL.md now documents that normal runs use the
  nf-core/sarek default `split_fastq = 50000000`, while `--demo` (`-profile test`)
  disables splitting (`0`). The wrapper faithfully passes the pipeline default rather
  than silently injecting `0`, so a normal run and a `--demo` run of the *same* tiny
  dataset can produce different variant sets (FASTQ split/merge perturbs calls on
  degenerate, test-scale inputs — a known nf-core artifact, not a wrapper bug). Pass
  `--split-fastq 0` to match `-profile test` when reproducing demo/test output; leave the
  default for production cohorts.
- **`--allow-remote-inputs` semantics clarified.** SKILL.md now states explicitly
  that the flag relaxes only the wrapper's own local-first preflight check: remote
  FASTQ/reference URIs are written into the normalized samplesheet/`params.yaml`
  verbatim and staged natively by Nextflow at run time. The wrapper does not
  download them, so remote inputs require outbound network access and are
  incompatible with `NXF_OFFLINE` (under which Nextflow's own nf-schema file-existence
  validation still runs and fails on the remote paths). Shared wording across the
  three wrappers.
- **"Output Structure" corrected.** The paragraph claimed the `reproducibility/`
  bundle held `report.md` / `result.json` and that execution logs lived under
  `reproducibility/logs/` with "no separate top-level `logs/` directory". The code
  writes `report.md`, `result.json`, and `logs/` at the output root (matching the rest
  of the same SKILL.md and the sibling wrappers). The paragraph now describes the real
  layout, and the output tree lists `check_result.json` at the root.

### Fixed

- **`checksums.sha256` no longer leaks the nested Nextflow work tree into the manifest.**
  The checksum generator excluded `work/`, `.nextflow/`, `reproducibility/` and `logs/`
  by matching only the *first* path component (`rel.parts[0]`). But the default work dir
  is `<output>/upstream/work`, so every ephemeral scratch file sat at `upstream/work/**`
  (first component `upstream`) and was hashed into the manifest — the bulk of its lines
  were transient task files that change every run, defeating the manifest's purpose and
  breaking `sha256sum -c` reproducibility. Both generators — `provenance._iter_checksum_paths`
  and the stdlib-only `remap_paths._regenerate_checksums` (kept in lockstep) — now exclude
  a file when **any** ancestor directory is an excluded tree, so the nested work dir is
  caught wherever it lives (including a custom `--work-dir` placed under the output root).
  The manifest now contains only the real `upstream/results/**` outputs. The existing test
  used an unrealistic flat `output/work/` layout (which the old first-component check
  happened to catch), masking the bug; new tests in `test_provenance.py` and
  `test_remap_paths.py` assert exclusion of the real `upstream/work/**` layout.
- **`--genome testdata.nf-core.sarek` is now rejected early without a matching
  `--igenomes-base`.** nf-core/sarek's tiny test genome resolves only under the
  test-datasets mirror (`conf/test.config` sets `igenomes_base` to
  `https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/`), never the
  default `s3://ngi-igenomes` catalogue — yet it was accepted like a normal iGenomes
  key, so a real run without that override failed late during nf-schema reference
  validation with no warning. `_check_genome_known` now special-cases it: accepted under
  `-profile test`/`--demo` or with a custom `--igenomes-base`, otherwise rejected at
  preflight (`INVALID_GENOME`) with a fix naming `--demo` and the exact test-datasets base.
- **The host resourceLimits cap now ships in the reproducibility bundle and replays.**
  Previously the cap was applied to the live run but stripped from `commands.sh`, so a
  from-scratch reproduction on the generating (memory-constrained) host re-aborted with
  `Process requirement exceeds available memory`. The cap is now written to
  `reproducibility/resource_limits.config` and `commands.sh` re-applies it via a
  `uname != Darwin` guard (complementary to the macOS-config guard), so a bundle
  reproduces the run on the machine that made it. The cap only ever clamps requests, so
  a larger replay host simply under-uses it.
- **Host memory auto-cap now also applies on Linux (not only macOS).** The host-scaled
  `process.resourceLimits` cap that prevents Nextflow's local executor from aborting a
  real run with `Process requirement exceeds available memory` was written only on
  macOS+docker (`_write_macos_docker_config` returned `None` off-Darwin), so a normal
  Linux docker run got no cap and failed whenever an nf-core default process request
  exceeded the host. A non-macOS docker run now writes a portable resourceLimits config
  (`.nextflow_resource_limits.config`) scaled to the machine — the smaller of physical
  RAM and the Docker `MemTotal`, minus headroom, no 15 GB macOS-VM ceiling — per
  nf-core's resourceLimits guidance. It is applied to the **live** run only and stripped
  from the portable `commands.sh` replay bundle so a single machine's RAM is never pinned
  into it. `--demo` and the macOS path are unchanged.
- **`commands.sh` replay is now idempotent (`-resume`).** The nextflow-direct replay
  script re-ran the whole pipeline from scratch against an already-completed output dir.
  It now emits a guard that adds `-resume` when a prior Nextflow session (`.nextflow/`)
  exists in the launch dir, so a first run starts fresh but a replay reuses the cache —
  matching the `nfcore-rnaseq-wrapper` bundle. Any `-resume` captured from the original
  run is stripped so the guard is the single source of truth.
- **IPv6/NAT64 hint now actually shows (scan `.nextflow.log`).** The `EXECUTION_FAILED`
  environment-hint scanner read only `logs/stdout.txt` and `logs/stderr.txt`. When
  Nextflow fails while parsing the config (before the pipeline starts), the underlying
  "Network is unreachable" cause is often recorded only in `.nextflow.log`, so the
  IPv6/NAT64 hint (`NXF_OPTS='-Djava.net.preferIPv6Addresses=true'`) never appeared.
  `.nextflow.log` (in the Nextflow launch cwd) is now scanned too. Shared verbatim with
  the sibling wrappers; covered by a new test.
- **Downstream handoff examples now use `python3`.** The generated
  `sarek_downstream_handoff.sh` listed cross-skill example commands (clinical-variant-reporter,
  clinical-trial-finder, omics-target-evidence-mapper, wes-clinical-report-en/es) as bare
  `python skills/<skill>/<script>.py`. On python3-only systems (modern macOS and many Linux
  distributions, PEP 394) that fails with `python: command not found`. They now use `python3`,
  matching the wrapper's other reproduction/handoff scripts. Covered by a new test.
- **User `-c` configs are now copied into the bundle and replayed portably.** A
  `--nextflow-config`/`-c` file living outside the output directory was rewritten to a
  `<EDIT_ME>` placeholder in `commands.sh` (the portable-argv rewriter cannot anchor an
  out-of-tree absolute path), so the reproduction script failed out-of-the-box unless
  the user hand-edited it. Each external config is now copied into
  `reproducibility/nextflow_configs/config_NN_<name>` before `commands.sh` is built, so
  the argv rewriter resolves it to `$SCRIPT_DIR/nextflow_configs/<name>` — no
  `<EDIT_ME>`, no host-specific absolute path. Configs already inside the bundle (e.g.
  `macos_docker.config`), non-absolute values, and remote URIs are left untouched.
  Matches the nfcore-scrnaseq-wrapper bundle.
- **Output root tolerates OS/VCS scratch files (macOS parity).** `_check_output_dir`
  previously rejected any root entry other than `reproducibility/`, so a stray
  `.DS_Store` / `.gitkeep` / `.gitignore` / `Thumbs.db` at `--output` raised
  `OUTPUT_DIR_NOT_EMPTY`. On macOS, Finder writes `.DS_Store` into any visited folder,
  so sarek would spuriously fail on directories the nfcore-rnaseq / nfcore-scrnaseq
  wrappers already accept. Sarek now shares the same `_IGNORED_ROOT_NAMES` tolerance
  set.
- **`check_result.json` moved to the output root.** `--check` wrote its summary to
  `reproducibility/check_result.json`, whereas nfcore-rnaseq / nfcore-scrnaseq write
  `<output>/check_result.json` — the check-mode parallel of `<output>/result.json`,
  which all three already place at the root. Sarek now writes `check_result.json` at
  the output root too; the new `_IGNORED_ROOT_NAMES` entry keeps a subsequent real run
  in the same directory from tripping `OUTPUT_DIR_NOT_EMPTY`, and it joins `report.md`
  / `result.json` in the checksum-manifest exclusion set (a wrapper summary is not a
  pipeline output).
- **Custom `--fasta` without `--genome` now disables iGenomes automatically.**
  nf-core/sarek 3.8.1 defaults `genome = 'GATK.GRCh38'` in `nextflow.config`, so a
  custom-reference run that supplied `--fasta` but left `--genome` unset would still
  load the full GATK.GRCh38 iGenomes configuration and abort while validating its
  ~20 `s3://ngi-igenomes/...` reference paths as local files — a wall of "does not
  exist" errors that never mentioned `--igenomes-ignore`. The wrapper now sets
  `igenomes_ignore=true` in that case (emitting a WARNING recorded in the report),
  matching the nf-core guidance for custom references and the equivalent automatic
  behaviour already in nfcore-scrnaseq-wrapper. An explicit `--genome` alongside
  `--fasta` (documented partial override) and `--demo`/`-profile test` runs are left
  untouched.
- **The BQSR preflight error now explains that baserecalibrator runs by default.**
  When a start step requires `--dbsnp`/`--known_indels` and neither is available,
  the `MISSING_REFERENCE` message now adds: "baserecalibrator runs by default in the
  Sarek mapping workflow, independently of the requested variant-calling --tools."
  A user who only requested `--tools haplotypecaller` previously had no way to know
  why base recalibration — and therefore the known-sites requirement — was active.
- **Detected Java/Nextflow versions preserve zero-padded components in reports.**
  The version string shown in `report.md` and `result.json` was reconstructed from
  the integer comparison tuple, so `int("04") = 4` turned `26.04.3` into `26.4.3` —
  not a real Nextflow release, and inconsistent with `manifest.json`, which kept the
  correct string. The check now records the version exactly as reported by the tool
  (via `_detected_version_string`) for display, using the integer tuple only for the
  minimum-version comparison. Parity with nfcore-scrnaseq-wrapper, which already did
  this.
- **Config-parse failures now point at `NXF_OFFLINE`.** On Nextflow 26.x the
  nf-core `nextflow.config`'s `includeConfig ... ? <url> : '/dev/null'` line fails
  to parse when the remote `nfcore_custom.config` cannot be fetched. The executor's
  `EXECUTION_FAILED` fix now detects `Unable to parse config file` /
  `ConfigParseException` and suggests `NXF_OFFLINE=true` for a fully local run (or
  confirming outbound HTTPS/DNS). Shared verbatim across the three wrappers.

- **Remote MultiQC config/logo URIs are preserved in `params.yaml`.** A remote
  `--multiqc-config`/`--multiqc-logo`/`--multiqc-methods-description` URL was
  resolved as a local path (collapsing the scheme, e.g. `https://host/x` →
  `<cwd>/https:/host/x`); those fields now pass URIs through verbatim via
  `_posix_or_uri`, parity with the genome reference fields (which already did).
  The now-redundant URI-unaware `_posix` helper was removed.
- **Output layout now matches nfcore-rnaseq/scrnaseq.** `report.md` and
  `result.json` are written at the output root (were under `reproducibility/`),
  and execution logs are written to `<output>/logs/` (were under
  `reproducibility/logs/`). A consumer now finds `<output>/result.json` and
  `<output>/report.md` for any of the three pipelines. The `reproducibility/`
  directory keeps the portable replay + provenance bundle (`commands.sh`,
  `params.yaml`, `manifest.json`, `checksums.sha256`, `environment.yml`,
  `remap_paths.py`, provenance JSON). The root `report.md`/`result.json` and the
  `logs/` tree are excluded from `checksums.sha256`. On failure, the error
  `result.json` marker is likewise written at the root.
- **Report `Profile`, `Java`, and `Nextflow` are populated.** They previously
  rendered as `-` because the runtime metadata was not threaded into the report
  generator. Preflight now surfaces the detected Java and Nextflow versions, and
  the composed Nextflow profile string is passed through, so all three appear in
  `report.md` and in `result.json` (`run.profile`/`run.java_version`/
  `run.nextflow_version`) as the sibling wrappers already do.
- **`--demo` reports the effective samples and tools.** The `Samples` count now
  falls back to the samples detected in the outputs when the local samplesheet is
  empty (the upstream `-profile test` supplies samples remotely), and `Tools`
  falls back to the tools observed in the parsed variant-calling/annotation
  outputs when none were requested — so a demo that runs Strelka no longer reports
  `Samples: 0` / `Tools: (none)`. `result.json` records both the requested
  `run.tools` and the effective `run.tools_effective` (+ `run.tools_from_outputs`).
- **`result.json` carries the shared `ok`/`status` contract.** A successful run is
  `ok: true` (existing `status` retained); a failed run is `ok: false` with
  `status: "error"` — a minimal discriminator shared with nfcore-rnaseq/scrnaseq.
- **`clawbio run sarek-pipeline` now forwards `-c`/`--config` Nextflow config
  files.** The launcher (`clawbio/cli.py`) accepted `-c`/`--config` but forwarded
  only `--nextflow-config`, silently dropping configs supplied with the short or
  `--config` spelling. All three spellings are now normalised and forwarded as
  `--nextflow-config` (the wrapper accepts them as aliases).
- **`--allow-remote-inputs` and `--allow-pipeline-version-override` are recorded as
  value-free in the launcher allowlist.** Both are `store_true` wrapper flags but
  were listed only among value-taking flags, so the extra-args filter could consume
  the following token as a spurious value. They are now in the without-values set.
- **`--demo` under `NXF_OFFLINE` now fails fast with a clear message.** Demo mode
  composes nf-core's upstream `-profile test`, whose FASTQs and references are remote
  GitHub URLs; on an offline/sandboxed host the nf-schema plugin previously aborted
  with a cryptic `does not exist`. Preflight now detects `NXF_OFFLINE` + the test
  profile and raises `DEMO_REQUIRES_NETWORK` with an actionable fix. Docs (SKILL.md,
  AGENTS.md) clarify that `--demo` downloads only nf-core public test data — no
  user/genetic data is uploaded — so it is compatible with the local-first guarantee.
- **macOS + Docker `/tmp` post-failure hint (demo-mode parity).** When a run fails
  with `--output` under `/tmp` on macOS, the executor now appends the actionable
  Colima / `.command.run: No such file or directory` hint to `EXECUTION_FAILED`,
  matching nfcore-scrnaseq (the preflight WARNING already existed). The
  `is_under_tmp` check is consolidated into `schemas.py` as the single source of
  truth shared by preflight and the executor.
- **nf-core-native snake_case flag spellings are now accepted.** nf-core's own
  parameters are snake_case (`--skip_tools`, `--fasta_fai`, …) while the wrapper
  exposes them as hyphenated flags (`--skip-tools`, `--fasta-fai`). A user copying
  an upstream nf-core command previously had the token silently dropped: the
  launcher's INT-001 allowlist filter matches exact tokens (all hyphenated), and
  the wrapper parser only registered the hyphen spelling. Now `clawbio run
  sarek-pipeline` canonicalises `_`↔`-` when matching the allowlist and forwards
  the wrapper's canonical hyphen spelling, and the wrapper parser registers the
  nf-core-native `--<param>` spelling as an alias of every schema passthrough flag
  — so `--skip_tools baserecalibrator` reaches the pipeline whether invoked via the
  launcher or the wrapper directly. The canonicalisation is scoped to the three
  nf-core pipeline skills; other skills keep exact-match filtering.

### Added

- **Environment post-failure hints on `EXECUTION_FAILED`.** When a run fails, the
  executor scans the captured Nextflow logs and appends an actionable hint if it
  finds a known environment signature — diagnosed from the actual error text, so no
  resource thresholds are invented. `Process requirement exceeds available memory`
  (nf-core's default request larger than the host) yields a hint to cap resources
  via a `-c` config using `process.resourceLimits`; `Network is unreachable` / a
  Java connection exception (common on IPv6-only / NAT64 hosts, where the JVM
  prefers IPv4) yields a hint to verify outbound DNS/HTTPS and to set
  `NXF_OPTS='-Djava.net.preferIPv6Addresses=true'`. Shared verbatim across the three
  wrappers.
- **`--allow-remote-inputs` opt-in (local-first by default).** Remote samplesheet
  inputs and user-supplied reference paths (`s3://`, `gs://`, `https://`, `ftp://`,
  …) are now rejected at preflight (`REMOTE_INPUT_NOT_ALLOWED`) unless the flag is
  passed, in which case a runtime warning names every path fetched over the
  network. The public iGenomes mirror base (`igenomes_base`) and the object-store
  `--work-dir` are intentionally remote and are not gated. Shared verbatim with
  `nfcore-scrnaseq/rnaseq`. (Flag surface 176 → 177.)
- **Control-flag parity with the sibling wrappers.** Added `--work-dir`
  (Nextflow work directory override, accepts a local path or an object-store URI
  for cloud executors; was hardcoded to `<output>/upstream/work`),
  `--allow-pipeline-version-override` plus a `_check_pipeline_version_supported`
  guard (a non-pinned `--pipeline-version` is now blocked unless explicitly
  overridden, so the 3.8.1-pinned validations can't silently run against a
  different version), and `-c`/`--config` as aliases of `--nextflow-config`.
  These bumped the flag surface to 176 (154 passthrough + 22 wrapper controls).
- **`--timeout-hours` is now configurable** (was hardcoded to 24h). Accepts a
  positive number of hours or `0` to disable the wall-clock cap for long
  HPC/cloud runs whose walltime the scheduler enforces — parity with
  `nfcore-scrnaseq/rnaseq`. The executor's `timeout_seconds` is now `int | None`.
  This bumped the flag surface 173 → 174 (154 passthrough + 20 wrapper controls).

### Changed

- **Demo coercions now emit `WARNING:` on stderr** (were `[demo]` on stdout):
  overriding a user-supplied flag is an advisory the user should notice, so it
  belongs on stderr — parity with the sibling wrappers.
- **Logging fully consistent with the siblings**: `_print` docstring unified and
  a dedicated `[provenance]` stage line added, so all three wrappers emit the
  identical stage-prefix set (`[preflight]`/`[execute]`/`[outputs]`/`[report]`/
  `[provenance]`/`[done]`/`[check]`/`[abort]`).
- **`SkillError` now exits `1`** (was `2`), matching `nfcore-scrnaseq/rnaseq`.
  Exit `2` is reserved by argparse for CLI-usage errors, so reusing it for
  validation/preflight failures made the two indistinguishable to machine
  consumers. Internal/unexpected errors still exit `1`.
- **Order-independent test suite.** `tests/conftest.py` now claims this skill's
  bare modules (`errors`, `schemas`, …) via a canonical-object cache — at
  collection and before each test — so `make test` (which collects all three
  wrapper suites together) no longer hits cross-skill module shadowing. The same
  block is shared verbatim across all three wrappers' conftests.

### Added

- **Cross-wrapper homogenization pass.** Centralized the import-isolation guard
  into `_isolated_imports.py` (replacing per-module copy-pasted
  `_purge_foreign_modules` and a buggy `if not in sys.path` path block that did
  not force this skill's dir to the front), renamed the inner runner `_run` →
  `_run_wrapper`, removed the dead `FASTQ_NOT_READABLE` error code (with the
  shared readability-policy NOTE), gave the README the common Scope/Out-of-Scope/
  Quick-Start skeleton, and added a machine-readable JSON error on **stderr**
  (in addition to the human box and `result.json`) so a structured error is
  available even when `result.json` cannot be written — matching
  `nfcore-scrnaseq-wrapper` and `nfcore-rnaseq-wrapper`.
- **Cross-OS reproducibility-bundle hardening.** The entire bundle is now
  byte-stable and portable across Linux, macOS, and Windows/WSL:
  - Every bundle text file (`commands.sh`, `report.md`, `result.json`,
    `params.yaml`, `checksums.sha256`, `manifest.json`, `environment.yml`, the
    JSON handoffs, the `.sh` handoff, the normalized samplesheet, and the macOS
    config) is now written with forced **LF** line endings via the shared
    `clawbio.common.textio.write_text_lf` (and a self-contained `_write_text_lf`
    in the standalone `remap_paths.py`). Previously they used platform-default
    newlines, so on Windows the bash `commands.sh` got a CRLF shebang
    (`/usr/bin/env bash\r` → bad interpreter) and `checksums.sha256` differed by
    OS. Matches `nfcore-scrnaseq-wrapper`.
  - `commands.sh` is now **self-contained and directly runnable** (`bash
    commands.sh`): it self-anchors from `BASH_SOURCE`, runs the pinned pipeline
    against the bundled `params.yaml`, requires no `CLAWBIO_REPO` and has no
    `<EDIT_ME>` in the runnable path, and **pins the Nextflow engine** via
    `export NXF_VER`. Re-running through ClawBio is now an optional, commented
    alternative. Previously a bare `bash commands.sh` aborted (`${CLAWBIO_REPO:?}`
    under `set -u`) and required editing an `<EDIT_ME>` output path.
  - `_portable_argv` now rewrites paths with resolved `Path.relative_to` instead
    of string prefixes, so a symlinked/unresolved output dir (e.g. macOS
    `/var`→`/private/var`) or `/` vs `\` separators no longer degrade valid paths
    to `<EDIT_ME>`.
  - The macOS-only `macos_docker.config` (Apple-Silicon `--platform`, VirtioFS
    `stageInMode`, host-scaled `resourceLimits`) is applied through a
    `uname`-gated `EXTRA_CONFIG` variable, so a macOS-generated bundle replays
    cleanly on Linux. The remap tool's samplesheet rewrite now forces
    `lineterminator="\n"` (csv defaults to CRLF).
- Preflight now validates an explicitly-provided `--genome` against the pinned
  42-key iGenomes catalogue (`SUPPORTED_IGENOMES_NAMES`, an exact mirror of
  `conf/igenomes.config`). The schema declares `genome` as a free string, so a
  typo'd/invalid key — e.g. bare `GRCh38` instead of `GATK.GRCh38` — was never
  caught and only failed late inside Nextflow with an opaque reference-resolution
  error. With the default iGenomes mirror an unknown key is now an `INVALID_GENOME`
  preflight error; with a custom `--igenomes-base` it softens to a warning (the
  mirror may define other keys). The absent-key (upstream default) and
  `null`/`igenomes_ignore` custom-FASTA paths are unaffected. This activates the
  previously-unused `SUPPORTED_IGENOMES_NAMES` constant and `INVALID_GENOME` code.
- `result.json` now carries a `samples_detected` field (the samples parsed from
  the actual upstream outputs) alongside the existing input-derived `samples`.
  In `--demo` the `test` profile supplies the samplesheet remotely, so `samples`
  is `[]` while the run still produces sample `test`; surfacing both keeps the
  machine-readable marker honest about what was actually processed.
- `remap_paths.py` now remaps reference/index paths in `params.yaml`
  (`--refs-old/--refs-new`), not just `commands.sh`. For nf-core/sarek the
  references (`--fasta`, `--dbsnp`, `--pon`, …) live in `params.yaml`, so this
  makes a bundle fully portable across machines/OSes in one command (data via
  `--old/--new`, references via `--refs-old/--refs-new`, then `--verify`). URIs,
  the `false` disable sentinel, and glob patterns are left untouched; `--verify`
  now also checks the params.yaml reference paths. The key list is kept in sync
  with `provenance.REFERENCE_PATH_PARAMS` by a parity test.

### Fixed

- **Direct invocation no longer fails with `ModuleNotFoundError: clawbio`.**
  CLAUDE.md documents `python skills/nfcore-sarek-wrapper/nfcore_sarek_wrapper.py
  --help` for the full flag surface, but Python only adds the *script's* directory
  to `sys.path`, not the repo root, so the top-level `from clawbio.common...` import
  failed unless the caller pre-set `PYTHONPATH`. The wrapper now bootstraps the repo
  root (`_SKILL_DIR.parent.parent`) onto `sys.path` before importing `clawbio`
  (a no-op under `clawbio.py run`, which already has it). Mirrors
  nfcore-scrnaseq-wrapper, which already did this. Guarded by
  `test_cli_direct_invocation.py` (subprocess `--help` from a clean env/cwd).
- **Replay instructions now match the self-contained `commands.sh`.**
  `remap_paths.py` (its docstring, the `--verify` success message, and the help
  epilog), the `report.md` reproducibility section, and `SKILL.md` still told
  users to replay with `CLAWBIO_REPO=… bash commands.sh` and claimed
  `commands.sh` "regenerates" the macOS Docker config. Since `commands.sh` is now
  self-contained (a bare `bash commands.sh`, no environment variable; the macOS
  config is a pre-written, `uname`-gated `-c reproducibility/macos_docker.config`),
  those hints were stale and contradicted the bundle itself and the
  `nfcore-scrnaseq-wrapper` replay hint. Corrected all five sites and added
  regression guards (`test_replay_hint_is_self_contained_no_required_clawbio_repo`,
  `test_report_md_describes_commands_sh_as_self_contained`). Also added an
  end-to-end remap roundtrip test that relocates FASTQs to a path containing a
  space and a non-ASCII character, then remaps and re-verifies — proving
  cross-machine portability (parity with `nfcore-scrnaseq-wrapper`).
- A missing `java` binary now raises `MISSING_JAVA` (previously the code
  contradicted its own "was not found" message by reporting the old
  `JAVA_TOO_OLD`), and a missing `nextflow` raises `MISSING_NEXTFLOW`. An
  unparseable Nextflow version now raises `NEXTFLOW_VERSION_TOO_OLD` rather than
  `NEXTFLOW_NOT_FOUND` (the binary is present; only the version gate fails).
  Matches the `nfcore-scrnaseq-wrapper` taxonomy and activates the
  previously-dead `MISSING_JAVA`/`MISSING_NEXTFLOW` codes.
- Renamed the `JAVA_TOO_OLD` error code to `JAVA_VERSION_TOO_OLD` so it is
  internally consistent with `NEXTFLOW_VERSION_TOO_OLD` and matches the sibling
  `nfcore-scrnaseq-wrapper` enum. No external consumers depended on the old code.
- Backend preflight now emits runtime-specific error codes (`MISSING_DOCKER`,
  `DOCKER_NOT_RUNNING`, `MISSING_PODMAN`, `PODMAN_NOT_RUNNING`, `MISSING_CONDA`,
  `MISSING_SINGULARITY`, `MISSING_HPC_RUNTIME`) instead of the single generic
  `BACKEND_UNAVAILABLE`, matching the sibling `nfcore-scrnaseq-wrapper` so both
  skills classify the same backend failure identically. `BACKEND_UNAVAILABLE`
  remains only as a fallback for any unmapped runtime. This activates the
  previously-dead specific codes; user-facing messages are unchanged.
- `--wes` help text no longer claims "(requires --intervals)". The preflight does
  not require an intervals file to be present with `--wes` (per the upstream docs a
  target BED is recommended, not mandatory); it only enforces BED format *when*
  `--intervals` is supplied. The help now reads "Whole-exome/panel mode (provide a
  target --intervals BED; enforced as BED when given)" to match actual behaviour.
- Consolidated the Sentieon emit-mode and `sentieon_dnascope_pcr_indel_model`
  enum checks onto their `schemas.py` constants (`SUPPORTED_SENTIEON_EMIT_MODE`,
  now the complete 10-value set including the `gvcf` pairings, and
  `SUPPORTED_GATK_PCR_INDEL_MODEL`) instead of inline regex/set literals. Removes
  duplicated allowed-value definitions and the previously-dead constants; behaviour
  is unchanged and now covered by explicit regression tests.
- Preflight now validates cheap, environment-independent input (step, aligner,
  enumerated params, required-tools-for-step, and unknown `--tools` tokens) BEFORE
  probing the container backend. Previously a typo'd `--step`/`--tools`/enum with the
  backend down surfaced `BACKEND_UNAVAILABLE` instead of the actionable input error.
  The unknown-tools check was extracted to `_check_tools_known` (called early and still
  self-guarded inside `_check_tools_against_pairing`); the backend probe still runs for
  valid input.
- The `mutect` profile is now rejected outside `--demo`. In 3.8.1 it only pulls in
  `conf/test_mutect2.config`, whose sole effect is forcing Mutect2's
  `--normal-sample normal` — a value valid only for the upstream test dataset. The
  wrapper previously exposed `--mutect-profile` (and the `mutect` token in `--profile`)
  as a general modifier with no guard, so a real somatic paired run whose normal was
  not named `normal` would be silently mislabelled. `_validate_wrapper_flags` now
  raises `INVALID_FLAG_COMBINATION` when `mutect` is composed without `--demo`; the
  CLI help and the `schemas.py` profile-token comment now describe it as test-data-only.
- The report's annotation table now links snpEff's stats HTML. In sarek 3.8.1
  the snpEff module runs `-csvStats <prefix>.csv` but passes no `-stats <prefix>`,
  so snpEff emits its summary under the default constant name `snpEff_summary.html`
  (in `reports/snpeff/<caller>/<sample>/`). The parser previously searched only for
  `<stem>_snpEff.html` / `<stem>_snpEff.ann.summary.html`, leaving the HTML column
  blank; it now also matches `snpEff_summary.html`, scoped to the per-sample report
  directory so it never leaks across samples. Verified on a live `--tools strelka,snpeff` run.
- The report's per-caller "Variant count" no longer counts all-sites gVCF /
  genome VCF records. When a caller emits both a called-variants VCF and a gVCF
  (Strelka `genome.vcf.gz`; HaplotypeCaller/DeepVariant/Sentieon `.g.vcf.gz`),
  `_best_variant_count` now deprioritises the gVCF — preferring a filtered call
  set, then the raw called VCF — so e.g. Strelka reports its 8 called variants
  instead of the 190 gVCF reference blocks. Falls back to the gVCF only when no
  other VCF is present.
- Execution logs (`stdout.txt`/`stderr.txt`) are now written under
  `reproducibility/logs/` instead of a stray top-level `logs/` directory, so the
  output root keeps to exactly two children (`upstream/`, `reproducibility/`) as
  documented. The whole `reproducibility/` tree is excluded from
  `checksums.sha256`, so logs no longer leak into the manifest.
- The macOS Docker compatibility config is now written to
  `reproducibility/macos_docker.config` (previously the output root).
- On error, the `result.json` marker is now written to
  `reproducibility/result.json` (falling back to the output root only if the
  bundle directory cannot be created), so the output root keeps to two children
  on failure paths too.
- `--demo` now also clears a user-supplied `--igenomes-ignore` (and any
  reference-path/genome keys passed through the `--extra-param` escape hatch).
  Previously `igenomes_ignore=true` (a boolean flag) and reference extras could
  leak into `params.yaml` in demo mode and break the upstream `test` profile,
  contradicting the "demo clears all reference flags" contract.
- Output-structure docs clarified: besides the two ClawBio children
  (`upstream/`, `reproducibility/`), Nextflow itself writes hidden `.nextflow/`
  and `.nextflow.log` in the launch dir (`= output_dir`, required so relative
  `input`/`outdir` resolve). Both are excluded from `checksums.sha256`. Verified
  against a live `--demo` run.
- Documentation flag-count corrected to a consistent **173** (154 Sarek
  passthrough + 19 wrapper controls) in README.md and CLAUDE.md; SKILL.md's
  "154 passthrough" was already correct. A regression test now locks these counts.
- A bare `--profile test` (no `--demo`, no `--input`) is recognised as a demo
  run before the input requirement is enforced, so it no longer fails with
  `MISSING_INPUT`. Profile-flag syncing now runs before wrapper-flag validation.
- `remap_paths.py._regenerate_checksums` now mirrors
  `provenance.py._iter_checksum_paths` exactly (excluding `work/`, `.nextflow/`,
  `reproducibility/`, `logs/`, and `.log` files), so `sha256sum -c` still
  verifies after a cross-machine remap. It previously hashed a divergent file
  set and referenced the removed top-level `logs/` and `provenance/` paths.

### Added

- Initial release of the wrapper around nf-core/sarek 3.8.1.
- 6-step orchestration: `mapping → markduplicates → prepare_recalibration →
  recalibrate → variant_calling → annotate`.
- Coverage for germline, tumor-only, and somatic paired analyses across the
  28 `--tools` tokens (variant/CNV/SV/MSI callers plus SnpEff / VEP / bcftools /
  merge annotation).
- Strict preflight (Java 17+, Nextflow >=25.10.2, profile composition,
  reference-path existence, annotation cache layout, flag compatibility,
  resume-manifest drift).
- 25 profiles supported (containers + arm64 + gpu + spark + mutect + 10 test
  variants).
- Reproducibility bundle: `params.yaml`, `samplesheet.valid.csv`,
  `commands.sh`, `manifest.json`, `checksums.sha256`, `environment.yml`,
  `pipeline_source.json`, `parameters.json`, `samplesheet.json`,
  `outputs.json`, `tool_versions.json`, `compatibility_policy.json`,
  `remap_paths.py`.
- Opt-in downstream handoff (`--run-downstream --downstream-skill <name>`)
  for `clinical-variant-reporter`, `clinical-trial-finder`,
  `omics-target-evidence-mapper`, `wes-clinical-report-en`,
  `wes-clinical-report-es`.
- The runtime floor now matches Sarek 3.8.1's official
  `manifest.nextflowVersion = '!>=25.10.2'` requirement.
- Provenance collects tool versions from the workflow-emitted
  `nf_core_sarek_software_mqc_versions.yml`, retaining compatibility with the
  generic `software_versions.yml` name shown in rendered output documentation.
- `varlociraptor` now requires a selected variant caller that supplies its
  candidate-variant input stream, as required by Sarek post-variant processing.
- Annotation routes now reject tools that `VCF_ANNOTATE_ALL` cannot execute at
  `--step annotate`, and require an upstream VCF-producing caller when
  annotators are selected from an earlier processing step.
- The preflight snapshot now carries direct CLI values for
  `only_paired_variant_calling`, `normalize_vcfs`, `concatenate_vcfs`, ASCAT
  purity/ploidy and `vep_custom_args`, so validations operate on the same
  effective parameters that are emitted to `params.yaml`.
- `varlociraptor` is now rejected with `filter_vcfs`, `normalize_vcfs`,
  `snv_consensus_calling`, or `concatenate_vcfs`, because upstream
  `POST_VARIANTCALLING` executes the Varlociraptor branch instead of those
  bcftools transformations when it is selected.
- `filter_vcfs` / `normalize_vcfs` now require at least one official
  `small_variantcallers` source, and `concatenate_vcfs` requires an effective
  germline VCF stream; requests that would publish no transformed output are
  rejected before execution.
- ClawBio integration: `python clawbio.py run sarek-pipeline ...` with full
  per-skill allowlist for 154 sarek passthrough flags + wrapper-only modifiers
  (174 long flags total).

### Verified aligned with nf-core/sarek 3.8.1

- All 25 profile names match `nextflow.config` verbatim.
- All 6 `--step` enum values match `nextflow_schema.json`.
- All 28 `--tools` and 15 `--skip-tools` enum tokens match the upstream regex.
- 5 `--aligner` choices match.
- 42 `--genome` (iGenomes) keys exact-mirror `conf/igenomes.config`.
- 47 reference-path params validated in preflight.
- Samplesheet schema fields (patient, sample, sex, status, lane, fastq_1/2,
  spring_1/2, bam, bai, cram, crai, table, vcf, variantcaller, contamination)
  match `assets/schema_input.json`.
- `--bcftools-filter-criteria` accepts free-form bcftools expressions
  (matches upstream default `-f PASS,.`).
- `--cf-window` accepts decimal values (upstream type `number`).
- `--freebayes-filter` typed as a string vcflib/vcffilter expression (upstream
  type `string`, default `"30"`) — accepts expressions like `QUAL > 20`.
- Output discovery mirrors all documented paths under
  `<outdir>/preprocessing/`, `<outdir>/variant_calling/`,
  `<outdir>/annotation/`, `<outdir>/reports/`, `<outdir>/pipeline_info/`,
  plus the documented `<outdir>/csv/` handoff (5 CSVs including
  `markduplicates_no_table.csv`), with legacy `<outdir>/preprocessing/csv/`
  fallback support.
- ASCAT auxiliary outputs (`metrics.txt`, `_tumourBAF.txt`,
  `_tumourLogR.txt`) captured.
- Mutect2 auxiliary outputs (`contamination.table`, `segmentation.table`,
  `pileups.table`, `artifactprior.tar.gz`, unfiltered VCF) captured.
- Manta tumor-only output (`tumor_sv.vcf.gz`) captured.
- MSIsensorPro output directory matched to `variant_calling/msisensorpro/`.
- Sentieon caller filename suffix (`.haplotyper`, `.dnascope`, `.tnscope`)
  resolved via `_VCF_SUFFIX`.
- ASCAT manual-purity override now enforces the official
  `ascat_purity`-requires-`ascat_ploidy` dependency before execution.
- Strelka pairing validation follows the executable `3.8.1` workflow
  (`STRELKA_GERMLINE` plus somatic Strelka); the usage-page caller matrix
  omits its executable germline branch.

### Documented gotchas

- ASCAT for WES emits the upstream warning recommending custom alleles/loci/loci_gc/loci_rt resources.
- Mutect2 without `--pon` (any mode) or `--germline_resource` produces unreliable calls (warning).
- GATK MarkDuplicates Spark cannot perform header/positional UMI deduplication
  (`--umi-in-read-header` or `--umi-location`).
- iGenomes auto-fills `genome: 'GATK.GRCh38'` unless `--igenomes-ignore` is set.
- Somatic mode requires same-`patient` rows with `status=0` AND `status=1`.
- `--demo` clears every reference flag; `--profile test*` / `test_full*`
  trigger demo cleanup automatically.
- `--resume` is rejected on manifest drift (step, aligner, tools, profile,
  arm/gpu/spark, params/samplesheet checksums).
- `arm64` profile implicitly enables Wave (needs outbound network).
- Sentieon callers/aligner need `SENTIEON_LICENSE` env var.

### Fixed (schema/output alignment audit)

- Samplesheet validation no longer rejects cross-row `sex` / `status` /
  `contamination` differences for one sample beyond Sarek's documented and
  executable checks; pairing inference now deduplicates by `(sample, status)`
  exactly as `samplesheet_to_channel/main.nf` does.
- Omitted `genome` now resolves to Sarek's official `GATK.GRCh38` default in
  preflight, while explicit `--genome null` / YAML-null custom-reference
  configurations still disable iGenomes.
- Azure `az://` resource and samplesheet URIs are now treated as remote paths,
  matching Sarek's official annotation-cache URI recognition; unsupported
  `azure://` is no longer accepted as a remote-path alias.
- `bbsplit` mapping runs now require either a pre-built index or FASTA list,
  and `ngscheckmate` runs require its SNP BED unless an official iGenomes
  bundle supplies it, matching the channels consumed by Sarek 3.8.1.
- iGenomes-backed warnings now report Sarek's effective default genome
  (`GATK.GRCh38`) instead of an omitted wrapper parameter.
- ASCAT preflight retains the `ascat_genome` requirement described by the
  official schema and consumed as `genomeVersion` by the ASCAT module, in
  addition to the archive checks in `samplesheet_to_channel`.
- Spark output compatibility now rejects only the executable conflict
  (`save_mapped` together with `save_output_as_bam`), and VCF filtering accepts
  the official omitted `bcftools_filter_criteria` default (`-f PASS,.`).
- Provenance, reports, output parsing and resume drift now record/compare the
  effective upstream defaults `step=mapping` and `aligner=bwa-mem` when those
  optional parameters are omitted from `params.yaml`.
- Partial iGenomes reference overrides remain allowed: Sarek's official usage
  guide documents replacing an individual reference file while keeping the
  selected `--genome` bundle. Passing `--fasta` no longer implicitly forces
  `--igenomes-ignore`; wholly custom runs still use the documented
  `--genome null --igenomes-ignore --fasta` form.
- The official `false` reference sentinel is preserved for first-class
  reference flags and is not incorrectly restored from iGenomes during
  preflight; `--extra-param` now overrides a supplied base `--params-file`.
- Mutect2 preflight now evaluates effective inherited PON/germline resources;
  `GATK.GRCh38` reports the bundled generic PON recommendation rather than
  incorrectly reporting that no PON exists.
- Downstream runs without `--input` now resolve Sarek's official CSV handoffs
  (`mapped.csv` through `variantcalled.csv`), validate them, and submit the
  normalized sheet as `input`; explicit `input_restart` is handled the same
  way because Sarek v3.8.1 overwrites that schema field internally.
- Steps and tools supplied by `--params-file` or `--extra-param` are resolved
  before input validation, so downstream handoff retrieval is not incorrectly
  rejected as an input-less `mapping` run.
- Remote `--input` samplesheets are now staged with `nextflow fs cp` before
  wrapper validation, matching the HTTPS samplesheets used by Sarek's official
  `test_full*` profiles while retaining cloud URI support through Nextflow.
- `freebayes_filter` retyped from integer to string (vcflib/vcffilter expression).
- `cf_ploidy` reclassified as a string param (schema type `string`, default `"2"`)
  so comma-separated ploidy lists like `"2,3,4"` are preserved.
- Strelka no longer forced to tumor mode — runs germline (single sample) and
  somatic (tumor/normal pair), matching upstream.
- `--genome null`/`none` sentinels no longer counted as an active iGenomes genome.
- `step`/`tools` supplied via `--params-file` now reach samplesheet validation
  (previously validated as `mapping`).
- `patient-sample-status-lane` uniqueness now enforced for every step (sarek
  validates it before step branching), not just `mapping`.
- `--snv-consensus-calling` maps the public `mpileup` tool token to its emitted
  `bcftools` VCF caller and includes `sentieon_tnscope`; only the internal
  `samtools` pileup stream used by ControlFREEC is excluded upstream.
- `--vep-phenotypes` permits VEP's automatic data download and only requires an
  index for a supplied gzipped phenotype file; `bcftools_columns` is validated
  as a path.
- LoFreq output discovery accepts both the rendered output-guide form
  `lofreq/<sample>/<tumorsample>.vcf.gz` and the executable
  `conf/modules/lofreq.config` form `<tumorsample>.lofreq.vcf.gz`; CNVkit files globbed by sample prefix so tumor/normal pairs and the
  plain `.cns` vs `.call.cns`/`.bintest.cns` variants are not missed/misclassified;
  Strelka germline `genome.vcf.gz` captured.
- Tool×mode validation is now evaluated against the SET of per-patient modes
  present, not a single global mode: mixed samplesheets (germline-only patients
  + tumor/normal pairs) are accepted, with each tool allowed when at least one
  patient matches its required mode (matches sarek's per-patient processing).
- `snpeff`/`merge` require `--snpeff_db` when an iGenomes genome is not active;
  with iGenomes selected, the database may be sourced from its genome config.
- SnpEff database validation also applies to `--download_cache` /
  `--build_only_index`, because the upstream cache download consumes
  `snpeff_db`.
- Interval validation follows the published `v3.8.1` parameter schema:
  `.bed` or `.interval_list` are accepted, with WES remaining BED-only.
- Output discovery now records the documented Parabricks preprocessing
  directory and the executable `sentieon_consensus` directory emitted when
  Sentieon consensus deduplication is enabled.
- Pair reporting follows Sarek's tumor/normal `cross()` behavior, and output
  discovery records the BAM/CRAM conversion directories published by its configs.
- MuSE result parsing now preserves its documented `.MuSE.txt` call table in
  addition to the compressed VCF.
- User-supplied `varlociraptor_scenario_*` inputs are resolved, validated, and
  included in provenance because Sarek consumes them with `Channel.fromPath`.
- Local `outdir_cache` values are converted to absolute paths before launch so
  Nextflow's output-directory working directory cannot redirect cache output.
- Relative resource paths loaded through `--params-file` are canonicalised
  before launch because this wrapper executes Nextflow from its output directory.
- `known_indels`/`known_indels_tbi` now honour their official
  `file-path-pattern` contract, including deterministic provenance for matches;
  existing reference/cache directories are no longer labelled as missing.
- Mapping-specific UMI incompatibility checks are no longer applied to restart
  runs from later steps, matching `samplesheet_to_channel/main.nf`.
- Joint germline runs now surface Sarek's non-blocking warning when VQSR
  resources or required intervals are absent.
- `--snv-consensus-calling` now requires `--normalize-vcfs`, matching the
  upstream workflow guard.
- Germline callers accept tumor/normal cohorts because Sarek routes the normal
  member through germline calling; they are rejected only when no normal
  sample is present.
- `mpileup` accepts tumor/normal cohorts for that same executable reason:
  Sarek calls its normal member through `BAM_VARIANT_CALLING_GERMLINE_ALL`.
- Sex checks for ASCAT, ControlFREEC and Varlociraptor mirror Sarek's
  `input_sample.map` guard and reject `NA` on every supplied sample row when
  any of those tools is selected.
- ASCAT validation now includes its required `ascat_genome`, while
  MSIsensorPro no longer requires a precomputed scan because Sarek generates
  one from the FASTA when absent; an iGenomes-supplied scan is now reported as
  inherited instead of incorrectly reported as newly generated.
- Mutect2 now warns about a missing `--pon` in ANY mode (not only tumor-only)
  and additionally warns about a missing `--germline_resource`, matching sarek.
- `indexcov` is rejected when combined with `--wes`: both variant-calling
  workflows guard it with `params.wes == false`, so accepting that combination
  would request an output Sarek never executes.
- `merge` annotation output suffix corrected to `_snpEff_VEP` (was `_merge`) and
  `bcfann` to `_BCF` (was `_bcf`); annotation suffixes are matched longest-first
  so a merge VCF is not mis-detected as a plain VEP output.
- Annotation outputs are now discovered recursively under
  `annotation/<variantcaller>/<sample_or_pair>/` (the real sarek 3.8.1 layout),
  not just one level below `annotation/`.
- SnpSift annotation output (`_snpSift.ann`) is now parsed and validated
  (`snpsift` added to the annotator set).
- VEP cache layout check corrected to `${vep_species}${suffix}/${vep_cache_version}_${vep_genome}`
  (where `${suffix}` is `_merged` or `_refseq` when requested), per the Sarek
  cache-initialisation workflow and usage documentation — no more
  false "cache layout" warnings.
- `lane` is validated against the schema pattern `^\S+$` (whitespace rejected).
- `--umi_read_structure` is validated against fgbio read-structure grammar at
  preflight (e.g. `8M2S+T`, `+T +T`), instead of failing late in Nextflow.
- DragMap aligner without `baserecalibrator` in `--skip_tools` now emits sarek's
  recommended warning.
- Enumerated pass-through params are validated at preflight: `--umi-location`,
  `--group-by-umi-strategy`, `--vep-out-format`, `--ascat-genome`,
  `--publish-dir-mode`, and `--use-gatk-spark` tokens.
- SnpEff QC reports discovered at `reports/snpeff/` (lowercase, per
  conf/modules/annotate.config) — the previous `reports/SnpEff/` lookup missed
  them on case-sensitive filesystems; the capitalised variant is still accepted
  as a fallback. Added `reports/sentieon_dedup/` to the QC report inventory.
- Preprocessing parsing now captures the `mapped/` CRAM/BAM index
  (`.sorted.cram.crai`) and the BQSR recalibration tables under
  `preprocessing/recal_table/<sample>/<sample>.recal.table`.
- Mapped BAM/CRAM discovery follows
  `conf/modules/{aligner,markduplicates}.config`: outputs are captured when
  `--save-mapped` is enabled and also when `--skip-tools markduplicates`
  makes mapped files the published deliverable, including the default CRAM
  route when `--save-output-as-bam` is omitted (except when
  `sentieon_dedup` takes over that route).
- Generic nf-core parameters supplied through `--extra-param` retain their
  official schema types: institutional string values such as
  `custom_config_version=false` are no longer serialized as booleans, while
  documented generic boolean fields remain validated booleans.
- `--only-paired-variant-calling` now rejects germline-only callers (and
  paired-only `mpileup`) when all available normals are matched and are
  therefore removed from Sarek's germline input channel; mixed cohorts that
  still supply a runnable germline or tumor-only input remain valid.
- Handoff-derived sample detection now uses the most advanced CSV available
  (`variantcalled.csv` back to `mapped.csv`) instead of accidentally preferring
  the earliest intermediate sheet.
- Reference output discovery is enabled for `--build-only-index` as well as
  `--save-reference`, matching the `PREPARE_GENOME` and `PREPARE_INTERVALS`
  publication conditions in Sarek.
- Reference discovery now captures non-empty files published directly under
  `reference/` as well as files inside named subdirectories, matching index
  publication paths in `conf/modules/prepare_genome.config`.
- Root-level `indexcov` artefacts are now catalogued as cohort outputs rather
  than inferred as a biological sample named `indexcov`, matching its official
  non-sample-subdirectory publication layout.
- `--build-only-index` follows Sarek's empty-samplesheet execution path:
  row/mode and caller-output checks are suppressed while upstream global
  resource guards (including BQSR on preprocessing start steps) and cache
  validation stay active.
- Resume drift detection now compares the complete effective params checksum
  and reference fingerprints already stored in the manifest, preventing cached
  work reuse after a scientifically material parameter or reference change.
- Manifest reference fingerprints record Sarek's official `false` reference
  sentinel as `<disabled>` instead of incorrectly labelling it as a missing
  local path.
- Local reference/index/cache directories are fingerprinted recursively from
  relative file names and contents, so `--resume` cannot reuse cached work
  after silent directory-content changes.
- `--input false` (sarek's no-samplesheet sentinel, used with
  `--build-only-index`; cache download may be layered on that run) is normalised to "no input" instead
  of being validated as a path literally named "false".
- Strelka is now rejected in tumor-only mode; Sarek's executable routing calls
  it in germline and somatic-paired workflows. `freebayes`/`cnvkit` remain
  available in every supported analysis mode.
- `spring_2` without `spring_1` is rejected (schema dependency `spring_2 -> spring_1`).
- Mapping-step BAM input now captures the optional `bai` (sarek emits
  `[meta, bam, bai]`).
- MarkDuplicatesSpark + UMI rule corrected: incompatible only with header/
  positional UMIs (`umi_in_read_header`/`umi_location`), NOT with
  `umi_read_structure` (fgbio consensus runs upstream).
- Mutect2 warns when the default GATK Panel-of-Normals (`1000g_pon.hg38.vcf.gz`)
  is used, matching sarek.
- VEP/merge annotation outputs in `json`/`tab` format (`*.ann.json.gz`,
  `*.ann.tab.gz`) are now discovered, not only `*.ann.vcf.gz`.
- `step_completed` recognises a SnpSift-only annotation run.
- Numeric params are bounds-checked before launch (schema minimums): `split_fastq`
  (0 or >=250), `umi_length`>=1, `umi_base_skip`>=0, `ascat_*`/`cf_*` mins, etc.
- Integrated `clawbio.py run sarek-pipeline --help` delegates to the wrapper's
  schema-derived parser, and shared parser flags are forwarded without being lost.
- `commands.sh` now replays the captured, path-remapped Nextflow argv directly,
  preserving custom configs, work directory and resume state.
- The macOS Docker compatibility config no longer replaces Sarek's official
  `docker`/`gpu`/`spark` run options with a platform-only value: AMD64
  emulation now retains UID/GID mapping and GPU access where the selected
  profile defines them.
- `--extra-param` now uses one effective value throughout samplesheet
  validation, preflight and `params.yaml`: exposed numeric/boolean Sarek
  parameters retain their schema types, final overrides cannot be shadowed by
  earlier dedicated flags, and unknown native pass-through keys cannot toggle
  wrapper-only controls such as `--demo` or `--resume`.
- `--extra-param` cannot replace wrapper-managed `input`, `input_restart`, or
  `outdir`, preventing a validated normalized samplesheet or the tracked
  results tree from being silently bypassed after preflight.
- Relative `igenomes_base` paths loaded through `--params-file` are now
  resolved before the wrapper changes Nextflow's launch directory, matching
  the already-normalised CLI path behaviour.
- Cache downloads accept supported remote `--outdir-cache` targets (`s3://`,
  `gs://`, `az://`) without misapplying local filesystem writability checks;
  local targets remain checked before launch.

### Test suite

- 258 focused pytest tests covering every module (preflight, samplesheet,
  params builder, command builder, executor, outputs parser, provenance,
  reporting, orchestrator, remap paths).
- Audit-driven regression tests for every issue caught during the alignment
  reviews against nf-core/sarek 3.8.1 official docs.
