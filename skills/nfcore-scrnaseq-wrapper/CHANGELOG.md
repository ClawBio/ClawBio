# nfcore-scrnaseq-wrapper — Changelog

All notable changes to the `nfcore-scrnaseq-wrapper` ClawBio skill are documented
here. The format roughly follows [Keep a Changelog](https://keepachangelog.com)
and the wrapper version is tracked in `SKILL.md` YAML frontmatter.

## [Unreleased] — 0.1.0

### Fixed

- **A lane-split sample is now counted once, not once per row.** The samplesheet summary
  counted raw CSV rows (`sample_count = len(normalized_rows)`) and appended to
  `sample_names` with no duplicate guard. nf-core/scrnaseq documents re-sequencing a sample
  across lanes by repeating the `sample` identifier on several rows — *"The `sample`
  identifiers have to be the same when you have re-sequenced the same sample more than
  once"*, and the pipeline *"concatenate[s] the raw reads before performing any downstream
  analysis"* — so rows > samples is normal and expected. A 2-sample sheet with one sample
  split over two lanes therefore reported **3** samples and listed the split sample twice.

  This was not purely cosmetic. The value is written into the provenance bundle
  (`inputs.json` → `sample_count`) and is *preferred* by `report.md` over the samples
  actually detected in the outputs, so the run's own audit trail claimed more samples than
  the run had — while the pipeline itself merged the lanes correctly, making the
  discrepancy invisible at execution time.

  `sample_names` is now deduplicated and `sample_count` is `len(sample_names)`, matching
  `nfcore-rnaseq-wrapper` and `nfcore-sarek-wrapper`, which both already did this;
  scrnaseq was the only outlier. Every lane row is still written to the normalized
  samplesheet unchanged, so the pipeline concatenates them exactly as before — only the
  reported/recorded count changes.
- **`remap_paths.py`: combining actions no longer discards one of them silently.**
  `main()` dispatched through a chain of early returns, so only the first matching flag
  ever ran. Here `--output-dir` came first, so `--output-dir <new> --verify` printed the
  self-relocation hint and then silently dropped `--verify`: the user asked for a
  readiness check, never got one, and saw exit code 0 implying all was well. (The
  rnaseq/sarek bundles had the mirror image — `--verify` ran first and silently discarded
  the `--output-dir` rewrite, then reported a stale bundle as "ready to replay".)

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

### Documentation

- **The `validation.*` / `validationSchemaIgnoreParams` warning is documented as upstream,
  with a regression guard.** A review reported the bundle's `commands.sh` emitting
  unrecognised-config warnings for `validation.defaultIgnoreParams` and
  `validation.monochromeLogs`, and attributed them to the wrapper. They are not the
  wrapper's: nf-core/scrnaseq 4.1.0's own `nextflow.config` declares
  `plugins { id 'nf-schema@2.5.1' }` and, directly beneath it, the matching
  `validation { defaultIgnoreParams = ["genomes"]; monochromeLogs = params.monochrome_logs }`
  scope, plus a `monochromeLogs = null` param that upstream itself annotates
  `// TODO temporary workaround a warning` (with a matching ignored lint check in its
  `.nf-core.yml`). nf-schema surfaces this as
  `WARN: The following invalid input values have been detected: * --validationSchemaIgnoreParams: genomes`.

  Verified against a real run: the wrapper's `params.yaml` contains none of these keys
  (only the audited CLI surface), the plugin loads correctly
  (`Plugins declared=[nf-schema@2.5.1]`), and the run's `.nextflow.log` contains zero
  "unknown config attribute" entries — the identical warning appears on a plain upstream
  `nextflow run nf-core/scrnaseq -r 4.1.0`. No wrapper change is warranted; the warning is
  harmless and belongs to the pipeline.

  Documented as a Gotcha in SKILL.md, and pinned by a new regression test
  (`test_params_never_carry_upstream_plugin_config_scopes`) asserting the wrapper never
  writes `validation` / `validationSchemaIgnoreParams` / `monochromeLogs` / `prov` into
  `params.yaml`. Injecting those to silence the warning would override the pipeline's own
  config through `-params-file` — the wrong fix, and the one a future reader is most likely
  to reach for.

### Fixed

- **Normalized replay samplesheets now match the nf-core/scrnaseq 4.1.0 input
  schema exactly.** Unsupported metadata headers are reported but omitted from
  `reproducibility/samplesheet.valid.csv`, preventing nf-schema's unidentified-
  header warning during replay. This specifically removes the non-schema
  `protocol` column; the effective protocol remains the documented global
  `--protocol` pipeline parameter.

### Added

- **Bundle `remap_paths.py` gains `--repair-bundle` (crash-recovery parity).** The
  shipped standalone `remap_paths.py` can now regenerate the three regenerable
  reproducibility-bundle files — `checksums.sha256` (recomputed from the current
  bundle contents), and post-hoc `manifest.json` / `environment.yml` stubs (marked
  `regenerated_post_hoc: true`) — if the wrapper crashed mid-post-processing and
  left them missing. This matches the `nfcore-rnaseq-wrapper` / `nfcore-sarek-wrapper`
  bundles, which already shipped `--repair-bundle`; scrnaseq was the only one
  without it. `checksums.sha256` is regenerated stdlib-only over **exactly the
  wrapper's original checksum allowlist** — the normalized samplesheet, `params.yaml`,
  `logs/stdout.txt` / `logs/stderr.txt`, the h5ad candidates (canonical
  `*/mtx_conversions/` matrices at both nesting depths, else a full-tree fallback),
  and the MultiQC report — reconstructed from the bundle on disk so a repaired manifest
  is byte-identical to one `provenance.write_reproducibility_checksums` would write. The
  `reproducibility/` provenance JSON tree, `commands.sh`, the `environment.yml` /
  `manifest.json` stubs, other `upstream/results` outputs, and non-`stdout/stderr` logs
  are deliberately excluded, matching the original manifest (labels are relative to the
  output directory so `sha256sum -c` passes). A self-contained `_write_text_lf`
  (mirroring the sibling bundles) keeps the regenerated files LF-only on every OS.
  Recompute is triggered **only when `checksums.sha256` itself is missing**, so an
  existing valid manifest is never silently overwritten when merely a stub is
  regenerated. `--output-dir` and the `--output`-based resume warning are deliberately
  **not** added: scrnaseq's `commands.sh` self-anchors to the bundle location via
  `BASH_SOURCE`, so the output path never needs patching.

### Added

- **`remap_paths.py --output-dir` accepted for CLI parity.** The rnaseq/sarek bundles
  expose `--output-dir` to rewrite a baked `--output`; the scrnaseq bundle self-relocates
  (its `commands.sh` self-anchors and `params.yaml` is output-relative), so it never
  needed one. `--output-dir` is now accepted as a documented no-op that explains the
  self-relocation, so a user reaching for it out of habit gets guidance instead of an
  argument error. SKILL.md documents the difference across all three wrappers.

### Documentation

- **`--allow-remote-inputs` semantics clarified.** SKILL.md now states explicitly
  that the flag relaxes only the wrapper's own local-first preflight check: remote
  FASTQ/reference URIs are written into the normalized samplesheet/`params.yaml`
  verbatim and staged natively by Nextflow at run time. The wrapper does not
  download them, so remote inputs require outbound network access and are
  incompatible with `NXF_OFFLINE` (under which Nextflow's own nf-schema file-existence
  validation still runs and fails on the remote paths). Shared wording across the
  three wrappers.

### Fixed

- **`--demo` no longer prints a misleading "requested" preset warning.** Running `--demo`
  without `--preset` fired `WARNING: --demo forces preset=star (requested: 'standard')`,
  implying the user had asked for `standard` when it was merely the parser default. The
  demo-override warning now consults the existing `preset_explicit` flag (set only by
  `_PresetAction` when `--preset` is actually supplied on the CLI), so it fires **only**
  when the user explicitly requested a non-`star` preset that `--demo` overrides. Purely
  cosmetic — the run result was always identical — but the message no longer fabricates a
  request the user never made.
- **Output-dir "not empty" check now ignores the whole `reproducibility/` bundle.** The
  check used a per-file allowlist (`_ALLOWED_REPRO_FILES` = samplesheet/params/macos
  config only), so an output directory containing a *complete* reproducibility bundle —
  with the wrapper's own `commands.sh`, `checksums.sha256`, `environment.yml`,
  `remap_paths.py`, `resource_limits.config` and provenance JSON — was rejected with
  `OUTPUT_DIR_NOT_EMPTY` on a non-resume run, while nfcore-sarek accepted it. The
  `reproducibility/` directory is entirely wrapper-generated (never user data), so it is
  now ignored wholesale (added to `_IGNORED_ROOT_NAMES`), matching nfcore-sarek and
  nfcore-rnaseq. Genuine result artifacts at the output root (`report.md`, `result.json`,
  `upstream/`, `logs/`) still block a non-resume re-run.
- **The host resourceLimits cap now ships in the reproducibility bundle and replays.**
  Previously the cap was applied to the live run but written outside the bundle, so a
  from-scratch reproduction on the generating host re-aborted at `STAR_GENOMEGENERATE`
  with `Process requirement exceeds available memory`. It is now written to
  `reproducibility/resource_limits.config` and `commands.sh` re-applies it via a
  `uname != Darwin` guard, so a bundle reproduces the run on the machine that made it.
  The cap only ever clamps requests, so a larger replay host simply under-uses it.
- **Host memory auto-cap now also applies on Linux (not only macOS).** The host-scaled
  `process.resourceLimits` cap that prevents Nextflow's local executor from aborting a
  real run with `Process requirement exceeds available memory` was written only on
  macOS+docker, so a normal Linux docker run got no cap and failed whenever a process's
  production request (e.g. STARsolo's index build) exceeded the host. A non-macOS docker
  run now writes a portable resourceLimits config (`.nextflow_resource_limits.config`)
  scaled to the machine — the smaller of physical RAM and the Docker `MemTotal`, minus
  headroom, no 15 GB macOS-VM ceiling — per nf-core's resourceLimits guidance. It is
  applied to the **live** run only; the portable `commands.sh` bundle (rebuilt from
  `params.yaml`) stays host-independent. `--demo` and the macOS path are unchanged.
- **`commands.sh` replay is now idempotent (`-resume`).** The nextflow-direct replay
  script re-ran the whole pipeline from scratch against an already-completed output dir.
  It now emits a guard that adds `-resume` when a prior Nextflow session (`.nextflow/`)
  exists in the launch dir, so a first run starts fresh but a replay reuses the cache —
  matching the `nfcore-rnaseq-wrapper` bundle. `build_nextflow_commands_sh` no longer
  takes a `resume` flag; the guard is the single source of truth.
- **IPv6/NAT64 hint now actually shows (scan `.nextflow.log`).** The `EXECUTION_FAILED`
  environment-hint scanner read only `logs/stdout.txt` and `logs/stderr.txt`. When
  Nextflow fails while parsing the config (before the pipeline starts), the underlying
  "Network is unreachable" cause is often recorded only in `.nextflow.log`, so the
  IPv6/NAT64 hint (`NXF_OPTS='-Djava.net.preferIPv6Addresses=true'`) never appeared.
  `.nextflow.log` (in the Nextflow launch cwd) is now scanned too. Shared verbatim with
  the sibling wrappers; covered by a new test.
- **Report "Next Steps" handoff is now portable.** The downstream handoff bullets in
  `report.md` baked an **absolute** path to `clawbio.py` (resolved on the machine that
  generated the report, e.g. `/mnt/.../ClawBio/clawbio.py`), so the suggested command did
  not exist on any other machine, and invoked it with a bare `python`. They now use the
  repository-relative `clawbio.py` (matching the convention documented across the repo)
  and `python3` (portable to python3-only systems, PEP 394). Covered by a new test.
- **`save_align_intermeds` no longer carries an invented `default: true` in the
  contract.** nf-core/scrnaseq 4.1.0's `nextflow_schema.json` declares no default for
  `save_align_intermeds` (an opt-in boolean, false when unset — same as
  `save_reference`); the `OFFICIAL_PARAMS` contract wrongly recorded `default: True`.
  The annotation was reference-only (never emitted to `params.yaml`, which the wrapper
  writes only for user-supplied flags, and never used for validation), so no run was
  affected — but it misrepresented the upstream schema. Corrected to `{"type":
  "boolean"}` and the contract test now asserts the upstream reality (no `default` key).
- **`params.yaml` now carries the explanatory header the sibling wrappers emit.**
  `serialize_params_yaml` wrote a bare `yaml.dump`, so the bundled `params.yaml`
  lacked the three-line comment header that nfcore-rnaseq / nfcore-sarek prepend
  (noting that `input`/`outdir` are relative to the Nextflow launch directory and
  pointing at `reproducibility/commands.sh` for replay). It now emits the same header,
  adapted to nf-core/scrnaseq. YAML comments do not affect parsing, the reference-path
  remapper, or the bundle checksum round-trip.
- **A `protocol` samplesheet column is no longer flagged as unrecognised.**
  nf-core/scrnaseq 4.1.0's `schema_input.json` does not define a `protocol` column,
  but the pipeline's own example samplesheet (`assets/samplesheet.csv`) ships one,
  so a user copying it was warned that `protocol` is an unrecognised column. It is
  now recognised (and preserved as before). The effective protocol is still taken
  from the `--protocol` flag — a global pipeline parameter, not a per-sample column
  — and the "explicit protocol required" preflight error now says so explicitly, so
  a user with a `protocol` column knows to pass `--protocol` as well.
- **Config-parse failures now point at `NXF_OFFLINE`.** On Nextflow 26.x the
  nf-core `nextflow.config`'s `includeConfig … ? <url> : '/dev/null'` line fails to
  parse when the remote `nfcore_custom.config` cannot be fetched. The executor's
  `EXECUTION_FAILED` fix now detects `Unable to parse config file` /
  `ConfigParseException` and suggests `NXF_OFFLINE=true` for a fully local run (or
  confirming outbound HTTPS/DNS). Shared verbatim across the three wrappers.
- **Remote reference URIs are no longer corrupted in `params.yaml`.** With
  `--allow-remote-inputs`, a remote `--fasta`/`--gtf`/`--transcript-fasta`/index
  URI (`https://`, `s3://`, `gs://`, `ftp://`, …) was resolved as a local path by
  `params_builder`, collapsing the scheme and anchoring it to the working
  directory (`https://host/x` → `<cwd>/https:/host/x`), which nf-core then rejected
  during parameter validation. Reference fields now preserve URIs verbatim (via a
  new `_posix_or_uri` helper, mirroring `igenomes_base` and the FASTQ handling) and
  only resolve genuinely local paths. Remote references are still gated behind
  `--allow-remote-inputs` at preflight and staged by Nextflow.
- **CellBender failures now point at `--skip-cellbender`.** `CELLBENDER_REMOVEBACKGROUND`
  estimates ambient RNA from the droplet-count distribution and errors on very small
  or test datasets (`IndexError: index -100 is out of bounds`). When it is the
  failing process, the `EXECUTION_FAILED` fix now explains that CellBender is optional
  and can be skipped with `--skip-cellbender` so the rest of the pipeline finishes.
  Diagnosed from the actual log signal, so it only appears when relevant.
- **`--demo` no longer reports `Samples: 0`.** The reported sample count trusted a
  local samplesheet count of `0` instead of falling back to the samples detected
  in the outputs. Under `--demo` the upstream `-profile test` supplies samples
  remotely (so the local count is 0) while the run still produces real samples;
  the count now uses a positive local count when present and otherwise falls back
  to `samples_detected`.
- **`result.json` carries the shared `ok`/`status` contract.** A successful run's
  envelope now includes `status: "ok"` and `ok: true`, and a failed run includes
  `status: "error"` alongside the existing `ok: false` — a minimal discriminator
  shared with nfcore-sarek/rnaseq (via the opt-in `status`/`ok` parameters on the
  shared `write_result_json` helper).

- **`clawbio run scrnaseq-pipeline` now forwards `--nextflow-config`.** The launcher
  (`clawbio/cli.py`) forwarded only `-c`/`--config`, silently dropping configs
  supplied via `--nextflow-config`. All three spellings are now a single repeatable
  option, normalised and forwarded as `--nextflow-config` (the wrapper accepts them
  as aliases).
- **Output directory inside the ClawBio source tree is now rejected** at preflight
  with the dedicated `OUTPUT_DIR_INSIDE_REPO` code, matching nfcore-rnaseq and
  nfcore-sarek. Previously this wrapper alone permitted writing multi-gigabyte
  pipeline outputs into the repository checkout.
- **`--demo` under `NXF_OFFLINE` now fails fast with a clear message.** Demo mode
  runs nf-core's upstream `-profile test`, whose FASTQs and references are remote
  GitHub URLs; on an offline/sandboxed host the nf-schema plugin previously aborted
  with a cryptic `does not exist`. Preflight now detects `NXF_OFFLINE` + demo and
  raises `DEMO_REQUIRES_NETWORK` with an actionable fix. Docs (SKILL.md, AGENTS.md)
  clarify that `--demo` downloads only nf-core public test data — no user/genetic
  data is uploaded — so it is compatible with the local-first guarantee.
- **nf-core-native snake_case flag spellings are now accepted via the launcher.**
  nf-core parameters are snake_case while the wrapper exposes them as hyphenated
  flags. A user copying an upstream nf-core command previously had the token
  silently dropped by the launcher's INT-001 allowlist filter, which matches exact
  (hyphenated) tokens. `clawbio run scrnaseq-pipeline` now canonicalises `_`↔`-`
  when matching the allowlist and forwards the wrapper's canonical hyphen spelling
  (which its parser already accepts). Scoped to the three nf-core pipeline skills;
  other skills keep exact-match filtering.
- **Environment post-failure hints on `EXECUTION_FAILED`.** When a run fails, the
  executor scans the captured Nextflow logs and appends an actionable hint on a
  known environment signature — diagnosed from the actual error text, so no
  resource thresholds are invented. `Process requirement exceeds available memory`
  yields a hint to cap resources via a `-c` config using `process.resourceLimits`;
  `Network is unreachable` / a Java connection exception (common on IPv6-only /
  NAT64 hosts, where the JVM prefers IPv4) yields a hint to verify outbound
  DNS/HTTPS and set `NXF_OPTS='-Djava.net.preferIPv6Addresses=true'`. Shared
  verbatim across the three wrappers.

### Changed

- **Remote input/reference URIs are gated behind `--allow-remote-inputs`
  (local-first by default).** Samplesheet FASTQs and reference/index paths must
  be local unless the flag is passed; otherwise preflight rejects remote URIs
  (`s3://`, `gs://`, `https://`, `ftp://`, …) with `REMOTE_INPUT_NOT_ALLOWED`, so
  the "keeps processing local / prevents accidental cloud access" guarantee is
  enforced by the code. With `--allow-remote-inputs` (a flag shared by all three
  nf-core wrappers), remote values pass through verbatim — only the basename is
  validated, existence is deferred to Nextflow staging — and preflight emits a
  runtime warning naming every path fetched over the network. Local paths are
  always existence-checked (`MISSING_FASTQ` / `MISSING_REFERENCE`). Resolved FASTQ
  values are `Path | str` so remote URIs are written unchanged to the normalized
  samplesheet and provenance (no `Path().as_posix()` `//`-collapse).
- **SKILL.md frontmatter aligned to the canonical `SKILL-TEMPLATE` schema** so the
  machine catalog (`skills/catalog.json`) and the generator never drift: the
  `trigger_keywords` list now lives under `metadata.openclaw` (was a bare
  `metadata.trigger_keywords` key) and a `metadata.demo_data` entry was added.
  Matches `nfcore-sarek-wrapper` and `nfcore-rnaseq-wrapper`.
- **Validation philosophy unified with the sibling wrappers.** The launcher-side
  `os.access(R_OK)` FASTQ readability pre-check was removed: the wrapper never
  reads FASTQ data — Nextflow does, often inside a root container under the
  default Docker profile — so a launcher-side readability probe can false-block
  valid runs. Existence and regular-file type are still validated; readability is
  deferred to Nextflow's staging in the true execution context. A path that
  exists but is not a regular file now raises `MISSING_FASTQ` (was
  `FASTQ_NOT_READABLE`). Mirrors `nfcore-rnaseq-wrapper`'s documented policy.

### Changed

- **Logging consistency**: added a dedicated `[provenance]` stage line so all
  three wrappers emit the identical stage-prefix set.
- **Config-flag parity**: `--nextflow-config` is now accepted as an alias of
  `-c`/`--config`, so the same custom-config flag names work across all three
  wrappers.

### Added

- **Order-independent test suite.** `tests/conftest.py` now claims this skill's
  bare modules (`errors`, `schemas`, …) via a canonical-object cache — at
  collection and before each test — so a single `pytest` session that also
  collects the sibling wrapper suites (e.g. `make test`) no longer hits
  cross-skill module shadowing. Shared verbatim across all three wrappers.
- **Converged CLI presentation with the sibling wrappers.** Added a startup
  banner (`--no-banner`), `-v/--verbose`, and sarek-style staged progress logs
  (`[preflight]`/`[execute]`/`[outputs]`/`[report]`/`[done]`) plus a
  human-readable boxed error on stdout — while keeping the machine-readable JSON
  error on stderr (a robust fallback that survives a failed `result.json` write)
  and `result.json` on disk as the machine contract. `--verbose`/`--no-banner`
  are registered in the `clawbio.py` runner allowlist.
- **Robust `main()` entrypoint.** `main()` now accepts an explicit `argv` list
  (unit-testable without mutating `sys.argv`) and returns exit code `130` on
  `KeyboardInterrupt` (the SIGINT convention) instead of dumping a traceback.
  Matches `nfcore-sarek-wrapper` and `nfcore-rnaseq-wrapper`.

### Removed

- Dead `ErrorCode.FASTQ_NOT_READABLE` constant (no longer referenced after the
  readability-precheck removal) and the now-unused `import os` in
  `samplesheet_builder.py`.
