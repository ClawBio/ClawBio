# nf-core wrapper shared modules

The three implemented nf-core skills use `clawbio.common.nfcore` for mechanical work while keeping scientific policy and orchestration in their own skill directories. Their CLI aliases, entry scripts, defaults, result formats, and pinned upstream versions are unchanged.

| Module | Responsibility |
|---|---|
| errors | The structured SkillError every wrapper raises, with JSON-safe failure details; each skill keeps its own ErrorCode constants |
| imports | Bare-module isolation with an explicit skill directory |
| paths | URI preservation, path containment and the single /tmp predicate |
| cli | `--timeout-hours` parsing (finite, non-negative), timeout seconds, `--work-dir` resolution and stable manifest spelling |
| commands | Nextflow argv and the display form, including the environment the release needs |
| parameters | params.yaml text: fixed header, key order, one value per line; strict mapping parsing and pinned scalar-contract validation that skips the pipeline's `defaultIgnoreParams`, as nf-schema does |
| execution | Logs, process-group lifecycle (timeout, Ctrl+C and SIGTERM all stop Nextflow), bounded log reads, ordered failure hints |
| launcher | POSIX root-launcher signal forwarding and bounded graceful shutdown before a wrapper exits or its runner deadline expires |
| failures | result.json, stdout box and stderr JSON for SkillError and unexpected errors |
| runtime | Executable, version, memory and Docker CPU probes; the `--demo` offline refusal |
| source | Bounded Git probes, required-file inspection and literal manifest-version parsing |
| artifacts | JSON/CSV reading and writing, config basenames, atomic prepared-manifest writing and the safe error result.json writer |
| resources | Local Docker eligibility and the standalone replay capacity guard, which reads the recorded `executor.cpus`/`executor.memory` lines; wrappers supply profiles and choose resource budgets |
| bundle_runtime | Standard-library relocation: shell and YAML quoting, remap/verify, repair, stubs |
| bundle_emit | Embed the relocation runtime in a single-file remap_paths.py |

Shared modules never import a skill's bare `schemas`, `errors`, or orchestrator. Pass values explicitly. Keep reference rules, samplesheet columns, presets, source fallback/dirty rules, resource calculations, manifests, checksum target selection, scientific output ranking and downstream handoff local. Preserve each skill's validation order and warning channel.

The compatibility audit builds on the shared package from PR #623. There is no
competing `clawbio/common/nfcore.py` module or second execution/relocation runtime.
Scalar contracts stay beside each pinned skill and record the SHA-256 of the pinned
`nextflow_schema.json` and `nextflow.config` and that config's `defaultIgnoreParams`;
the shared parameter validator receives the selected contract explicitly. Prepared-manifest payloads and resource
ceilings remain pipeline policy; shared helpers only write the supplied checkpoint
or render a replay guard. Combined tests also check that collecting wrapper tests
preserves the real public `clawbio` package instead of substituting a partial stub.

Each skill's `executor.py` and `command_builder.py` are thin adapters that supply pipeline policy: timeout advice, scRNA's CellBender hint, profile composition, and `schemas.NEXTFLOW_ENV`. Pipeline pinning remains in each skill's schemas and compatibility JSON rather than a global registry.

## Root launcher shutdown

The root launcher uses `launcher.run_wrapper` for POSIX skills with a declared
`timeout_hours_flag`. SIGINT and SIGTERM reach the wrapper before the launcher
waits for cleanup; an explicit root timeout also sends SIGTERM before reporting
its timeout result. The wrapper then stops its detached Nextflow session using
the shared executor. A 30-second launcher grace period covers the executor's
10-second termination interval and output draining. Repeated signals cannot
interrupt cleanup. SIGINT remains a catchable `KeyboardInterrupt` for Python API
callers, SIGTERM exits with 143, and prior handlers are restored. Worker-thread
calls do not install signal handlers. Native Windows retains its existing path;
Nextflow's supported Windows route is WSL2/Linux.

## Pinned releases and Nextflow

Nextflow 26.04 parses with its strict v2 syntax by default. nf-core/sarek 3.8.1 and nf-core/scrnaseq 4.1.0 predate it and fail to parse there, so their `NEXTFLOW_ENV` sets `NXF_SYNTAX_PARSER=v1`; nf-core/rnaseq 3.26.0 needs nothing. The executor applies it over the caller's environment, the recorded command shows it, and the direct-replay `commands.sh` exports it. Never make a release run by editing Nextflow's cached clone: provenance would still name the release while different code ran.

## Portable bundles

`bundle_runtime.py` imports only Python's standard library. Source remappers import its functions inside one `SHARED NFCORE RUNTIME` marker block. Bundle generation replaces that block with the runtime source and binds `_nfcore_*` aliases for every name in `bundle_emit._EXPORTS`; adapters call the runtime only through those aliases. The generated script has no package dependency and remains executable with `python -I -S` on another machine.

Because the runtime is pasted into the adapter's namespace, the emitter refuses an adapter that binds a name the runtime defines (before the block) or calls (after it), unless it is the same standard-library import. Do not copy the source adapter into a bundle directly. Use `write_standalone_remapper`, or `render_standalone_remapper` after scRNA's canonical reference-key regeneration. Missing or duplicate runtime markers fail before writing.

Every value the remappers rewrite must read back exactly. `SHELL_WORD` and `unquote_shell_word` are the inverse of `format_output_value` and of `shlex.quote`, which writes the original commands.sh; `yaml_scalar_value` and `format_yaml_scalar` follow PyYAML, which writes params.yaml. Property tests pin both round trips with bash and `yaml.safe_load` as oracles. params.yaml is written unwrapped so a line-based edit always sees a whole value.

CSV backup naming, reference locations, the checksum policy (`checksums_cover_stubs`), the environment-stub name and remap CLI text remain local. RNA's commands.sh is assembled and written atomically, so an in-place replay cannot truncate the script bash is reading; its resume guard escapes the output directory for bash double quotes. Sarek/scRNA replay Nextflow directly and resolve their output directory from the bundle location. Old bundles do not require regeneration.

## Validation

Run the shared tests and all three wrapper suites together. CI does this on Python 3.11, 3.12 and 3.13, including changes confined to common code. Loading-order tests cover local schema ownership and shared exception identity. Tests retain each wrapper's different parameter fingerprint serialization and artifact conventions; scRNA's resume also accepts fingerprints of params.yaml written before it stopped wrapping long values.

The catalog derives CI evidence from every explicit skill argument in workflow pytest commands. Regenerate it when changing CI coverage, including a mixed-wrapper test step, even if SKILL.md metadata is unchanged. Its workflow reader ignores comments and quoted command examples.

Bundle tests exercise previews, backups, prefix boundaries, Unicode and line-separator characters, shell metacharacters, action order and standalone execution without site packages. Filesystem property tests disable timing deadlines because a busy server or network volume must not cause timing-only failures.

For distribution changes, build a wheel from the sdist, install outside the checkout, unset PYTHONPATH, exercise all pipeline --help routes, and generate a standalone remapper through the installed package. The emitter needs its runtime source in the wheel; the existing package build includes this Python resource.

Before shipping a refactor, run each pinned public demo through its wrapper against pristine pipeline clones (a fresh `NXF_HOME`, warmed sequentially before parallel runs), validate artifacts and checksums, resume, confirm drift is refused, replay in place, and replay a relocated bundle. Use separate output roots. Preserve the generating revision, runtime versions, commands, exit status, and evidence logs. An environment failure needs diagnosis before changes to wrapper logic.
