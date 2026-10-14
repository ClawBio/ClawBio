# nf-core wrapper shared modules

The three implemented nf-core skills use `clawbio.common.nfcore` for mechanical work while keeping scientific policy and orchestration in their existing skill directories. Their CLI aliases, entry scripts, defaults, result formats, and pinned upstream versions are unchanged.

| Module | Responsibility |
|---|---|
| errors | Shared structured SkillError; each skill retains its own ErrorCode constants |
| imports | Bare-module isolation with an explicit skill directory |
| paths | URI preservation and path containment |
| commands | Nextflow argv and display quoting; callers compose profiles |
| parameters | Exact YAML headers, ordering and file writing |
| execution | Logs, subprocess lifecycle and shared environment hints; caller-supplied failure policy |
| runtime | Raw executable, version and memory evidence |
| source | Bounded Git probes and required-file inspection |
| artifacts | JSON/CSV writing and caller-named config copies |
| bundle_runtime | Standard-library relocation, quoting, ordered actions and compatible repair mechanics |
| bundle_emit | Embed the relocation runtime in a single-file remap_paths.py |

Shared modules never import a skill's bare `schemas`, `errors`, or orchestrator. Pass values explicitly. Keep reference rules, samplesheet columns, presets, source fallback/dirty rules, resource calculations, manifests, checksum target selection, scientific output ranking and downstream handoff local. Preserve each skill's validation order and warning channel.

The local `errors.py`, `_isolated_imports.py`, `command_builder.py` and `executor.py` remain compatibility adapters. Sarek/RNA and scRNA have different interruption/Windows cleanup policies; these are intentional adapter inputs. Pipeline pinning remains in the existing schemas and compatibility JSON, rather than a new global registry.

## Portable bundles

`bundle_runtime.py` imports only Python's standard library. Source remappers import its functions inside one `SHARED NFCORE RUNTIME` marker block. Bundle generation replaces that block with the runtime source and binds private aliases before the local adapter definitions. The generated script has no package dependency and remains executable with `python -I -S` on another machine.

Do not copy the source adapter into a bundle directly. Use `write_standalone_remapper`, or render after scRNA's canonical reference-key regeneration. Missing or duplicate runtime markers fail before writing. Keep future imports in the adapter header. CSV backup naming, reference locations, repair/checksum policies and remap CLI text remain local.

RNA's commands.sh is still assembled and written atomically, so an in-place replay cannot truncate the script bash is reading. Sarek/scRNA retain their direct Nextflow replay mode. Old bundles do not require regeneration.

## Validation

Run the shared tests and all three wrapper suites together. CI does this on Python 3.11, 3.12 and 3.13, including changes confined to common code. Six loading-order tests cover local schema ownership and shared exception identity. Tests retain each wrapper's different parameter fingerprint serialization and artifact conventions.

The catalog derives CI evidence from every explicit skill argument in workflow pytest commands. Regenerate it when changing CI coverage, including a mixed-wrapper test step, even if SKILL.md metadata is unchanged. Its workflow reader ignores comments and quoted command examples.

Bundle tests exercise previews, backups, prefix boundaries, Unicode, literal shell characters, action order and standalone execution without site packages. Filesystem property tests disable timing deadlines because a busy server or network volume must not cause timing-only failures.

For distribution changes, build a wheel from the sdist, install outside the checkout, unset PYTHONPATH, exercise all pipeline --help routes, and generate a standalone remapper through the installed package. The emitter needs its runtime source in the wheel; the existing package build includes this Python resource.

Before shipping a refactor, run each pinned public demo through its wrapper, validate artifacts and checksums, replay in place, and replay a relocated bundle. Use separate output roots. Preserve the generating revision, runtime versions, commands, exit status, and evidence logs. An environment failure needs diagnosis before changes to wrapper logic.
