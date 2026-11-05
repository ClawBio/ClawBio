# Using the `clawbio-flow` CLI

This is the **driver guide** for any ClawBio 2.0 skill. Each skill's `SKILL.md`
links here for the lifecycle and invariants; the SKILL.md itself only describes
*that* skill's DAG, inputs, and outputs.

If you have run a skill before, `clawbio-flow --help` and `clawbio-flow status` are
usually enough. Read this once if you haven't.

## The contract

> **The CLI orchestrates. You execute. Neither side improvises.**

A skill is a DAG declared in `pipeline.yaml`. The CLI tracks state in
`.clawbio-run/state.json` and refuses any transition that would skip a prereq,
fabricate an output, or end a run without a reproducibility bundle.

You drive the run with five verbs.

## The five verbs

```
start  →  next  →  open  →  done  →  bundle
              ↑_____________|
                 (loop until next prints nothing)
```

| Verb | What it does | Refuses if |
| --- | --- | --- |
| `start <skill-dir> --input k=v ...` | Loads `pipeline.yaml`, validates it, writes `.clawbio-run/state.json`. | A run already exists; required inputs missing; unknown input given; YAML invalid. |
| `next` | Prints the id of the next ready step (prereqs done, required emits present). | — (prints nothing when no step is ready). |
| `open <step-id>` | Renders `{placeholders}` in `command_template`, prints the command, marks the step `in_progress`. | Step is blocked, already in progress, or already done. |
| `done <step-id> --output k=path ...` | Hashes each declared output (sha256), records absolute paths, recomputes downstream blocks. | Step isn't in progress; output keys don't match the step's `emits` exactly; an output file is missing. |
| `skip <step-id> --reason "..."` | Marks an `optional: true` step `skipped` instead of run; recorded on the events log. | Step isn't `optional`; step isn't `pending` (blocked steps must be unblocked first; done/in-progress/already-skipped steps can't be re-skipped). |
| `bundle` | Emits `commands.sh` (execution order), `checksums.sha256`, `environment.yml`. | Any step isn't `done` or `skipped`. |

`status` is the seventh verb but doesn't change state — it prints a
self-sufficient human-readable summary, written so an agent that just woke up
from `/compact` or `/clear` can fully re-orient from `clawbio-flow status` alone.

## Worked example

The repo ships [`tests/workflow/fixtures/variant_annotation/`](../../tests/workflow/fixtures/variant_annotation/) (placeholder `cp`/`echo`
commands — no bcftools/VEP required), so this runs anywhere.

```bash
# 1. Initialize
clawbio-flow start tests/workflow/fixtures/variant_annotation \
  --input vcf_path=/tmp/x.vcf \
  --input reference=/tmp/ref.fa

# 2. Ask what's ready
clawbio-flow next                # -> normalize

# 3. Open it — prints the rendered command, marks in_progress
clawbio-flow open normalize
# (the CLI prints, e.g.:
#   cp /tmp/x.vcf /…/.clawbio-run/work/normalize/normalized.vcf
#  you — the agent — run that command yourself)

# 4. Report outputs back
clawbio-flow done normalize --output normalized_vcf=.clawbio-run/work/normalize/normalized.vcf

# 5. Loop
clawbio-flow next                # -> annotate
clawbio-flow open annotate
# (run the printed command)
clawbio-flow done annotate --output annotated_vcf=.clawbio-run/work/annotate/annotated.vcf

# 6. Seal the run
clawbio-flow bundle              # writes commands.sh, checksums.sha256, environment.yml
```

`bundle` is mandatory. If you skip it, the Claude Code `Stop` hook
([hooks/stop-bundle.json](../../hooks/stop-bundle.json)) refuses to let the session end.

## Invariants worth knowing

- **`open` is a refusal point, not a permission slip.** If `clawbio-flow open X`
  prints "blocked by …", trust it — do not run X's command anyway. The block
  exists because a prereq's output is missing or unhashed.
- **`done` is the only way outputs enter the run.** The CLI never inspects the
  filesystem to *discover* outputs; you must declare each emit explicitly with
  `--output emit=path`. Missing files or mismatched keys fail loudly.
- **Paths in `state.json` are absolute.** Relative `--output` paths get
  resolved against the cwd at the time of `done`, then stored absolute so the
  bundle is cwd-independent.
- **The event log, not the step list, defines bundle order.** `commands.sh`
  reflects the order steps actually completed, not the YAML declaration order
  — important when sibling steps run in parallel.
- **`pipeline.yaml` is validated once, at `start`.** Unknown placeholders,
  forward prereqs, duplicate ids, undeclared `requires` keys all fail then.
  Once `start` succeeds, the shape is trusted.
- **`/compact` and `/clear` do not destroy progress.** Everything load-bearing
  lives in `.clawbio-run/state.json`. After a context wipe, `clawbio-flow status`
  alone re-orients an agent.

## Common refusals and what they mean

| Message | Cause | Fix |
| --- | --- | --- |
| `blocked by [stepX, emit:foo]` | A prereq isn't `done` or its emit isn't hashed | Finish stepX with `done --output foo=…` |
| `already done` | Re-running a completed step | Don't. Re-runs aren't a feature — they corrupt reproducibility. If you need to redo, start a fresh run. |
| `missing required inputs: [vcf_path]` | `start` called without an input declared `required: true` | Pass `--input vcf_path=…` |
| `unknown input 'bogus'` | `--input` key not declared in `pipeline.yaml: inputs:` | Drop it, or add it to the skill's inputs. |
| `not in emits: bonus` | `done --output` key not in the step's `emits` list | Use only the keys the step declares. |
| `not a file: /path` | `done` got a path that doesn't exist | Run the rendered command first; verify the output. |
| `refusing to overwrite .clawbio-run/` | `start` called inside an existing run | `rm -rf .clawbio-run/` to restart from scratch, or `cd` elsewhere. |
| `bundle: stepX is not done` | `bundle` before the run is complete | Finish the remaining steps. |
| `stepX is not optional; cannot skip` | `clawbio-flow skip` on a step without `optional: true` | Run it with `open`/`done` like any other step. |
| `stepX is not pending; cannot skip` | `clawbio-flow skip` on a blocked, done, in-progress, or already-skipped step | Unblock it first, or don't re-skip. |

## What lives where

```
<skill-dir>/
  pipeline.yaml        # the DAG declaration — validated at start time
  SKILL.md             # this skill's DAG, inputs, outputs
  scripts/             # whatever the command_templates invoke
  data/                # bundled fixtures, if any

<cwd>/.clawbio-run/
  state.json           # the source of truth — read this if anything is unclear
  events.log           # append-only transition log
  work/<step-id>/      # per-step output directory (the {out} placeholder)
  commands.sh          # emitted by `bundle`
  checksums.sha256     # emitted by `bundle`
  environment.yml      # emitted by `bundle`
```

## See also

- [pipeline-schema.md](pipeline-schema.md) — authoritative `pipeline.yaml` reference (writing skills, not running them).
- [pipeline.schema.json](pipeline.schema.json) — machine-readable draft-07 JSON Schema.
- [hooks/stop-bundle.json](../../hooks/stop-bundle.json) — Claude Code hook that enforces `bundle` before session end.
