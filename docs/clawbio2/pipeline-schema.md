# `pipeline.yaml` schema (v1)

A ClawBio 2.0 skill ships a `pipeline.yaml` next to its `SKILL.md`. The `clawbio-flow`
CLI loads this file at `clawbio-flow start <skill-path>` to seed the run state.

A machine-readable JSON Schema (draft-07) is published at
[`pipeline.schema.json`](pipeline.schema.json) for editor integration and
external tooling. It pins **shape** (required fields, types, allowed
properties); **cross-step semantics** (unknown prereqs, forward references,
unresolved `{var}` placeholders, `requires` keys not in a prior `emits`) are
enforced by the python loader in [`clawbio/workflow/pipeline.py`](../../clawbio/workflow/pipeline.py)
and cannot be expressed in JSON Schema. Both checks run at `clawbio-flow start`.

## Top-level fields

| Field | Type | Required | Notes |
|---|---|---|---|
| `skill` | string | yes | Slug identifying the skill. |
| `version` | int | yes | Schema version. Only `1` accepted today. |
| `inputs` | list | no | Inputs the user/agent supplies at `clawbio-flow start`. |
| `steps` | list | yes | Ordered DAG of pipeline steps. Must be non-empty. |

## `inputs[]`

| Field | Type | Required | Notes |
|---|---|---|---|
| `name` | string | yes | Key used by `--input name=value`. |
| `required` | bool | no (default `true`) | If `true`, `start` fails when missing. |

## `steps[]`

| Field | Type | Required | Notes |
|---|---|---|---|
| `id` | string | yes | Unique step slug. |
| `summary` | string | yes | One-line description. |
| `command_template` | string | yes | Shell snippet with `{var}` placeholders. |
| `requires_inputs` | list[string] | no | Inputs that must be present at run start. |
| `prereqs` | list[string] | no | Step IDs that must be `done` before this step is unblocked. |
| `requires` | mapping | no | Emit-key gates. Each key must be in some prior step's `emits`. |
| `emits` | list[string] | no | Output keys this step produces. Matches `clawbio-flow done --output`. |
| `optional` | bool | no (default `false`) | If `true`, this step may be skipped with `clawbio-flow skip <id> --reason ...` instead of run. |

### Optional steps

A step marked `optional: true` may be skipped instead of executed, via
`clawbio-flow skip <step-id> --reason "..."`. `skip` is only valid when the step is
`pending` (its own prereqs are met but it hasn't been opened) — a blocked step
must be unblocked first. `skipped` is a terminal state alongside `done`:
`clawbio-flow bundle` accepts a run where every step is `done` *or* `skipped`, and
`commands.sh` lists skipped steps in a trailing `# skipped: <id> (reason: ...)`
comment block instead of a command invocation.

Downstream steps: a skipped step satisfies plain `prereqs` (DAG ordering — it
ran its course) but never satisfies `requires: {emit: present}` gates (it
produced no output). A downstream step gated on a skipped step's emit blocks
forever; the agent must skip it too. This is deliberate — see the "strict by default" tradeoff: a skipped step produces no output, so gates on its emits stay closed.

## Placeholders in `command_template`

`{var}` is substituted with, in order of priority:

1. `{out}` — the run's output directory (`.clawbio-run/out`).
2. `{skill_dir}` — the absolute path of the skill (the directory containing `pipeline.yaml`). Use this to invoke scripts/data bundled with the skill.
3. A declared input (from `inputs[]`).
4. An emitted output from any prior step.

Validation at load time fails if any placeholder is unresolved.

## Forbidden in v1

- Conditionals, loops, `if:`, `when:` clauses
- Parallel/fan-out semantics
- Nested pipelines
- Templating beyond simple `{var}` substitution

If a step needs branching logic, write a Python script and make that script the
step body.

## Example

See [`tests/workflow/fixtures/variant_annotation/pipeline.yaml`](../../tests/workflow/fixtures/variant_annotation/pipeline.yaml).
