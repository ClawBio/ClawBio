# Skill Intent Descriptors

ClawBio chat adapters use `clawbio.skill_intents.plan_skill_intent()` to map
natural-language requests into deterministic `clawbio.py run ...` plans. A
skill can publish optional routing metadata beside its `SKILL.md` without
adding chat-platform-specific code.

## Location

Place one JSON file in the skill directory:

- `skills/<skill>/INTENTS.json` preferred
- `skills/<skill>/skill_intents.json` supported alias

RoboTerri and the Discord bot discover descriptors both for registered skills
and by scanning `skills/*/INTENTS.json`, so descriptor-only symlinked or copied
skill directories can publish routing metadata before they are added to the
main `SKILLS` registry. Descriptor-only skills are not exposed in executable
chat tool enums and route matches return `needs_registration` until every
`skill_run.skill` target is registered in `clawbio.py`'s `SKILLS` dictionary.
If neither file exists, chat adapters fall back to the legacy skill/mode
behavior. Demo mode is only planned when the raw user text explicitly asks for a
demo, example, synthetic data, or sample data, or when the user confirms an
already proposed demo with text such as "yes" or "go ahead".

## Schema

```json
{
  "schema": "clawbio.skill_intents.v1",
  "skill": "gentle-cloning",
  "aliases": ["gentle", "cloning"],
  "routes": [
    {
      "intent_id": "runtime_version",
      "description": "Check the installed runtime version for this skill backend.",
      "trigger_terms": ["version", "runtime version", "installed version", "status"],
      "demo_policy": "never_unless_explicit",
      "plan": [
        {
          "kind": "skill_run",
          "skill": "gentle-cloning",
          "input": "examples/request_runtime_version.json"
        }
      ]
    }
  ]
}
```

Required top-level fields:

- `schema`: must be `clawbio.skill_intents.v1`
- `skill`: ClawBio skill name used with `python clawbio.py run <skill>`
- `routes`: list of intent routes

Optional top-level fields:

- `aliases`: skill names or terms users may type

Route fields:

- `intent_id`: stable machine-readable route identifier
- `description`: human-readable explanation for logs and review
- `trigger_terms`: words or phrases matched deterministically against raw user text
- `aliases`: extra route trigger terms
- `demo_policy`: use `never_unless_explicit` by default; use `only_when_explicit` for routes that should only match demo requests
- `requires_confirmation`: optional route-level confirmation gate
- `plan`: one or more execution steps

Plan step fields:

- `kind`: currently only `skill_run`
- `skill`: optional override; defaults to the descriptor `skill`
- `input`: optional request JSON or data file path, resolved relative to the skill directory
- `demo`: optional boolean; `true` only becomes `--demo` when the user explicitly asks for demo/example/synthetic/sample data
- `args`: optional array of literal CLI arguments; no shell is used
- `output`: optional path, resolved relative to the skill directory
- `confirmation`: optional object, for example `{"required": true, "reason": "Writes cached backend state."}`

The descriptor is data only. It cannot run arbitrary code, define shell
commands, or call chat-platform APIs.

## Planner Output

`plan_skill_intent(...)` returns a structured `SkillExecutionPlan` with:

- `status`: `planned`, `needs_input`, `needs_confirmation`, or `needs_registration`
- `raw_user_text` and `raw_user_text_sha256`
- `skill`, `intent_id`, `confidence`, and `reason`
- `matched_route`: route id, matched terms, score, and demo policy
- `executions`: command argv arrays and output/input paths

Adapters log the raw text hash or preview, selected skill, selected intent,
matched route, final command(s), and output bundle path(s).

## GENtle Descriptor Guidance

GENtle should add `skills/gentle-cloning/INTENTS.json` using the schema above.
For runtime/status requests, provide a route whose `trigger_terms` include
`version`, `runtime version`, `installed version`, and `status`, and whose plan
uses an `input` JSON request such as `examples/request_runtime_version.json`.
Multi-step workflows should list each `skill_run` step in order and mark
mutating or expensive steps with `confirmation.required: true`.
