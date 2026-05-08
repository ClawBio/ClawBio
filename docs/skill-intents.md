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
skill directories can publish routing metadata. A descriptor-only skill becomes
executable when it provides a safe local Python entrypoint, either through
top-level `entrypoint`/`script`, `execution.entrypoint`, or a conventional file
such as `<skill_name_with_underscores>.py` in the skill directory. Descriptors
without an executable entrypoint are still discoverable, but route matches
return `needs_registration`. If neither file exists, chat adapters fall back to
the legacy skill/mode behavior. Demo mode is only planned when the raw user text
explicitly asks for a demo, example, synthetic data, or sample data, or when the
user confirms an already proposed demo with text such as "yes" or "go ahead".

Descriptor paths are confined to the resolved skill directory. `entrypoint`,
`script`, plan-step `input`, and plan-step `output` may be absolute only if they
resolve inside that directory; relative paths containing `..` that escape the
skill directory are rejected.

## Schema

```json
{
  "schema": "clawbio.skill_intents.v1",
  "skill": "gentle-cloning",
  "entrypoint": "gentle_cloning.py",
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
- `entrypoint` or `script`: local Python skill script, resolved relative to the skill directory
- `execution.entrypoint`: nested form of `entrypoint`

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
- `input_template`: optional request JSON object filled from extracted slots and materialized as a temporary input file
- `slots`: optional slot extraction specs for `input_template`; supported fields include `pattern`, `choices`, `aliases`, `default`, and `required`
- `demo`: optional boolean; `true` only becomes `--demo` when the user explicitly asks for demo/example/synthetic/sample data
- `args`: optional array of literal CLI arguments; no shell is used
- `output`: optional path, resolved relative to the skill directory
- `confirmation`: optional object, for example `{"required": true, "reason": "Writes cached backend state."}`

The descriptor is data only. The planner never runs arbitrary code, invokes a
shell, or calls chat-platform APIs from metadata. A descriptor may still create
ordinary request JSON fields such as `shell_line` for a registered skill that
validates and interprets those fields itself.

## Validation And Safety

Descriptors are validated before use:

- `schema` must equal `clawbio.skill_intents.v1`
- `skill`, `intent_id`, step `skill`, and slot names must match conservative identifier patterns
- `routes`, `plan`, `aliases`, `trigger_terms`, and `args` have item and length caps
- slot regexes must compile and must stay short
- `input`, `output`, `entrypoint`, and `script` must stay within the skill directory
- descriptor prompt text is sanitized as untrusted labels and capped before entering Roboterri's system prompt

Plan-step `args` are deliberately narrow. Descriptor metadata cannot grant
itself new CLI privileges. Extra args are emitted only when the target skill is
already registered with a static ClawBio `allowed_extra_flags` allow-list, and
blocked core/sensitive flags such as `--input`, `--output`, `--profile`, and
`--demo` are never accepted through descriptor `args`. Sensitive path-oriented
flags and absolute, traversal, or path-separator values are also rejected there;
use `input` or `input_template` for request data instead.

## Chat Adapter Behavior

Executable descriptor-backed skills are added to the chat tool schema together
with their aliases, intent ids, and leading trigger terms. When a user names a
descriptor skill or alias and asks a matching question, such as a version,
status, runtime, guide, isoforms, 2D gel, or descriptor-specific analysis
request, RoboTerri and the Discord bot should call `clawbio` before answering.
If the same name also refers to public upstream software, answers should keep
the two meanings separate: public/latest upstream information is not the same
as the locally installed ClawBio runtime or rewrite reported by the descriptor
route.

Parameterized request example:

```json
{
  "kind": "skill_run",
  "skill": "gentle-cloning",
  "input_template": {
    "mode": "gene-protein-2d-gel",
    "gene_symbol": "{gene_symbol}",
    "species": "{species}",
    "source": "{source}"
  },
  "slots": {
    "gene_symbol": {"pattern": "\\b([A-Z][A-Z0-9]{2,15})\\b"},
    "species": {
      "aliases": {"human": "homo_sapiens", "homo sapiens": "homo_sapiens"},
      "default": "homo_sapiens"
    },
    "source": {"choices": ["ensembl", "refseq", "uniprot"], "default": "ensembl"}
  }
}
```

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
For parameterized protein-gel requests, use an `input_template` with slots for
`gene_symbol`, `species`, and `source`, as shown above. The generated temporary
JSON request is passed to `gentle-cloning` via `--input`.

GENtle can route Telegram guide requests with a shell request payload:

```json
{
  "intent_id": "isoforms_guide",
  "description": "Show a GENtle isoforms guide for a gene.",
  "trigger_terms": ["isoforms guide", "isoforms", "guide"],
  "plan": [
    {
      "kind": "skill_run",
      "skill": "gentle-cloning",
      "input_template": {
        "mode": "shell",
        "shell_line": "services guide --channel telegram --section isoforms --gene {gene_symbol}"
      },
      "slots": {
        "gene_symbol": {"pattern": "\\b([A-Z][A-Z0-9]{2,15})\\b"}
      }
    }
  ]
}
```

And a 2D protein gel request with this payload shape:

```json
{
  "mode": "gene-protein-2d-gel",
  "source": "ensembl",
  "species": "homo_sapiens",
  "gene_symbol": "{gene_symbol}",
  "state_path": ".gentle_state.json",
  "ladders": ["Protein Ladder 10-100 kDa"],
  "timeout_secs": 600
}
```
