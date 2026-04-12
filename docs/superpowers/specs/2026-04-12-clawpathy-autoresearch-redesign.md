# Clawpathy-Autoresearch Redesign

**Date:** 2026-04-12
**Status:** Draft for review

## Goal

A task-agnostic framework that iteratively improves a `SKILL.md` artefact against a task-defined scorer. The framework does not assume what "improvement" means — the task does.

## Why a redesign

The previous version conflated three different loops (paper-reproduction, method-discovery, generalisation) and added scaffolding for each (critic, held-out split, prior-art, auditors). On a real test, none of the scaffolding earned its keep: the executor short-circuited by fetching the paper's Table 1 and returning the answer verbatim. The "anti-memorisation" penalty scanned `SKILL.md` for leaked answers but could not stop the executor from looking them up at runtime. The loop scored 0.000 without exercising the method it was supposed to evolve.

Root causes:

1. No per-task sandbox — the executor had the same tools on every task, including web access to source papers.
2. Framework baked in scoring assumptions (provenance, memorisation) that only made sense for one kind of task.
3. Held-out split pattern imported from ML without a training/test distinction that mapped onto the task.
4. Critic and auditor subagents never fired in practice (LLM backend not wired).

## Design

### Task layout

```
tasks/<task-name>/
  task.md                  # prompt describing what SKILL.md should do
  sandbox.yaml             # tools, data, network, timeout
  scorer.py                # score(skill_output, task_dir) -> (float, dict)
  data/                    # inputs the executor may read
  ground_truth.json        # optional, scorer reads if it wants
```

### Workspace layout

```
workspace/<run-id>/
  skill/SKILL.md           # artefact under improvement
  history.jsonl            # one line per iteration
  snapshots/<iter>/        # revert points
```

### Loop engine

`autoresearch.py` runs one loop:

1. Load task. Initialise empty `SKILL.md`.
2. For each iteration up to `max_iterations`:
   1. Dispatch **proposer** subagent (see below). Receive replacement `SKILL.md`.
   2. Snapshot current `SKILL.md`. Write the proposal.
   3. Dispatch **executor** subagent with the sandbox contract. Receive JSON output.
   4. Run `scorer(output, task_dir)`. Compare to best-so-far.
   5. If score improved, keep. Else, revert from snapshot.
   6. Append to `history.jsonl`.
3. Stop on `max_iterations`, `early_stop_n` consecutive non-improvements, or `target_score` hit.

No critic. No prior-art check. No held-out split. These can be added later as opt-in modules if a real task shows they help.

### Sandbox

`sandbox.yaml` format:

```yaml
allowed_tools: [Read, Bash, Grep]
network:
  allowed: false           # or list of hostnames
data_dir: ./data
python_packages: [pandas, numpy]
timeout_seconds: 600
```

Enforcement is prompt-level plus tool-whitelist in the Agent dispatch call. The loop parses the executor subagent's tool-use log and flags any call outside the whitelist. Not hermetic — a determined subagent could work around the prompt — but enough that accidental and opportunistic violations leave a trail. Tasks needing hard isolation run their executor in a container inside the scorer.

### Proposer contract

Dispatched fresh each iteration. Input:

- Task prompt (verbatim `task.md`)
- Current `SKILL.md`
- Last iteration's score breakdown (if any)
- Last three history entries, each one line

Output: a full replacement `SKILL.md`. Not a diff. The loop owns snapshots for revert.

### Executor contract

Dispatched fresh each iteration. Input:

- `SKILL.md` content (verbatim)
- `sandbox.yaml` content (verbatim, framed as a contract)
- Tool whitelist (enforced in the Agent dispatch call)
- Instruction: follow `SKILL.md`, produce the output format `SKILL.md` specifies, return JSON only.

The executor does not see `task.md`, `scorer.py`, or `ground_truth.json`. It sees only what the skill tells it to do.

### Scorer contract

Task author writes:

```python
def score(skill_output: dict, task_dir: Path) -> tuple[float, dict]:
    """Lower is better. Returns (total, breakdown_dict)."""
```

- `total` is a float, any range, lower is better.
- `breakdown` is a flat dict of named sub-errors the proposer can learn from.
- Scorer may read anything in `task_dir`.
- Scorer is trusted code. The framework does not sandbox it.

The framework provides no scoring primitives. Provenance checks, memorisation penalties, and the like live in the task's scorer if the task wants them.

### History format

`history.jsonl`, one line per iteration:

```json
{
  "iter": 3,
  "score": 0.025,
  "kept": true,
  "breakdown": {"rsid_error": 0.0, "gene_error": 0.5},
  "proposer_summary": "widened gene annotation window to 500kb"
}
```

### Stopping conditions

- `max_iterations` reached
- `early_stop_n` consecutive non-improvements
- `score <= target_score` (optional, task-defined in `task.md` front-matter)

## What gets deleted

- `skills/clawpathy-autoresearch/critic.py` and `tests/test_critic.py`
- `prior-art-prompt.md`, `ground-truth-auditor-prompt.md`, `reproduction-auditor-prompt.md`, `clarifier-prompt.md`
- `heldout_ground_truth` field on `Workspace` and its load path
- `enable_critic` wiring in `autoresearch.py`
- Any scoring primitives in the framework itself

## What survives

- Workspace dataclass and loader (flattened — no dev/heldout)
- Snapshot and revert
- Proposer subagent dispatch pattern (already refactored away from `claude -p`)
- History logging

## Testing

- Unit: workspace load/save, history append, snapshot/revert, scorer-contract shape check.
- Integration: a mock "hello world" task where the scorer returns 0 if `skill_output["text"] == "hello"`. Real subagent dispatch, trivial executor. Loop converges in one or two iterations.
- No mock-LLM tests. If the LLM isn't real, the test tells you nothing.

## Out of scope

- Containerised sandboxes
- Multi-proposer ensembles
- Critic / reviewer agents
- Held-out generalisation splits
- Built-in scoring primitives (provenance, memorisation, etc)
- Any task-specific logic in framework code

Each of these is addable later when a real task produces evidence the framework needs it.
