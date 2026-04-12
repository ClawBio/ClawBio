# Clawpathy-Autoresearch Redesign Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Rebuild clawpathy-autoresearch as a task-agnostic loop that iteratively improves a SKILL.md artefact against a task-defined scorer, with per-task sandbox enforcement.

**Architecture:** A thin loop engine (proposer → executor → scorer → keep/revert) plus a pluggable dispatcher callback so the outer agent (or a future SDK integration) drives real subagent calls. No critic, no held-out split, no built-in scoring primitives. Tasks own everything task-specific.

**Tech Stack:** Python 3.11+, pytest, PyYAML (for sandbox.yaml), stdlib json/subprocess. Follows ClawBio TDD rules.

**Spec:** `docs/superpowers/specs/2026-04-12-clawpathy-autoresearch-redesign.md`

---

## File Map

**New:**
- `skills/clawpathy-autoresearch/sandbox.py` — loads and validates `sandbox.yaml`
- `skills/clawpathy-autoresearch/dispatcher.py` — pluggable subagent-call callback interface
- `skills/clawpathy-autoresearch/executor.py` — builds executor prompt, parses JSON return
- `skills/clawpathy-autoresearch/loop.py` — the iteration loop engine
- `skills/clawpathy-autoresearch/tasks/hello-world/` — integration-test task
- `skills/clawpathy-autoresearch/tests/test_sandbox.py`
- `skills/clawpathy-autoresearch/tests/test_dispatcher.py`
- `skills/clawpathy-autoresearch/tests/test_executor.py`
- `skills/clawpathy-autoresearch/tests/test_loop.py`

**Modified:**
- `skills/clawpathy-autoresearch/workspace.py` — flatten (remove heldout), add sandbox_path
- `skills/clawpathy-autoresearch/proposer.py` — rewrite around dispatcher
- `skills/clawpathy-autoresearch/autoresearch.py` — rewrite as CLI entry calling `loop.run`

**Deleted:**
- `skills/clawpathy-autoresearch/critic.py`
- `skills/clawpathy-autoresearch/tests/test_critic.py`
- `skills/clawpathy-autoresearch/prior-art-prompt.md`
- `skills/clawpathy-autoresearch/ground-truth-auditor-prompt.md`
- `skills/clawpathy-autoresearch/reproduction-auditor-prompt.md`
- `skills/clawpathy-autoresearch/clarifier-prompt.md`
- `skills/clawpathy-autoresearch/final_judge.py` (audit-era artefact)

---

### Task 1: Delete obsolete files

**Files:**
- Delete: `skills/clawpathy-autoresearch/critic.py`
- Delete: `skills/clawpathy-autoresearch/tests/test_critic.py`
- Delete: `skills/clawpathy-autoresearch/prior-art-prompt.md`
- Delete: `skills/clawpathy-autoresearch/ground-truth-auditor-prompt.md`
- Delete: `skills/clawpathy-autoresearch/reproduction-auditor-prompt.md`
- Delete: `skills/clawpathy-autoresearch/clarifier-prompt.md`
- Delete: `skills/clawpathy-autoresearch/final_judge.py`

- [ ] **Step 1: Delete the files**

```bash
rm skills/clawpathy-autoresearch/critic.py
rm skills/clawpathy-autoresearch/tests/test_critic.py
rm skills/clawpathy-autoresearch/prior-art-prompt.md
rm skills/clawpathy-autoresearch/ground-truth-auditor-prompt.md
rm skills/clawpathy-autoresearch/reproduction-auditor-prompt.md
rm skills/clawpathy-autoresearch/clarifier-prompt.md
rm skills/clawpathy-autoresearch/final_judge.py
```

- [ ] **Step 2: Remove references from autoresearch.py**

Open `skills/clawpathy-autoresearch/autoresearch.py`, remove any `from .critic import ...`, `from .final_judge import ...`, `self.enable_critic`, and any block that calls `review_edit` or the deleted prompt files. Leave a minimal stub `def run(*args, **kwargs): raise NotImplementedError("rewritten in Task 8")` so the module still imports.

```python
"""Stub. The loop engine is rewritten in Task 8."""
def run(*args, **kwargs):
    raise NotImplementedError("rewritten in Task 8")
```

- [ ] **Step 3: Run existing test suite to find other breakage**

Run: `uv run pytest skills/clawpathy-autoresearch/tests/ -x`
Expected: any test that imported `critic` or `final_judge` fails to collect. Note which, delete those tests in Step 4 if they only tested deleted modules.

- [ ] **Step 4: Remove any remaining tests that reference deleted modules**

For each failing collection, open the test file. If the file exists only to test a deleted module, delete it. If it tests surviving code but imports a deleted module, remove just the offending import and test function.

- [ ] **Step 5: Rerun suite**

Run: `uv run pytest skills/clawpathy-autoresearch/tests/ -v`
Expected: collection succeeds, remaining tests pass or fail on real (non-import) reasons.

- [ ] **Step 6: Commit**

```bash
git add -u skills/clawpathy-autoresearch/
git commit -m "refactor(autoresearch): delete obsolete critic/auditor/judge code"
```

---

### Task 2: Flatten Workspace (remove heldout, add sandbox_path)

**Files:**
- Modify: `skills/clawpathy-autoresearch/workspace.py`
- Modify: `skills/clawpathy-autoresearch/tests/test_workspace.py`

- [ ] **Step 1: Write failing test for flat ground truth with sandbox_path**

Add to `skills/clawpathy-autoresearch/tests/test_workspace.py`:

```python
def test_workspace_has_sandbox_path(tmp_path):
    ws_dir = tmp_path / "ws"
    (ws_dir / "skill").mkdir(parents=True)
    (ws_dir / "skill" / "SKILL.md").write_text("# skill\n")
    (ws_dir / "task.json").write_text('{"name":"t","max_iterations":5,"early_stop_n":2}')
    (ws_dir / "ground_truth.json").write_text('{"targets":[]}')
    (ws_dir / "scorer.py").write_text("def score(o, d): return (0.0, {})\n")
    (ws_dir / "sandbox.yaml").write_text("allowed_tools: [Read]\n")
    from skills.clawpathy_autoresearch.workspace import load_workspace
    ws = load_workspace(ws_dir)
    assert ws.sandbox_path == ws_dir / "sandbox.yaml"
    assert ws.ground_truth == {"targets": []}
    assert not hasattr(ws, "heldout_ground_truth")
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest skills/clawpathy-autoresearch/tests/test_workspace.py::test_workspace_has_sandbox_path -v`
Expected: FAIL — `heldout_ground_truth` still on Workspace, or `sandbox_path` missing.

- [ ] **Step 3: Rewrite workspace.py**

Replace entire file with:

```python
"""Workspace loader for clawpathy-autoresearch.

A workspace is a self-contained directory with task config, ground truth,
scorer, sandbox spec, and the SKILL.md being optimised.
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any


@dataclass
class Workspace:
    name: str
    description: str
    max_iterations: int
    early_stop_n: int
    target_score: float | None
    ground_truth: dict[str, Any]
    workspace_dir: Path
    skill_dir: Path
    scorer_path: Path
    sandbox_path: Path
    data_dir: Path


def validate_workspace(workspace_dir: Path) -> list[str]:
    workspace_dir = Path(workspace_dir)
    errors: list[str] = []
    if not workspace_dir.exists():
        return [f"Workspace directory does not exist: {workspace_dir}"]
    required = [
        ("task.json", workspace_dir / "task.json"),
        ("ground_truth.json", workspace_dir / "ground_truth.json"),
        ("scorer.py", workspace_dir / "scorer.py"),
        ("sandbox.yaml", workspace_dir / "sandbox.yaml"),
        ("skill/SKILL.md", workspace_dir / "skill" / "SKILL.md"),
    ]
    for label, path in required:
        if not path.exists():
            errors.append(f"Missing required file: {label}")
    return errors


def load_workspace(workspace_dir: Path) -> Workspace:
    workspace_dir = Path(workspace_dir)
    if not workspace_dir.exists():
        raise FileNotFoundError(f"Workspace not found: {workspace_dir}")
    errors = validate_workspace(workspace_dir)
    if errors:
        raise FileNotFoundError(f"Invalid workspace: {'; '.join(errors)}")
    task_data = json.loads((workspace_dir / "task.json").read_text())
    ground_truth = json.loads((workspace_dir / "ground_truth.json").read_text())
    return Workspace(
        name=task_data["name"],
        description=task_data.get("description", ""),
        max_iterations=task_data.get("max_iterations", 20),
        early_stop_n=task_data.get("early_stop_n", 5),
        target_score=task_data.get("target_score"),
        ground_truth=ground_truth,
        workspace_dir=workspace_dir,
        skill_dir=workspace_dir / "skill",
        scorer_path=workspace_dir / "scorer.py",
        sandbox_path=workspace_dir / "sandbox.yaml",
        data_dir=workspace_dir / "data",
    )
```

- [ ] **Step 4: Delete dev/heldout split tests**

In `tests/test_workspace.py`, delete any test named `test_ground_truth_split_dev_heldout` or similar that asserts on `heldout_ground_truth`.

- [ ] **Step 5: Run tests**

Run: `uv run pytest skills/clawpathy-autoresearch/tests/test_workspace.py -v`
Expected: PASS on all remaining tests, including the new one.

- [ ] **Step 6: Commit**

```bash
git add skills/clawpathy-autoresearch/workspace.py skills/clawpathy-autoresearch/tests/test_workspace.py
git commit -m "refactor(autoresearch): flatten workspace, add sandbox_path"
```

---

### Task 3: Sandbox loader

**Files:**
- Create: `skills/clawpathy-autoresearch/sandbox.py`
- Create: `skills/clawpathy-autoresearch/tests/test_sandbox.py`

- [ ] **Step 1: Write failing test**

Create `skills/clawpathy-autoresearch/tests/test_sandbox.py`:

```python
from pathlib import Path
import pytest
from skills.clawpathy_autoresearch.sandbox import Sandbox, load_sandbox


def test_load_sandbox_full(tmp_path: Path):
    path = tmp_path / "sandbox.yaml"
    path.write_text(
        "allowed_tools: [Read, Bash]\n"
        "network:\n  allowed: false\n"
        "data_dir: ./data\n"
        "python_packages: [pandas]\n"
        "timeout_seconds: 300\n"
    )
    sb = load_sandbox(path)
    assert sb.allowed_tools == ["Read", "Bash"]
    assert sb.network_allowed is False
    assert sb.network_hosts == []
    assert sb.data_dir == "./data"
    assert sb.python_packages == ["pandas"]
    assert sb.timeout_seconds == 300


def test_load_sandbox_network_hostlist(tmp_path: Path):
    path = tmp_path / "sandbox.yaml"
    path.write_text(
        "allowed_tools: [WebFetch]\n"
        "network:\n  allowed: [ebi.ac.uk, ensembl.org]\n"
    )
    sb = load_sandbox(path)
    assert sb.network_allowed is True
    assert sb.network_hosts == ["ebi.ac.uk", "ensembl.org"]


def test_load_sandbox_defaults(tmp_path: Path):
    path = tmp_path / "sandbox.yaml"
    path.write_text("allowed_tools: [Read]\n")
    sb = load_sandbox(path)
    assert sb.network_allowed is False
    assert sb.timeout_seconds == 600
    assert sb.python_packages == []


def test_missing_allowed_tools_raises(tmp_path: Path):
    path = tmp_path / "sandbox.yaml"
    path.write_text("network:\n  allowed: false\n")
    with pytest.raises(ValueError, match="allowed_tools"):
        load_sandbox(path)


def test_sandbox_as_contract_text():
    sb = Sandbox(
        allowed_tools=["Read"],
        network_allowed=False,
        network_hosts=[],
        data_dir="./data",
        python_packages=["pandas"],
        timeout_seconds=300,
    )
    text = sb.as_contract_text()
    assert "allowed_tools" in text
    assert "Read" in text
    assert "network" in text
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest skills/clawpathy-autoresearch/tests/test_sandbox.py -v`
Expected: FAIL — `sandbox` module does not exist.

- [ ] **Step 3: Implement sandbox.py**

Create `skills/clawpathy-autoresearch/sandbox.py`:

```python
"""Sandbox spec loader for clawpathy-autoresearch.

Tasks declare executor constraints in sandbox.yaml: which tools the executor
may call, whether it has network, which hosts it may reach, what packages
are guaranteed importable, and a timeout.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

import yaml


@dataclass
class Sandbox:
    allowed_tools: list[str]
    network_allowed: bool = False
    network_hosts: list[str] = field(default_factory=list)
    data_dir: str = "./data"
    python_packages: list[str] = field(default_factory=list)
    timeout_seconds: int = 600

    def as_contract_text(self) -> str:
        """Render as a YAML-like block suitable for inclusion in a prompt."""
        return (
            f"allowed_tools: {self.allowed_tools}\n"
            f"network:\n"
            f"  allowed: {self.network_allowed}\n"
            f"  hosts: {self.network_hosts}\n"
            f"data_dir: {self.data_dir}\n"
            f"python_packages: {self.python_packages}\n"
            f"timeout_seconds: {self.timeout_seconds}\n"
        )


def load_sandbox(path: Path) -> Sandbox:
    raw = yaml.safe_load(Path(path).read_text()) or {}
    if "allowed_tools" not in raw or not raw["allowed_tools"]:
        raise ValueError(f"sandbox.yaml missing 'allowed_tools': {path}")
    net = raw.get("network") or {}
    allowed = net.get("allowed", False)
    if isinstance(allowed, list):
        network_allowed = True
        network_hosts = list(allowed)
    else:
        network_allowed = bool(allowed)
        network_hosts = []
    return Sandbox(
        allowed_tools=list(raw["allowed_tools"]),
        network_allowed=network_allowed,
        network_hosts=network_hosts,
        data_dir=raw.get("data_dir", "./data"),
        python_packages=list(raw.get("python_packages") or []),
        timeout_seconds=int(raw.get("timeout_seconds", 600)),
    )
```

- [ ] **Step 4: Install PyYAML if missing**

Run: `uv add pyyaml` (or add to project deps however this repo tracks them).

- [ ] **Step 5: Run tests**

Run: `uv run pytest skills/clawpathy-autoresearch/tests/test_sandbox.py -v`
Expected: all 5 tests PASS.

- [ ] **Step 6: Commit**

```bash
git add skills/clawpathy-autoresearch/sandbox.py skills/clawpathy-autoresearch/tests/test_sandbox.py
git commit -m "feat(autoresearch): sandbox spec loader"
```

---

### Task 4: Dispatcher interface

**Files:**
- Create: `skills/clawpathy-autoresearch/dispatcher.py`
- Create: `skills/clawpathy-autoresearch/tests/test_dispatcher.py`

The dispatcher is the seam between the loop and the real world. The loop calls `dispatcher.dispatch(role, prompt, allowed_tools)` and gets back a string. The default dispatcher raises with guidance. Tests inject a callable. The outer agent can provide its own.

- [ ] **Step 1: Write failing test**

Create `skills/clawpathy-autoresearch/tests/test_dispatcher.py`:

```python
import pytest
from skills.clawpathy_autoresearch.dispatcher import (
    Dispatcher,
    CallableDispatcher,
    DispatchRequest,
)


def test_callable_dispatcher_passes_through():
    calls = []

    def fake(req: DispatchRequest) -> str:
        calls.append(req)
        return "ok-response"

    d = CallableDispatcher(fake)
    out = d.dispatch(DispatchRequest(role="proposer", prompt="hi", allowed_tools=["Read"]))
    assert out == "ok-response"
    assert calls[0].role == "proposer"
    assert calls[0].prompt == "hi"
    assert calls[0].allowed_tools == ["Read"]


def test_base_dispatcher_raises_with_guidance():
    d = Dispatcher()
    with pytest.raises(RuntimeError, match="outer agent"):
        d.dispatch(DispatchRequest(role="proposer", prompt="x", allowed_tools=[]))
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest skills/clawpathy-autoresearch/tests/test_dispatcher.py -v`
Expected: FAIL — module does not exist.

- [ ] **Step 3: Implement dispatcher.py**

Create `skills/clawpathy-autoresearch/dispatcher.py`:

```python
"""Pluggable subagent dispatcher.

The loop does not call LLMs directly. It asks a Dispatcher to execute a
role (proposer or executor) with a prompt and a tool whitelist, and
receives the response as text. The outer agent provides a concrete
dispatcher; tests inject a callable.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable


@dataclass
class DispatchRequest:
    role: str              # "proposer" or "executor"
    prompt: str
    allowed_tools: list[str]


class Dispatcher:
    """Base dispatcher. Override `dispatch` in subclasses or use CallableDispatcher."""

    def dispatch(self, request: DispatchRequest) -> str:
        raise RuntimeError(
            "No dispatcher configured. Provide a CallableDispatcher "
            "whose callback routes requests to an outer agent's Agent-tool "
            "invocation, or subclass Dispatcher."
        )


class CallableDispatcher(Dispatcher):
    """Adapter that wraps any callable taking a DispatchRequest."""

    def __init__(self, func: Callable[[DispatchRequest], str]):
        self._func = func

    def dispatch(self, request: DispatchRequest) -> str:
        return self._func(request)
```

- [ ] **Step 4: Run tests**

Run: `uv run pytest skills/clawpathy-autoresearch/tests/test_dispatcher.py -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add skills/clawpathy-autoresearch/dispatcher.py skills/clawpathy-autoresearch/tests/test_dispatcher.py
git commit -m "feat(autoresearch): pluggable dispatcher interface"
```

---

### Task 5: Proposer

**Files:**
- Modify: `skills/clawpathy-autoresearch/proposer.py` (rewrite)
- Create/modify: `skills/clawpathy-autoresearch/tests/test_proposer.py`

- [ ] **Step 1: Write failing test**

Overwrite `skills/clawpathy-autoresearch/tests/test_proposer.py`:

```python
from skills.clawpathy_autoresearch.dispatcher import CallableDispatcher, DispatchRequest
from skills.clawpathy_autoresearch.proposer import propose_skill


def test_proposer_returns_replacement_skill():
    captured = {}

    def fake(req: DispatchRequest) -> str:
        captured["prompt"] = req.prompt
        captured["role"] = req.role
        return "```markdown\n# NEW SKILL\nStep 1: do thing\n```"

    d = CallableDispatcher(fake)
    out = propose_skill(
        dispatcher=d,
        task_prompt="make the thing",
        current_skill="# OLD\n",
        last_score=0.5,
        last_breakdown={"err_a": 0.5},
        recent_history=[{"iter": 1, "score": 0.6, "kept": False}],
    )
    assert out == "# NEW SKILL\nStep 1: do thing"
    assert captured["role"] == "proposer"
    assert "make the thing" in captured["prompt"]
    assert "# OLD" in captured["prompt"]
    assert "err_a" in captured["prompt"]


def test_proposer_no_history_first_iteration():
    def fake(req):
        return "# SKILL\ncontent"
    d = CallableDispatcher(fake)
    out = propose_skill(
        dispatcher=d, task_prompt="t", current_skill="",
        last_score=None, last_breakdown=None, recent_history=[],
    )
    assert out == "# SKILL\ncontent"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest skills/clawpathy-autoresearch/tests/test_proposer.py -v`
Expected: FAIL — either old API or extraction logic missing.

- [ ] **Step 3: Rewrite proposer.py**

Replace entire file with:

```python
"""Proposer: asks a subagent to produce a replacement SKILL.md."""
from __future__ import annotations

import json
import re
from typing import Any

from .dispatcher import Dispatcher, DispatchRequest


PROPOSER_PROMPT = """\
You are a methodology author. Propose an improved SKILL.md for the task below.

Task:
{task_prompt}

Current SKILL.md:
```
{current_skill}
```

Last score: {last_score_line}
Last score breakdown: {last_breakdown_line}
Recent history (most recent last):
{history_block}

Return a full replacement SKILL.md as a markdown code block. Do not include \
commentary outside the code block. Describe methodology only — never paste \
ground-truth answers, rsIDs, specific gene names, or target values into the \
skill. The executor that runs this skill will have sandbox-restricted tools.
"""


def _extract_markdown(text: str) -> str:
    m = re.search(r"```(?:markdown|md)?\s*\n(.*?)\n```", text, re.DOTALL)
    if m:
        return m.group(1).strip()
    return text.strip()


def propose_skill(
    dispatcher: Dispatcher,
    task_prompt: str,
    current_skill: str,
    last_score: float | None,
    last_breakdown: dict[str, Any] | None,
    recent_history: list[dict[str, Any]],
) -> str:
    last_score_line = "n/a (first iteration)" if last_score is None else f"{last_score:.4f}"
    last_breakdown_line = "n/a" if not last_breakdown else json.dumps(last_breakdown)
    history_lines = [
        f"  iter={h.get('iter')} score={h.get('score')} kept={h.get('kept')}"
        for h in recent_history
    ]
    history_block = "\n".join(history_lines) if history_lines else "  (no prior iterations)"

    prompt = PROPOSER_PROMPT.format(
        task_prompt=task_prompt,
        current_skill=current_skill or "(empty)",
        last_score_line=last_score_line,
        last_breakdown_line=last_breakdown_line,
        history_block=history_block,
    )
    response = dispatcher.dispatch(DispatchRequest(
        role="proposer",
        prompt=prompt,
        allowed_tools=[],
    ))
    return _extract_markdown(response)
```

- [ ] **Step 4: Run tests**

Run: `uv run pytest skills/clawpathy-autoresearch/tests/test_proposer.py -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add skills/clawpathy-autoresearch/proposer.py skills/clawpathy-autoresearch/tests/test_proposer.py
git commit -m "feat(autoresearch): rewrite proposer around dispatcher"
```

---

### Task 6: Executor

**Files:**
- Create: `skills/clawpathy-autoresearch/executor.py`
- Create: `skills/clawpathy-autoresearch/tests/test_executor.py`

The executor dispatches a subagent with the skill content and sandbox contract, parses the returned JSON, and returns it. Tool whitelist from sandbox is passed into `DispatchRequest.allowed_tools`.

- [ ] **Step 1: Write failing test**

Create `skills/clawpathy-autoresearch/tests/test_executor.py`:

```python
import pytest
from skills.clawpathy_autoresearch.dispatcher import CallableDispatcher, DispatchRequest
from skills.clawpathy_autoresearch.sandbox import Sandbox
from skills.clawpathy_autoresearch.executor import execute_skill


def _sandbox():
    return Sandbox(
        allowed_tools=["Read", "Bash"],
        network_allowed=False,
        network_hosts=[],
        data_dir="./data",
        python_packages=["pandas"],
        timeout_seconds=300,
    )


def test_executor_passes_whitelist_and_returns_json():
    captured = {}

    def fake(req: DispatchRequest) -> str:
        captured["req"] = req
        return '```json\n{"result": 42}\n```'

    d = CallableDispatcher(fake)
    out = execute_skill(
        dispatcher=d,
        skill_content="# SKILL\nDo X",
        sandbox=_sandbox(),
        data_dir_abs="/tmp/data",
    )
    assert out == {"result": 42}
    assert captured["req"].role == "executor"
    assert captured["req"].allowed_tools == ["Read", "Bash"]
    assert "# SKILL" in captured["req"].prompt
    assert "allowed_tools" in captured["req"].prompt  # contract included


def test_executor_raises_on_bad_json():
    d = CallableDispatcher(lambda req: "not json at all")
    with pytest.raises(ValueError, match="JSON"):
        execute_skill(
            dispatcher=d,
            skill_content="# SKILL",
            sandbox=_sandbox(),
            data_dir_abs="/tmp/data",
        )


def test_executor_accepts_raw_json_no_fence():
    d = CallableDispatcher(lambda req: '{"a": 1}')
    out = execute_skill(
        dispatcher=d, skill_content="# SKILL", sandbox=_sandbox(), data_dir_abs="/tmp/data",
    )
    assert out == {"a": 1}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest skills/clawpathy-autoresearch/tests/test_executor.py -v`
Expected: FAIL — module does not exist.

- [ ] **Step 3: Implement executor.py**

Create `skills/clawpathy-autoresearch/executor.py`:

```python
"""Executor: runs SKILL.md via a subagent under a sandbox contract."""
from __future__ import annotations

import json
import re
from typing import Any

from .dispatcher import Dispatcher, DispatchRequest
from .sandbox import Sandbox


EXECUTOR_PROMPT = """\
You are executing the methodology described in the SKILL.md below. Follow it \
exactly. Do not consult external information beyond what the skill tells you \
to fetch. Return only the JSON output the skill specifies.

Sandbox contract (enforced):
```yaml
{sandbox_contract}
```
Your task data directory (absolute path): {data_dir_abs}

SKILL.md:
```
{skill_content}
```

Return ONLY a JSON object, optionally wrapped in a ```json code fence. No \
prose.
"""


def _extract_json(text: str) -> dict[str, Any]:
    m = re.search(r"```json\s*\n(.*?)\n```", text, re.DOTALL)
    candidate = m.group(1) if m else text
    candidate = candidate.strip()
    start = candidate.find("{")
    end = candidate.rfind("}")
    if start == -1 or end == -1:
        raise ValueError(f"Executor response contained no JSON object: {text[:200]}")
    try:
        return json.loads(candidate[start:end + 1])
    except json.JSONDecodeError as exc:
        raise ValueError(f"Executor response JSON was malformed: {exc}") from exc


def execute_skill(
    dispatcher: Dispatcher,
    skill_content: str,
    sandbox: Sandbox,
    data_dir_abs: str,
) -> dict[str, Any]:
    prompt = EXECUTOR_PROMPT.format(
        sandbox_contract=sandbox.as_contract_text(),
        data_dir_abs=data_dir_abs,
        skill_content=skill_content,
    )
    response = dispatcher.dispatch(DispatchRequest(
        role="executor",
        prompt=prompt,
        allowed_tools=list(sandbox.allowed_tools),
    ))
    return _extract_json(response)
```

- [ ] **Step 4: Run tests**

Run: `uv run pytest skills/clawpathy-autoresearch/tests/test_executor.py -v`
Expected: all 3 tests PASS.

- [ ] **Step 5: Commit**

```bash
git add skills/clawpathy-autoresearch/executor.py skills/clawpathy-autoresearch/tests/test_executor.py
git commit -m "feat(autoresearch): executor with sandbox contract and JSON parse"
```

---

### Task 7: Loop engine

**Files:**
- Create: `skills/clawpathy-autoresearch/loop.py`
- Create: `skills/clawpathy-autoresearch/tests/test_loop.py`

The loop orchestrates proposer, executor, scorer, snapshot, history. Scorer is loaded via `importlib.util.spec_from_file_location` since it lives in the task directory.

- [ ] **Step 1: Write failing test**

Create `skills/clawpathy-autoresearch/tests/test_loop.py`:

```python
import json
from pathlib import Path

from skills.clawpathy_autoresearch.dispatcher import CallableDispatcher, DispatchRequest
from skills.clawpathy_autoresearch.loop import run_loop


def _build_task(tmp_path: Path, target_text: str) -> Path:
    ws = tmp_path / "ws"
    (ws / "skill").mkdir(parents=True)
    (ws / "skill" / "SKILL.md").write_text("")
    (ws / "task.json").write_text(json.dumps({
        "name": "hello",
        "max_iterations": 3,
        "early_stop_n": 2,
        "target_score": 0.0,
    }))
    (ws / "ground_truth.json").write_text(json.dumps({"target": target_text}))
    (ws / "sandbox.yaml").write_text("allowed_tools: [Read]\n")
    (ws / "scorer.py").write_text(
        "import json\n"
        "from pathlib import Path\n"
        "def score(output, task_dir):\n"
        "    gt = json.loads((Path(task_dir)/'ground_truth.json').read_text())\n"
        "    return (0.0 if output.get('text')==gt['target'] else 1.0, {})\n"
    )
    return ws


def test_loop_converges_on_trivial_task(tmp_path):
    ws = _build_task(tmp_path, target_text="hello")
    calls = {"n": 0}

    def fake(req: DispatchRequest) -> str:
        calls["n"] += 1
        if req.role == "proposer":
            # Always propose the same minimal skill
            return '```markdown\n# SKILL\nReturn {"text": "hello"}\n```'
        # executor
        return '```json\n{"text": "hello"}\n```'

    d = CallableDispatcher(fake)
    result = run_loop(workspace_dir=ws, dispatcher=d)
    assert result["best_score"] == 0.0
    assert result["iterations_run"] >= 1
    history = [json.loads(l) for l in (ws / "history.jsonl").read_text().splitlines()]
    assert history[0]["kept"] is True
    assert (ws / "skill" / "SKILL.md").read_text().startswith("# SKILL")


def test_loop_reverts_when_score_worsens(tmp_path):
    ws = _build_task(tmp_path, target_text="hello")
    # First proposal hits target, second proposal breaks it.
    proposals = iter([
        '```markdown\n# good\nreturn hello\n```',
        '```markdown\n# bad\nreturn wrong\n```',
    ])
    outputs = iter([
        '{"text": "hello"}',
        '{"text": "goodbye"}',
    ])

    def fake(req):
        return next(proposals) if req.role == "proposer" else next(outputs)

    d = CallableDispatcher(fake)
    result = run_loop(workspace_dir=ws, dispatcher=d)
    # best score should be 0.0 and final SKILL.md should be the good one
    assert result["best_score"] == 0.0
    assert "good" in (ws / "skill" / "SKILL.md").read_text()
    history = [json.loads(l) for l in (ws / "history.jsonl").read_text().splitlines()]
    assert history[0]["kept"] is True
    assert history[1]["kept"] is False


def test_loop_stops_on_early_stop(tmp_path):
    ws = _build_task(tmp_path, target_text="hello")
    # Every proposal/execution misses — 2 non-improvements in a row triggers stop.
    def fake(req):
        if req.role == "proposer":
            return '```markdown\nbad\n```'
        return '{"text": "wrong"}'
    d = CallableDispatcher(fake)
    result = run_loop(workspace_dir=ws, dispatcher=d)
    # early_stop_n = 2 in fixture, max_iterations = 3
    assert result["iterations_run"] <= 3
    assert result["best_score"] == 1.0
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest skills/clawpathy-autoresearch/tests/test_loop.py -v`
Expected: FAIL — `loop` module does not exist.

- [ ] **Step 3: Implement loop.py**

Create `skills/clawpathy-autoresearch/loop.py`:

```python
"""Iteration loop engine for clawpathy-autoresearch."""
from __future__ import annotations

import importlib.util
import json
import math
import shutil
from pathlib import Path
from typing import Any

from .dispatcher import Dispatcher
from .executor import execute_skill
from .proposer import propose_skill
from .sandbox import load_sandbox
from .workspace import load_workspace


def _load_scorer(scorer_path: Path):
    spec = importlib.util.spec_from_file_location("_autoresearch_scorer", scorer_path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not load scorer: {scorer_path}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    if not hasattr(mod, "score"):
        raise AttributeError(f"Scorer at {scorer_path} has no score() function")
    return mod.score


def _snapshot(skill_md: Path, snap_dir: Path, iter_num: int) -> None:
    snap_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(skill_md, snap_dir / f"iter-{iter_num:03d}.md")


def _revert(skill_md: Path, snap_dir: Path, iter_num: int) -> None:
    snap = snap_dir / f"iter-{iter_num:03d}.md"
    shutil.copy2(snap, skill_md)


def run_loop(workspace_dir: Path, dispatcher: Dispatcher) -> dict[str, Any]:
    ws = load_workspace(Path(workspace_dir))
    sandbox = load_sandbox(ws.sandbox_path)
    scorer = _load_scorer(ws.scorer_path)

    skill_md = ws.skill_dir / "SKILL.md"
    history_path = ws.workspace_dir / "history.jsonl"
    snap_dir = ws.workspace_dir / "snapshots"
    history_path.write_text("")  # reset each run

    task_prompt = ws.description or ws.name
    best_score = math.inf
    last_score: float | None = None
    last_breakdown: dict[str, Any] | None = None
    recent_history: list[dict[str, Any]] = []
    consecutive_non_improvements = 0
    iterations_run = 0

    for i in range(1, ws.max_iterations + 1):
        iterations_run = i
        _snapshot(skill_md, snap_dir, i)

        proposal = propose_skill(
            dispatcher=dispatcher,
            task_prompt=task_prompt,
            current_skill=skill_md.read_text(),
            last_score=last_score,
            last_breakdown=last_breakdown,
            recent_history=recent_history[-3:],
        )
        skill_md.write_text(proposal)

        try:
            output = execute_skill(
                dispatcher=dispatcher,
                skill_content=proposal,
                sandbox=sandbox,
                data_dir_abs=str(ws.data_dir.resolve()),
            )
            total, breakdown = scorer(output, ws.workspace_dir)
        except Exception as exc:
            total = math.inf
            breakdown = {"error": str(exc)}

        kept = total < best_score
        if kept:
            best_score = total
            consecutive_non_improvements = 0
        else:
            _revert(skill_md, snap_dir, i)
            consecutive_non_improvements += 1

        entry = {
            "iter": i,
            "score": total if total != math.inf else None,
            "kept": kept,
            "breakdown": breakdown,
        }
        recent_history.append(entry)
        last_score = total if total != math.inf else None
        last_breakdown = breakdown
        with history_path.open("a") as f:
            f.write(json.dumps(entry) + "\n")

        if ws.target_score is not None and best_score <= ws.target_score:
            break
        if consecutive_non_improvements >= ws.early_stop_n:
            break

    return {
        "best_score": best_score if best_score != math.inf else None,
        "iterations_run": iterations_run,
        "history_path": str(history_path),
    }
```

- [ ] **Step 4: Run tests**

Run: `uv run pytest skills/clawpathy-autoresearch/tests/test_loop.py -v`
Expected: all 3 tests PASS.

- [ ] **Step 5: Commit**

```bash
git add skills/clawpathy-autoresearch/loop.py skills/clawpathy-autoresearch/tests/test_loop.py
git commit -m "feat(autoresearch): loop engine with snapshot/revert and history"
```

---

### Task 8: Rewrite autoresearch.py CLI entry

**Files:**
- Modify: `skills/clawpathy-autoresearch/autoresearch.py`
- Create: `skills/clawpathy-autoresearch/tests/test_cli.py`

The CLI does not include a real dispatcher — running the loop autonomously requires an LLM backend the outer agent provides. The CLI supports two modes: `--dry-run` (mock dispatcher, used by tests) and a default that exits with a helpful message if no dispatcher is wired.

- [ ] **Step 1: Write failing test**

Create `skills/clawpathy-autoresearch/tests/test_cli.py`:

```python
import json
import subprocess
import sys
from pathlib import Path


def test_cli_dry_run_on_hello_task(tmp_path):
    ws = tmp_path / "ws"
    (ws / "skill").mkdir(parents=True)
    (ws / "skill" / "SKILL.md").write_text("")
    (ws / "task.json").write_text(json.dumps({
        "name": "hello", "max_iterations": 2, "early_stop_n": 2, "target_score": 0.0,
    }))
    (ws / "ground_truth.json").write_text(json.dumps({"target": "hello"}))
    (ws / "sandbox.yaml").write_text("allowed_tools: [Read]\n")
    (ws / "scorer.py").write_text(
        "import json\nfrom pathlib import Path\n"
        "def score(o, d):\n"
        "    gt = json.loads((Path(d)/'ground_truth.json').read_text())\n"
        "    return (0.0 if o.get('text')==gt['target'] else 1.0, {})\n"
    )
    result = subprocess.run(
        [sys.executable, "-m", "skills.clawpathy_autoresearch.autoresearch",
         "--workspace", str(ws), "--dry-run"],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    assert "best_score" in result.stdout
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest skills/clawpathy-autoresearch/tests/test_cli.py -v`
Expected: FAIL — current autoresearch.py is a stub.

- [ ] **Step 3: Rewrite autoresearch.py**

Overwrite `skills/clawpathy-autoresearch/autoresearch.py`:

```python
"""CLI entry for clawpathy-autoresearch.

Usage:
    python -m skills.clawpathy_autoresearch.autoresearch \\
        --workspace path/to/workspace [--dry-run]

The loop requires a subagent dispatcher. The CLI does not provide one —
the outer agent (or a future SDK integration) must wire it. With --dry-run,
a mock dispatcher is used that returns the current SKILL.md for both
proposer and executor roles; useful for validating workspace layout.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .dispatcher import CallableDispatcher, DispatchRequest
from .loop import run_loop


def _dry_run_dispatcher() -> CallableDispatcher:
    """Returns a skeletal output. Produces malformed JSON for executor, which
    the loop catches and treats as a non-improvement. Useful for smoke-testing
    workspace layout without an LLM."""
    def fake(req: DispatchRequest) -> str:
        if req.role == "proposer":
            return "```markdown\n# stub skill\n```"
        return '```json\n{}\n```'
    return CallableDispatcher(fake)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workspace", required=True, type=Path)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)

    if args.dry_run:
        dispatcher = _dry_run_dispatcher()
    else:
        print(
            "error: no dispatcher configured. Use --dry-run, or invoke "
            "run_loop() directly with a CallableDispatcher provided by the "
            "outer agent.",
            file=sys.stderr,
        )
        return 2

    result = run_loop(workspace_dir=args.workspace, dispatcher=dispatcher)
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

- [ ] **Step 4: Run tests**

Run: `uv run pytest skills/clawpathy-autoresearch/tests/test_cli.py -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add skills/clawpathy-autoresearch/autoresearch.py skills/clawpathy-autoresearch/tests/test_cli.py
git commit -m "feat(autoresearch): CLI entry with dry-run mode"
```

---

### Task 9: Hello-world task + full-suite sanity

**Files:**
- Create: `skills/clawpathy-autoresearch/tasks/hello-world/task.json`
- Create: `skills/clawpathy-autoresearch/tasks/hello-world/ground_truth.json`
- Create: `skills/clawpathy-autoresearch/tasks/hello-world/sandbox.yaml`
- Create: `skills/clawpathy-autoresearch/tasks/hello-world/scorer.py`
- Create: `skills/clawpathy-autoresearch/tasks/hello-world/skill/SKILL.md`
- Create: `skills/clawpathy-autoresearch/tasks/hello-world/README.md`

- [ ] **Step 1: Create task.json**

```bash
mkdir -p skills/clawpathy-autoresearch/tasks/hello-world/skill
```

`skills/clawpathy-autoresearch/tasks/hello-world/task.json`:

```json
{
  "name": "hello-world",
  "description": "Return a JSON object whose 'text' field equals the string 'hello'.",
  "max_iterations": 5,
  "early_stop_n": 2,
  "target_score": 0.0
}
```

- [ ] **Step 2: Create ground_truth.json, sandbox.yaml, empty SKILL.md**

`skills/clawpathy-autoresearch/tasks/hello-world/ground_truth.json`:

```json
{"target": "hello"}
```

`skills/clawpathy-autoresearch/tasks/hello-world/sandbox.yaml`:

```yaml
allowed_tools: [Read]
network:
  allowed: false
python_packages: []
timeout_seconds: 60
```

`skills/clawpathy-autoresearch/tasks/hello-world/skill/SKILL.md`:

```
```
(An empty file.)

- [ ] **Step 3: Create scorer.py**

`skills/clawpathy-autoresearch/tasks/hello-world/scorer.py`:

```python
"""Scorer for the hello-world integration task."""
from __future__ import annotations

import json
from pathlib import Path


def score(skill_output: dict, task_dir: Path) -> tuple[float, dict]:
    gt = json.loads((Path(task_dir) / "ground_truth.json").read_text())
    got = skill_output.get("text")
    if got == gt["target"]:
        return 0.0, {"text_match": 0.0}
    return 1.0, {"text_match": 1.0, "got": got, "want": gt["target"]}
```

- [ ] **Step 4: Create README**

`skills/clawpathy-autoresearch/tasks/hello-world/README.md`:

```markdown
# hello-world

Smallest possible autoresearch task. The executor must return
`{"text": "hello"}`. Used as an integration smoke-test for the loop.
Run via the outer agent providing a dispatcher, or via
`autoresearch.py --workspace tasks/hello-world --dry-run`.
```

- [ ] **Step 5: Run full autoresearch test suite**

Run: `uv run pytest skills/clawpathy-autoresearch/tests/ -v`
Expected: all tests in sandbox, dispatcher, proposer, executor, loop, cli, workspace pass.

- [ ] **Step 6: Run dry-run CLI against the new task**

Run: `uv run python -m skills.clawpathy_autoresearch.autoresearch --workspace skills/clawpathy-autoresearch/tasks/hello-world --dry-run`
Expected: prints a JSON result with `best_score` (likely 1.0, since stub dispatcher won't produce "hello") and exit code 0. Confirms the workspace is well-formed.

- [ ] **Step 7: Commit**

```bash
git add skills/clawpathy-autoresearch/tasks/hello-world/
git commit -m "feat(autoresearch): hello-world integration task"
```

---

## Post-plan notes

- The outer agent runs a real loop by constructing a `CallableDispatcher` whose callback invokes the Agent tool with the appropriate subagent prompt and tool whitelist, then calls `run_loop(workspace_dir, dispatcher)`.
- Future work (not in this plan): Claude Agent SDK integration so the loop runs headless; a tool-log parser that detects sandbox violations; optional opt-in modules (critic, held-out eval) when a real task shows they help.
