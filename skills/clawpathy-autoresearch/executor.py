"""Executor: runs SKILL.md via a subagent under a sandbox contract."""
from __future__ import annotations

import json
import re
from typing import Any

from .dispatcher import Dispatcher, DispatchRequest
from .sandbox import Sandbox


EXECUTOR_PROMPT = """\
You are executing the methodology described in the SKILL.md below. Execute \
every numbered step of the workflow, in order, using the shell and file \
tools available. Do not skip steps. Do not substitute your own estimates \
for values the skill says to compute. Do not consult external information \
beyond what the skill tells you to fetch.

Sandbox contract (enforced):
```yaml
{sandbox_contract}
```
Your task data directory (absolute path): {data_dir_abs}

SKILL.md:
```
{skill_content}
```

After completing every step, return ONLY the JSON object the skill \
specifies, wrapped in a ```json code fence. No prose before or after. \
If a step fails, fix it and retry before responding.
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
