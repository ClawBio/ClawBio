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
