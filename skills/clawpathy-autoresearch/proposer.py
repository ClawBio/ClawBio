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

Last score (lower is better): {last_score_line}
Last score breakdown: {last_breakdown_line}
Recent history (most recent last):
{history_block}

Write a **minimal, imperative, executor-facing playbook**. Constraints:
- No YAML frontmatter. No Trigger/Scope/Gotchas/Safety sections.
- One short `## Workflow` section with numbered steps.
- Each step must be either (a) an exact bash command to run, (b) an exact \
file to read, or (c) an exact rule for assembling output JSON. No \
conditionals the executor has to resolve at runtime unless strictly needed.
- End with the exact output JSON shape the executor must return.
- The executor is an LLM with shell access but no memory of this task; it \
will read this SKILL.md cold and act. Optimise for "a smart intern follows \
this in one pass without asking questions".
- Never paste ground-truth answers, rsIDs, gene names, or target values \
into the skill — the executor must derive them by running the methodology.

Return the replacement SKILL.md as a single markdown code block. No \
commentary outside the code block.
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
