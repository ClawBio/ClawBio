"""Proposer: asks an LLM to rewrite the SKILL.md based on last judge feedback."""
from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

from .dispatcher import Dispatcher, DispatchRequest


def _extract_markdown(text: str) -> str:
    # Greedy — match the outermost fenced block. Non-greedy was the old bug.
    m = re.search(r"```(?:markdown|md)?\s*\n(.*)\n```", text, re.DOTALL)
    return (m.group(1) if m else text).strip()


def propose_skill(
    dispatcher: Dispatcher,
    prompt_template: str,
    task_description: str,
    rubric: str,
    current_skill: str,
    last_score: float | None,
    last_judgement: dict[str, Any] | None,
    recent_history: list[dict[str, Any]],
    model: str = "sonnet",
) -> str:
    history_block = "\n".join(
        f"  iter={h.get('iter')} score={h.get('score')} kept={h.get('kept')}"
        for h in recent_history
    ) or "  (no prior iterations)"
    prompt = prompt_template.format(
        task_description=task_description,
        rubric=rubric,
        current_skill=current_skill or "(empty)",
        last_score="n/a (first iteration)" if last_score is None else f"{last_score:.4f}",
        last_judgement=json.dumps(last_judgement, indent=2) if last_judgement else "n/a",
        history_block=history_block,
    )
    response = dispatcher.dispatch(DispatchRequest(
        role="proposer", prompt=prompt, model=model, allowed_tools=[],
    ))
    return _extract_markdown(response)
