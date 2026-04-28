"""LLM judge. Reads SKILL.md, executor transcript, outputs, rubric — emits score."""
from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

from .dispatcher import Dispatcher, DispatchRequest


MAX_TRANSCRIPT_CHARS = 40_000
MAX_FILE_CHARS = 8_000


def _truncate(text: str, n: int) -> str:
    if len(text) <= n:
        return text
    half = n // 2
    return text[:half] + f"\n\n... [truncated {len(text) - n} chars] ...\n\n" + text[-half:]


def _summarise_output_dir(output_dir: Path) -> str:
    if not output_dir.is_dir():
        return "(no output directory)"
    parts = []
    for f in sorted(output_dir.iterdir()):
        if not f.is_file():
            continue
        size = f.stat().st_size
        header = f"### {f.name} ({size} bytes)"
        if f.suffix.lower() in (".json", ".txt", ".md", ".tsv", ".csv", ".log"):
            try:
                parts.append(header + "\n```\n" + _truncate(f.read_text(errors="replace"), MAX_FILE_CHARS) + "\n```")
                continue
            except Exception:
                pass
        parts.append(header + "  (binary / unreadable)")
    return "\n\n".join(parts) if parts else "(empty output directory)"


def _extract_json(text: str) -> dict[str, Any]:
    # Prefer fenced json block; else first {...} span.
    m = re.search(r"```json\s*\n(.*?)\n```", text, re.DOTALL)
    if m:
        return json.loads(m.group(1))
    m = re.search(r"\{.*\}", text, re.DOTALL)
    if m:
        return json.loads(m.group(0))
    raise ValueError(f"no JSON in judge response: {text[:500]}")


def judge_run(
    dispatcher: Dispatcher,
    prompt_template: str,
    task_description: str,
    rubric: str,
    skill_content: str,
    transcript: str,
    output_dir: Path,
    reference_dir: Path,
    model: str = "opus",
) -> tuple[float, dict[str, Any]]:
    reference_summary = _summarise_output_dir(reference_dir) if reference_dir.is_dir() else "(no reference artifacts)"
    prompt = prompt_template.format(
        task_description=task_description,
        rubric=rubric,
        skill=skill_content,
        transcript=_truncate(transcript, MAX_TRANSCRIPT_CHARS),
        outputs=_summarise_output_dir(output_dir),
        reference=reference_summary,
    )
    response = dispatcher.dispatch(DispatchRequest(
        role="judge", prompt=prompt, model=model, allowed_tools=[],
    ))
    verdict = _extract_json(response)
    # Convention: score in [0,1], LOWER is better (0 = perfect).
    score = float(verdict.get("score"))
    if not 0.0 <= score <= 1.0:
        raise ValueError(f"judge score out of range [0,1]: {score}")
    return score, verdict
