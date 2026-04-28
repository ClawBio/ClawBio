"""Executor: runs the current SKILL.md against the task, captures full transcript."""
from __future__ import annotations

from pathlib import Path

from .dispatcher import Dispatcher, DispatchRequest


DEFAULT_TOOLS = ["Bash", "Read", "Write", "Edit", "Glob", "Grep", "WebFetch"]


def execute_skill(
    dispatcher: Dispatcher,
    prompt_template: str,
    task_description: str,
    skill_content: str,
    workspace_root: Path,
    model: str = "sonnet",
    allowed_tools: list[str] | None = None,
) -> str:
    """Return the raw transcript from the executor subagent.

    The executor runs with shell access inside workspace_root. The judge
    will read this transcript plus whatever the executor wrote to output/.
    """
    prompt = prompt_template.format(
        task_description=task_description,
        skill=skill_content,
        workspace=str(workspace_root.resolve()),
    )
    return dispatcher.dispatch(DispatchRequest(
        role="executor",
        prompt=prompt,
        model=model,
        allowed_tools=allowed_tools or DEFAULT_TOOLS,
        cwd=str(workspace_root.resolve()),
    ))
