"""Shared helpers for structured skill-suggested follow-up actions.

Skills may return `suggested_actions[]` in result.json. The chat adapters use
this module to render those actions, remember them briefly, resolve natural
follow-up replies, and execute only the stored nested request payload.
"""

from __future__ import annotations

import difflib
import json
import re
import tempfile
import time
from datetime import datetime
from pathlib import Path
from typing import Any, Callable

from clawbio.contract_alerts import normalise_contract_alerts, render_contract_alerts  # noqa: F401  (render_contract_alerts re-exported for bot/roboterri*.py)

DEFAULT_PENDING_ACTION_TTL_SECONDS = 30 * 60

# These are stand-alone confirmation/navigation phrases. They are intentionally
# separate from numbered choices so "yes" can stay ambiguous when several
# pending actions exist.
_AFFIRM_ONLY = {
    "yes",
    "y",
    "yeah",
    "yep",
    "ok",
    "okay",
    "sure",
    "please do",
    "go ahead",
    "do it",
    "run it",
    "please show me",
    "show me",
}

_CANCEL_WORDS = {
    "cancel",
    "no",
    "nope",
    "never mind",
    "nevermind",
    "stop",
    "abort",
}

_INDEX_CONFIRM_PREFIXES = (
    "yes",
    "run",
    "do",
    "execute",
    "show",
    "please run",
    "please do",
    "please show me",
    "show me",
)
_INDEX_CONFIRM_PATTERN = re.compile(
    r"\s*(?:" + "|".join(re.escape(prefix) for prefix in _INDEX_CONFIRM_PREFIXES) + r")?\s*(\d+)\s*",
    flags=re.IGNORECASE,
)


def _normalize_text(text: str) -> str:
    text = text.strip().lower()
    return re.sub(r"\s+", " ", text)


def _token_count(text: str) -> int:
    return len([part for part in text.split(" ") if part])


def looks_like_action_followup(reply_text: str) -> bool:
    """Return True for short replies that only make sense after an offer."""
    normalized = _normalize_text(reply_text)
    if normalized in _AFFIRM_ONLY or normalized in _CANCEL_WORDS:
        return True
    return _INDEX_CONFIRM_PATTERN.fullmatch(reply_text.strip()) is not None


def extract_action_offer(result: dict[str, Any]) -> list[dict[str, Any]] | None:
    """Return executable suggested actions from a skill result.

    A displayed offer must be executable later. For that reason a valid action
    needs a human label, a stable id, and a structured nested request. The chat
    adapters deliberately ignore shell-only or prose-only suggestions.
    """
    actions = result.get("suggested_actions")
    if not isinstance(actions, list):
        return None
    cleaned: list[dict[str, Any]] = []
    for action in actions:
        if not isinstance(action, dict):
            continue
        action_id = action.get("action_id")
        label = action.get("label")
        request = action.get("request")
        if not isinstance(action_id, str) or not action_id.strip():
            continue
        if not isinstance(label, str) or not label.strip():
            continue
        if not isinstance(request, dict):
            continue
        cleaned.append(action)
    return cleaned or None


def extract_chat_summary_lines(result: dict[str, Any]) -> list[str]:
    """Return concise skill-authored lines intended for chat display."""
    lines = result.get("chat_summary_lines")
    if not isinstance(lines, list):
        return []
    return [str(line).strip() for line in lines if str(line).strip()]


def load_bundle_fields(output_dir: str | Path | None) -> dict[str, Any]:
    """Load structured display/action fields from a skill output bundle."""
    if output_dir is None:
        return {}

    out_dir = Path(output_dir)
    if not out_dir.exists():
        return {}

    fields: dict[str, Any] = {}
    result_json_path = out_dir / "result.json"
    if result_json_path.exists():
        fields["result_json_path"] = str(result_json_path)
        try:
            payload = json.loads(result_json_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            payload = None
        if isinstance(payload, dict):
            fields["skill_result_json"] = payload
            # Field meanings:
            # - chat_summary_lines: concise text written by the skill for chat
            # - preferred_artifacts: files the UI should prioritize displaying
            # - suggested_actions: deterministic follow-up requests to offer
            # - workflow_state: skill-emitted state identity/lifecycle metadata
            # - contract_alerts: structured contract/path discrepancy alerts
            # - report_md: full markdown report embedded in result.json
            for key in (
                "chat_summary_lines",
                "preferred_artifacts",
                "suggested_actions",
                "workflow_state",
                "contract_alerts",
                "report_md",
            ):
                if key in payload:
                    fields[key] = (
                        normalise_contract_alerts(payload[key])
                        if key == "contract_alerts"
                        else payload[key]
                    )

    if "report_md" not in fields:
        for pattern in ("report.md", "*_report.md", "*.md"):
            for md_file in sorted(out_dir.glob(pattern)):
                if md_file.name.startswith("."):
                    continue
                try:
                    fields["report_md"] = md_file.read_text(encoding="utf-8")
                except OSError:
                    continue
                return fields
    return fields


def render_workflow_state_header(workflow_state: dict[str, Any] | None) -> str:
    """Render the optional cross-skill workflow lifecycle header."""
    if not isinstance(workflow_state, dict):
        return ""
    lifecycle = workflow_state.get("lifecycle")
    if not isinstance(lifecycle, str) or not lifecycle.strip():
        return ""
    lifecycle = lifecycle.strip()
    state_label = workflow_state.get("state_label")
    if isinstance(state_label, str) and state_label.strip():
        return f"State: {lifecycle} — {state_label.strip()}"
    return f"State: {lifecycle}"


def render_action_offer(
    actions: list[dict[str, Any]],
    *,
    workflow_state: dict[str, Any] | None = None,
) -> str:
    """Render actions as numbered chat choices."""
    header = render_workflow_state_header(workflow_state)
    if not actions:
        return header
    lines = []
    if header:
        lines.append(header)
    lines.append("I can do the next step for you:")
    for idx, action in enumerate(actions, start=1):
        suffix = ""
        estimate = action.get("estimate")
        if isinstance(estimate, str) and estimate.strip():
            suffix = f" ({estimate.strip()})"
        elif action.get("requires_confirmation") is False:
            suffix = " (safe refresh)"
        lines.append(f"{idx}. {action.get('label', f'Action {idx}')}{suffix}")
    if len(actions) == 1:
        lines.append("Reply with `yes`, `1`, the action name, or `cancel`.")
    else:
        choices = ", ".join(f"`{i}`" for i in range(1, len(actions) + 1))
        lines.append(f"Reply with {choices}, `yes 1`, the action name, or `cancel`.")
    return "\n".join(lines)


def is_cancel_reply(reply_text: str) -> bool:
    """Return True when a reply clearly declines the pending offer."""
    return _normalize_text(reply_text) in _CANCEL_WORDS


def _match_by_index(reply_text: str, actions: list[dict[str, Any]]) -> dict[str, Any] | None:
    m = _INDEX_CONFIRM_PATTERN.fullmatch(reply_text.strip())
    if not m:
        return None
    idx = int(m.group(1)) - 1
    if 0 <= idx < len(actions):
        return actions[idx]
    return None


def _match_by_label(reply_text: str, actions: list[dict[str, Any]]) -> dict[str, Any] | None:
    normalized = _normalize_text(reply_text)
    exact: list[dict[str, Any]] = []
    loose: list[dict[str, Any]] = []
    for action in actions:
        label = _normalize_text(str(action.get("label", "")))
        action_id = _normalize_text(str(action.get("action_id", "")))
        if normalized in {label, action_id}:
            exact.append(action)
        if (
            normalized
            and _token_count(normalized) > 1
            and _token_count(label) > 1
            and (normalized in label or label in normalized)
        ):
            loose.append(action)
    if len(exact) == 1:
        return exact[0]
    if len(loose) == 1:
        return loose[0]

    label_map = {
        _normalize_text(str(action.get("label", ""))): action for action in actions
    }
    close = difflib.get_close_matches(normalized, list(label_map), n=2, cutoff=0.86)
    if len(close) == 1:
        return label_map[close[0]]
    return None


def parse_action_reply(reply_text: str, actions: list[dict[str, Any]]) -> dict[str, Any]:
    """Classify a follow-up reply against stored suggested actions.

    The returned `confirmed` flag means the wording was an explicit approval
    ("yes", "go ahead", "please show me"). A numeric or label selection may
    still need a second confirmation when the action requires confirmation.
    """
    if not actions:
        return {"kind": "none"}
    normalized = _normalize_text(reply_text)
    if not normalized:
        return {"kind": "none"}
    if is_cancel_reply(normalized):
        return {"kind": "cancel"}

    if normalized in _AFFIRM_ONLY:
        if len(actions) == 1:
            return {"kind": "matched", "action": actions[0], "confirmed": True}
        return {"kind": "ambiguous"}

    index_match = _match_by_index(reply_text, actions)
    if index_match is not None:
        confirmed = re.fullmatch(r"\s*\d+\s*", reply_text) is None
        return {"kind": "matched", "action": index_match, "confirmed": confirmed}

    label_match = _match_by_label(reply_text, actions)
    if label_match is not None:
        return {"kind": "matched", "action": label_match, "confirmed": False}

    return {"kind": "none"}


def make_pending_action_entry(
    *,
    skill: str,
    actions: list[dict[str, Any]],
    source_summary: list[str] | None = None,
    source_output_dir: str | None = None,
    timestamp: float | None = None,
    selected_action: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Create the per-chat offer state stored by adapters."""
    return {
        "skill": skill,
        "timestamp": float(timestamp if timestamp is not None else time.time()),
        "actions": actions,
        "source_summary": list(source_summary or []),
        "source_output_dir": source_output_dir,
        "selected_action": selected_action,
    }


def is_pending_action_expired(
    entry: dict[str, Any],
    *,
    now_ts: float | None = None,
    ttl_seconds: int = DEFAULT_PENDING_ACTION_TTL_SECONDS,
) -> bool:
    """Return True when a stored offer is too old to accept safely."""
    try:
        ts = float(entry.get("timestamp", 0))
    except (TypeError, ValueError):
        return True
    now = float(now_ts if now_ts is not None else time.time())
    return now - ts > ttl_seconds


def choice_list_text(actions: list[dict[str, Any]]) -> str:
    choices = [str(i) for i in range(1, len(actions) + 1)]
    if not choices:
        return ""
    if len(choices) == 1:
        return choices[0]
    if len(choices) == 2:
        return f"{choices[0]} or {choices[1]}"
    return ", ".join(choices[:-1]) + f", or {choices[-1]}"


def execute_stored_action(
    pending_entry: dict[str, Any],
    action: dict[str, Any],
    *,
    runner: Callable[..., dict[str, Any]],
    output_root: str | Path,
) -> dict[str, Any]:
    """Execute a selected stored action through the normal ClawBio runner.

    The nested request is written to a temporary JSON file and passed as
    `--input`. `shell_line` is never executed from this path.
    """
    request = action.get("request")
    if not isinstance(request, dict):
        raise ValueError("stored action is missing a structured request payload")
    skill = str(pending_entry.get("skill") or "").strip()
    if not skill:
        raise ValueError("pending action entry is missing the source skill")

    with tempfile.TemporaryDirectory(prefix="clawbio_action_") as tmpdir:
        request_path = Path(tmpdir) / "request.json"
        request_path.write_text(
            json.dumps(request, indent=2, ensure_ascii=True) + "\n",
            encoding="utf-8",
        )

        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_dir = Path(output_root) / f"{skill}_{ts}"
        runner_kwargs: dict[str, Any] = {
            "skill_name": skill,
            "input_path": str(request_path),
            "output_dir": str(output_dir),
            "demo": False,
        }
        timeout_secs = action.get("timeout_secs")
        if isinstance(timeout_secs, int) and timeout_secs > 0:
            runner_kwargs["timeout"] = timeout_secs
        return runner(**runner_kwargs)
