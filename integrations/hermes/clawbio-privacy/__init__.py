"""Hermes Agent plugin: keep genomic file contents out of the model's context.

A thin adapter over ClawBio's genomic-data policy (clawbio/privacy.py): four hooks and one
middleware, stdlib only. The policy is loaded from the ClawBio checkout: this file's resolved path
(plugin symlinked from integrations/hermes/), else CLAWBIO_DIR (copied plugin). Not found = fail closed.

- pre_gateway_dispatch   Hermes pastes small text uploads into the prompt. Swap a
                         genomic one for its saved path and a --sealed instruction.
- pre_tool_call          Block tools that would read or send a genomic file or ClawBio's output
                         (reports), except a single `clawbio.py run ... --sealed` command.
- transform_tool_result  Withhold tool output that looks like genotype data, and
                         remember a sealed run's output_dir.
- transform_llm_output   Drop every MEDIA: tag the model wrote for ClawBio output, from any turn,
                         and attach this turn's report (report.md, nutrigx_report.md) with a clean MEDIA: line.
- llm_request            Scrub genomic text from every model request. Hermes queues a
  (middleware)           message sent while a turn is running and replays it without
                         pre_gateway_dispatch, so this catches those uploads.

Hermes runs pre_gateway_dispatch, transform_tool_result and middleware fail-open, so
each catches its own errors and withholds rather than passes content through.
"""

from __future__ import annotations

import importlib.util
import json
import os
import re
from functools import cache
from pathlib import Path
from types import ModuleType
from typing import Any

FILE_TOOLS = {"read_file", "search_files", "patch", "write_file"}
SHELL_TOOLS = {"terminal", "execute_code"}

# Same text as clawbio.privacy.WITHHELD; defined here too because it is needed when the policy is missing.
WITHHELD = "[withheld by clawbio-privacy: genomic data is not shown to the model]"
_NO_POLICY = (f"{WITHHELD} clawbio-privacy cannot find ClawBio's policy (clawbio/privacy.py), so every file "
              "and command is blocked. Symlink this plugin from a ClawBio checkout, or set CLAWBIO_DIR.")
_OUTPUT_DIR_RE = re.compile(r'output_dir\\?"\s*:\s*\\?"([^"\\]+)')
# A MEDIA: tag plus any markdown or spoiler wrapping the model put around it.
_MEDIA_RE = re.compile(r"[|`*_]*MEDIA:\s*[|`*_]*([^\s|`*]+)[|`*_]*")

# ponytail: in-memory per session; a gateway restart between the run and the reply drops the attachment
_PENDING_RUNS: dict[str, list[str]] = {}


def _hermes_home() -> Path | None:
    """The profile serving this turn; Hermes sets it per routed profile."""
    try:
        from hermes_constants import get_hermes_home  # type: ignore[import-not-found]
        return Path(get_hermes_home())
    except Exception:
        return None


def _env_file_value(path: Path, key: str) -> str | None:
    try:
        lines = path.read_text().splitlines()
    except OSError:
        return None
    for line in lines:
        k, sep, v = line.strip().partition("=")
        if sep and k.strip() == key:
            return v.strip().strip("'\"") or None
    return None


def _clawbio_dirs() -> list[str]:
    """CLAWBIO_DIR from the routed profile's .env first, then the process env.

    A multiplexed gateway never merges a profile's .env into os.environ.
    """
    home = _hermes_home()
    dirs = [_env_file_value(home / ".env", "CLAWBIO_DIR") if home else None, os.environ.get("CLAWBIO_DIR")]
    return list(dict.fromkeys(d.rstrip("/") for d in dirs if d))


@cache
def _load(path: str) -> ModuleType:
    spec = importlib.util.spec_from_file_location("clawbio_privacy_policy", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _policy() -> ModuleType:
    """clawbio/privacy.py from this plugin's checkout (symlinked install), else from CLAWBIO_DIR (copied install)."""
    here = Path(__file__).resolve().parents
    roots = ([here[3]] if len(here) > 3 else []) + [Path(d).expanduser() for d in _clawbio_dirs()]
    for root in roots:
        candidate = root / "clawbio" / "privacy.py"
        if candidate.is_file():
            return _load(str(candidate))
    raise LookupError(_NO_POLICY)


def _block(message: str) -> dict:
    return {"action": "block", "message": message}


def on_pre_gateway_dispatch(event: Any = None, **_: Any) -> dict | None:
    text = getattr(event, "text", None) or ""
    paths = list(getattr(event, "media_urls", None) or [])
    if not text or not paths:
        return None
    try:
        new_text = _policy().replace_uploads(text, paths, _clawbio_dirs())
        return None if new_text == text else {"action": "rewrite", "text": new_text}
    except Exception:
        # ponytail: drops the whole message on a bug; a precise fallback is not worth the risk
        return {"action": "rewrite", "text": WITHHELD}


def on_pre_tool_call(tool_name: str = "", args: Any = None, **_: Any) -> dict | None:
    args = args if isinstance(args, dict) else {}
    try:
        policy, dirs = _policy(), _clawbio_dirs()
        if tool_name in FILE_TOOLS:
            path = str(args.get("path") or "")
            message = policy.check_path(path, dirs) if path else None
        elif tool_name in SHELL_TOOLS:
            message = policy.check_command(str(args.get("command") or args.get("code") or ""), dirs)
        elif tool_name.startswith("connectors__"):
            # Hermes 0.22 connector tools upload any local path they are given; connectors__execute batches them.
            message = next((m for v in _strings(args) if (m := policy.check_path(v, dirs))), None)
        else:
            message = None
    except LookupError:
        message = _NO_POLICY
    except Exception:
        # An older clawbio/privacy.py without these checks, or a bug in one: block rather than let the tool run.
        message = f"{WITHHELD} (privacy check failed; deploy this plugin and its ClawBio checkout together)"
    return _block(message) if message else None


def _strings(obj: Any):
    """Every string in a nested tool-argument structure."""
    if isinstance(obj, str):
        yield obj
    elif isinstance(obj, dict):
        for v in obj.values():
            yield from _strings(v)
    elif isinstance(obj, list):
        for v in obj:
            yield from _strings(v)


def on_transform_tool_result(tool_name: str = "", result: Any = None, session_id: str = "", **_: Any) -> str | None:
    if tool_name not in FILE_TOOLS | SHELL_TOOLS or not isinstance(result, str):
        return None
    try:
        if tool_name in SHELL_TOOLS and "sealed" in result:
            _PENDING_RUNS.setdefault(session_id, []).extend(_OUTPUT_DIR_RE.findall(result))
        return WITHHELD if _policy().is_genomic_text(_tool_text(result)) else None
    except Exception:
        return WITHHELD


def _tool_text(result: str) -> str:
    """The text inside Hermes's JSON wrapper ({"output": ...}, {"content": ...}): escaped, its newlines hide headers."""
    try:
        data = json.loads(result)
    except ValueError:
        return result
    return "\n".join(v for v in data.values() if isinstance(v, str)) if isinstance(data, dict) else result


def on_transform_llm_output(response_text: str = "", session_id: str = "", **_: Any) -> str | None:
    runs = _PENDING_RUNS.pop(session_id, [])
    dirs = _clawbio_dirs()
    try:
        under_output = _policy().is_under_output
    except Exception:  # no usable policy: no runs were allowed, so nothing to attach or drop
        return None
    # Most skills write report.md; some prefix it (nutrigx_report.md).
    reports = [str(md) for d in runs for md in sorted(Path(d).glob("*report.md"))[:1]]

    # Remove every MEDIA: tag the model wrote for ClawBio output, from this run or an earlier one.
    def drop(m: re.Match) -> str:
        path = m.group(1)
        return "" if under_output(path, dirs) or any(path.startswith(r) for r in runs) else m.group(0)

    text = _MEDIA_RE.sub(drop, response_text)
    if text == response_text and not reports:
        return None
    # Attach this turn's report.md files, each on its own line. Hermes always delivers a bare MEDIA: line.
    return text.rstrip() + "".join(f"\nMEDIA:{md}" for md in reports)


def _walk(obj: Any, fix) -> Any:
    """Apply ``fix`` to every message text: ``content`` strings and content parts' ``text``."""
    if isinstance(obj, dict):
        return {k: fix(v) if k in ("content", "text") and isinstance(v, str) else _walk(v, fix) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_walk(v, fix) for v in obj]
    return obj


def on_llm_request(request: Any = None, **_: Any) -> dict | None:
    if not isinstance(request, dict):
        return None
    try:
        dirs = _clawbio_dirs()
        policy = _policy()
        scrubbed = _walk(request, lambda text: policy.scrub_text(text, dirs))
    except Exception:
        # ponytail: withholds every message on a bug, so the turn is useless but nothing leaks
        scrubbed = _walk(request, lambda _text: WITHHELD)
    if scrubbed == request:
        return None
    return {"request": scrubbed, "source": "clawbio-privacy", "reason": "genomic data withheld"}


def register(ctx) -> None:
    ctx.register_hook("pre_gateway_dispatch", on_pre_gateway_dispatch)
    ctx.register_hook("pre_tool_call", on_pre_tool_call)
    ctx.register_hook("transform_tool_result", on_transform_tool_result)
    ctx.register_hook("transform_llm_output", on_transform_llm_output)
    ctx.register_middleware("llm_request", on_llm_request)
