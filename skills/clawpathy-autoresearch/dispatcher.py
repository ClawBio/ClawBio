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
    role: str
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


class ClaudeCLIDispatcher(Dispatcher):
    """Dispatches each request as a headless `claude -p` call.

    Each dispatch runs in a fresh Claude CLI process with --bare so no
    project CLAUDE.md or auto-memory leaks into the subagent. The tool
    whitelist is passed via --allowed-tools and a permissive mode so the
    subagent can act without interactive approval.
    """

    def __init__(self, model: str = "haiku", extra_dirs: list[str] | None = None,
                 timeout_seconds: int = 900):
        self.model = model
        self.extra_dirs = extra_dirs or []
        self.timeout_seconds = timeout_seconds

    def dispatch(self, request: DispatchRequest) -> str:
        import os
        import subprocess
        cmd = [
            "claude", "-p", request.prompt,
            "--model", self.model,
            "--output-format", "text",
            "--permission-mode", "bypassPermissions",
            "--no-session-persistence",
        ]
        if request.allowed_tools:
            cmd += ["--allowed-tools", ",".join(request.allowed_tools)]
        for d in self.extra_dirs:
            cmd += ["--add-dir", d]
        env = os.environ.copy()
        # Stale/invalid ANTHROPIC_API_KEY shadows OAuth; strip it so the CLI
        # falls back to keychain credentials.
        env.pop("ANTHROPIC_API_KEY", None)
        proc = subprocess.run(
            cmd, capture_output=True, text=True,
            timeout=self.timeout_seconds, env=env,
        )
        if proc.returncode != 0:
            raise RuntimeError(
                f"claude -p failed (code {proc.returncode}): {proc.stderr[-500:]}"
            )
        return proc.stdout
