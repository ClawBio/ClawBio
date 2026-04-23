"""Thin `claude -p` wrapper. Each call is an independent, stateless subagent."""
from __future__ import annotations

import os
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path


@dataclass
class DispatchRequest:
    role: str
    prompt: str
    model: str = "sonnet"
    allowed_tools: list[str] | None = None
    cwd: str | None = None
    timeout: int = 28800


class Dispatcher:
    def dispatch(self, req: DispatchRequest) -> str:
        raise NotImplementedError


class ClaudeCLIDispatcher(Dispatcher):
    def __init__(self, default_model: str = "sonnet", binary: str = "claude"):
        self.default_model = default_model
        self.binary = binary

    def dispatch(self, req: DispatchRequest) -> str:
        cmd = [
            self.binary, "-p",
            "--model", req.model or self.default_model,
            "--output-format", "text",
            "--permission-mode", "bypassPermissions",
            "--no-session-persistence",
        ]
        if req.allowed_tools:
            cmd += ["--allowed-tools", ",".join(req.allowed_tools)]
        debug_dir = Path(os.environ.get("CLAWPATHY_DEBUG_DIR", "/tmp/clawpathy_debug"))
        debug_dir.mkdir(parents=True, exist_ok=True)
        tag = f"{req.role}-{int(time.time())}"
        proc = subprocess.run(
            cmd, input=req.prompt, capture_output=True, text=True,
            cwd=req.cwd, timeout=req.timeout,
        )
        (debug_dir / f"{tag}.argv").write_text(repr(cmd) + f"\nprompt_len={len(req.prompt)}\nreturncode={proc.returncode}\n")
        (debug_dir / f"{tag}.stdout").write_text(proc.stdout or "")
        (debug_dir / f"{tag}.stderr").write_text(proc.stderr or "")
        if proc.returncode != 0:
            raise RuntimeError(
                f"[{req.role}] claude -p failed rc={proc.returncode} "
                f"stderr_len={len(proc.stderr or '')} stdout_len={len(proc.stdout or '')} "
                f"debug={debug_dir}/{tag}.* | stderr_tail: {(proc.stderr or '')[-2000:]!r} "
                f"stdout_tail: {(proc.stdout or '')[-500:]!r}"
            )
        return proc.stdout
