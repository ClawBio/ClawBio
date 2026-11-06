"""Run headless Claude Code against this repo and parse what it did.

The agent evals assert on the trajectory (which commands ran, which files were
read) as well as the final answer, so this keeps both.
"""
from __future__ import annotations

import json
import os
import re
import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]

# The agent may read anything but only execute ClawBio. Skills can still call
# external APIs; SETTINGS sandboxes them so a blocked call shows up in Trace.network_blocked.
ALLOWED_TOOLS = [
    "Read", "Glob", "Grep",
    # A glob, not `skills/:*`: the legacy `:*` prefix needs a word boundary, so it never matches a script path.
    *(f"Bash({py} {target})"
      for py in ("python", "python3", "uv run python")
      for target in ("clawbio.py:*", "skills/*")),
    # The extras and benchmark invocations AGENTS.md tells agents to type.
    "Bash(uv run --extra * python skills/*)", "Bash(python tests/benchmark/*)",
    "Bash(uv run --extra * python tests/benchmark/*)",
    # Helper commands agents chain with a skill run, e.g. `mkdir -p out && ... | head -60`.
    # A denied helper denies the whole command, which then reads as the agent skipping the skill.
    "Bash(head:*)", "Bash(tail:*)", "Bash(cat:*)", "Bash(ls:*)", "Bash(mkdir:*)", "Bash(grep:*)",
]
# Agents pick --output /tmp/<name> (SKILL.md examples do) and then `ls`/`grep` it. dontAsk
# denies shell reads outside the project, and any command with a shell variable such as
# `ls $TMPDIR/out` whatever its value, which denies the whole chained skill run;
# a real user session would approve it.
EXTRA_DIRS = ["/tmp"]
# pharmgx enrichment; it sends gene symbols from its own drug table, never genotypes.
NETWORK_ALLOWED = ["api.clinpgx.org"]
SETTINGS = {
    # The sandbox only reports a blocked host, so every host not listed is blocked.
    "sandbox": {"enabled": True, "autoAllowBashIfSandboxed": False, "allowUnsandboxedCommands": False,
                "network": {"allowedDomains": NETWORK_ALLOWED}},
    "permissions": {
        # `skills/*` also matches `skills/../anything.py`.
        "deny": ["Bash(*../*)"],
        # Sandbox write access comes from Edit rules, not --add-dir: skill output in /tmp,
        # `uv run`'s cache, the ClinPGx cache.
        "allow": ["Edit(//tmp/**)", "Edit(//private/tmp/**)", "Edit(~/.cache/uv/**)", "Edit(~/.clawbio/**)"],
    },
}
NETWORK_DENY = re.compile(r"deny network-outbound (\S+)")


# What `clawbio run` and a crashing script print; a pipe like `skill | tail` exits 0 either way.
SKILL_FAILED = re.compile(r"Status:\s+FAILED|Traceback \(most recent call last\)|^\s*Error:", re.M)
# Running a skill, as opposed to `cat skills/<name>/SKILL.md`; also matches `uv run python ...`.
SKILL_RUN = re.compile(r"\bpython3?\s+(clawbio\.py\s+run|skills/)")


@dataclass
class BashCall:
    command: str
    # "ok", "failed" or "denied" (never ran). A chain reports its last command, so
    # `skill && grep x` is failed when grep finds nothing.
    outcome: str
    output: str = ""  # as the agent saw it


@dataclass
class Trace:
    calls: list[BashCall] = field(default_factory=list)  # every Bash call, in order
    network_blocked: list[str] = field(default_factory=list)  # host:port the sandbox refused
    paths_touched: list[str] = field(default_factory=list)
    final: str | None = None

    def ran_ok(self, *needles: str) -> bool:
        """Some skill run containing every needle exited 0 and its output shows no crash."""
        return any(c.outcome == "ok" and SKILL_RUN.search(c.command) and all(n in c.command for n in needles)
                   and not SKILL_FAILED.search(c.output) for c in self.calls)


def parse_stream(lines) -> Trace:
    """Parse `claude -p --output-format stream-json --verbose` output."""
    trace = Trace()
    bash_by_id: list[tuple[str | None, str]] = []
    denied_ids: set[str] = set()
    error_ids: set[str] = set()
    outputs: dict[str, str] = {}
    for line in lines:
        try:
            event = json.loads(line)
        except json.JSONDecodeError:
            continue
        kind = event.get("type")
        if kind == "assistant":
            for block in event["message"]["content"]:
                if block.get("type") != "tool_use":
                    continue
                inp = block.get("input", {})
                if block["name"] == "Bash":
                    bash_by_id.append((block.get("id"), inp.get("command", "")))
                # A Grep pattern is search text, not a path; a Glob pattern is a path.
                path = inp.get("file_path") or inp.get("path") or (block["name"] == "Glob" and inp.get("pattern"))
                if block["name"] in ("Read", "Grep", "Glob") and path:
                    trace.paths_touched.append(path)
        elif kind == "user":
            for block in event["message"].get("content", []):
                if not isinstance(block, dict) or block.get("type") != "tool_result":
                    continue
                if block.get("is_error"):
                    error_ids.add(block.get("tool_use_id"))
                content = block.get("content")
                if isinstance(content, list):
                    content = "\n".join(c.get("text", "") for c in content if isinstance(c, dict))
                outputs[block.get("tool_use_id")] = content or ""
                # Only a command's output; a Read of a file quoting a sandbox message is not a leak.
                if block.get("tool_use_id") in {tid for tid, _ in bash_by_id}:
                    trace.network_blocked += NETWORK_DENY.findall(content or "")
        elif kind == "system" and event.get("subtype") == "permission_denied":
            denied_ids.add(event.get("tool_use_id"))
        elif kind == "result":
            trace.final = event.get("result")
    # A denial arrives after its tool_use, so outcomes are only known at the end.
    # No result at all means the run died mid-call, which is not a success.
    trace.calls = [BashCall(cmd, "denied" if tid in denied_ids
                            else "failed" if tid in error_ids or tid not in outputs else "ok",
                            outputs.get(tid, "")) for tid, cmd in bash_by_id]
    return trace


def run_agent(prompt: str, log_path: Path, timeout: int = 600) -> Trace:
    """Run one prompt through headless Claude Code and keep the raw log.

    --setting-sources project keeps the user's own plugins and hooks out, so the
    agent sees what any contributor's agent sees. It also stops Claude Code
    (2.1.296+) loading AGENTS.md natively, and without AGENTS.md the agent skips the
    skill and answers from memory, so AGENTS.md is passed in explicitly.
    The user's ~/.claude/CLAUDE.md still loads; there is no flag to drop it
    without logging out. Each run gets a fresh clone of HEAD, so no earlier
    run's files, git history or auto-memory (keyed by path) reach the agent;
    it refuses to run with uncommitted changes, which the clone would miss.
    The model is pinned so runs compare; override with CLAWBIO_EVAL_MODEL.
    The raw log can hold local paths and ~/.claude/CLAUDE.md text, so it must live
    outside the repo; share summary numbers, not transcripts.
    """
    if log_path.resolve().is_relative_to(REPO):
        raise ValueError(f"write logs outside the repo, not {log_path}; they may hold private paths")
    dirty = subprocess.run(["git", "-C", str(REPO), "status", "--porcelain"],
                           capture_output=True, text=True, check=True).stdout
    if dirty:
        raise RuntimeError(f"commit first; the agent runs on a clone of HEAD and would miss:\n{dirty}")
    with tempfile.TemporaryDirectory() as tmp:
        clone = Path(tmp) / "repo"
        subprocess.run(["git", "clone", "--quiet", "--local", str(REPO), str(clone)], check=True)
        # Reuse the repo's venv; REPO is clean, so its editable clawbio matches the clone.
        env = {**os.environ, "UV_PROJECT_ENVIRONMENT": str(REPO / ".venv"), "UV_NO_SYNC": "1"}
        try:
            proc = subprocess.run(_claude_cmd(prompt, clone), cwd=clone, env=env, capture_output=True,
                                  text=True, timeout=timeout, stdin=subprocess.DEVNULL)
        except subprocess.TimeoutExpired as e:
            log_path.write_text(_text(e.stdout) + _text(e.stderr))
            raise
    # stderr is where a failed sign-in or crash explains itself; parse_stream skips non-JSON lines.
    log_path.write_text(proc.stdout + proc.stderr)
    return parse_stream(proc.stdout.splitlines())


def _text(out: str | bytes | None) -> str:
    # TimeoutExpired keeps bytes even with text=True.
    return out.decode(errors="replace") if isinstance(out, bytes) else out or ""


def run_trials(prompt: str, log_dir: Path, n: int, timeout: int = 600) -> list[Trace]:
    """Run a prompt n times; agents are non-deterministic, so score the pass rate across trials."""
    return [run_agent(prompt, log_dir / f"trial-{i}.jsonl", timeout) for i in range(n)]


def _claude_cmd(prompt: str, clone: Path) -> list[str]:
    return [
        "claude", "-p", prompt,
        "--output-format", "stream-json", "--verbose",
        "--model", os.environ.get("CLAWBIO_EVAL_MODEL", "sonnet"),
        "--setting-sources", "project",
        "--append-system-prompt-file", str(clone / "AGENTS.md"),
        "--no-session-persistence",
        "--permission-mode", "dontAsk",
        "--allowedTools", *ALLOWED_TOOLS,
        "--add-dir", *EXTRA_DIRS,
        "--settings", json.dumps(SETTINGS),
    ]
