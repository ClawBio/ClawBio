#!/usr/bin/env python3
"""
run_compat.py - Check that a ClawBio skill works in third-party coding agents.

For each agent CLI, the harness builds a throwaway project that contains only
one ClawBio skill folder (copied into that agent's skills directory, exactly as
the README tells users to do), then runs two headless checks:

  discover  Ask the agent, without reading files, which skills it has.
            Passes if the skill name appears in its answer.
  execute   Ask for the task in plain words, without naming the skill.
            Passes if the agent routed to the skill and wrote the expected
            output file (a pharmacogenomic report with real gene calls).

Agents whose CLI is not on PATH are reported as "not installed", not failed.
Each agent needs its own API key in the environment; the harness never reads
or prints keys.

Usage:
    python scripts/agent_compat/run_compat.py                    # all agents
    python scripts/agent_compat/run_compat.py --agent claude-code
    python scripts/agent_compat/run_compat.py --skills-dir .agents/skills \
        --agent claude-code                                      # try another path
    python scripts/agent_compat/run_compat.py --json results.json
"""
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass, field
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]

SKILL = "pharmgx-reporter"
DEMO_INPUT = "demo_patient.txt"
EXPECTED_OUTPUT = "out/report.md"
EXPECTED_MARKERS = ("CYP2C19", "CYP2D6")

DISCOVER_PROMPT = (
    "Without reading or listing any files, name the skills available to you "
    "in this session. Names only, one per line."
)
EXECUTE_PROMPT = (
    "I have 23andMe raw data in {input}. Make me a pharmacogenomic report "
    "and save it in ./out. Do not write your own analysis code; use the "
    "tools you have."
)


@dataclass
class Agent:
    name: str
    binary: str
    # Project-relative directory the agent loads skills from. Where a vendor
    # documents several, this is the one we want to verify first; override
    # with --skills-dir to test another.
    skills_dir: str
    # argv after the binary; "{prompt}" is substituted.
    argv: list[str]
    notes: str = ""


AGENTS: dict[str, Agent] = {
    a.name: a
    for a in [
        Agent(
            "claude-code", "claude", ".claude/skills",
            ["-p", "{prompt}", "--allowedTools", "Bash,Read,Write,Glob,Grep,Skill"],
            notes="Project .agents/skills was not discovered in 2.1.328.",
        ),
        Agent(
            "codex-cli", "codex", ".agents/skills",
            ["exec", "--skip-git-repo-check", "--full-auto", "{prompt}"],
        ),
        Agent(
            "gemini-cli", "gemini", ".gemini/skills",
            ["-p", "{prompt}", "--yolo"],
            notes="Skills path unverified; try --skills-dir .agents/skills too.",
        ),
        Agent(
            "opencode", "opencode", ".opencode/skills",
            ["run", "{prompt}"],
            notes="Skills path unverified; try --skills-dir .agents/skills too.",
        ),
    ]
}


@dataclass
class CheckResult:
    status: str  # pass | fail | error | not installed
    seconds: float = 0.0
    detail: str = ""


@dataclass
class AgentResult:
    agent: str
    skills_dir: str
    version: str = ""
    checks: dict[str, CheckResult] = field(default_factory=dict)


def _run(argv: list[str], cwd: Path, timeout: int) -> tuple[int, str, float]:
    start = time.monotonic()
    try:
        proc = subprocess.run(
            argv, cwd=cwd, capture_output=True, text=True, timeout=timeout,
            stdin=subprocess.DEVNULL,
        )
        out = proc.stdout + proc.stderr
        code = proc.returncode
    except subprocess.TimeoutExpired as exc:
        out = f"timeout after {timeout}s\n{exc.stdout or ''}"
        code = -1
    return code, out, time.monotonic() - start


def _make_project(base: Path, skills_dir: str) -> Path:
    project = Path(tempfile.mkdtemp(prefix="compat-", dir=base))
    dest = project / skills_dir / SKILL
    shutil.copytree(
        PROJECT_ROOT / "skills" / SKILL, dest,
        ignore=shutil.ignore_patterns("tests", "evals", "demo_report", "__pycache__"),
    )
    shutil.copy(dest / DEMO_INPUT, project / "patient.txt")
    return project


def _argv(agent: Agent, prompt: str) -> list[str]:
    return [agent.binary] + [a.replace("{prompt}", prompt) for a in agent.argv]


def run_agent(agent: Agent, skills_dir: str, base: Path, timeout: int,
              keep: bool) -> AgentResult:
    result = AgentResult(agent.name, skills_dir)
    if shutil.which(agent.binary) is None:
        for check in ("discover", "execute"):
            result.checks[check] = CheckResult("not installed")
        return result

    _, ver, _ = _run([agent.binary, "--version"], base, 30)
    result.version = ver.strip().splitlines()[0] if ver.strip() else ""

    # discover
    project = _make_project(base, skills_dir)
    code, out, secs = _run(_argv(agent, DISCOVER_PROMPT), project, timeout)
    if code != 0:
        result.checks["discover"] = CheckResult("error", secs, out[-500:])
    elif SKILL in out:
        result.checks["discover"] = CheckResult("pass", secs)
    else:
        result.checks["discover"] = CheckResult("fail", secs, out[-500:])

    # execute (fresh project so discovery output cannot leak in)
    project2 = _make_project(base, skills_dir)
    prompt = EXECUTE_PROMPT.format(input="patient.txt")
    code, out, secs = _run(_argv(agent, prompt), project2, timeout)
    report = project2 / EXPECTED_OUTPUT
    if report.exists():
        text = report.read_text(errors="replace")
        missing = [m for m in EXPECTED_MARKERS if m not in text]
        if missing:
            result.checks["execute"] = CheckResult(
                "fail", secs, f"report written but missing {missing}")
        else:
            result.checks["execute"] = CheckResult("pass", secs)
    else:
        status = "error" if code != 0 else "fail"
        result.checks["execute"] = CheckResult(
            status, secs, f"no {EXPECTED_OUTPUT}\n{out[-800:]}")

    if not keep:
        shutil.rmtree(project, ignore_errors=True)
        shutil.rmtree(project2, ignore_errors=True)
    return result


def to_markdown(results: list[AgentResult]) -> str:
    lines = [
        f"## Agent compatibility: `{SKILL}`",
        "",
        "| Agent | Version | Skills dir | Discover | Execute |",
        "|---|---|---|---|---|",
    ]
    for r in results:
        d, e = r.checks["discover"], r.checks["execute"]
        lines.append(
            f"| {r.agent} | {r.version or '-'} | `{r.skills_dir}` | "
            f"{d.status} | {e.status} |"
        )
    details = [
        f"### {r.agent} / {name}\n```\n{c.detail.strip()}\n```"
        for r in results for name, c in r.checks.items()
        if c.status in ("fail", "error") and c.detail
    ]
    if details:
        lines += ["", "<details><summary>Failure output</summary>", ""]
        lines += details + ["", "</details>"]
    return "\n".join(lines) + "\n"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--agent", action="append", choices=sorted(AGENTS),
                    help="Agent to test (repeatable). Default: all.")
    ap.add_argument("--skills-dir", help="Override the agent's skills dir.")
    ap.add_argument("--timeout", type=int, default=600)
    ap.add_argument("--json", type=Path, help="Write raw results here.")
    ap.add_argument("--keep", action="store_true",
                    help="Keep the temp projects for inspection.")
    args = ap.parse_args(argv)

    names = args.agent or sorted(AGENTS)
    base = Path(tempfile.mkdtemp(prefix="clawbio-agent-compat-"))
    results = []
    for name in names:
        agent = AGENTS[name]
        skills_dir = args.skills_dir or agent.skills_dir
        print(f"[{name}] testing with {skills_dir} ...", file=sys.stderr)
        results.append(run_agent(agent, skills_dir, base, args.timeout, args.keep))

    print(to_markdown(results))
    if args.json:
        args.json.write_text(json.dumps(
            [{"agent": r.agent, "version": r.version, "skills_dir": r.skills_dir,
              "checks": {k: vars(v) for k, v in r.checks.items()}}
             for r in results], indent=2))
    if args.keep:
        print(f"Temp projects kept in {base}", file=sys.stderr)
    else:
        shutil.rmtree(base, ignore_errors=True)

    tested = [c for r in results for c in r.checks.values()
              if c.status != "not installed"]
    return 0 if tested and all(c.status == "pass" for c in tested) else 1


if __name__ == "__main__":
    sys.exit(main())
