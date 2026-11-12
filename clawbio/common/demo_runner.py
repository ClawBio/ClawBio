"""Run a skill's catalog demo and read its outputs in a form two runs can compare.

Shared by scripts/nightly_demo_sweep.py (does every demo still run?) and the
per-skill demo snapshot tests (did a code change move a demo number?), so both
exercise exactly the demo that skills/catalog.json advertises.
"""

from __future__ import annotations

import json
import re
import shlex
import subprocess
import sys
import time
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
CATALOG = PROJECT_ROOT / "skills" / "catalog.json"
DEFAULT_TIMEOUT = 180

# Fields that differ between two runs of the same demo (measured across every
# CI demo on 2026-11-05). Anything not listed here must match exactly.
_TIMESTAMP = re.compile(r"\d{4}-\d{2}-\d{2}[T ]\d{2}:\d{2}(:\d{2})?([.,]\d+)?(Z|[+-]\d{2}:?\d{2})?")


def with_output(argv: list[str], output: Path) -> list[str]:
    """Return argv with any --output removed and `--output <output>` appended."""
    kept, skip_value = [], False
    for i, arg in enumerate(argv):
        if skip_value:
            skip_value = False
            continue
        if arg == "--output":
            skip_value = i + 1 < len(argv) and not argv[i + 1].startswith("--")
            continue
        if arg.startswith("--output="):
            continue
        kept.append(arg)
    return kept + ["--output", str(output)]


def _demo_command(name: str) -> str:
    for skill in json.loads(CATALOG.read_text())["skills"]:
        if skill["name"] == name and skill.get("demo_command"):
            return skill["demo_command"]
    raise LookupError(f"{name} has no demo_command in {CATALOG.name}")


def run_demo(name: str, output_dir: Path, *, timeout: int = DEFAULT_TIMEOUT) -> dict:
    """Run the catalog demo for skill `name` with `--output output_dir`. Never raises.

    Returns {name, command, exit_code, elapsed, stderr_tail, passed}. exit_code is
    -1 on timeout and -2 when the demo could not be started (unknown skill,
    unparseable command); stderr_tail then carries the reason.
    """
    t0 = time.time()
    cmd = ""
    try:
        cmd = _demo_command(name)
        # demo_command is built from skill filenames, which contributors
        # choose, so it is split into argv and never handed to a shell.
        argv = shlex.split(cmd)
        if argv and argv[0] == "python":
            argv[0] = sys.executable
        # Every skill must write only under --output (AGENTS.md), so every demo
        # gets one. A skill that rejects it fails here, which is the point.
        argv = with_output(argv, Path(output_dir))
        cmd = shlex.join(argv)
        proc = subprocess.run(
            argv,
            cwd=str(PROJECT_ROOT),
            capture_output=True,
            text=True,
            timeout=timeout,
        )
        exit_code, stderr_tail = proc.returncode, proc.stderr[-500:] if proc.stderr else ""
    except subprocess.TimeoutExpired:
        exit_code, stderr_tail = -1, f"TIMEOUT after {timeout}s"
    except Exception as e:
        exit_code, stderr_tail = -2, str(e)
    return {
        "name": name,
        "command": cmd,
        "exit_code": exit_code,
        "elapsed": round(time.time() - t0, 1),
        "stderr_tail": stderr_tail,
        "passed": exit_code == 0,
    }


MAX_SNAPSHOT_BYTES = 50_000  # a larger text output would swamp an inline snapshot
_NOT_SNAPSHOTTED_DIRS = {"reproducibility"}  # replay paths and checksums of timestamped files
_NOT_SNAPSHOTTED_SUFFIXES = {".log"}


def demo_outputs(name: str, output_dir: Path, *, timeout: int = DEFAULT_TIMEOUT) -> dict[str, str]:
    """Run skill `name`'s catalog demo under output_dir and return {relative path: masked text} of its outputs.

    Covers every UTF-8 text file up to MAX_SNAPSHOT_BYTES, except reproducibility/
    and logs. Binary files (figures, PDFs, gzip) are left out: their bytes change
    with library versions. Raises RuntimeError with the demo's stderr if it fails.
    """
    # A fixed folder name: reports that echo `--output <name>` then read the same in every run.
    output_dir = Path(output_dir) / "demo_output"
    result = run_demo(name, output_dir, timeout=timeout)
    if not result["passed"]:
        raise RuntimeError(f"{name} demo failed (exit {result['exit_code']}): {result['stderr_tail']}")
    outputs = {}
    for path in sorted(output_dir.rglob("*")):
        rel = path.relative_to(output_dir)
        if (not path.is_file() or rel.parts[0] in _NOT_SNAPSHOTTED_DIRS
                or path.suffix in _NOT_SNAPSHOTTED_SUFFIXES or path.stat().st_size > MAX_SNAPSHOT_BYTES):
            continue
        try:
            text = path.read_text(encoding="utf-8")
        except UnicodeDecodeError:
            continue
        outputs[rel.as_posix()] = mask(text, output_dir)
    return outputs


def mask(text: str, output_dir: Path) -> str:
    """Replace the output directory and timestamps so two runs of one demo compare equal."""
    text = text.replace(str(output_dir), "<OUTPUT_DIR>")
    return _TIMESTAMP.sub("<TIMESTAMP>", text)
