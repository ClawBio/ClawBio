"""Run or submit a generated job script — only when explicitly asked.

Shared by the skills that emit runnable scripts (`sra-fetch`, and the
`download-script` command on the archive-fetch skills).

**These functions execute code and move data over the network.** They are never
reached by default: a skill generates its script and stops. `--run` and
`--submit` are opt-in flags, and the skill's SKILL.md instructs the agent to
show the user what the job will do and ask before either. One approval covers
one invocation; it does not carry to a later run.

Preflight exists so a job fails here, in a second, with a readable message,
rather than twenty minutes into a compute-node allocation.
"""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path


class ExecutionRefused(RuntimeError):
    """Raised instead of running anything when preflight finds a problem."""


def preflight(
    scripts: list[Path | str],
    *,
    submit: bool = False,
    require: tuple[str, ...] = (),
    write_targets: list[Path | str] | None = None,
) -> list[str]:
    """Return a list of human-readable problems; empty means good to go."""
    problems: list[str] = []

    for script in scripts:
        path = Path(script)
        if not path.is_file():
            problems.append(f"script not found: {path}")
        elif not os.access(path, os.R_OK):
            problems.append(f"script is not readable: {path}")

    for tool in require:
        if shutil.which(tool) is None:
            problems.append(
                f"{tool} is not on PATH; install it or choose another tool")

    if submit and shutil.which("sbatch") is None:
        problems.append(
            "sbatch is not on PATH; this host is not a SLURM submit node")

    for target in write_targets or []:
        path = Path(target)
        probe = path if path.is_dir() else path.parent
        if not os.access(probe, os.W_OK):
            problems.append(f"target directory is not writable: {probe}")

    return problems


def _refuse_unless_clean(problems: list[str], action: str) -> None:
    if problems:
        raise ExecutionRefused(
            f"refusing to {action}:\n  - " + "\n  - ".join(problems))


def run_scripts(
    scripts: list[Path | str],
    *,
    require: tuple[str, ...] = (),
    write_targets: list[Path | str] | None = None,
    cwd: Path | str | None = None,
) -> list[subprocess.CompletedProcess]:
    """Execute each script locally, in order, stopping at the first failure.

    Stopping matters: the sra-fetch pair is prefetch-then-dump, and dumping over
    an empty prefetch directory produces a confusing success rather than an
    honest failure.
    """
    _refuse_unless_clean(
        preflight(scripts, require=require, write_targets=write_targets), "run")

    results: list[subprocess.CompletedProcess] = []
    for script in scripts:
        proc = subprocess.run(["bash", str(script)], cwd=str(cwd) if cwd else None)
        results.append(proc)
        if proc.returncode != 0:
            break
    return results


def submit_scripts(
    scripts: list[Path | str],
    *,
    write_targets: list[Path | str] | None = None,
    cwd: Path | str | None = None,
) -> list[str]:
    """sbatch each script, chaining them with `--dependency=afterok:<prev>`.

    Returns the job ids in submission order.
    """
    _refuse_unless_clean(
        preflight(scripts, submit=True, write_targets=write_targets), "submit")

    job_ids: list[str] = []
    for script in scripts:
        cmd = ["sbatch", "--parsable"]
        if job_ids:
            cmd.append(f"--dependency=afterok:{job_ids[-1]}")
        cmd.append(str(script))
        proc = subprocess.run(cmd, capture_output=True, text=True,
                              cwd=str(cwd) if cwd else None)
        if proc.returncode != 0:
            raise ExecutionRefused(
                f"sbatch rejected {Path(script).name}: "
                f"{(proc.stderr or proc.stdout or '').strip()}")
        job_ids.append((proc.stdout or "").strip().split(";")[0])
    return job_ids
