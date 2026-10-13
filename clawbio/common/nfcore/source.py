"""Bounded Git and file-presence probes; callers decide source selection."""

import subprocess
from pathlib import Path


def git_stdout(path: Path, args: list[str], *, timeout_seconds: float = 10) -> str:
    try:
        proc = subprocess.run(
            ["git", *args],
            cwd=str(path),
            capture_output=True,
            text=True,
            check=True,
            timeout=timeout_seconds,
        )
    except (
        subprocess.CalledProcessError,
        FileNotFoundError,
        OSError,
        subprocess.TimeoutExpired,
    ):
        return ""
    return proc.stdout.strip()


def missing_pipeline_files(
    path: Path, required_files: list[str] | tuple[str, ...]
) -> list[str]:
    return [name for name in required_files if not (path / name).exists()]
