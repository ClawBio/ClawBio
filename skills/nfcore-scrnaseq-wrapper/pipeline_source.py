from __future__ import annotations

import subprocess
import sys
from pathlib import Path

_SKILL_DIR = Path(__file__).resolve().parent
if str(_SKILL_DIR) not in sys.path:
    sys.path.insert(0, str(_SKILL_DIR))

from errors import ErrorCode, SkillError
from schemas import DEFAULT_LOCAL_PIPELINE_DIR, DEFAULT_REMOTE_PIPELINE, PIPELINE_REQUIRED_FILES


_GIT_TIMEOUT = 10


def _git_stdout(path: Path, args: list[str]) -> str:
    try:
        proc = subprocess.run(
            ["git", *args],
            cwd=str(path),
            capture_output=True,
            text=True,
            check=True,
            timeout=_GIT_TIMEOUT,
        )
    except (subprocess.CalledProcessError, FileNotFoundError, OSError, subprocess.TimeoutExpired):
        return ""
    return proc.stdout.strip()


def resolve_pipeline_source(
    *,
    requested_version: str,
    local_pipeline_dir: Path | None = None,
) -> dict[str, str | bool]:
    local_dir = (local_pipeline_dir or DEFAULT_LOCAL_PIPELINE_DIR).resolve()
    if local_dir.exists():
        missing = [name for name in PIPELINE_REQUIRED_FILES if not (local_dir / name).exists()]
        if missing:
            raise SkillError(
                stage="preflight",
                error_code=ErrorCode.PIPELINE_SOURCE_INVALID,
                message="Local scrnaseq checkout is missing required files.",
                fix="Ensure the sibling `scrnaseq` repository is present and complete.",
                details={"path": str(local_dir), "missing_files": missing},
            )

        branch = _git_stdout(local_dir, ["branch", "--show-current"])
        commit = _git_stdout(local_dir, ["rev-parse", "HEAD"])
        # No --untracked-files=no: untracked files are intentionally included so that
        # locally-added config overrides or data files are reflected in dirty=True.
        # git respects .gitignore, so build artifacts (work/, .nextflow/) stay invisible
        # as long as the pipeline checkout carries a proper .gitignore.
        status = _git_stdout(local_dir, ["status", "--porcelain"])
        return {
            "source_kind": "local_checkout",
            "source_ref": str(local_dir),
            "resolved_version": commit or requested_version,
            "branch": branch,
            "dirty": bool(status),
        }

    return {
        "source_kind": "remote_repo",
        "source_ref": DEFAULT_REMOTE_PIPELINE,
        "resolved_version": requested_version,
        "branch": "",
        "dirty": False,
    }
