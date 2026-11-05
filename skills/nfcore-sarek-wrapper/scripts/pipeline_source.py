# ruff: noqa: E402
# Ported from ClawBio 1.x nfcore-sarek-wrapper/pipeline_source.py. The 1.x
# `_isolated_imports.purge_foreign_bare_modules` dance guarded against a
# sibling wrapper skill (nfcore-rnaseq-wrapper, nfcore-scrnaseq-wrapper)
# shadowing this file's bare `errors`/`schemas` names in a shared pytest
# session. ClawBio 2.0 runs each step as its own `python3` subprocess (see
# pipeline.yaml command_templates), so that cross-skill shadowing risk does
# not apply here — dropped in favor of a plain sys.path insert.
from __future__ import annotations

from pathlib import Path
import subprocess
import sys

_SCRIPT_DIR = Path(__file__).resolve().parent
if str(_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_DIR))

from errors import ErrorCode, SkillError
from schemas import DEFAULT_LOCAL_PIPELINE_DIR, DEFAULT_REMOTE_PIPELINE, PIPELINE_REQUIRED_FILES


_GIT_TIMEOUT = 10


def _path_has_whitespace(path: Path) -> bool:
    return any(" " in part for part in path.parts)


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
    if local_dir.exists() and _path_has_whitespace(local_dir):
        print(
            f"WARNING: Local sarek checkout at '{local_dir}' contains whitespace in its path. "
            "Docker on macOS cannot reliably execute scripts from paths with spaces (errno 35). "
            "Falling back to the remote nf-core/sarek pipeline.",
            file=sys.stderr,
        )
    elif local_dir.exists():
        missing = [name for name in PIPELINE_REQUIRED_FILES if not (local_dir / name).exists()]
        if missing:
            raise SkillError(
                stage="preflight",
                error_code=ErrorCode.PIPELINE_SOURCE_INVALID,
                message="Local sarek checkout is missing required files.",
                fix="Ensure the sibling `sarek` repository is present and complete.",
                details={"path": str(local_dir), "missing_files": missing},
            )

        branch = _git_stdout(local_dir, ["branch", "--show-current"])
        commit = _git_stdout(local_dir, ["rev-parse", "HEAD"])
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
