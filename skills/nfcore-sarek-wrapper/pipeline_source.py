# ruff: noqa: E402
from __future__ import annotations

from pathlib import Path
import sys

from clawbio.common.nfcore import source as _source

_SKILL_DIR = Path(__file__).resolve().parent
if str(_SKILL_DIR) in sys.path:
    sys.path.remove(str(_SKILL_DIR))
sys.path.insert(0, str(_SKILL_DIR))


sys.modules.pop("_isolated_imports", None)
from _isolated_imports import purge_foreign_bare_modules

purge_foreign_bare_modules("errors", "schemas")

from errors import ErrorCode, SkillError
from schemas import DEFAULT_LOCAL_PIPELINE_DIR, DEFAULT_REMOTE_PIPELINE, PIPELINE_REQUIRED_FILES


_GIT_TIMEOUT = 10


def _path_has_whitespace(path: Path) -> bool:
    return any(char.isspace() for char in path.as_posix())


def _git_stdout(path: Path, args: list[str]) -> str:
    return _source.git_stdout(path, args, timeout_seconds=_GIT_TIMEOUT)


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
        missing = _source.missing_pipeline_files(local_dir, PIPELINE_REQUIRED_FILES)
        if missing:
            raise SkillError(
                stage="preflight",
                error_code=ErrorCode.PIPELINE_SOURCE_INVALID,
                message="Local sarek checkout is missing required files.",
                fix="Ensure the sibling `sarek` repository is present and complete.",
                details={"path": str(local_dir), "missing_files": missing},
            )

        from clawbio.common.nfcore.source import read_manifest_version

        branch = _git_stdout(local_dir, ["branch", "--show-current"])
        commit = _git_stdout(local_dir, ["rev-parse", "HEAD"])
        status = _git_stdout(local_dir, ["status", "--porcelain"])
        return {
            "source_kind": "local_checkout",
            "source_ref": str(local_dir),
            "resolved_version": commit or requested_version,
            "manifest_version": read_manifest_version(local_dir / "nextflow.config"),
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
