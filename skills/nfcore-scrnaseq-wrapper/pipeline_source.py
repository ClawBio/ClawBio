from __future__ import annotations

import sys
from pathlib import Path

from clawbio.common.nfcore import source as _source

_SKILL_DIR = Path(__file__).resolve().parent
if str(_SKILL_DIR) in sys.path:
    sys.path.remove(str(_SKILL_DIR))
sys.path.insert(0, str(_SKILL_DIR))
sys.modules.pop("_isolated_imports", None)
from _isolated_imports import purge_foreign_bare_modules

purge_foreign_bare_modules("errors", "schemas")

from errors import ErrorCode, SkillError
from schemas import (
    DEFAULT_LOCAL_PIPELINE_DIR,
    DEFAULT_REMOTE_PIPELINE,
    PIPELINE_REQUIRED_FILES,
)


_GIT_TIMEOUT = 10


def _path_has_whitespace(path: Path) -> bool:
    return any(" " in part for part in path.parts)


def _git_stdout(path: Path, args: list[str]) -> str:
    return _source.git_stdout(path, args, timeout_seconds=_GIT_TIMEOUT)


def _remote_pipeline_source(requested_version: str) -> dict[str, object]:
    return {
        "source_kind": "remote_repo",
        "source_ref": DEFAULT_REMOTE_PIPELINE,
        "resolved_version": requested_version,
        "branch": "",
        "dirty": False,
    }


def _validate_local_checkout_version(
    local_dir: Path,
    requested_version: str,
) -> tuple[str, str, str] | None:
    tag = _git_stdout(local_dir, ["describe", "--tags", "--exact-match"])
    commit = _git_stdout(local_dir, ["rev-parse", "HEAD"])

    if tag == requested_version and commit:
        return commit, commit, tag
    if commit and commit == requested_version:
        return commit, commit, tag
    if not tag and not commit:
        return None

    raise SkillError(
        stage="preflight",
        error_code=ErrorCode.PIPELINE_SOURCE_INVALID,
        message="Local scrnaseq checkout does not match the requested pipeline version.",
        fix=(
            "Checkout the requested nf-core/scrnaseq tag or commit, "
            "or remove the local checkout to use the remote pipeline."
        ),
        details={
            "path": str(local_dir),
            "requested_version": requested_version,
            "detected_tag": tag,
            "detected_commit": commit,
        },
    )


def resolve_pipeline_source(
    *,
    requested_version: str,
    local_pipeline_dir: Path | None = None,
    allow_dirty: bool = False,
    require_local: bool = False,
) -> dict[str, object]:
    local_dir = (local_pipeline_dir or DEFAULT_LOCAL_PIPELINE_DIR).resolve()
    if local_dir.exists() and _path_has_whitespace(local_dir):
        if require_local:
            raise SkillError(
                stage="preflight",
                error_code=ErrorCode.PIPELINE_SOURCE_INVALID,
                message="Local scrnaseq checkout path contains whitespace and cannot be used safely.",
                fix="Move the nf-core/scrnaseq checkout to a whitespace-free sibling path, or remove --require-local-pipeline.",
                details={"path": str(local_dir), "require_local": True},
            )
        print(
            f"WARNING: Local scrnaseq checkout at '{local_dir}' contains whitespace in its path. "
            "Docker on macOS cannot reliably execute scripts from paths with spaces (errno 35). "
            "Falling back to the remote nf-core/scrnaseq pipeline.",
            file=sys.stderr,
        )
    elif local_dir.exists():
        missing = _source.missing_pipeline_files(local_dir, PIPELINE_REQUIRED_FILES)
        if missing:
            raise SkillError(
                stage="preflight",
                error_code=ErrorCode.PIPELINE_SOURCE_INVALID,
                message="Local scrnaseq checkout is missing required files.",
                fix="Ensure the sibling `scrnaseq` repository is present and complete.",
                details={"path": str(local_dir), "missing_files": missing},
            )

        version_metadata = _validate_local_checkout_version(
            local_dir, requested_version
        )
        if version_metadata is None:
            print(
                f"WARNING: Local scrnaseq checkout at '{local_dir}' has no verifiable git metadata. "
                "Falling back to the remote nf-core/scrnaseq pipeline.",
                file=sys.stderr,
            )
            if require_local:
                raise SkillError(
                    stage="preflight",
                    error_code=ErrorCode.PIPELINE_SOURCE_INVALID,
                    message="A local scrnaseq checkout was required but its git metadata is not verifiable.",
                    fix="Use a git checkout at the requested nf-core/scrnaseq tag/commit, or remove --require-local-pipeline.",
                    details={"path": str(local_dir), "require_local": True},
                )
            return _remote_pipeline_source(requested_version)

        resolved_version, commit, tag = version_metadata
        branch = _git_stdout(local_dir, ["branch", "--show-current"])
        # No --untracked-files=no: untracked files are intentionally included so that
        # locally-added config overrides or data files are reflected in dirty=True.
        # git respects .gitignore, so build artifacts (work/, .nextflow/) stay invisible
        # as long as the pipeline checkout carries a proper .gitignore.
        status = _git_stdout(local_dir, ["status", "--porcelain"])
        if status and not allow_dirty:
            raise SkillError(
                stage="preflight",
                error_code=ErrorCode.PIPELINE_SOURCE_INVALID,
                message="Local scrnaseq checkout has uncommitted or untracked changes.",
                fix=(
                    "Commit, stash, or remove local changes before running this production wrapper. "
                    "Use --allow-dirty-pipeline only for intentional local pipeline development."
                ),
                details={
                    "path": str(local_dir),
                    "requested_version": requested_version,
                    "detected_tag": tag,
                    "detected_commit": commit,
                    "dirty": True,
                    "status": status,
                },
            )
        return {
            "source_kind": "local_checkout",
            "source_ref": str(local_dir),
            "resolved_version": resolved_version,
            "commit": commit,
            "tag": tag,
            "branch": branch,
            "dirty": bool(status),
        }

    if require_local:
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.PIPELINE_SOURCE_INVALID,
            message="A local scrnaseq checkout was required but no checkout was found.",
            fix="Place nf-core/scrnaseq as a sibling checkout next to ClawBio, or remove --require-local-pipeline.",
            details={"path": str(local_dir), "require_local": True},
        )
    return _remote_pipeline_source(requested_version)
