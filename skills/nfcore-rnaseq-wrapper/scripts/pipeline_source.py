from __future__ import annotations

from pathlib import Path
import re
import subprocess
import sys

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

# Matches `version = '3.26.0'` inside the manifest block of nextflow.config.
# nextflow.config also carries `nextflowVersion`; the trailing-boundary group
# (`\b`) prevents `nextflowVersion` from matching the bare `version` key.
_MANIFEST_VERSION_RE = re.compile(r"(?<![A-Za-z])version\s*=\s*['\"]([^'\"]+)['\"]")


def _path_has_whitespace(path: Path) -> bool:
    return any(" " in part for part in path.parts)


def _read_manifest_version(local_dir: Path) -> str:
    """Best-effort parse of manifest.version from a local nextflow.config.

    Returns "" when the file is unreadable or the key is absent — callers must
    treat an empty string as "unknown" and not as a version mismatch.
    """
    config_path = local_dir / "nextflow.config"
    try:
        text = config_path.read_text(encoding="utf-8")
    except (OSError, UnicodeError):
        return ""
    match = _MANIFEST_VERSION_RE.search(text)
    return match.group(1).strip() if match else ""


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
    # Diagnostic recorded on the remote fallback when a present local checkout is
    # rejected, so provenance/upstream.json (build_upstream_payload) can explain why
    # the sibling checkout was not used instead of silently going remote.
    local_rejection: dict[str, str] = {}
    if local_dir.exists() and _path_has_whitespace(local_dir):
        reason = (
            "whitespace in checkout path — Docker on macOS cannot reliably execute "
            "scripts from paths with spaces (errno 35)"
        )
        print(
            f"WARNING: Local rnaseq checkout at '{local_dir}' contains whitespace in its path. "
            "Docker on macOS cannot reliably execute scripts from paths with spaces (errno 35). "
            "Falling back to the remote nf-core/rnaseq pipeline.",
            file=sys.stderr,
        )
        local_rejection = {"local_attempted": str(local_dir), "local_rejected_reason": reason}
    elif local_dir.exists():
        missing = [name for name in PIPELINE_REQUIRED_FILES if not (local_dir / name).exists()]
        if missing:
            raise SkillError(
                stage="preflight",
                error_code=ErrorCode.PIPELINE_SOURCE_INVALID,
                message="Local rnaseq checkout is missing required files.",
                fix="Ensure the sibling `rnaseq` repository is present and complete.",
                details={"path": str(local_dir), "missing_files": missing},
            )

        branch = _git_stdout(local_dir, ["branch", "--show-current"])
        commit = _git_stdout(local_dir, ["rev-parse", "HEAD"])
        status = _git_stdout(local_dir, ["status", "--porcelain"])
        return {
            "source_kind": "local_checkout",
            "source_ref": str(local_dir),
            "resolved_version": commit or requested_version,
            # manifest.version of the checkout itself — used to enforce the
            # pinned-version contract on local checkouts (audit F-01). "" = unknown.
            "manifest_version": _read_manifest_version(local_dir),
            "branch": branch,
            "dirty": bool(status),
        }

    return {
        "source_kind": "remote_repo",
        "source_ref": DEFAULT_REMOTE_PIPELINE,
        "resolved_version": requested_version,
        "manifest_version": requested_version,
        "branch": "",
        "dirty": False,
        **local_rejection,
    }
