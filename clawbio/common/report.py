"""Common report generation helpers for ClawBio skills."""

from __future__ import annotations

import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from clawbio.common.audit import write as _audit_write
from clawbio.common.checksums import sha256_file

_DEFAULT_AUDIT_LOG = Path.home() / ".clawbio" / "audit.jsonl"

DISCLAIMER = (
    "ClawBio is a research and educational tool. It is not a medical device "
    "and does not provide clinical diagnoses. Consult a healthcare "
    "professional before making any medical decisions."
)


def write_audit_log(
    skill: str,
    version: str,
    input_checksum: str,
    output_dir: str,
    log_path: Path | str = _DEFAULT_AUDIT_LOG,
) -> None:
    """Append a skill_run record to the central audit log."""
    _audit_write(
        "skill_run",
        skill=skill,
        version=version,
        input_checksum=input_checksum,
        output_dir=output_dir,
        log_path=log_path,
    )


def generate_report_header(
    title: str,
    skill_name: str,
    skill_version: str = "",
    input_files: list[Path] | None = None,
    extra_metadata: dict[str, str] | None = None,
) -> str:
    """Generate the standard markdown report header.

    Args:
        title: Report title.
        skill_name: Name of the skill that generated the report.
        input_files: List of input file paths (checksums computed automatically).
        extra_metadata: Additional key-value pairs to include in the header.
    """
    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")

    checksums = []
    if input_files:
        for f in input_files:
            f = Path(f)
            if f.exists():
                checksums.append(f"- `{f.name}`: `{sha256_file(f)}`")
            else:
                checksums.append(f"- `{f.name}`: (file not found)")

    lines = [
        f"# {title}",
        "",
        f"**Date**: {now}",
        f"**Skill**: {skill_name}",
        f"**Version**: {skill_version}",
    ]
    if extra_metadata:
        for key, val in extra_metadata.items():
            lines.append(f"**{key}**: {val}")
    if checksums:
        lines.append("**Input files**:")
        lines.extend(checksums)
    lines.extend(["", "---", ""])

    return "\n".join(lines)


def generate_report_footer() -> str:
    """Generate the standard markdown report footer with disclaimer."""
    return f"""
---

## Disclaimer

*{DISCLAIMER}*
"""


def version_from_skill_md(text: str) -> str:
    """Return `metadata.version` from a SKILL.md's YAML frontmatter.

    Raises ValueError when there is no frontmatter or no version: a missing
    version is an error, never a silent `0.1.0`.
    """
    import yaml  # lazy, like scripts/generate_catalog.py

    match = re.match(r"---\n(.*?)\n---\n", text, re.S)
    front = yaml.safe_load(match.group(1)) if match else None
    version = ((front or {}).get("metadata") or {}).get("version")
    if not version:
        raise ValueError("SKILL.md frontmatter has no metadata.version")
    return str(version)


def skill_version(path: str | Path) -> str:
    """Version of the skill that owns `path`, read from the nearest SKILL.md above it.

    Call it as `skill_version(__file__)`; SKILL.md is the only place a skill's
    version is written.
    """
    for folder in Path(path).resolve().parents:
        skill_md = folder / "SKILL.md"
        if skill_md.is_file():
            return version_from_skill_md(skill_md.read_text(encoding="utf-8"))
    raise FileNotFoundError(f"no SKILL.md above {path}")


def write_result_json(
    output_dir: str | Path,
    skill: str,
    version: str,
    summary: dict[str, Any],
    data: dict[str, Any],
    input_checksum: str = "",
    datasets: dict[str, str] | None = None,
    status: str | None = None,
    ok: bool | None = None,
) -> Path:
    """Write the standardized result.json envelope alongside report.md.

    Args:
        output_dir: Directory to write result.json to.
        skill: Skill name (e.g., "pharmgx").
        version: Skill version string.
        summary: High-level summary dict (skill-specific).
        data: Full result data dict (skill-specific).
        input_checksum: SHA-256 hex digest of the input file.
        status: Optional run status (e.g. "ok"/"error"). When provided it is
            written as a top-level key. Opt-in and backward-compatible: callers
            that omit it get the historical envelope unchanged.
        ok: Optional boolean success discriminator, written as a top-level key
            when provided. Together with ``status`` this is the minimal shared
            contract across the nf-core pipeline wrappers.

    Returns:
        Path to the written result.json file.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    envelope: dict[str, Any] = {}
    if status is not None:
        envelope["status"] = status
    if ok is not None:
        envelope["ok"] = ok
    envelope.update(
        {
            "skill": skill,
            "version": version,
            "completed_at": datetime.now(timezone.utc).isoformat(),
            "input_checksum": f"sha256:{input_checksum}" if input_checksum else "",
            "datasets": datasets or {},
            "summary": summary,
            "data": data,
        }
    )

    result_path = output_dir / "result.json"
    result_path.write_text(json.dumps(envelope, indent=2, default=str))
    return result_path
