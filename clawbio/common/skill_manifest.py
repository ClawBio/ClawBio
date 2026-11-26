"""What one skill folder declares in its SKILL.md frontmatter.

The only reader of SKILL.md frontmatter in ClawBio, apart from hatch_build.py's
licence reader (the build environment has no PyYAML). It reads; it never
validates: `agentskills validate` and the conformance tests own that.

Every field except the four top-level keys the agentskills spec allows (name,
description, license, metadata) is read from `metadata:` only.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

_FRONTMATTER = re.compile(r"^---\s*\n(.*?)\n---", re.S)


class SkillManifestError(ValueError):
    """A SKILL.md is missing, has no frontmatter, or lacks a value a caller asked for."""


@dataclass(frozen=True)
class SkillManifest:
    path: Path
    name: str = ""
    description: str = ""
    license: str = ""
    model_license: str = ""
    data_license: str = ""
    tags: tuple[str, ...] = ()
    dependencies: Any = field(default_factory=list)
    openclaw: dict = field(default_factory=dict)
    _version: str | None = None

    @property
    def version(self) -> str:
        """`metadata.version`. Raises when absent: never a silent default."""
        if not self._version:
            raise SkillManifestError(f"{self.path} frontmatter has no metadata.version")
        return self._version


def _find_skill_md(path: Path) -> Path:
    path = Path(path).resolve()
    for folder in (path, *path.parents):
        if (folder / "SKILL.md").is_file():
            return folder / "SKILL.md"
    raise SkillManifestError(f"no SKILL.md at or above {path}")


def read_skill_manifest(path: str | Path) -> SkillManifest:
    """Read the SKILL.md for a skill folder, or for any file inside one."""
    import yaml  # lazy; no regex fallback (ReDoS on crafted frontmatter)

    skill_md = _find_skill_md(Path(path))
    try:
        match = _FRONTMATTER.match(skill_md.read_text(encoding="utf-8"))
        raw = yaml.safe_load(match.group(1)) if match else None
    except (yaml.YAMLError, UnicodeDecodeError) as exc:
        raise SkillManifestError(f"{skill_md} frontmatter is not readable YAML: {exc}") from exc
    if not isinstance(raw, dict):
        raise SkillManifestError(f"{skill_md} has no YAML frontmatter mapping")

    meta = raw.get("metadata") if isinstance(raw.get("metadata"), dict) else {}
    openclaw = meta.get("openclaw") if isinstance(meta.get("openclaw"), dict) else {}
    version = meta.get("version")
    return SkillManifest(
        path=skill_md,
        name=str(raw.get("name") or ""),
        description=str(raw.get("description") or ""),
        license=str(raw.get("license") or ""),
        model_license=str(meta.get("model_license") or ""),
        data_license=str(meta.get("data_license") or ""),
        tags=tuple(str(t) for t in meta.get("tags") or ()),
        dependencies=meta.get("dependencies") or [],
        openclaw=openclaw,
        _version=str(version) if version else None,
    )
