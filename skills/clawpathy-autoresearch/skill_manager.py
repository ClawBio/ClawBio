"""Skill manager for clawpathy-autoresearch.

Reads, modifies, creates, snapshots, and restores SKILL.md files
in a skills directory. All changes are reversible via snapshot/restore.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class SkillSnapshot:
    """Frozen state of all SKILL.md files for rollback."""

    skills: dict[str, str]  # skill_name -> SKILL.md content


class SkillManager:
    """Manages SKILL.md files in a skills directory."""

    def __init__(self, skills_dir: Path) -> None:
        self.skills_dir = Path(skills_dir)

    def list_skills(self) -> list[str]:
        """Return names of all skills (directories containing SKILL.md)."""
        skills = []
        if not self.skills_dir.exists():
            return skills
        for d in sorted(self.skills_dir.iterdir()):
            if d.is_dir() and (d / "SKILL.md").exists():
                skills.append(d.name)
        return skills

    def read_skill(self, name: str) -> str:
        """Read the SKILL.md content for a skill."""
        path = self.skills_dir / name / "SKILL.md"
        if not path.exists():
            raise FileNotFoundError(f"Skill not found: {name}")
        return path.read_text()

    def write_skill(self, name: str, content: str) -> None:
        """Overwrite a skill's SKILL.md with new content."""
        path = self.skills_dir / name / "SKILL.md"
        if not path.parent.exists():
            raise FileNotFoundError(f"Skill directory not found: {name}")
        path.write_text(content)

    def create_skill(self, name: str, content: str) -> None:
        """Create a new skill directory with a SKILL.md.

        Raises FileExistsError if the skill already exists.
        """
        skill_dir = self.skills_dir / name
        if skill_dir.exists():
            raise FileExistsError(f"Skill already exists: {name}")
        skill_dir.mkdir(parents=True)
        (skill_dir / "SKILL.md").write_text(content)

    def snapshot(self) -> SkillSnapshot:
        """Capture the current state of all SKILL.md files."""
        skills = {}
        for name in self.list_skills():
            skills[name] = self.read_skill(name)
        return SkillSnapshot(skills=skills)

    def restore(self, snapshot: SkillSnapshot) -> None:
        """Restore all SKILL.md files to a previous snapshot state.

        Skills created after the snapshot are removed.
        """
        current_skills = set(self.list_skills())
        snapshot_skills = set(snapshot.skills.keys())

        # Restore existing skills
        for name, content in snapshot.skills.items():
            path = self.skills_dir / name / "SKILL.md"
            if path.parent.exists():
                path.write_text(content)

        # Remove skills created after snapshot
        for name in current_skills - snapshot_skills:
            skill_dir = self.skills_dir / name
            if skill_dir.exists():
                for f in skill_dir.iterdir():
                    f.unlink()
                skill_dir.rmdir()

    def diff_from_snapshot(self, snapshot: SkillSnapshot) -> dict[str, str]:
        """Return a dict of skill_name -> 'modified' | 'created' | 'deleted'."""
        current_skills = set(self.list_skills())
        snapshot_skills = set(snapshot.skills.keys())
        diff = {}

        for name in current_skills & snapshot_skills:
            if self.read_skill(name) != snapshot.skills[name]:
                diff[name] = "modified"

        for name in current_skills - snapshot_skills:
            diff[name] = "created"

        for name in snapshot_skills - current_skills:
            diff[name] = "deleted"

        return diff
