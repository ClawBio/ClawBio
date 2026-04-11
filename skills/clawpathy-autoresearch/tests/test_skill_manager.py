"""Tests for skill manager — read, modify, create SKILL.md files."""
from __future__ import annotations

from pathlib import Path

import pytest

from skills.clawpathy_autoresearch.skill_manager import (
    SkillManager,
    SkillSnapshot,
)


@pytest.fixture
def skills_dir(tmp_path: Path) -> Path:
    """Create a fake skills directory with one skill."""
    skill_dir = tmp_path / "gwas-lookup"
    skill_dir.mkdir()
    skill_md = skill_dir / "SKILL.md"
    skill_md.write_text(
        "---\nname: gwas-lookup\nversion: 0.1.0\n---\n\n"
        "# GWAS Lookup\n\n## Workflow\n\n"
        "1. Resolve rsID\n2. Query databases\n3. Generate report\n"
    )
    return tmp_path


def test_list_skills(skills_dir: Path):
    mgr = SkillManager(skills_dir)
    skills = mgr.list_skills()
    assert "gwas-lookup" in skills


def test_read_skill(skills_dir: Path):
    mgr = SkillManager(skills_dir)
    content = mgr.read_skill("gwas-lookup")
    assert "# GWAS Lookup" in content
    assert "## Workflow" in content


def test_snapshot_and_restore(skills_dir: Path):
    mgr = SkillManager(skills_dir)
    snapshot = mgr.snapshot()
    assert isinstance(snapshot, SkillSnapshot)
    assert "gwas-lookup" in snapshot.skills

    # Modify skill
    mgr.write_skill("gwas-lookup", "# Modified content")
    assert "Modified" in mgr.read_skill("gwas-lookup")

    # Restore
    mgr.restore(snapshot)
    assert "# GWAS Lookup" in mgr.read_skill("gwas-lookup")


def test_create_new_skill(skills_dir: Path):
    mgr = SkillManager(skills_dir)
    mgr.create_skill("new-skill", "---\nname: new-skill\n---\n\n# New Skill\n")
    assert "new-skill" in mgr.list_skills()
    assert "# New Skill" in mgr.read_skill("new-skill")


def test_create_skill_overwrites_false_by_default(skills_dir: Path):
    mgr = SkillManager(skills_dir)
    with pytest.raises(FileExistsError):
        mgr.create_skill("gwas-lookup", "# Overwrite attempt")


def test_diff_from_snapshot(skills_dir: Path):
    mgr = SkillManager(skills_dir)
    snapshot = mgr.snapshot()
    mgr.write_skill("gwas-lookup", "# Modified GWAS Lookup\n\nNew workflow.")
    diff = mgr.diff_from_snapshot(snapshot)
    assert len(diff) > 0
    assert "gwas-lookup" in diff
