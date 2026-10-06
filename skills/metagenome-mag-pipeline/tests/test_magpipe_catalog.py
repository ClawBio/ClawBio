"""The skills/catalog.json entry must mirror SKILL.md and the CLI registration.

The catalog is the machine-readable index an agent routes on, so it must not
drift from the skill's own description/version/tags/trigger_keywords, and the
CLI alias the catalog advertises must be the alias `clawbio.py run` accepts.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import yaml

SKILL_DIR = Path(__file__).resolve().parent.parent
REPO_ROOT = SKILL_DIR.parent.parent
CATALOG = SKILL_DIR.parent / "catalog.json"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from clawbio.cli import SKILLS

ALIAS = "mag-pipeline"
FOLDER = "metagenome-mag-pipeline"


def _frontmatter() -> dict:
    text = (SKILL_DIR / "SKILL.md").read_text(encoding="utf-8")
    assert text.startswith("---"), "SKILL.md must open with a YAML frontmatter block"
    _, block, _ = text.split("---", 2)
    return yaml.safe_load(block)


def _entry() -> dict:
    catalog = json.loads(CATALOG.read_text(encoding="utf-8"))
    items = catalog if isinstance(catalog, list) else catalog.get("skills", [])
    return next(entry for entry in items if entry.get("name") == FOLDER)


def test_catalog_entry_mirrors_skill_frontmatter() -> None:
    frontmatter = _frontmatter()
    metadata = frontmatter.get("metadata", {})
    entry = _entry()
    assert entry["description"] == frontmatter["description"]
    assert entry["version"] == metadata["version"]
    assert entry["tags"] == metadata["tags"]
    assert entry["trigger_keywords"] == metadata["openclaw"]["trigger_keywords"]
    assert entry["license"] == frontmatter["license"]


def test_catalog_entry_advertises_the_cli_alias() -> None:
    entry = _entry()
    assert entry["cli_alias"] == ALIAS
    assert entry["has_script"] is True
    assert entry["has_tests"] is True
    assert entry["demo_command"] == f"python clawbio.py run {ALIAS} --demo"


def test_alias_is_registered_with_the_same_script() -> None:
    info = SKILLS[ALIAS]
    assert Path(info["script"]) == SKILL_DIR / "metagenome_mag_pipeline.py"
    assert info["demo_args"] == ["--demo"]
    assert info["accepts_genotypes"] is False
    assert Path(info["script"]).is_file()


def test_chaining_partners_are_registered_skills() -> None:
    entry = _entry()
    assert entry["chaining_partners"], "the skill chains with claw-metagenomics and others"
    for partner in entry["chaining_partners"]:
        assert (SKILL_DIR.parent / partner / "SKILL.md").is_file(), (
            f"{partner} is advertised as a chaining partner but has no SKILL.md"
        )


def test_catalog_counts_include_this_skill() -> None:
    catalog = json.loads(CATALOG.read_text(encoding="utf-8"))
    assert catalog["skill_count"] == len(catalog["skills"])
    assert any(entry["name"] == FOLDER for entry in catalog["skills"])
