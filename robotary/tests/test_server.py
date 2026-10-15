"""Tests for Robotary server."""
import sys
from pathlib import Path

import pytest

pytest.importorskip("fastapi", reason="needs fastapi (uv sync --extra robotary)")

# Add robotary to path
ROBOTARY_DIR = Path(__file__).resolve().parent.parent
PROJECT_ROOT = ROBOTARY_DIR.parent
sys.path.insert(0, str(ROBOTARY_DIR))


def test_skill_catalog_loads():
    """Skill catalog should discover at least 7 skills from SKILL.md files."""
    from server import build_skill_catalog

    catalog = build_skill_catalog()
    assert len(catalog) >= 7
    # Core skills must be present
    for skill in ["pharmgx-reporter", "nutrigx", "genome-compare",
                  "gwas-prs", "clinpgx", "gwas-lookup", "profile-report"]:
        assert skill in catalog, f"{skill} not in catalog"
        assert len(catalog[skill]) > 0  # has a description


def test_skill_catalog_descriptions_are_strings():
    from server import build_skill_catalog

    catalog = build_skill_catalog()
    for name, desc in catalog.items():
        assert isinstance(desc, str)
        assert len(desc) <= 200
