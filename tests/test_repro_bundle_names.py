"""Bare-name imports collide in sys.modules, so bundle module names must be unique across skills."""

from pathlib import Path

SKILLS = Path(__file__).resolve().parents[1] / "skills"


def test_repro_bundle_module_names_are_unique_across_skills():
    names = [p.name for p in SKILLS.glob("*/*repro_bundle*.py")]
    assert len(names) == len(set(names))
