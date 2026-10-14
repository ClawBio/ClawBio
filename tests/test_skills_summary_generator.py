"""Inventory tests use a synthetic tree, not the changing repository registry."""

from __future__ import annotations

import importlib.util
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest
from hypothesis import given, strategies as st

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = REPO_ROOT / "scripts" / "generate_skills_summary.py"


def _load_generator():
    spec = importlib.util.spec_from_file_location("generate_skills_summary", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def gen(tmp_path, monkeypatch):
    module = _load_generator()
    skills_dir = tmp_path / "skills"
    skills_dir.mkdir()
    for name, description, domain in [
        ("alpha", "Alpha description.", "quality"),
        ("beta", "Beta description.", "quality"),
        ("gamma", "Gamma description.", None),
    ]:
        folder = skills_dir / name
        folder.mkdir()
        metadata = f"metadata:\n  domain: {domain}\n" if domain else ""
        (folder / "SKILL.md").write_text(
            f"---\nname: {name}\ndescription: {description}\n{metadata}---\n"
        )
        if name != "gamma":
            (folder / f"{name}.py").write_text("# synthetic implementation\n")
    registry = {
        "alpha": {
            "script": skills_dir / "alpha" / "alpha.py",
            "description": "Old registry description.",
        }
    }
    from clawbio import cli
    monkeypatch.setattr(cli, "SKILLS", registry)
    monkeypatch.setattr(cli, "_skill_md_only_directories",
                        lambda: [skills_dir / "beta", skills_dir / "gamma"])
    monkeypatch.setattr(module, "CLAWBIO_DIR", tmp_path)
    monkeypatch.setattr(module, "SKILLS_DIR", skills_dir)
    monkeypatch.setattr(module, "DOMAINS_PATH", tmp_path / "skill_domains.json")
    monkeypatch.setattr(module, "OUTPUT_PATH", tmp_path / "SKILLS_SUMMARY.md")
    module.DOMAINS_PATH.write_text(json.dumps({
        "domains": [
            {"id": "analysis", "title": "Analysis", "members": ["alpha"]},
            {"id": "quality", "title": "Quality", "members": []},
        ],
        "curated_findings": {},
    }))
    return module


@pytest.fixture
def skills(gen):
    return gen.collect_skills()


@pytest.fixture
def taxonomy(gen):
    return gen.load_taxonomy()


def test_unclassified_skill_warns_and_is_rendered(gen, taxonomy, skills, capsys):
    resolved = gen.validate_coverage(taxonomy, skills)
    content = gen.render(resolved, skills)
    assert "Uncategorised" in content and "`gamma`" in content
    assert "warning:" in capsys.readouterr().err
    assert sum(len(d["members"]) for d in resolved["domains"]) == len(skills)


def test_stale_taxonomy_entry_warns_and_is_dropped(gen, taxonomy, skills, capsys):
    taxonomy["domains"][0]["members"].append("retired-skill")
    resolved = gen.validate_coverage(taxonomy, skills)
    assert "retired-skill" not in {m for d in resolved["domains"] for m in d["members"]}
    assert "retired-skill" in capsys.readouterr().err
    assert gen.render(resolved, skills).startswith("# ClawBio")


def test_double_assignment_is_rejected(gen, taxonomy, skills):
    taxonomy["domains"][1]["members"].append("alpha")
    with pytest.raises(SystemExit) as exc:
        gen.validate_coverage(taxonomy, skills)
    assert exc.value.code == 1


def test_metadata_domain_is_used_without_a_json_entry(gen, taxonomy, skills):
    assert skills["beta"]["domain"] == "quality"
    resolved = gen.validate_coverage(taxonomy, skills)
    assert "beta" in next(d for d in resolved["domains"] if d["id"] == "quality")["members"]


def test_json_membership_can_override_coarse_frontmatter(gen, taxonomy, skills):
    resolved = gen.validate_coverage(taxonomy, skills)
    assert "alpha" in next(d for d in resolved["domains"] if d["id"] == "analysis")["members"]
    assert "alpha" not in next(d for d in resolved["domains"] if d["id"] == "quality")["members"]


def test_new_frontmatter_domain_is_preserved(gen, taxonomy, skills):
    skills["gamma"]["domain"] = "new-domain"
    resolved = gen.validate_coverage(taxonomy, skills)
    assert any(d["id"] == "new-domain" and "gamma" in d["members"] for d in resolved["domains"])


def test_registered_description_comes_from_skill_frontmatter(skills):
    assert skills["alpha"]["description"] == "Alpha description."


def test_tiers_match_synthetic_registry(gen, skills):
    assert {n for n, s in skills.items() if s["tier"] == gen.TIER_REGISTERED} == {"alpha"}
    assert {n for n, s in skills.items() if s["tier"] == gen.TIER_AGENT} == {"beta", "gamma"}


def test_alias_folder_collision_is_rejected(gen, monkeypatch):
    from clawbio import cli
    monkeypatch.setattr(cli, "_skill_md_only_directories", lambda: [gen.SKILLS_DIR / "alpha"])
    with pytest.raises(SystemExit):
        gen.assert_flat_namespace({})


def test_spec_only_classification_matches_python_indicator(skills):
    assert skills["beta"]["has_script"] and not skills["beta"]["spec_only"]
    assert not skills["gamma"]["has_script"] and skills["gamma"]["spec_only"]


def test_every_skill_appears_exactly_once_in_domain_entries(gen, taxonomy, skills):
    resolved = gen.validate_coverage(taxonomy, skills)
    members = [n for d in resolved["domains"] for n in d["members"]]
    assert sorted(members) == sorted(skills)
    content = gen.render(resolved, skills)
    for name in skills:
        assert f"`{name}`" in content


def test_missing_scripts_are_flagged(gen, taxonomy, skills):
    skills["alpha"]["script_exists"] = False
    assert "script missing: `alpha.py`" in gen.render(taxonomy, skills)


def test_render_has_a_golden_output(gen, taxonomy, skills):
    content = gen.render(taxonomy, skills)
    expected = (REPO_ROOT / "tests/fixtures/skills_summary_expected.txt").read_text().rstrip("\n") + "\n"
    assert content == expected


def test_catalog_regeneration_cannot_change_rendered_inventory(gen, taxonomy, skills):
    before = gen.render(taxonomy, skills)
    catalog = gen.SKILLS_DIR / "catalog.json"
    catalog.write_text('{"skill_count": 123, "skills": []}')
    after = gen.render(taxonomy, skills)
    assert before == after
    assert "indexes 123" not in after


@given(st.permutations(("alpha", "beta", "gamma")))
def test_input_order_does_not_change_output(order):
    gen = _load_generator()
    skills = {
        n: {"name": n, "tier": gen.TIER_AGENT, "description": "", "domain": "quality",
            "folder": n, "has_script": False, "spec_only": True,
            "script_exists": None, "script_name": None}
        for n in order
    }
    taxonomy = {"domains": [{"id": "quality", "title": "Quality", "members": []}]}
    ordered = {n: skills[n] for n in sorted(skills)}
    assert gen.render(taxonomy, skills) == gen.render(taxonomy, ordered)


def test_cli_check_on_synthetic_tree(gen, monkeypatch):
    content = gen.render(gen.load_taxonomy(), gen.collect_skills())
    gen.OUTPUT_PATH.write_text(content)
    monkeypatch.setattr(sys, "argv", ["generate_skills_summary", "--check"])
    gen.main()


def test_cli_check_stale_synthetic_file_exits_two(gen, monkeypatch):
    gen.OUTPUT_PATH.write_text("stale\n")
    monkeypatch.setattr(sys, "argv", ["generate_skills_summary", "--check"])
    with pytest.raises(SystemExit) as exc:
        gen.main()
    assert exc.value.code == 2


def test_cli_stdout_does_not_modify_synthetic_file(gen, monkeypatch, capsys):
    gen.OUTPUT_PATH.write_text("preserve me\n")
    monkeypatch.setattr(sys, "argv", ["generate_skills_summary", "--output", "-"])
    gen.main()
    assert capsys.readouterr().out.startswith("# ClawBio Skill Library")
    assert gen.OUTPUT_PATH.read_text() == "preserve me\n"


@pytest.mark.parametrize("docstring", [None, "", "one line"])
def test_cli_description_handles_missing_docstrings(gen, monkeypatch, docstring):
    monkeypatch.setattr(gen, "__doc__", docstring)
    monkeypatch.setattr(sys, "argv", ["generate_skills_summary", "--help"])
    with pytest.raises(SystemExit) as exc:
        gen.main()
    assert exc.value.code == 0


def test_cli_runs_with_optimized_docstrings():
    result = subprocess.run([sys.executable, "-OO", str(SCRIPT_PATH), "--help"],
                            capture_output=True, text=True)
    assert result.returncode == 0, result.stderr


@pytest.mark.integration
@pytest.mark.skipif(os.environ.get("CLAWBIO_CHECK_SKILLS_SUMMARY") != "1",
                    reason="Live inventory freshness is opt-in, not an unrelated-PR gate")
def test_check_mode_reports_committed_file_is_current():
    result = subprocess.run([sys.executable, str(SCRIPT_PATH), "--check"],
                            cwd=REPO_ROOT, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
