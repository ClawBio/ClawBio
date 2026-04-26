"""Tests for the public importable clawbio package API."""

from __future__ import annotations

from pathlib import Path


FIXTURES = Path(__file__).resolve().parent / "fixtures"


def test_package_exports_documented_symbols():
    from clawbio import __all__, list_skills, run_skill, upload_profile

    assert "run_skill" in __all__
    assert "list_skills" in __all__
    assert "upload_profile" in __all__
    assert callable(run_skill)
    assert callable(list_skills)
    assert callable(upload_profile)


def test_package_exports_same_runner_implementation():
    from clawbio import list_skills, run_skill, upload_profile
    from clawbio import runner as package_runner

    root_runner = package_runner._load_root_runner()

    assert list_skills is root_runner.list_skills
    assert run_skill is root_runner.run_skill
    assert upload_profile is root_runner.upload_profile


def test_list_skills_returns_registry_with_known_skill():
    from clawbio import list_skills

    skills = list_skills()

    assert isinstance(skills, dict)
    assert "pharmgx" in skills


def test_run_skill_unknown_skill_returns_structured_failure():
    from clawbio import run_skill

    result = run_skill("definitely-not-a-real-skill")

    assert result["success"] is False
    assert result["skill"] == "definitely-not-a-real-skill"
    assert "Unknown skill" in result["stderr"]


def test_upload_profile_works_with_fixture(tmp_path, monkeypatch):
    from clawbio import upload_profile
    from clawbio import runner as package_runner

    root_runner = package_runner._load_root_runner()
    monkeypatch.setattr(root_runner, "PROFILES_DIR", tmp_path)

    result = upload_profile(str(FIXTURES / "mock_23andme.txt"), patient_id="pkg-api")

    assert result["success"] is True
    assert result["patient_id"] == "pkg-api"
    assert Path(result["profile_path"]).exists()
    assert Path(result["profile_path"]).parent == tmp_path
