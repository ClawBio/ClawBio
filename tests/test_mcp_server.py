"""Tests for the ClawBio MCP server core.

The core is deliberately independent of the `mcp` SDK so that it can be tested
without the optional dependency installed. The SDK binding in
`clawbio.mcp_server.serve` is a thin wrapper over these functions.
"""

import json

import pytest

from clawbio import mcp_server


class TestCatalog:
    def test_catalog_skills_loads_every_skill(self):
        skills = mcp_server.catalog_skills()
        assert isinstance(skills, list)
        # The catalog is the authoritative count; it should agree with itself.
        raw = json.loads((mcp_server.CATALOG_PATH).read_text())
        assert len(skills) == raw["skill_count"]

    def test_catalog_entries_have_the_fields_the_tools_expose(self):
        for entry in mcp_server.catalog_skills():
            assert entry["name"]
            assert "description" in entry
            assert "status" in entry


class TestListSkills:
    def test_empty_query_returns_all_skills(self):
        assert len(mcp_server.list_skills()) == len(mcp_server.catalog_skills())

    def test_query_filters_by_name(self):
        results = mcp_server.list_skills("pharmgx")
        assert results
        assert any("pharmgx" in r["name"] for r in results)

    def test_query_is_case_insensitive(self):
        assert mcp_server.list_skills("PHARMGX") == mcp_server.list_skills("pharmgx")

    def test_query_matches_description_and_tags(self):
        results = mcp_server.list_skills("pharmacogenomic")
        assert results, "expected a description/tag match for 'pharmacogenomic'"

    def test_unmatched_query_returns_empty_list(self):
        assert mcp_server.list_skills("zzz-no-such-skill-zzz") == []

    def test_entries_declare_whether_they_can_actually_be_run(self):
        """Only CLI-registered skills are executable; agents must not guess."""
        from clawbio.cli import SKILLS

        for entry in mcp_server.list_skills():
            expected = entry["name"] in SKILLS or (entry.get("cli_alias") or "") in SKILLS
            assert entry["runnable"] is expected, entry["name"]

    def test_some_skills_are_runnable_and_some_are_not(self):
        flags = {e["runnable"] for e in mcp_server.list_skills()}
        assert flags == {True, False}, "expected a mix of runnable and spec-only skills"

    def test_results_are_compact(self):
        """Listing 94 skills must not blow up an agent's context."""
        for entry in mcp_server.list_skills():
            assert set(entry) <= {
                "name",
                "cli_alias",
                "description",
                "status",
                "maturity_tier",
                "has_demo",
                "demo_command",
                "tags",
                "runnable",
            }


class TestRunnableIdentifier:
    """Regression: the CLI registry is keyed by cli_alias, not catalog name."""

    def test_run_dispatches_using_the_cli_registry_key(self, monkeypatch):
        from clawbio.cli import SKILLS

        seen = {}
        monkeypatch.setattr(
            mcp_server,
            "_cli_run_skill",
            lambda skill_name, **kw: seen.setdefault("name", skill_name)
            and None or {"skill": skill_name, "success": True, "stdout": "", "files": []},
        )
        mcp_server.run_skill("pharmgx-reporter", demo=True)
        assert seen["name"] in SKILLS, f"{seen['name']!r} is not a runnable CLI skill"

    def test_non_runnable_skill_is_refused_with_a_clear_message(self):
        from clawbio.cli import SKILLS

        spec_only = next(
            e for e in mcp_server.list_skills() if not e["runnable"]
        )
        with pytest.raises(mcp_server.SkillNotRunnableError):
            mcp_server.run_skill(spec_only["name"], demo=True)


class TestDescribeSkill:
    def test_returns_spec_text_and_metadata(self):
        info = mcp_server.describe_skill("pharmgx-reporter")
        assert info["name"] == "pharmgx-reporter"
        assert "maturity_tier" in info
        assert info["spec"], "SKILL.md contents should be returned"

    def test_unknown_skill_raises(self):
        with pytest.raises(mcp_server.SkillNotFoundError):
            mcp_server.describe_skill("zzz-no-such-skill-zzz")

    def test_accepts_cli_alias(self):
        info = mcp_server.describe_skill("pharmgx")
        assert info["name"] == "pharmgx-reporter"


class TestRunSkillSafety:
    def test_demo_run_is_allowed_by_default(self, monkeypatch):
        captured = {}

        def fake_run_skill(skill_name, **kwargs):
            captured.update({"skill_name": skill_name, **kwargs})
            return {"skill": skill_name, "success": True, "stdout": "ok", "files": []}

        monkeypatch.setattr(mcp_server, "_cli_run_skill", fake_run_skill)
        result = mcp_server.run_skill("pharmgx", demo=True)
        assert result["success"] is True
        assert captured["demo"] is True

    def test_local_input_paths_are_refused_unless_opted_in(self, monkeypatch):
        monkeypatch.delenv(mcp_server.ALLOW_LOCAL_FILES_ENV, raising=False)
        with pytest.raises(mcp_server.LocalFileAccessDenied):
            mcp_server.run_skill("pharmgx", input_path="/tmp/my_genome.vcf")

    def test_local_input_paths_allowed_when_opted_in(self, monkeypatch):
        monkeypatch.setenv(mcp_server.ALLOW_LOCAL_FILES_ENV, "1")
        monkeypatch.setattr(
            mcp_server,
            "_cli_run_skill",
            lambda skill_name, **kw: {"skill": skill_name, "success": True, "stdout": "", "files": []},
        )
        result = mcp_server.run_skill("pharmgx", input_path="/tmp/my_genome.vcf")
        assert result["success"] is True

    def test_unknown_skill_raises_before_execution(self):
        with pytest.raises(mcp_server.SkillNotFoundError):
            mcp_server.run_skill("zzz-no-such-skill-zzz", demo=True)

    def test_stdout_is_truncated(self, monkeypatch):
        monkeypatch.setattr(
            mcp_server,
            "_cli_run_skill",
            lambda skill_name, **kw: {
                "skill": skill_name,
                "success": True,
                "stdout": "x" * (mcp_server.MAX_OUTPUT_CHARS + 5000),
                "files": [],
            },
        )
        result = mcp_server.run_skill("pharmgx", demo=True)
        assert len(result["stdout"]) <= mcp_server.MAX_OUTPUT_CHARS + 200
        assert "truncated" in result["stdout"]
