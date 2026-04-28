"""
Tests for the ClawBio Skill Builder.

Run:
    pytest skills/skill-builder/tests/ -v
or via ClawBio runner:
    python -m pytest skills/skill-builder/tests/ -v
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

# Make the skill importable regardless of working directory
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from skill_builder import (
    collect_spec_interactive,
    generate_catalog_entry,
    generate_skill_md,
    generate_skill_py,
    generate_test_py,
    load_spec,
    patch_clawbio_py,
    run_demo,
    scaffold_skill,
    update_catalog_json,
    validate_skill_md,
    validate_spec,
    REQUIRED_SECTIONS,
    DEMO_SPEC_PATH,
    # interactive / draft helpers
    _Jump,
    _ask,
    _ask_list,
    _clear_draft,
    _find_resume_cursor,
    _load_draft,
    _NUM_FIELDS,
    _prompt_field,
    _run_review_loop,
    _save_draft,
    DRAFT_SPEC_PATH,
)
import skill_builder as _sb  # used for monkeypatching module-level attributes


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

MINIMAL_SPEC = {
    "name": "test-skill",
    "description": "A minimal test skill for unit testing",
    "author": "Test Author",
    "version": "0.1.0",
    "domain": "genomics",
    "cli_alias": "testskill",
    "tags": ["test"],
    "trigger_keywords": ["test skill", "run test"],
    "capabilities": ["Do something", "Do something else"],
    "input_formats": [],
    "dependencies": {"required": [], "optional": []},
    "chaining_partners": [],
}


@pytest.fixture
def spec() -> dict:
    """Return a minimal valid skill spec."""
    return dict(MINIMAL_SPEC)


@pytest.fixture
def tmp_output(tmp_path: Path) -> Path:
    """Provide a temporary output directory per test."""
    return tmp_path / "output"


@pytest.fixture
def demo_spec() -> dict:
    """Load the real demo_spec.json bundled with the skill."""
    assert DEMO_SPEC_PATH.exists(), f"demo_spec.json not found at {DEMO_SPEC_PATH}"
    return json.loads(DEMO_SPEC_PATH.read_text(encoding="utf-8"))


# ---------------------------------------------------------------------------
# validate_spec
# ---------------------------------------------------------------------------

class TestValidateSpec:
    def test_valid_minimal_spec(self, spec: dict) -> None:
        assert validate_spec(spec) == []

    def test_missing_name(self, spec: dict) -> None:
        del spec["name"]
        errors = validate_spec(spec)
        assert any("name" in e for e in errors)

    def test_missing_description(self, spec: dict) -> None:
        del spec["description"]
        errors = validate_spec(spec)
        assert any("description" in e for e in errors)

    def test_missing_author(self, spec: dict) -> None:
        del spec["author"]
        errors = validate_spec(spec)
        assert any("author" in e for e in errors)

    def test_invalid_name_uppercase(self, spec: dict) -> None:
        spec["name"] = "My-Skill"
        errors = validate_spec(spec)
        assert any("name" in e for e in errors)

    def test_invalid_name_spaces(self, spec: dict) -> None:
        spec["name"] = "my skill"
        errors = validate_spec(spec)
        assert any("name" in e for e in errors)

    def test_valid_name_with_digits(self, spec: dict) -> None:
        spec["name"] = "skill-v2-beta"
        assert validate_spec(spec) == []


# ---------------------------------------------------------------------------
# generate_skill_md
# ---------------------------------------------------------------------------

class TestGenerateSkillMd:
    def test_generates_string(self, spec: dict) -> None:
        md = generate_skill_md(spec)
        assert isinstance(md, str)
        assert len(md) > 500

    def test_yaml_frontmatter_present(self, spec: dict) -> None:
        md = generate_skill_md(spec)
        assert md.startswith("---"), "SKILL.md must start with YAML frontmatter"

    def test_required_frontmatter_fields(self, spec: dict) -> None:
        md = generate_skill_md(spec)
        for field in ["name:", "description:", "version:", "author:", "license:", "trigger_keywords:"]:
            assert field in md, f"Missing frontmatter field: {field}"

    def test_all_required_sections_present(self, spec: dict) -> None:
        import re
        md = generate_skill_md(spec)
        for section_name, pattern in REQUIRED_SECTIONS:
            assert re.search(pattern, md, re.MULTILINE), f"Missing section: {section_name}"

    def test_skill_name_in_content(self, spec: dict) -> None:
        md = generate_skill_md(spec)
        assert spec["name"] in md

    def test_trigger_keywords_included(self, spec: dict) -> None:
        md = generate_skill_md(spec)
        for kw in spec["trigger_keywords"]:
            assert kw in md, f"trigger_keyword '{kw}' not found in SKILL.md"

    def test_output_structure_section_correct(self, spec: dict) -> None:
        md = generate_skill_md(spec)
        assert "report.md" in md
        assert "result.json" in md
        assert "commands.sh" in md
        assert "environment.yml" in md

    def test_os_values_use_process_platform_names(self, spec: dict) -> None:
        """
        ``scripts/lint_skills.py`` (run in CI) rejects ``macos``/``windows`` and
        requires Node's ``process.platform`` names. Ensure the generator emits
        values that pass that lint so no scaffolded skill is DOA.
        """
        md = generate_skill_md(spec)
        assert "os: [darwin" in md, "os field must use 'darwin' (not 'macos')"
        assert "macos" not in md.split("os:")[1].split("\n")[0]
        assert "windows" not in md.split("os:")[1].split("\n")[0]


# ---------------------------------------------------------------------------
# generate_skill_py
# ---------------------------------------------------------------------------

class TestGenerateSkillPy:
    def test_generates_string(self, spec: dict) -> None:
        py = generate_skill_py(spec)
        assert isinstance(py, str)
        assert len(py) > 200

    def test_has_shebang(self, spec: dict) -> None:
        py = generate_skill_py(spec)
        assert py.startswith("#!/usr/bin/env python3")

    def test_has_argparse(self, spec: dict) -> None:
        py = generate_skill_py(spec)
        assert "argparse" in py

    def test_has_demo_flag(self, spec: dict) -> None:
        py = generate_skill_py(spec)
        assert "--demo" in py

    def test_has_input_output_flags(self, spec: dict) -> None:
        py = generate_skill_py(spec)
        assert "--input" in py
        assert "--output" in py

    def test_has_run_function(self, spec: dict) -> None:
        py = generate_skill_py(spec)
        assert "def run(" in py

    def test_has_run_demo_function(self, spec: dict) -> None:
        py = generate_skill_py(spec)
        assert "def run_demo(" in py

    def test_has_main_guard(self, spec: dict) -> None:
        py = generate_skill_py(spec)
        assert 'if __name__ == "__main__"' in py

    def test_output_bundle_code_present(self, spec: dict) -> None:
        py = generate_skill_py(spec)
        assert "report.md" in py
        assert "result.json" in py
        assert "commands.sh" in py
        assert "environment.yml" in py


# ---------------------------------------------------------------------------
# generate_test_py
# ---------------------------------------------------------------------------

class TestGenerateTestPy:
    def test_generates_string(self, spec: dict) -> None:
        tpy = generate_test_py(spec)
        assert isinstance(tpy, str)
        assert len(tpy) > 100

    def test_has_four_standard_tests(self, spec: dict) -> None:
        tpy = generate_test_py(spec)
        assert "def test_demo_runs" in tpy
        assert "def test_report_generated" in tpy
        assert "def test_result_json_valid" in tpy
        assert "def test_reproducibility_bundle" in tpy

    def test_imports_run_and_run_demo(self, spec: dict) -> None:
        tpy = generate_test_py(spec)
        assert "from test_skill import run, run_demo" in tpy or "run, run_demo" in tpy


# ---------------------------------------------------------------------------
# generate_catalog_entry
# ---------------------------------------------------------------------------

class TestGenerateCatalogEntry:
    def test_entry_has_required_keys(self, spec: dict) -> None:
        entry = generate_catalog_entry(spec)
        for key in ["name", "cli_alias", "description", "version", "has_script", "has_tests", "has_demo"]:
            assert key in entry, f"Missing key in catalog entry: {key}"

    def test_name_matches_spec(self, spec: dict) -> None:
        entry = generate_catalog_entry(spec)
        assert entry["name"] == spec["name"]

    def test_cli_alias_matches_spec(self, spec: dict) -> None:
        entry = generate_catalog_entry(spec)
        assert entry["cli_alias"] == spec["cli_alias"]

    def test_flags_true(self, spec: dict) -> None:
        entry = generate_catalog_entry(spec)
        assert entry["has_script"] is True
        assert entry["has_tests"] is True
        assert entry["has_demo"] is True


# ---------------------------------------------------------------------------
# update_catalog_json
# ---------------------------------------------------------------------------

class TestUpdateCatalogJson:
    def _make_catalog(self, tmp_path: Path, skills: list[dict] | None = None) -> Path:
        p = tmp_path / "catalog.json"
        p.write_text(
            json.dumps({"skill_count": len(skills or []), "skills": skills or []}),
            encoding="utf-8",
        )
        return p

    def test_adds_new_skill(self, spec: dict, tmp_path: Path) -> None:
        catalog = self._make_catalog(tmp_path)
        entry = generate_catalog_entry(spec)
        result = update_catalog_json(catalog, entry)
        assert result is True
        data = json.loads(catalog.read_text())
        assert any(s["name"] == spec["name"] for s in data["skills"])

    def test_skill_count_updated(self, spec: dict, tmp_path: Path) -> None:
        catalog = self._make_catalog(tmp_path)
        entry = generate_catalog_entry(spec)
        update_catalog_json(catalog, entry)
        data = json.loads(catalog.read_text())
        assert data["skill_count"] == 1

    def test_duplicate_skipped(self, spec: dict, tmp_path: Path) -> None:
        entry = generate_catalog_entry(spec)
        catalog = self._make_catalog(tmp_path, skills=[entry])
        result = update_catalog_json(catalog, entry)
        assert result is False
        data = json.loads(catalog.read_text())
        assert data["skill_count"] == 1

    def test_invalid_json_returns_false(self, spec: dict, tmp_path: Path) -> None:
        catalog = tmp_path / "catalog.json"
        catalog.write_text("NOT JSON", encoding="utf-8")
        entry = generate_catalog_entry(spec)
        result = update_catalog_json(catalog, entry)
        assert result is False

    def test_preserves_unicode_in_other_skills(self, spec: dict, tmp_path: Path) -> None:
        """
        Regression guard: updating catalog.json must not re-encode non-ASCII
        characters already present in unrelated entries. Earlier versions used
        ``json.dumps`` with the default ``ensure_ascii=True``, so em-dashes in
        every other skill's description were escaped to ``\\u2014`` on every run.
        """
        pre_existing = {
            "name": "existing-skill",
            "description": "Analyses — with em-dash — intact",
        }
        catalog = self._make_catalog(tmp_path, skills=[pre_existing])
        entry = generate_catalog_entry(spec)
        update_catalog_json(catalog, entry)
        raw = catalog.read_text(encoding="utf-8")
        assert "—" in raw, "em-dash should survive the round-trip verbatim"
        assert "\\u2014" not in raw, "em-dash must not be escaped to \\u2014"


class TestGenerateCatalogEntryDependencies:
    def test_dependencies_copied_from_spec(self) -> None:
        spec = dict(MINIMAL_SPEC)
        spec["dependencies"] = {"required": ["pysam>=0.22", "biopython"], "optional": []}
        entry = generate_catalog_entry(spec)
        assert entry["dependencies"] == ["pysam>=0.22", "biopython"]

    def test_dependencies_defaults_to_empty(self) -> None:
        spec = dict(MINIMAL_SPEC)
        spec.pop("dependencies", None)
        entry = generate_catalog_entry(spec)
        assert entry["dependencies"] == []


# ---------------------------------------------------------------------------
# patch_clawbio_py
# ---------------------------------------------------------------------------

class TestPatchClawbioPy:
    """
    Fixture mirrors the real ``clawbio.py`` layout: the anchor comment sits
    *after* the ``}`` that closes ``SKILLS``, not inside the dict. An earlier
    version of this test built a fixture where the anchor lived inside the
    dict, which masked a bug that produced invalid Python in the real file.
    """

    _ANCHOR = "# Skills that run in the full-profile pipeline (order matters)"

    def _make_clawbio(self, tmp_path: Path, already_has_alias: bool = False) -> Path:
        p = tmp_path / "clawbio.py"
        extra = (
            '    "testskill": {\n'
            '        "script": SKILLS_DIR / "x" / "x.py",\n'
            "    },\n"
            if already_has_alias
            else ""
        )
        content = (
            "from pathlib import Path\n"
            "SKILLS_DIR = Path('skills')\n"
            "\n"
            "SKILLS = {\n"
            '    "existing": {\n'
            '        "script": SKILLS_DIR / "existing" / "existing.py",\n'
            '        "demo_args": ["--demo"],\n'
            '        "description": "Pre-existing entry",\n'
            "        \"allowed_extra_flags\": set(),\n"
            '        "accepts_genotypes": False,\n'
            "    },\n"
            f"{extra}"
            "}\n"
            "\n"
            f"{self._ANCHOR}\n"
            'FULL_PROFILE_PIPELINE = ["existing"]\n'
        )
        p.write_text(content, encoding="utf-8")
        return p

    def test_inserts_entry(self, spec: dict, tmp_path: Path) -> None:
        clawbio = self._make_clawbio(tmp_path)
        result = patch_clawbio_py(clawbio, spec)
        assert result is True
        source = clawbio.read_text()
        assert '"testskill"' in source

    def test_anchor_still_present(self, spec: dict, tmp_path: Path) -> None:
        clawbio = self._make_clawbio(tmp_path)
        patch_clawbio_py(clawbio, spec)
        assert self._ANCHOR in clawbio.read_text()

    def test_duplicate_alias_skipped(self, spec: dict, tmp_path: Path) -> None:
        clawbio = self._make_clawbio(tmp_path, already_has_alias=True)
        result = patch_clawbio_py(clawbio, spec)
        assert result is False

    def test_missing_anchor_returns_false(self, spec: dict, tmp_path: Path) -> None:
        clawbio = tmp_path / "clawbio.py"
        clawbio.write_text("SKILLS = {}\n", encoding="utf-8")
        result = patch_clawbio_py(clawbio, spec)
        assert result is False

    def test_patched_source_is_valid_python(self, spec: dict, tmp_path: Path) -> None:
        """
        Regression guard for the bug where the new entry was inserted *after*
        the closing ``}`` of ``SKILLS``, producing ``IndentationError`` when
        the file was imported.
        """
        import py_compile

        clawbio = self._make_clawbio(tmp_path)
        assert patch_clawbio_py(clawbio, spec) is True
        # py_compile raises PyCompileError on any syntax/indentation error
        py_compile.compile(str(clawbio), doraise=True)

    def test_entry_inserted_inside_dict(self, spec: dict, tmp_path: Path) -> None:
        """The new entry must appear between ``SKILLS = {`` and its close ``}``."""
        clawbio = self._make_clawbio(tmp_path)
        patch_clawbio_py(clawbio, spec)
        source = clawbio.read_text()
        open_idx  = source.index("SKILLS = {")
        # The dict close is the first `^}$` after the open.
        import re as _re
        close_match = _re.search(r"^\}$", source[open_idx:], _re.MULTILINE)
        assert close_match is not None
        close_idx = open_idx + close_match.start()
        alias_idx = source.index('"testskill"')
        assert open_idx < alias_idx < close_idx, (
            f"entry at {alias_idx} is not inside the SKILLS dict "
            f"(open={open_idx}, close={close_idx})"
        )

    def test_patch_on_real_clawbio_py(self, spec: dict, tmp_path: Path) -> None:
        """
        End-to-end guard: run the patch against a copy of the repo's actual
        ``clawbio.py`` and confirm it compiles. This would have caught the
        shipped bug because no synthetic fixture matches the real layout
        exactly.
        """
        import py_compile
        from shutil import copyfile

        real = Path(__file__).resolve().parents[3] / "clawbio.py"
        if not real.exists():
            pytest.skip("real clawbio.py not found — running outside the repo")
        target = tmp_path / "clawbio.py"
        copyfile(real, target)
        # Pick an alias we know isn't already used.
        spec["cli_alias"] = "__skillbuilder_probe__"
        assert patch_clawbio_py(target, spec) is True
        py_compile.compile(str(target), doraise=True)


# ---------------------------------------------------------------------------
# validate_skill_md
# ---------------------------------------------------------------------------

class TestValidateSkillMd:
    def test_generated_md_passes_all_checks(self, spec: dict, tmp_path: Path) -> None:
        md_path = tmp_path / "SKILL.md"
        md_path.write_text(generate_skill_md(spec), encoding="utf-8")
        results = validate_skill_md(md_path)
        failures = [(name, detail) for name, ok, detail in results if not ok]
        assert failures == [], f"Validation failures: {failures}"

    def test_empty_file_fails_all(self, tmp_path: Path) -> None:
        md_path = tmp_path / "SKILL.md"
        md_path.write_text("", encoding="utf-8")
        results = validate_skill_md(md_path)
        failures = [name for name, ok, _ in results if not ok]
        assert len(failures) > 0

    def test_missing_file_returns_error(self, tmp_path: Path) -> None:
        md_path = tmp_path / "NONEXISTENT.md"
        results = validate_skill_md(md_path)
        assert len(results) == 1
        _, ok, _ = results[0]
        assert ok is False


# ---------------------------------------------------------------------------
# scaffold_skill (integration)
# ---------------------------------------------------------------------------

class TestScaffoldSkill:
    def test_creates_expected_files(self, spec: dict, tmp_output: Path) -> None:
        manifest = scaffold_skill(spec, tmp_output, dry_run=False, repo_root=None)
        skill_dir = tmp_output / spec["name"]
        assert (skill_dir / "SKILL.md").exists()
        assert (skill_dir / f"{spec['name'].replace('-','_')}.py").exists()
        assert (skill_dir / "tests" / f"test_{spec['name'].replace('-','_')}.py").exists()
        assert (skill_dir / "examples" / "example_spec.json").exists()

    def test_manifest_has_expected_keys(self, spec: dict, tmp_output: Path) -> None:
        manifest = scaffold_skill(spec, tmp_output, dry_run=False, repo_root=None)
        for key in ["skill", "output_dir", "files_written", "files_skipped", "generated_at"]:
            assert key in manifest, f"Manifest missing key: {key}"

    def test_dry_run_does_not_write_files(self, spec: dict, tmp_output: Path) -> None:
        scaffold_skill(spec, tmp_output, dry_run=True, repo_root=None)
        skill_dir = tmp_output / spec["name"]
        assert not skill_dir.exists(), "dry-run should not write any files"

    def test_idempotent_second_run_skips(self, spec: dict, tmp_output: Path) -> None:
        scaffold_skill(spec, tmp_output, dry_run=False, repo_root=None)
        manifest2 = scaffold_skill(spec, tmp_output, dry_run=False, repo_root=None)
        assert len(manifest2["files_skipped"]) > 0, "Second run should skip already-existing files"

    def test_output_outside_repo_does_not_touch_registry(
        self, spec: dict, tmp_path: Path
    ) -> None:
        """
        Regression guard: running ``--demo --output /tmp/foo`` from inside a
        clone of the repo used to patch ``clawbio.py`` and ``catalog.json`` to
        point at ``/tmp/foo``, leaving the registry pointing at files that
        didn't exist inside the repo.

        Build a minimal fake "repo" with both sentinel files and a ``skills/``
        dir, then scaffold into a *sibling* directory. The registry files must
        stay untouched.
        """
        fake_repo = tmp_path / "repo"
        (fake_repo / "skills").mkdir(parents=True)
        (fake_repo / "scripts").mkdir()
        (fake_repo / "clawbio.py").write_text(
            "SKILLS = {\n}\n\n# Skills that run in the full-profile pipeline\n",
            encoding="utf-8",
        )
        catalog_original = '{"skill_count": 0, "skills": []}\n'
        (fake_repo / "skills" / "catalog.json").write_text(
            catalog_original, encoding="utf-8"
        )

        elsewhere = tmp_path / "somewhere_else"
        manifest = scaffold_skill(
            spec, elsewhere, dry_run=False, repo_root=fake_repo
        )
        assert manifest["catalog_updated"] is False
        assert manifest["clawbio_updated"] is False
        assert (fake_repo / "skills" / "catalog.json").read_text(
            encoding="utf-8"
        ) == catalog_original


# ---------------------------------------------------------------------------
# load_spec
# ---------------------------------------------------------------------------

class TestLoadSpec:
    def test_load_json(self, spec: dict, tmp_path: Path) -> None:
        p = tmp_path / "spec.json"
        p.write_text(json.dumps(spec), encoding="utf-8")
        loaded = load_spec(p)
        assert loaded["name"] == spec["name"]

    def test_load_invalid_json_exits(self, tmp_path: Path) -> None:
        p = tmp_path / "spec.json"
        p.write_text("{bad json", encoding="utf-8")
        with pytest.raises(SystemExit):
            load_spec(p)


# ---------------------------------------------------------------------------
# Demo spec sanity checks
# ---------------------------------------------------------------------------

class TestDemoSpec:
    def test_demo_spec_exists(self) -> None:
        assert DEMO_SPEC_PATH.exists(), f"demo_spec.json not found at {DEMO_SPEC_PATH}"

    def test_demo_spec_is_valid_json(self) -> None:
        data = json.loads(DEMO_SPEC_PATH.read_text(encoding="utf-8"))
        assert isinstance(data, dict)

    def test_demo_spec_passes_validation(self, demo_spec: dict) -> None:
        errors = validate_spec(demo_spec)
        assert errors == [], f"demo_spec.json validation errors: {errors}"

    def test_demo_spec_has_trigger_keywords(self, demo_spec: dict) -> None:
        assert len(demo_spec.get("trigger_keywords", [])) > 0


# ---------------------------------------------------------------------------
# run_demo (end-to-end)
# ---------------------------------------------------------------------------

class TestRunDemo:
    def test_demo_completes(self, tmp_output: Path) -> None:
        """run_demo() should complete without raising an exception."""
        run_demo(tmp_output, dry_run=False)

    def test_demo_creates_skill_directory(self, tmp_output: Path) -> None:
        run_demo(tmp_output, dry_run=False)
        skill_name = json.loads(DEMO_SPEC_PATH.read_text())["name"]
        assert (tmp_output / skill_name).is_dir()

    def test_demo_creates_report_md(self, tmp_output: Path) -> None:
        run_demo(tmp_output, dry_run=False)
        assert (tmp_output / "report.md").exists()

    def test_demo_creates_result_json(self, tmp_output: Path) -> None:
        run_demo(tmp_output, dry_run=False)
        result_file = tmp_output / "result.json"
        assert result_file.exists()
        data = json.loads(result_file.read_text())
        assert isinstance(data, dict)

    def test_demo_creates_reproducibility_bundle(self, tmp_output: Path) -> None:
        run_demo(tmp_output, dry_run=False)
        repro = tmp_output / "reproducibility"
        assert (repro / "commands.sh").exists()

    def test_demo_dry_run_no_files(self, tmp_output: Path) -> None:
        run_demo(tmp_output, dry_run=True)
        skill_name = json.loads(DEMO_SPEC_PATH.read_text())["name"]
        assert not (tmp_output / skill_name).exists(), "dry-run should not create files"


# ---------------------------------------------------------------------------
# _Jump
# ---------------------------------------------------------------------------

class TestJump:
    def test_raise_and_catch(self) -> None:
        with pytest.raises(_Jump):
            raise _Jump(3)

    def test_target_attribute(self) -> None:
        j = _Jump(5)
        assert j.target == 5

    def test_zero_target(self) -> None:
        j = _Jump(0)
        assert j.target == 0

    def test_is_exception_subclass(self) -> None:
        assert issubclass(_Jump, Exception)


# ---------------------------------------------------------------------------
# _save_draft / _load_draft / _clear_draft
# ---------------------------------------------------------------------------

class TestDraft:
    @pytest.fixture(autouse=True)
    def _patch_draft_path(self, tmp_path: Path, monkeypatch) -> None:
        """Redirect DRAFT_SPEC_PATH to a temp file for every test in this class."""
        monkeypatch.setattr(_sb, "DRAFT_SPEC_PATH", tmp_path / "draft.json")

    def test_save_and_load_roundtrip(self, spec: dict) -> None:
        _save_draft(spec)
        loaded = _load_draft()
        assert loaded == spec

    def test_load_returns_none_when_no_file(self) -> None:
        result = _load_draft()
        assert result is None

    def test_save_creates_file(self, spec: dict) -> None:
        _save_draft(spec)
        assert _sb.DRAFT_SPEC_PATH.exists()

    def test_clear_removes_file(self, spec: dict) -> None:
        _save_draft(spec)
        assert _sb.DRAFT_SPEC_PATH.exists()
        _clear_draft()
        assert not _sb.DRAFT_SPEC_PATH.exists()

    def test_clear_is_safe_when_no_file(self) -> None:
        _clear_draft()  # must not raise

    def test_saved_draft_is_valid_json(self, spec: dict) -> None:
        _save_draft(spec)
        raw = _sb.DRAFT_SPEC_PATH.read_text(encoding="utf-8")
        assert json.loads(raw) == spec


# ---------------------------------------------------------------------------
# _find_resume_cursor
# ---------------------------------------------------------------------------

class TestFindResumeCursor:
    def test_empty_spec_returns_zero(self) -> None:
        assert _find_resume_cursor({}) == 0

    def test_name_missing_returns_zero(self) -> None:
        assert _find_resume_cursor({"description": "foo"}) == 0

    def test_name_filled_description_missing_returns_one(self) -> None:
        assert _find_resume_cursor({"name": "foo"}) == 1

    def test_name_description_filled_author_missing_returns_two(self) -> None:
        assert _find_resume_cursor({"name": "foo", "description": "bar"}) == 2

    def test_all_str_required_filled_no_lists_returns_five(self) -> None:
        # name(0), description(1), author(2) filled → next required check is capabilities(5, list)
        s = {"name": "foo", "description": "bar", "author": "baz"}
        assert _find_resume_cursor(s) == 5

    def test_capabilities_empty_returns_five(self) -> None:
        s = {"name": "foo", "description": "bar", "author": "baz", "capabilities": []}
        assert _find_resume_cursor(s) == 5

    def test_all_required_fields_filled_returns_num_fields(self) -> None:
        s = {
            "name": "foo",
            "description": "bar",
            "author": "baz",
            "capabilities": ["cap1"],
            "trigger_keywords": ["kw"],
            "tags": ["t"],
        }
        assert _find_resume_cursor(s) == _NUM_FIELDS

    def test_num_fields_is_ten(self) -> None:
        assert _NUM_FIELDS == 10


# ---------------------------------------------------------------------------
# _ask
# ---------------------------------------------------------------------------

class TestAsk:
    def test_returns_typed_value(self, monkeypatch) -> None:
        monkeypatch.setattr("builtins.input", lambda _: "hello")
        assert _ask("Label") == "hello"

    def test_default_used_on_empty_input(self, monkeypatch) -> None:
        monkeypatch.setattr("builtins.input", lambda _: "")
        assert _ask("Label", default="fallback") == "fallback"

    def test_required_rejects_empty_then_accepts_value(self, monkeypatch) -> None:
        responses = iter(["", "valid"])
        monkeypatch.setattr("builtins.input", lambda _: next(responses))
        assert _ask("Label", required=True) == "valid"

    def test_back_raises_jump_with_target_zero(self, monkeypatch) -> None:
        # field_num=1 → max(0, 1-2) = 0
        monkeypatch.setattr("builtins.input", lambda _: "!back")
        with pytest.raises(_Jump) as exc:
            _ask("Label", field_num=1)
        assert exc.value.target == 0

    def test_numeric_jump_raises_jump(self, monkeypatch) -> None:
        # !3 → target = 3-1 = 2
        monkeypatch.setattr("builtins.input", lambda _: "!3")
        with pytest.raises(_Jump) as exc:
            _ask("Label", field_num=1)
        assert exc.value.target == 2

    def test_out_of_range_jump_loops_then_returns(self, monkeypatch) -> None:
        # !99 is out of range → error printed, loop continues; next returns "ok"
        responses = iter(["!99", "ok"])
        monkeypatch.setattr("builtins.input", lambda _: next(responses))
        assert _ask("Label", field_num=1) == "ok"

    def test_unknown_bang_command_loops(self, monkeypatch) -> None:
        responses = iter(["!notacommand", "ok"])
        monkeypatch.setattr("builtins.input", lambda _: next(responses))
        assert _ask("Label") == "ok"


# ---------------------------------------------------------------------------
# _ask_list
# ---------------------------------------------------------------------------

class TestAskList:
    def test_collects_items_until_blank(self, monkeypatch) -> None:
        responses = iter(["alpha", "beta", ""])
        monkeypatch.setattr("builtins.input", lambda _: next(responses))
        result = _ask_list("Items", field_num=1)
        assert result == ["alpha", "beta"]

    def test_empty_input_returns_empty_list(self, monkeypatch) -> None:
        monkeypatch.setattr("builtins.input", lambda _: "")
        assert _ask_list("Items") == []

    def test_back_raises_jump(self, monkeypatch) -> None:
        # field_num=1 → max(0, 1-2) = 0
        monkeypatch.setattr("builtins.input", lambda _: "!back")
        with pytest.raises(_Jump) as exc:
            _ask_list("Items", field_num=1)
        assert exc.value.target == 0

    def test_numeric_jump_raises_jump(self, monkeypatch) -> None:
        # !4 → target = 4-1 = 3
        monkeypatch.setattr("builtins.input", lambda _: "!4")
        with pytest.raises(_Jump) as exc:
            _ask_list("Items", field_num=2)
        assert exc.value.target == 3

    def test_single_item_list(self, monkeypatch) -> None:
        responses = iter(["only-item", ""])
        monkeypatch.setattr("builtins.input", lambda _: next(responses))
        assert _ask_list("Items") == ["only-item"]


# ---------------------------------------------------------------------------
# _prompt_field
# ---------------------------------------------------------------------------

class TestPromptField:
    def test_deps_raw_parsed_into_dependencies_required(self, spec: dict, monkeypatch) -> None:
        # Index 8 = _deps_raw (str kind) → stored into spec["dependencies"]["required"]
        monkeypatch.setattr("builtins.input", lambda _: "pysam>=0.22, biopython")
        result = _prompt_field(spec, 8)
        assert result["dependencies"]["required"] == ["pysam>=0.22", "biopython"]

    def test_deps_raw_empty_gives_empty_required_list(self, spec: dict, monkeypatch) -> None:
        monkeypatch.setattr("builtins.input", lambda _: "")
        result = _prompt_field(spec, 8)
        assert result["dependencies"]["required"] == []

    def test_chaining_partners_parsed_into_list(self, spec: dict, monkeypatch) -> None:
        # Index 9 = chaining_partners (str kind) → stored as list
        monkeypatch.setattr("builtins.input", lambda _: "fastqc, multiqc")
        result = _prompt_field(spec, 9)
        assert result["chaining_partners"] == ["fastqc", "multiqc"]

    def test_chaining_partners_empty_gives_empty_list(self, spec: dict, monkeypatch) -> None:
        monkeypatch.setattr("builtins.input", lambda _: "")
        result = _prompt_field(spec, 9)
        assert result["chaining_partners"] == []

    def test_deps_raw_preserves_optional_deps(self, spec: dict, monkeypatch) -> None:
        spec["dependencies"] = {"required": [], "optional": ["matplotlib"]}
        monkeypatch.setattr("builtins.input", lambda _: "numpy")
        result = _prompt_field(spec, 8)
        assert result["dependencies"]["optional"] == ["matplotlib"]
        assert result["dependencies"]["required"] == ["numpy"]


# ---------------------------------------------------------------------------
# _run_review_loop
# ---------------------------------------------------------------------------

class TestRunReviewLoop:
    @pytest.fixture(autouse=True)
    def _patch_draft_path(self, tmp_path: Path, monkeypatch) -> None:
        monkeypatch.setattr(_sb, "DRAFT_SPEC_PATH", tmp_path / "draft.json")

    def test_y_confirms_and_returns_spec(self, spec: dict, monkeypatch) -> None:
        monkeypatch.setattr("builtins.input", lambda _: "Y")
        result = _run_review_loop(spec)
        assert result["name"] == spec["name"]

    def test_empty_input_confirms(self, spec: dict, monkeypatch) -> None:
        monkeypatch.setattr("builtins.input", lambda _: "")
        result = _run_review_loop(spec)
        assert result["name"] == spec["name"]

    def test_yes_string_confirms(self, spec: dict, monkeypatch) -> None:
        monkeypatch.setattr("builtins.input", lambda _: "yes")
        result = _run_review_loop(spec)
        assert result["name"] == spec["name"]

    def test_n_aborts_with_sys_exit_zero(self, spec: dict, monkeypatch) -> None:
        monkeypatch.setattr("builtins.input", lambda _: "n")
        with pytest.raises(SystemExit) as exc_info:
            _run_review_loop(spec)
        assert exc_info.value.code == 0

    def test_abort_keyword_exits(self, spec: dict, monkeypatch) -> None:
        monkeypatch.setattr("builtins.input", lambda _: "abort")
        with pytest.raises(SystemExit) as exc_info:
            _run_review_loop(spec)
        assert exc_info.value.code == 0

    def test_invalid_choice_loops_then_confirm(self, spec: dict, monkeypatch) -> None:
        responses = iter(["what?", "Y"])
        monkeypatch.setattr("builtins.input", lambda _: next(responses))
        result = _run_review_loop(spec)
        assert result["name"] == spec["name"]

    def test_field_number_re_prompts_then_confirm(self, spec: dict, monkeypatch) -> None:
        # "1" → edit field 0 (name) → type new name → "Y" → return updated spec
        responses = iter(["1", "new-name", "Y"])
        monkeypatch.setattr("builtins.input", lambda _: next(responses))
        result = _run_review_loop(spec)
        assert result["name"] == "new-name"


# ---------------------------------------------------------------------------
# scaffold_skill — edit mode (ask_overwrite=True)
# ---------------------------------------------------------------------------

class TestScaffoldSkillEditMode:
    def test_overwrite_yes_updates_all_files(self, spec: dict, tmp_path: Path, monkeypatch) -> None:
        output = tmp_path / "out"
        # First run: create files fresh
        scaffold_skill(spec, output, dry_run=False, repo_root=None)
        # Second run with ask_overwrite=True, answer "y" to all overwrite prompts
        monkeypatch.setattr("builtins.input", lambda _: "y")
        manifest = scaffold_skill(spec, output, dry_run=False, repo_root=None, ask_overwrite=True)
        # SKILL.md, .py, test_.py all "updated"; example_spec.json force-refreshed → 4 files
        assert len(manifest["files_skipped"]) == 0
        assert any("SKILL.md" in f for f in manifest["files_written"])

    def test_overwrite_no_skips_existing_files(self, spec: dict, tmp_path: Path, monkeypatch) -> None:
        output = tmp_path / "out"
        scaffold_skill(spec, output, dry_run=False, repo_root=None)
        # Answer "n" to all overwrite prompts
        monkeypatch.setattr("builtins.input", lambda _: "n")
        manifest = scaffold_skill(spec, output, dry_run=False, repo_root=None, ask_overwrite=True)
        # SKILL.md, .py, test_.py all skipped
        assert len(manifest["files_skipped"]) == 3
        # example_spec.json is always force-refreshed (no prompt, always written)
        assert any("example_spec.json" in f for f in manifest["files_written"])

    def test_example_spec_content_is_current_spec(self, spec: dict, tmp_path: Path, monkeypatch) -> None:
        output = tmp_path / "out"
        scaffold_skill(spec, output, dry_run=False, repo_root=None)
        # Mutate the spec and run in edit mode
        spec["description"] = "Updated description for edit test"
        monkeypatch.setattr("builtins.input", lambda _: "n")
        scaffold_skill(spec, output, dry_run=False, repo_root=None, ask_overwrite=True)
        example_path = output / spec["name"] / "examples" / "example_spec.json"
        content = json.loads(example_path.read_text(encoding="utf-8"))
        assert content["description"] == "Updated description for edit test"


# ---------------------------------------------------------------------------
# Agent mode tests
# ---------------------------------------------------------------------------

# Paths used by all agent tests — use sys.executable so the path is always
# valid on whatever OS/environment pytest is running under.
_PYTHON = sys.executable
_SKILL_BUILDER_PY = str(Path(__file__).resolve().parent.parent / "skill_builder.py")

_AGENT_SPEC = {
    "name": "agent-test-skill",
    "description": "A skill created by the agent mode test suite",
    "author": "Test Author",
    "version": "0.1.0",
    "domain": "genomics",
    "cli_alias": "agenttest",
    "tags": ["agent", "test"],
    "trigger_keywords": ["agent test", "run agent"],
    "capabilities": ["Do something via agent", "Do something else via agent"],
    "input_formats": [],
    "dependencies": {"required": [], "optional": []},
    "chaining_partners": [],
}


@pytest.fixture
def agent_spec_file(tmp_path: Path) -> Path:
    """Write the agent spec to a temp JSON file and return the path."""
    spec_path = tmp_path / "agent_spec.json"
    spec_path.write_text(json.dumps(_AGENT_SPEC), encoding="utf-8")
    return spec_path


class TestAgentMode:
    """End-to-end tests for --agent mode via subprocess (full CLI path)."""

    def _run(self, args: list[str], *, input_text: str | None = None, tmp_path: Path) -> subprocess.CompletedProcess:
        """Helper: run skill_builder.py with the given extra args and capture output."""
        cmd = [_PYTHON, _SKILL_BUILDER_PY] + args + ["--output", str(tmp_path / "out")]
        return subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            input=input_text,
        )

    def test_agent_valid_spec_exits_zero(self, agent_spec_file: Path, tmp_path: Path) -> None:
        result = self._run(["--agent", "--input", str(agent_spec_file)], tmp_path=tmp_path)
        assert result.returncode == 0, f"Expected exit 0, got {result.returncode}.\nstderr: {result.stderr}"

    def test_agent_stdout_is_valid_json(self, agent_spec_file: Path, tmp_path: Path) -> None:
        result = self._run(["--agent", "--input", str(agent_spec_file)], tmp_path=tmp_path)
        try:
            data = json.loads(result.stdout)
            assert isinstance(data, dict)
        except json.JSONDecodeError as exc:
            pytest.fail(f"stdout was not valid JSON: {exc}\nstdout: {result.stdout!r}")

    def test_agent_status_ok(self, agent_spec_file: Path, tmp_path: Path) -> None:
        result = self._run(["--agent", "--input", str(agent_spec_file)], tmp_path=tmp_path)
        data = json.loads(result.stdout)
        assert data["status"] == "ok", f"Expected status 'ok', got: {data.get('status')!r}"

    def test_agent_manifest_keys(self, agent_spec_file: Path, tmp_path: Path) -> None:
        result = self._run(["--agent", "--input", str(agent_spec_file)], tmp_path=tmp_path)
        data = json.loads(result.stdout)
        required_keys = {"skill", "output_dir", "files_written", "files_skipped", "validation", "generated_at"}
        missing = required_keys - data.keys()
        assert not missing, f"Manifest missing keys: {missing}"

    def test_agent_validation_all_pass(self, agent_spec_file: Path, tmp_path: Path) -> None:
        result = self._run(["--agent", "--input", str(agent_spec_file)], tmp_path=tmp_path)
        data = json.loads(result.stdout)
        val = data["validation"]
        assert val["failed"] == 0, f"Expected 0 failures, got: {val['failures']}"
        assert val["passed"] == 13, f"Expected 13 passed, got: {val['passed']}"

    def test_agent_invalid_spec_exits_nonzero(self, tmp_path: Path) -> None:
        # Spec missing required 'name' field
        bad_spec = {k: v for k, v in _AGENT_SPEC.items() if k != "name"}
        bad_spec_file = tmp_path / "bad_spec.json"
        bad_spec_file.write_text(json.dumps(bad_spec), encoding="utf-8")
        result = self._run(["--agent", "--input", str(bad_spec_file)], tmp_path=tmp_path)
        assert result.returncode != 0, "Expected non-zero exit for invalid spec"
        data = json.loads(result.stdout)
        assert data["status"] == "error", f"Expected status 'error', got: {data.get('status')!r}"
        assert len(data["errors"]) > 0, "Expected at least one error message"

    def test_agent_dry_run_includes_generated_content(self, agent_spec_file: Path, tmp_path: Path) -> None:
        result = self._run(["--agent", "--dry-run", "--input", str(agent_spec_file)], tmp_path=tmp_path)
        assert result.returncode == 0, f"Exit {result.returncode}.\nstderr: {result.stderr}"
        data = json.loads(result.stdout)
        assert data.get("dry_run") is True, "Expected dry_run=true in manifest"
        assert "generated_content" in data, "Expected 'generated_content' key in dry-run manifest"
        script_name = _AGENT_SPEC["name"].replace("-", "_")
        expected_keys = {
            "SKILL.md",
            f"{script_name}.py",
            f"tests/test_{script_name}.py",
            "examples/example_spec.json",
        }
        missing = expected_keys - data["generated_content"].keys()
        assert not missing, f"generated_content missing keys: {missing}"
        # No files should be written to disk in dry-run mode
        skill_dir = tmp_path / "out" / _AGENT_SPEC["name"]
        assert not skill_dir.exists(), f"skill_dir should not exist in dry-run, but found: {skill_dir}"

    def test_agent_stdin_pipe(self, tmp_path: Path) -> None:
        spec_json = json.dumps(_AGENT_SPEC)
        result = self._run(["--agent", "--input", "-"], input_text=spec_json, tmp_path=tmp_path)
        assert result.returncode == 0, f"Exit {result.returncode}.\nstderr: {result.stderr}"
        data = json.loads(result.stdout)
        assert data["status"] == "ok", f"Expected status 'ok', got: {data.get('status')!r}"
