import importlib.util
import json
from pathlib import Path


def _load_generate_catalog_module():
    module_path = Path(__file__).resolve().parents[1] / "scripts" / "generate_catalog.py"
    spec = importlib.util.spec_from_file_location("generate_catalog", module_path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_load_skills_registry_reads_package_cli_registry():
    generate_catalog = _load_generate_catalog_module()

    registered_aliases = generate_catalog.load_skills_registry()

    assert "fastreer" in registered_aliases
    assert "analyze-fasta" in registered_aliases
    assert "pharmgx" in registered_aliases


def test_build_catalog_populates_cli_aliases_from_package_registry():
    generate_catalog = _load_generate_catalog_module()

    catalog = {entry["name"]: entry for entry in generate_catalog.build_catalog()}

    assert catalog["fastreer"]["cli_alias"] == "fastreer"
    assert catalog["fastreer"]["demo_command"] == "python clawbio.py run fastreer --demo"
    assert catalog["analyze-fasta"]["cli_alias"] == "analyze-fasta"
    assert catalog["analyze-fasta"]["demo_command"] == "python clawbio.py run analyze-fasta --demo"


def test_build_catalog_adds_objective_maturity_tiers():
    generate_catalog = _load_generate_catalog_module()

    catalog = {entry["name"]: entry for entry in generate_catalog.build_catalog()}

    assert catalog["pharmgx-reporter"]["maturity_tier"] == "ci-validated"
    assert catalog["pharmgx-reporter"]["maturity_evidence"] == {
        "has_skill_md": True,
        "has_script": True,
        "has_tests": True,
        "has_demo": True,
        "cli_registered": True,
        "ci_tested": True,
        "benchmark_validated": False,
    }

    assert catalog["fastreer"]["maturity_tier"] == "cli-registered"
    assert catalog["fastreer"]["maturity_evidence"]["cli_registered"] is True
    assert catalog["fastreer"]["maturity_evidence"]["ci_tested"] is False

    assert catalog["claw-semantic-sim"]["maturity_tier"] == "spec-only"
    assert catalog["claw-semantic-sim"]["maturity_evidence"]["has_script"] is False


def test_checked_in_catalog_is_current():
    generate_catalog = _load_generate_catalog_module()
    root = Path(__file__).resolve().parents[1]
    catalog_path = root / "skills" / "catalog.json"

    checked_in = json.loads(catalog_path.read_text(encoding="utf-8"))
    generated_skills = generate_catalog.build_catalog()

    assert checked_in["skill_count"] == len(generated_skills)
    assert checked_in["skills"] == generated_skills


def test_fallback_demo_script_selection_is_deterministic():
    generate_catalog = _load_generate_catalog_module()
    root = Path(__file__).resolve().parents[1]

    assert (
        generate_catalog.select_demo_script(
            root / "skills" / "clinical-trial-finder",
            "clinical-trial-finder",
        ).name
        == "clinical_trial_finder.py"
    )
    assert (
        generate_catalog.select_demo_script(
            root / "skills" / "turingdb-graph",
            "turingdb-graph",
        ).name
        == "turingdb_graph.py"
    )


# --- Upstream artefact licences are declared separately from the wrapper's own
# licence. Putting MIT on a wrapper does not make non-commercial weights
# commercially usable, so the catalog must carry the artefact terms as their
# own fields rather than letting `license` imply anything about them.

def test_normalize_carries_model_and_data_licences():
    generate_catalog = _load_generate_catalog_module()

    normalized = generate_catalog.normalize_skill_metadata(
        {
            "name": "demo",
            "license": "MIT",
            "model_license": "cc-by-nc-sa-4.0",
            "data_license": "CC0-1.0",
        }
    )

    assert normalized["license"] == "MIT"
    assert normalized["model_license"] == "cc-by-nc-sa-4.0"
    assert normalized["data_license"] == "CC0-1.0"


def test_normalize_reads_artefact_licences_nested_under_metadata():
    """Where new skills have to put them.

    The agentskills spec allows six top-level frontmatter keys and
    `agentskills validate` fails a skill on any other, so `model_license`
    beside `license` is not a legal placement. `metadata` is the free-form key.
    """
    generate_catalog = _load_generate_catalog_module()

    normalized = generate_catalog.normalize_skill_metadata(
        {
            "name": "demo",
            "license": "MIT",
            "metadata": {
                "model_license": "cc-by-nc-sa-4.0",
                "data_license": "CC0-1.0",
            },
        }
    )

    assert normalized["license"] == "MIT"
    assert normalized["model_license"] == "cc-by-nc-sa-4.0"
    assert normalized["data_license"] == "CC0-1.0"


def test_normalize_defaults_artefact_licences_to_empty():
    """A skill with no third-party artefact declares nothing; absence is not
    a claim that the artefact is permissive."""
    generate_catalog = _load_generate_catalog_module()

    normalized = generate_catalog.normalize_skill_metadata(
        {"name": "demo", "license": "MIT"}
    )

    assert normalized["model_license"] == ""
    assert normalized["data_license"] == ""


def test_frontmatter_parser_reads_artefact_licences():
    """Both fields survive frontmatter parsing end to end."""
    generate_catalog = _load_generate_catalog_module()

    raw = "\n".join(
        [
            "---",
            "name: demo",
            "description: A demo skill",
            "license: MIT",
            "model_license: cc-by-nc-sa-4.0",
            "data_license: CC0-1.0",
            "---",
            "",
            "## Trigger",
        ]
    )

    parsed = generate_catalog.parse_yaml_frontmatter(raw)

    assert parsed["model_license"] == "cc-by-nc-sa-4.0"
    assert parsed["data_license"] == "CC0-1.0"


def test_frontmatter_parser_reads_nested_artefact_licences():
    """The no-yaml fallback path has to find them under `metadata` too."""
    generate_catalog = _load_generate_catalog_module()

    raw = "\n".join(
        [
            "---",
            "name: demo",
            "description: A demo skill",
            "license: MIT",
            "metadata:",
            "  model_license: cc-by-nc-sa-4.0",
            "  data_license: CC0-1.0",
            "---",
            "",
            "## Trigger",
        ]
    )

    normalized = generate_catalog.normalize_skill_metadata(
        generate_catalog.parse_yaml_frontmatter(raw)
    )

    # `license` still resolves to the wrapper licence, not to `model_license`.
    assert normalized["license"] == "MIT"
    assert normalized["model_license"] == "cc-by-nc-sa-4.0"
    assert normalized["data_license"] == "CC0-1.0"


# ---------------------------------------------------------------------------
# Issue #359: the generator derived every maturity signal from Python file
# patterns, so a skill written in R or shell could never register above
# spec-only however much it shipped. claw-amplicon-qc (#329) ships a
# 1,274-line R script, a 17-assertion shell harness and a working --demo, and
# generated as planned / spec-only / all-false / demo_command null.
# ---------------------------------------------------------------------------


def _skill(tmp_path, name, files, frontmatter="name: x\n"):
    """Build a throwaway skill directory. Returns its path."""
    d = tmp_path / name
    (d / "tests").mkdir(parents=True)
    (d / "SKILL.md").write_text(f"---\n{frontmatter}---\n\n# {name}\n")
    for rel, body in files.items():
        p = d / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(body)
    return d


class TestEntryPointDetection:
    """`has_script` must not mean `has_python`."""

    def test_python_entry_point_is_found(self, tmp_path):
        gc = _load_generate_catalog_module()
        d = _skill(tmp_path, "demo-skill", {"demo_skill.py": "print(1)\n"})
        found = gc.detect_entry_point(d, "demo-skill")
        assert found is not None and found.name == "demo_skill.py"

    def test_r_entry_point_is_found(self, tmp_path):
        gc = _load_generate_catalog_module()
        d = _skill(tmp_path, "amplicon-qc", {"amplicon_qc.R": "cat('hi')\n"})
        found = gc.detect_entry_point(d, "amplicon-qc")
        assert found is not None and found.name == "amplicon_qc.R"

    def test_shell_entry_point_is_found(self, tmp_path):
        gc = _load_generate_catalog_module()
        d = _skill(tmp_path, "pipeline", {"pipeline.sh": "#!/bin/sh\n"})
        found = gc.detect_entry_point(d, "pipeline")
        assert found is not None and found.name == "pipeline.sh"

    def test_a_skill_with_no_code_has_no_entry_point(self, tmp_path):
        """claw-semantic-sim: SKILL.md only. This is the #353 case."""
        gc = _load_generate_catalog_module()
        d = _skill(tmp_path, "spec-only-skill", {})
        assert gc.detect_entry_point(d, "spec-only-skill") is None

    def test_python_wins_when_a_skill_ships_both(self, tmp_path):
        """Deterministic precedence, so the catalog cannot flap between runs."""
        gc = _load_generate_catalog_module()
        d = _skill(tmp_path, "both", {"both.py": "x\n", "both.R": "y\n"})
        assert gc.detect_entry_point(d, "both").suffix == ".py"

    def test_entry_point_in_a_subdirectory_is_found(self, tmp_path):
        """xena-tcga-gene-query keeps its script in `scripts/`. The old rglob
        found it; a root-only glob silently demoted the skill."""
        gc = _load_generate_catalog_module()
        d = _skill(tmp_path, "xena-thing", {"scripts/query_tcga_api.py": "x\n"})
        found = gc.detect_entry_point(d, "xena-thing")
        assert found is not None and found.name == "query_tcga_api.py"

    def test_a_root_entry_point_beats_a_nested_one(self, tmp_path):
        gc = _load_generate_catalog_module()
        d = _skill(tmp_path, "s", {"s.py": "x\n", "scripts/helper.py": "y\n"})
        assert gc.detect_entry_point(d, "s").name == "s.py"

    def test_a_nested_test_file_is_never_an_entry_point(self, tmp_path):
        gc = _load_generate_catalog_module()
        d = _skill(tmp_path, "s", {"tests/helpers/mock_thing.py": "x\n"})
        assert gc.detect_entry_point(d, "s") is None

    def test_tests_and_helpers_are_not_mistaken_for_entry_points(self, tmp_path):
        gc = _load_generate_catalog_module()
        d = _skill(tmp_path, "helper-only", {
            "tests/test_helper_only.py": "def test_x(): pass\n",
            "__init__.py": "",
            "api.py": "x = 1\n",
        })
        assert gc.detect_entry_point(d, "helper-only") is None


class TestTestDetection:
    """`has_tests` must recognise a shell harness, and must not be satisfied
    by an empty file. The zero-byte case granted the `tested` tier to
    my-awesome-skill (#345), whose test file is the canonical empty blob."""

    def test_pytest_suite_counts(self, tmp_path):
        gc = _load_generate_catalog_module()
        d = _skill(tmp_path, "s", {"tests/test_s.py": "def test_a(): pass\n"})
        assert gc.detect_tests(d) is True

    def test_shell_harness_counts(self, tmp_path):
        gc = _load_generate_catalog_module()
        d = _skill(tmp_path, "s", {"tests/run_test.sh": "#!/bin/bash\nexit 0\n"})
        assert gc.detect_tests(d) is True

    def test_r_test_counts(self, tmp_path):
        gc = _load_generate_catalog_module()
        d = _skill(tmp_path, "s", {"tests/test_s.R": "stopifnot(TRUE)\n"})
        assert gc.detect_tests(d) is True

    def test_an_empty_test_file_does_not_count(self, tmp_path):
        gc = _load_generate_catalog_module()
        d = _skill(tmp_path, "s", {"tests/test_s.py": ""})
        assert gc.detect_tests(d) is False

    def test_a_whitespace_only_test_file_does_not_count(self, tmp_path):
        gc = _load_generate_catalog_module()
        d = _skill(tmp_path, "s", {"tests/test_s.py": "\n\n   \n"})
        assert gc.detect_tests(d) is False

    def test_no_tests_directory(self, tmp_path):
        gc = _load_generate_catalog_module()
        d = tmp_path / "bare"
        d.mkdir()
        assert gc.detect_tests(d) is False


class TestDemoCommandUsesTheRightInterpreter:
    """A `.R` entry point invoked with `python` is a command that cannot run."""

    def test_python_skill_gets_python(self, tmp_path):
        gc = _load_generate_catalog_module()
        d = _skill(tmp_path, "s", {"s.py": "x\n"})
        assert gc.demo_invocation(gc.detect_entry_point(d, "s")).startswith("python ")

    def test_r_skill_gets_rscript(self, tmp_path):
        gc = _load_generate_catalog_module()
        d = _skill(tmp_path, "amplicon-qc", {"amplicon_qc.R": "x\n"})
        cmd = gc.demo_invocation(gc.detect_entry_point(d, "amplicon-qc"))
        assert cmd.startswith("Rscript ")
        assert "python" not in cmd

    def test_shell_skill_gets_bash(self, tmp_path):
        gc = _load_generate_catalog_module()
        d = _skill(tmp_path, "s", {"s.sh": "x\n"})
        assert gc.demo_invocation(gc.detect_entry_point(d, "s")).startswith("bash ")


class TestStatusAndEvidenceCannotContradict:
    """MVP_FOLDERS is hand-maintained and its own comment defines membership
    as "have working Python". claw-semantic-sim was in it with no code at all,
    so `status: mvp` sat beside `maturity_tier: spec-only` in one object."""

    def test_no_mvp_skill_lacks_an_entry_point(self):
        gc = _load_generate_catalog_module()
        catalog = {e["name"]: e for e in gc.build_catalog()}
        offenders = [
            name for name in gc.MVP_FOLDERS
            if name in catalog and not catalog[name]["has_script"]
        ]
        assert offenders == [], (
            f"MVP_FOLDERS claims these ship code, and they do not: {offenders}"
        )

    def test_status_mvp_implies_scripted_or_better(self):
        gc = _load_generate_catalog_module()
        for entry in gc.build_catalog():
            if entry["status"] == "mvp":
                assert entry["maturity_tier"] != "spec-only", (
                    f"{entry['name']} is advertised mvp but has no evidence of code"
                )
