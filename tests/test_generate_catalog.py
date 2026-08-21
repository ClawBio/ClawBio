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
