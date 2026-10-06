"""SKILL.md and schema conformance for metagenome-mag-pipeline.

These are the checks AGENTS.md's conformance table turns into CI: frontmatter
parses, `name` matches the folder, version is semver, the trigger list is loud
and does not steal `claw-metagenomics` traffic, every required section exists,
the gotchas are populated, the exact disclaimer is present, the file stays under
500 lines, and every error code the interface contracts name exists.
"""

from __future__ import annotations

import importlib
import importlib.util
import re
import sys
from pathlib import Path

import pytest
import yaml

SKILL_DIR = Path(__file__).resolve().parents[1]
if str(SKILL_DIR) not in sys.path:
    sys.path.insert(0, str(SKILL_DIR))

SKILL_MD = SKILL_DIR / "SKILL.md"


def load(name: str):
    """Import `magpipe_<name>` from the skill directory."""
    return importlib.import_module(f"magpipe_{name}")

# Sections AGENTS.md's conformance table requires, plus the template's own.
REQUIRED_SECTIONS = [
    "## Trigger",
    "## Why This Exists",
    "## Core Capabilities",
    "## Scope",
    "## Input Formats",
    "## Workflow",
    "## CLI Reference",
    "## Demo",
    "## Algorithm / Methodology",
    "## Example Output",
    "## Output Structure",
    "## Dependencies",
    "## Gotchas",
    "## Safety",
    "## Agent Boundary",
    "## Integration with Bio Orchestrator",
    "## Maintenance",
    "## Citations",
]

# Bare words that belong to skills/claw-metagenomics (read-based Kraken2/RGI/HUMAnN3).
FORBIDDEN_KEYWORDS = ["metagenomics", "microbiome", "Kraken2", "RGI", "CARD", "HUMAnN3"]

# Every code named in the plan's interface contract.
CONTRACT_ERROR_CODES = [
    "MISSING_INPUT",
    "INVALID_SAMPLESHEET",
    "MISSING_FASTQ",
    "REMOTE_INPUT_NOT_ALLOWED",
    "INVALID_PRESET",
    "INVALID_CONFIG",
    "MISSING_DATABASE",
    "ASSEMBLY_REQUIRES_PAIRED",
    "MISSING_SNAKEMAKE",
    "SNAKEMAKE_VERSION_TOO_OLD",
    "MISSING_CONDA",
    "PIPELINE_SOURCE_INVALID",
    "PIPELINE_FETCH_FAILED",
    "DEMO_REQUIRES_NETWORK",
    "OUTPUT_DIR_NOT_EMPTY",
    "OUTPUT_DIR_NOT_WRITABLE",
    "EXECUTION_FAILED",
    "EXECUTION_TIMEOUT",
    "EXPECTED_OUTPUTS_NOT_FOUND",
    "UNEXPECTED_ERROR",
]


@pytest.fixture(scope="module")
def text() -> str:
    return SKILL_MD.read_text(encoding="utf-8")


@pytest.fixture(scope="module")
def frontmatter(text: str) -> dict:
    match = re.match(r"^---\n(.*?)\n---\n", text, re.DOTALL)
    assert match, "SKILL.md has no YAML frontmatter block"
    return yaml.safe_load(match.group(1))


class TestFrontmatter:
    def test_frontmatter_parses(self, frontmatter: dict) -> None:
        assert isinstance(frontmatter, dict)

    def test_name_matches_folder(self, frontmatter: dict) -> None:
        assert frontmatter["name"] == SKILL_DIR.name

    def test_version_is_semver(self, frontmatter: dict) -> None:
        assert re.fullmatch(r"\d+\.\d+\.\d+", str(frontmatter["metadata"]["version"]))

    def test_author_present(self, frontmatter: dict) -> None:
        assert str(frontmatter["metadata"]["author"]).strip()

    def test_description_is_one_line_and_specific(self, frontmatter: dict) -> None:
        description = str(frontmatter["description"]).strip()
        assert description and "\n" not in description
        assert "metagenom" in description.lower()

    def test_license_is_mit(self, frontmatter: dict) -> None:
        assert frontmatter["license"] == "MIT"

    def test_linux_only(self, frontmatter: dict) -> None:
        assert frontmatter["metadata"]["openclaw"]["os"] == ["linux"]

    def test_inputs_and_outputs_declared(self, frontmatter: dict) -> None:
        metadata = frontmatter["metadata"]
        assert metadata["inputs"][0]["name"] == "input_path"
        for entry in metadata["inputs"]:
            assert entry["format"], f"input {entry['name']} declares no format"
            assert isinstance(entry["required"], bool)
        assert {entry["name"] for entry in metadata["outputs"]} >= {"report", "result"}
        assert {entry["format"][0] for entry in metadata["outputs"]} == {"md", "json"}

    def test_has_demo_data_and_cli_endpoint(self, frontmatter: dict) -> None:
        assert frontmatter["metadata"]["demo_data"][0]["path"]
        assert "metagenome_mag_pipeline.py" in frontmatter["metadata"]["endpoints"]["cli"]

    def test_requires_snakemake_and_conda(self, frontmatter: dict) -> None:
        bins = frontmatter["metadata"]["openclaw"]["requires"]["bins"]
        assert "snakemake" in bins and "conda" in bins


@pytest.fixture(scope="module")
def keywords(frontmatter: dict) -> list[str]:
    return list(frontmatter["metadata"]["openclaw"]["trigger_keywords"])


class TestTriggerKeywords:
    def test_at_least_three(self, keywords: list[str]) -> None:
        assert len(keywords) >= 3

    def test_names_the_upstream_tools(self, keywords: list[str]) -> None:
        joined = " ".join(keywords).lower()
        for expected in ("mag", "metaspades", "metaphlan", "checkm2", "gtdb-tk", "das tool"):
            assert expected in joined, f"no trigger keyword mentions {expected!r}"

    @pytest.mark.parametrize("forbidden", FORBIDDEN_KEYWORDS)
    def test_does_not_steal_claw_metagenomics_traffic(
        self, keywords: list[str], forbidden: str
    ) -> None:
        exact = [k for k in keywords if k.strip().lower() == forbidden.lower()]
        assert not exact, (
            f"{forbidden!r} is a bare keyword and belongs to skills/claw-metagenomics; "
            "keep this skill's triggers about assembly, binning and MAG recovery"
        )


class TestSections:
    @pytest.mark.parametrize("section", REQUIRED_SECTIONS)
    def test_section_present(self, text: str, section: str) -> None:
        assert section in text, f"SKILL.md is missing {section}"

    def test_workflow_is_numbered(self, text: str) -> None:
        body = text.split("## Workflow", 1)[1].split("\n## ", 1)[0]
        steps = [line for line in body.splitlines() if re.match(r"^\s*\d+\.\s+\S", line)]
        assert len(steps) >= 5, "## Workflow must be numbered steps, not prose"

    def test_gotchas_at_least_three(self, text: str) -> None:
        body = text.split("## Gotchas", 1)[1].split("\n## ", 1)[0]
        gotchas = [line for line in body.splitlines() if line.strip().startswith("- ")]
        assert len(gotchas) >= 3

    def test_gotchas_cover_the_five_documented_mistakes(self, text: str) -> None:
        body = text.split("## Gotchas", 1)[1].split("\n## ", 1)[0].lower()
        for phrase in ("relative abundance", "not zero", "high-quality mag", "synthetic", "database"):
            assert phrase in body, f"## Gotchas does not cover {phrase!r}"

    def test_gotchas_cover_the_real_run_lessons(self, text: str) -> None:
        body = text.split("## Gotchas", 1)[1].split("\n## ", 1)[0].lower()
        for phrase in ("two binners", "35x", "gtdb-tk database"):
            assert phrase in body, f"## Gotchas does not cover {phrase!r}"

    def test_example_output_shows_a_rendered_report(self, text: str) -> None:
        body = text.split("## Example Output", 1)[1]
        rendered = re.search(r"```markdown\n(.*?)```", body, re.DOTALL)
        assert rendered, "## Example Output must show a rendered report, not a description"
        assert "| " in rendered.group(1), "the rendered report must show real tables"

    def test_output_structure_lists_the_reproducibility_bundle(self, text: str) -> None:
        body = text.split("## Output Structure", 1)[1].split("\n## ", 1)[0]
        assert "reproducibility/" in body
        assert "commands.sh" in body and "environment.yml" in body and "checksums.sha256" in body


class TestContent:
    def test_disclaimer_matches_the_schema_constant(self, text: str) -> None:
        # Markdown wraps the disclaimer across lines; compare on collapsed space.
        collapsed = " ".join(text.split())
        assert " ".join(load("schemas").DISCLAIMER.split()) in collapsed

    def test_under_500_lines(self, text: str) -> None:
        assert len(text.splitlines()) < 500

    def test_documents_the_pinned_ref(self, text: str) -> None:
        assert load("schemas").PINNED_REF in text

    def test_documents_the_pinned_tag(self, text: str) -> None:
        assert load("schemas").PINNED_TAG in text

    def test_documents_every_cli_flag(self, text: str) -> None:
        body = text.split("## CLI Reference", 1)[1].split("\n## ", 1)[0]
        spec = importlib.util.spec_from_file_location(
            "metagenome_mag_pipeline", SKILL_DIR / "metagenome_mag_pipeline.py"
        )
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        missing = [
            flag
            for action in module.build_parser()._actions
            for flag in action.option_strings
            if flag not in ("-h", "--help") and flag not in body
        ]
        assert not missing, f"## CLI Reference does not document {missing}"


class TestErrorCodes:
    @pytest.mark.parametrize("code", CONTRACT_ERROR_CODES)
    def test_error_code_exists(self, code: str) -> None:
        assert getattr(load("errors").ErrorCode, code) == code


class TestDemoAndExamples:
    def test_demo_readme_exists(self) -> None:
        assert (SKILL_DIR / "demo" / "README.md").is_file()

    def test_example_samplesheet_is_upstreams(self) -> None:
        path = SKILL_DIR / "examples" / "samplesheet.example.tsv"
        header = path.read_text().splitlines()[0]
        assert header.split("\t") == list(load("schemas").SAMPLESHEET_COLUMNS)

    def test_tests_directory_populated(self) -> None:
        tests = sorted((SKILL_DIR / "tests").glob("test_*.py"))
        assert len(tests) >= 1
