"""Smoke tests for bgpt-mcp skill.

Validates that the SKILL.md is well-formed and that the skill metadata
is parseable without hitting the network.
"""

import os
import yaml
import json
import pytest

SKILL_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SKILL_MD = os.path.join(SKILL_DIR, "SKILL.md")


def _load_frontmatter(path):
    """Extract YAML frontmatter from a SKILL.md file."""
    with open(path, "r") as f:
        content = f.read()
    if not content.startswith("---"):
        raise ValueError("SKILL.md missing YAML frontmatter")
    parts = content.split("---", 2)
    if len(parts) < 3:
        raise ValueError("SKILL.md frontmatter not properly closed")
    return yaml.safe_load(parts[1])


class TestSkillMetadata:
    """Verify SKILL.md frontmatter is complete and valid."""

    @pytest.fixture(autouse=True)
    def load_meta(self):
        self.meta = _load_frontmatter(SKILL_MD)

    def test_frontmatter_loads(self):
        assert self.meta is not None

    def test_required_top_level_fields(self):
        for field in ("name", "description", "version", "author", "domain", "license"):
            assert field in self.meta, f"Missing required field: {field}"

    def test_name_matches_directory(self):
        assert self.meta["name"] == "bgpt-mcp"

    def test_openclaw_block_exists(self):
        assert "metadata" in self.meta
        assert "openclaw" in self.meta["metadata"]

    def test_trigger_keywords_present(self):
        keywords = self.meta["metadata"]["openclaw"]["trigger_keywords"]
        assert len(keywords) >= 3, "Need at least 3 trigger keywords"

    def test_os_values_are_node_platform(self):
        valid = {"darwin", "linux", "win32"}
        for os_val in self.meta["metadata"]["openclaw"]["os"]:
            assert os_val in valid, f"Invalid os value: {os_val}"

    def test_emoji_present(self):
        assert self.meta["metadata"]["openclaw"]["emoji"]

    def test_homepage_is_url(self):
        hp = self.meta["metadata"]["openclaw"]["homepage"]
        assert hp.startswith("https://"), f"Homepage should be HTTPS URL: {hp}"

    def test_inputs_defined(self):
        assert "inputs" in self.meta
        assert len(self.meta["inputs"]) >= 1

    def test_query_input_required(self):
        query_inputs = [i for i in self.meta["inputs"] if i["name"] == "query"]
        assert len(query_inputs) == 1
        assert query_inputs[0]["required"] is True

    def test_endpoints_defined(self):
        assert "endpoints" in self.meta
        eps = self.meta["endpoints"]
        assert "mcp_sse" in eps or "mcp_stream" in eps


class TestSkillSections:
    """Verify SKILL.md contains all required sections."""

    @pytest.fixture(autouse=True)
    def load_content(self):
        with open(SKILL_MD, "r") as f:
            self.content = f.read()

    @pytest.mark.parametrize("section", [
        "## Trigger",
        "## Why This Exists",
        "## Core Capabilities",
        "## Workflow",
        "## CLI Reference",
        "## Demo",
        "## Output Structure",
        "## Dependencies",
        "## Gotchas",
        "## Safety",
        "## Citations",
    ])
    def test_required_section_exists(self, section):
        assert section in self.content, f"Missing required section: {section}"

    def test_cli_reference_has_demo_flag(self):
        cli_start = self.content.index("## CLI Reference")
        next_section = self.content.index("\n## ", cli_start + 1)
        cli_section = self.content[cli_start:next_section]
        assert "--demo" in cli_section, "CLI Reference must include --demo example"

    def test_cli_reference_has_clawbio_runner(self):
        cli_start = self.content.index("## CLI Reference")
        next_section = self.content.index("\n## ", cli_start + 1)
        cli_section = self.content[cli_start:next_section]
        assert "clawbio.py run" in cli_section, "CLI Reference must include clawbio.py runner"


class TestDemoMode:
    """Verify demo can be described without network access."""

    @pytest.fixture(autouse=True)
    def load_meta(self):
        self.meta = _load_frontmatter(SKILL_MD)

    def test_demo_query_is_deterministic(self):
        """The demo section should reference a fixed example query."""
        with open(SKILL_MD, "r") as f:
            content = f.read()
        demo_start = content.index("## Demo")
        next_section = content.index("\n## ", demo_start + 1)
        demo_section = content[demo_start:next_section]
        assert "CAR-T" in demo_section or "CRISPR" in demo_section, \
            "Demo section should contain a concrete example query"

    def test_mcp_endpoints_are_https(self):
        """MCP endpoints must use HTTPS."""
        for key, url in self.meta["endpoints"].items():
            assert url.startswith("https://"), f"Endpoint {key} must be HTTPS: {url}"

    def test_version_is_semver(self):
        """Version should follow semver format."""
        version = self.meta["version"]
        parts = str(version).split(".")
        assert len(parts) >= 2, f"Version should be semver: {version}"
