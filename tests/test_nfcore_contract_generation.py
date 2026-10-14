from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_scalar_contracts_match_pinned_sources_when_available():
    """The full upstream schema comparison is explicit, local and repeatable."""
    import os, pytest

    upstream = os.environ.get("CLAWBIO_NFCORE_UPSTREAM_DIR")
    if not upstream:
        pytest.skip("set CLAWBIO_NFCORE_UPSTREAM_DIR to downloaded pinned source trees")
    result = subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts/generate_nfcore_scalar_contracts.py"),
            "--upstream-dir",
            upstream,
            "--check",
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr


def test_bulk_wrapper_python_metadata_matches_package_floor():
    import json, yaml

    skill = ROOT / "skills/nfcore-rnaseq-wrapper"
    frontmatter = yaml.safe_load((skill / "SKILL.md").read_text().split("---", 2)[1])
    pin = json.loads((skill / "reproducibility/pinned_versions.json").read_text())
    assert frontmatter["metadata"]["dependencies"]["python"] == ">=3.11"
    assert pin["minimum_versions"]["python"] == "3.11"


def test_pinned_contracts_include_official_patterns_and_enums():
    import json

    for name in ("rnaseq", "scrnaseq", "sarek"):
        contract = json.loads(
            (
                ROOT
                / "skills"
                / f"nfcore-{name}-wrapper"
                / "reproducibility/upstream_scalar_contract.json"
            ).read_text()
        )
        assert "pattern" in contract["parameters"]["fasta"]
        assert "enum" in contract["parameters"]["aligner"]
