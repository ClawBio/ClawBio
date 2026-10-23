"""Skill metadata that must match the pinned nf-core releases.

Set CLAWBIO_NFCORE_UPSTREAM_DIR to the downloaded release trees to run the
comparisons that need them.
"""

from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_bulk_wrapper_python_metadata_matches_package_floor():
    import json, yaml

    skill = ROOT / "skills/nfcore-rnaseq-wrapper"
    frontmatter = yaml.safe_load((skill / "SKILL.md").read_text().split("---", 2)[1])
    pin = json.loads((skill / "reproducibility/pinned_versions.json").read_text())
    assert frontmatter["metadata"]["dependencies"]["python"] == ">=3.11"
    assert pin["minimum_versions"]["python"] == "3.11"


def test_demo_resource_limits_match_the_pinned_test_profiles_when_available():
    """The demo cap must equal each pinned release's own -profile test resourceLimits."""
    import importlib.util
    import os
    import re

    import pytest

    upstream = os.environ.get("CLAWBIO_NFCORE_UPSTREAM_DIR")
    if not upstream:
        pytest.skip("set CLAWBIO_NFCORE_UPSTREAM_DIR to downloaded pinned source trees")
    for name in ("rnaseq", "sarek", "scrnaseq"):
        skill = ROOT / "skills" / f"nfcore-{name}-wrapper"
        spec = importlib.util.spec_from_file_location(f"_{name}_schemas", skill / "schemas.py")
        schemas = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(schemas)
        source = Path(upstream) / f"{name}-{schemas.DEFAULT_PIPELINE_VERSION}"
        test_config = (source / "conf" / "test.config").read_text()
        if "resourceLimits" not in test_config:
            # rnaseq declares them in the `test` profile block of nextflow.config.
            text = (source / "nextflow.config").read_text()
            test_config = text[text.index("    test {"):]
        found = re.search(r"resourceLimits\s*=\s*\[\s*cpus:\s*(\d+),\s*memory:\s*'(\d+)\.GB'", test_config)
        assert found, name
        assert schemas.DEMO_RESOURCE_LIMITS == {"cpus": int(found.group(1)), "memory_gb": int(found.group(2))}, name
