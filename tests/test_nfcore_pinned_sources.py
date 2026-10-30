"""Skill metadata that must match the pinned nf-core releases.

The comparisons read exact copies of the upstream files they need, vendored in
tests/fixtures/nfcore_pinned/ with their source commit and sha256 in
SOURCES.json, so they run in CI. Set CLAWBIO_NFCORE_UPSTREAM_DIR to downloaded
release trees (<pipeline>-<version>/) to read those instead and to check the
vendored copies for drift.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import os
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / "tests" / "fixtures" / "nfcore_pinned"
SOURCES = json.loads((FIXTURES / "SOURCES.json").read_text(encoding="utf-8"))
PIPELINES = ("rnaseq", "sarek", "scrnaseq")


def _schemas(name: str):
    skill = ROOT / "skills" / f"nfcore-{name}-wrapper"
    spec = importlib.util.spec_from_file_location(f"_{name}_pinned_schemas", skill / "schemas.py")
    schemas = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(schemas)
    return schemas


def _pinned_tree(name: str) -> Path:
    """The release a wrapper pins: the downloaded tree when CLAWBIO_NFCORE_UPSTREAM_DIR is set, else the vendored copy."""
    release = f"{name}-{_schemas(name).DEFAULT_PIPELINE_VERSION}"
    upstream = os.environ.get("CLAWBIO_NFCORE_UPSTREAM_DIR")
    return Path(upstream) / release if upstream else FIXTURES / release


def test_bulk_wrapper_python_metadata_matches_package_floor():
    import yaml

    skill = ROOT / "skills/nfcore-rnaseq-wrapper"
    frontmatter = yaml.safe_load((skill / "SKILL.md").read_text().split("---", 2)[1])
    pin = json.loads((skill / "reproducibility/pinned_versions.json").read_text())
    assert frontmatter["metadata"]["dependencies"]["python"] == ">=3.11"
    assert pin["minimum_versions"]["python"] == "3.11"


def test_vendored_releases_are_the_pinned_versions():
    """A wrapper that pins a new release needs its files re-vendored, or the checks compare the old one."""
    pinned = {f"{name}-{_schemas(name).DEFAULT_PIPELINE_VERSION}" for name in PIPELINES}
    assert set(SOURCES["releases"]) == pinned


def test_vendored_files_match_their_recorded_sha256():
    on_disk = {
        path.relative_to(FIXTURES).as_posix()
        for path in FIXTURES.rglob("*") if path.is_file() and path.name != "SOURCES.json"
    }
    recorded = {f"{release}/{relative}" for release, source in SOURCES["releases"].items() for relative in source["files"]}
    assert on_disk == recorded
    for release, source in SOURCES["releases"].items():
        for relative, digest in source["files"].items():
            assert hashlib.sha256((FIXTURES / release / relative).read_bytes()).hexdigest() == digest, f"{release}/{relative}"


def test_vendored_files_match_the_release_trees_when_available():
    upstream = os.environ.get("CLAWBIO_NFCORE_UPSTREAM_DIR")
    if not upstream:
        pytest.skip("set CLAWBIO_NFCORE_UPSTREAM_DIR to downloaded pinned source trees")
    for release, source in SOURCES["releases"].items():
        for relative in source["files"]:
            vendored = (FIXTURES / release / relative).read_bytes()
            assert vendored == (Path(upstream) / release / relative).read_bytes(), f"{release}/{relative}"


def test_demo_resource_limits_match_the_pinned_test_profiles():
    """The demo cap must equal each pinned release's own -profile test resourceLimits."""
    for name in PIPELINES:
        source = _pinned_tree(name)
        test_config = (source / "conf" / "test.config").read_text()
        if "resourceLimits" not in test_config:
            # rnaseq declares them in the `test` profile block of nextflow.config.
            text = (source / "nextflow.config").read_text()
            test_config = text[text.index("    test {"):]
        found = re.search(r"resourceLimits\s*=\s*\[\s*cpus:\s*(\d+),\s*memory:\s*'(\d+)\.GB'", test_config)
        assert found, name
        assert _schemas(name).DEMO_RESOURCE_LIMITS == {"cpus": int(found.group(1)), "memory_gb": int(found.group(2))}, name


def test_recorded_minimum_nextflow_matches_each_release_manifest():
    """The wrappers leave the minimum Nextflow to the release's own
    manifest.nextflowVersion (`!` makes Nextflow refuse an older engine); the
    minimum each skill records in pinned_versions.json must be that same one."""
    for name in PIPELINES:
        skill = ROOT / "skills" / f"nfcore-{name}-wrapper"
        config = (_pinned_tree(name) / "nextflow.config").read_text()
        gate = re.search(r"nextflowVersion\s*=\s*'!>=([\d.]+)'", config)
        assert gate, f"{name}: manifest.nextflowVersion is not a strict '!>=' gate"
        manifest = tuple(int(part) for part in gate.group(1).split("."))
        pin = json.loads((skill / "reproducibility/pinned_versions.json").read_text())
        recorded = pin.get("minimum_versions", {}).get("nextflow") or pin["runtime"]["nextflow_min_version"]
        recorded = tuple(int(part) for part in (recorded.split(".") if isinstance(recorded, str) else recorded))
        assert recorded == manifest, name
