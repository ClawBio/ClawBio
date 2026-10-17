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


# nf-schema skips each pipeline's validation.defaultIgnoreParams, read from the
# pinned nextflow.config; the contracts must skip the same keys.
_PINNED_IGNORED = {
    "rnaseq": ["genomes", "igenomes_base"],
    "sarek": ["genomes", "freebayes_filter", "vep_cache_version"],
    "scrnaseq": ["genomes"],
}


def test_pinned_contracts_record_upstream_ignored_parameters():
    import json

    for name, ignored in _PINNED_IGNORED.items():
        contract = json.loads(
            (ROOT / "skills" / f"nfcore-{name}-wrapper" / "reproducibility/upstream_scalar_contract.json").read_text()
        )
        assert contract["ignored_parameters"] == ignored, name
        assert len(contract["config_sha256"]) == 64, name


def _generator():
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "generate_nfcore_scalar_contracts", ROOT / "scripts/generate_nfcore_scalar_contracts.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_generator_reads_default_ignore_params_from_nextflow_config(tmp_path):
    import hashlib
    import json

    source = tmp_path / "sarek-3.8.1"
    source.mkdir()
    (source / "nextflow_schema.json").write_text(json.dumps(
        {"$defs": {"g": {"properties": {"freebayes_filter": {"type": "string", "default": "x"}}}}}
    ))
    config = (
        "validation {\n"
        "    defaultIgnoreParams = [\n        \"genomes\", 'freebayes_filter',\n        \"vep_cache_version\" ]\n"
        "    monochromeLogs = params.monochrome_logs\n}\n"
    )
    (source / "nextflow.config").write_text(config)
    contract = json.loads(_generator().generate_contract("sarek", "3.8.1", source))
    assert contract["ignored_parameters"] == ["genomes", "freebayes_filter", "vep_cache_version"]
    assert contract["config_sha256"] == hashlib.sha256(config.encode()).hexdigest()
    assert contract["config_url"].endswith("/nf-core/sarek/3.8.1/nextflow.config")
    assert contract["parameters"] == {"freebayes_filter": {"type": "string"}}


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
