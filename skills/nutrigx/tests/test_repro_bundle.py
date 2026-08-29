"""
test_repro_bundle.py — Tests for the NutriGx reproducibility bundle.

The bundle must reuse the shared clawbio.common reproducibility layer
(write_commands_sh, write_environment_yml, write_checksums, sha256_file)
and write artefacts into <output_dir>/reproducibility/ following the
ClawBio-wide convention.
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))
_PROJECT_ROOT = Path(__file__).resolve().parents[3]
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

import repro_bundle
from clawbio.common.checksums import sha256_file as common_sha256_file

SKILL_DIR = Path(__file__).parent.parent
SYNTHETIC = Path(__file__).parent / "synthetic_patient.csv"
PANEL = SKILL_DIR / "data" / "snp_panel.json"


def make_bundle(tmp_path):
    output_dir = tmp_path / "out"
    output_dir.mkdir()
    (output_dir / "nutrigx_report.md").write_text("# report\n")
    repro_bundle.create_reproducibility_bundle(
        input_file=str(SYNTHETIC),
        output_dir=str(output_dir),
        panel_path=str(PANEL),
        args={"input": str(SYNTHETIC), "output": str(output_dir)},
    )
    return output_dir


def test_sha256_delegates_to_common_layer():
    assert repro_bundle.sha256_file is common_sha256_file


def test_commands_sh_in_reproducibility_dir(tmp_path):
    output_dir = make_bundle(tmp_path)
    cmd_path = output_dir / "reproducibility" / "commands.sh"
    assert cmd_path.exists()
    content = cmd_path.read_text()
    assert content.startswith("#!/usr/bin/env bash")
    assert "nutrigx.py" in content
    assert "--input" in content
    assert cmd_path.stat().st_mode & 0o111, "commands.sh must be executable"


def test_environment_yml_in_reproducibility_dir(tmp_path):
    output_dir = make_bundle(tmp_path)
    env_path = output_dir / "reproducibility" / "environment.yml"
    assert env_path.exists()
    content = env_path.read_text()
    assert "name: nutrigx" in content
    assert "python=3.11" in content
    assert "clawbio==0.1.0" in content


def test_checksums_are_outputs_only_and_resolve_from_output_dir(tmp_path):
    output_dir = make_bundle(tmp_path)
    checksum_path = output_dir / "reproducibility" / "checksums.sha256"
    assert checksum_path.exists()
    lines = checksum_path.read_text().strip().splitlines()
    assert lines, "manifest must not be empty"
    for line in lines:
        digest, label = line.split("  ", 1)
        target = output_dir / label
        assert target.exists(), f"label {label!r} must resolve from output_dir"
        assert common_sha256_file(target) == digest
    labels = [line.split("  ", 1)[1] for line in lines]
    assert "nutrigx_report.md" in labels
    # Inputs are attested in provenance.json, not the sha256sum manifest.
    assert SYNTHETIC.name not in labels
    assert PANEL.name not in labels


def test_commands_sh_verify_step_resolves(tmp_path):
    output_dir = make_bundle(tmp_path)
    content = (output_dir / "reproducibility" / "commands.sh").read_text()
    assert "sha256sum -c reproducibility/checksums.sha256" in content
    assert "sha256sum -c checksums.sha256" not in content


def test_missing_report_raises(tmp_path):
    output_dir = tmp_path / "out"
    output_dir.mkdir()  # no nutrigx_report.md on disk
    try:
        repro_bundle.create_reproducibility_bundle(
            input_file=str(SYNTHETIC),
            output_dir=str(output_dir),
            panel_path=str(PANEL),
            args={},
        )
    except FileNotFoundError as exc:
        assert "nutrigx_report.md" in str(exc)
    else:
        raise AssertionError("bundle built without the report it must attest")


def test_provenance_json(tmp_path):
    output_dir = make_bundle(tmp_path)
    prov_path = output_dir / "reproducibility" / "provenance.json"
    assert prov_path.exists()
    prov = json.loads(prov_path.read_text())
    assert prov["tool"] == "ClawBio NutriGx Advisor"
    assert prov["input_file"] == SYNTHETIC.name
    assert prov["input_sha256"] == common_sha256_file(SYNTHETIC)
    assert prov["panel_sha256"] == common_sha256_file(PANEL)
    assert "timestamp" in prov
    assert "args" in prov
    assert b"\r" not in prov_path.read_bytes(), "provenance.json must be LF-only"


def test_version_is_consistent():
    assert repro_bundle.VERSION == "0.2.0"
    skill_md = (SKILL_DIR / "SKILL.md").read_text()
    assert "version: 0.2.0" in skill_md
    nutrigx_src = (SKILL_DIR / "nutrigx.py").read_text()
    assert 'version="0.1.0"' not in nutrigx_src
    assert "version=VERSION" in nutrigx_src


def test_report_footer_points_at_reproducibility_dir():
    import generate_report as gr

    src = Path(gr.__file__).read_text()
    assert "`reproducibility/` subdirectory" in src
