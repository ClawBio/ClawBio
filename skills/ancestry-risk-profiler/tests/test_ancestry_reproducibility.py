"""Observable reproducibility contract for ancestry-risk-profiler CLI runs."""

import hashlib
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

SKILL_DIR = Path(__file__).resolve().parents[1]
PROJECT_ROOT = SKILL_DIR.parents[1]
DATA_DIR = SKILL_DIR / "data"
sys.path.insert(0, str(SKILL_DIR))

import ancestry_risk_profiler as arp  # noqa: E402


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _assert_output_checksums(output_dir: Path) -> set[str]:
    lines = (output_dir / "reproducibility" / "checksums.sha256").read_text().splitlines()
    assert lines
    labels = set()
    for line in lines:
        digest, label = line.split("  ", 1)
        target = (output_dir / label).resolve(strict=True)
        assert target.is_relative_to(output_dir.resolve())
        assert digest == _sha256(target)
        labels.add(label)
    return labels


def test_demo_writes_portable_bundle_with_source_and_output_hashes(tmp_path):
    output_dir = tmp_path / "demo output"
    arp.main(["--demo", "--ancestry", "SAS", "--output", str(output_dir)])

    repro = output_dir / "reproducibility"
    assert (repro / "commands.sh").is_file()
    assert (repro / "environment.yml").is_file()
    assert (repro / "inputs.json").is_file()
    source_hashes = json.loads((repro / "inputs.json").read_text())
    assert source_hashes == {
        "input_sha256": _sha256(DATA_DIR / "demo_patient_south_asian.txt"),
        "aisnp_panel_sha256": _sha256(DATA_DIR / "aisnp_panel.csv"),
        "associations_sha256": _sha256(DATA_DIR / "ancestry_risk_associations.json"),
    }

    commands = (repro / "commands.sh").read_text()
    assert "--demo" in commands
    assert "--ancestry" in commands and "SAS" in commands
    assert "--input" not in commands
    assert "$OUTPUT_DIR" in commands
    assert str(output_dir) not in commands
    environment = (repro / "environment.yml").read_text()
    assert f"python={sys.version_info.major}.{sys.version_info.minor}" in environment
    assert "matplotlib" in environment

    labels = _assert_output_checksums(output_dir)
    assert {
        "ancestry_risk_report.md", "ancestry_risk_result.json",
        "reproducibility/commands.sh", "reproducibility/environment.yml",
        "reproducibility/inputs.json",
    } <= labels
    chart = output_dir / "figures" / "aes_chart.png"
    assert ("figures/aes_chart.png" in labels) == chart.exists()


def test_repo_local_input_uses_the_portable_root_anchor(tmp_path):
    output_dir = tmp_path / "repo input output"
    arp.main(["--input", str(DATA_DIR / "demo_patient_south_asian.txt"),
              "--ancestry", "SAS", "--output", str(output_dir)])
    commands = (output_dir / "reproducibility" / "commands.sh").read_text()
    assert '"$CLAWBIO_ROOT"/skills/ancestry-risk-profiler/data/demo_patient_south_asian.txt' in commands
    assert "--demo" not in commands
    assert _assert_output_checksums(output_dir)


def test_input_replay_preserves_override_and_quotes_the_genotype_path(tmp_path):
    input_file = tmp_path / "sample $(touch REPLAY_MARKER) O'Brien.txt"
    shutil.copyfile(DATA_DIR / "demo_patient_south_asian.txt", input_file)
    output_dir = tmp_path / "first output"
    arp.main(["--input", str(input_file), "--ancestry", "eur", "--output", str(output_dir)])

    commands = output_dir / "reproducibility" / "commands.sh"
    assert commands.is_file()
    commands_text = commands.read_text()
    assert "--input" in commands_text
    assert "--demo" not in commands_text
    assert "--ancestry" in commands_text and "eur" in commands_text
    assert json.loads((output_dir / "reproducibility" / "inputs.json").read_text())[
        "input_sha256"
    ] == _sha256(input_file)

    replay_dir = tmp_path / "replay output"
    replay_script = replay_dir / "reproducibility" / "commands.sh"
    replay_script.parent.mkdir(parents=True)
    shutil.copyfile(commands, replay_script)
    env = dict(os.environ)
    env.pop("CLAWBIO_OTLP_ENDPOINT", None)
    env.update({
        "CLAWBIO_ROOT": str(PROJECT_ROOT),
        "PYTHON": sys.executable,
        "CLAWBIO_AUDIT_LOG": str(tmp_path / "audit.jsonl"),
    })
    proc = subprocess.run(
        ["bash", str(replay_script)], cwd=tmp_path, env=env,
        capture_output=True, text=True, timeout=45, check=False,
    )
    assert proc.returncode == 0, proc.stderr
    assert not (tmp_path / "REPLAY_MARKER").exists()
    result = json.loads((replay_dir / "ancestry_risk_result.json").read_text())
    assert result["inferred_ancestry"] == "EUR"
    assert result["confidence"] == "user-supplied"
    assert result["overridden"] is True
    assert _assert_output_checksums(replay_dir)


def test_empty_risk_run_does_not_checksum_a_stale_chart(tmp_path):
    input_file = tmp_path / "no-scored-variants.txt"
    input_file.write_text(
        "# rsid\tchromosome\tposition\tgenotype\nrs9999999\t1\t100\tAG\n",
        encoding="utf-8",
    )
    output_dir = tmp_path / "existing output"
    chart = output_dir / "figures" / "aes_chart.png"
    chart.parent.mkdir(parents=True)
    chart.write_bytes(b"chart from a previous run")

    arp.main(["--input", str(input_file), "--ancestry", "EUR", "--output", str(output_dir)])

    result = json.loads((output_dir / "ancestry_risk_result.json").read_text())
    assert result["risks"] == []
    labels = _assert_output_checksums(output_dir)
    assert "figures/aes_chart.png" not in labels


def test_unavailable_chart_dependency_does_not_checksum_a_stale_chart(tmp_path, monkeypatch):
    output_dir = tmp_path / "old output"
    chart = output_dir / "figures" / "aes_chart.png"
    chart.parent.mkdir(parents=True)
    chart.write_bytes(b"chart from a previous run")
    monkeypatch.setitem(sys.modules, "matplotlib", None)

    arp.main(["--demo", "--ancestry", "SAS", "--output", str(output_dir)])

    result = json.loads((output_dir / "ancestry_risk_result.json").read_text())
    assert result["risks"]
    labels = _assert_output_checksums(output_dir)
    assert "figures/aes_chart.png" not in labels


def test_insufficient_coverage_does_not_write_a_replay_bundle(tmp_path):
    output_dir = tmp_path / "abstained"
    with pytest.raises(SystemExit) as exc_info:
        arp.main(["--demo", "--output", str(output_dir)])
    assert exc_info.value.code == 1
    assert not (output_dir / "reproducibility").exists()
