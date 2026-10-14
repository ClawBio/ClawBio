"""Profile-report reproducibility through the real CLI and generated replay."""

import hashlib
import json
import os
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import pytest

SKILL_DIR = Path(__file__).resolve().parents[1]
PROJECT_ROOT = SKILL_DIR.parents[1]


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _run_cli(args: list[str], tmp_path: Path) -> subprocess.CompletedProcess:
    env = dict(os.environ)
    env.pop("CLAWBIO_OTLP_ENDPOINT", None)
    env["CLAWBIO_AUDIT_LOG"] = str(tmp_path / "audit.jsonl")
    proc = subprocess.run(
        [sys.executable, str(SKILL_DIR / "profile_report.py"), *args],
        cwd=PROJECT_ROOT, env=env, capture_output=True, text=True, timeout=30, check=False,
    )
    assert proc.returncode == 0, proc.stderr
    return proc


def _assert_output_checksums(output_dir: Path) -> None:
    lines = (output_dir / "reproducibility" / "checksums.sha256").read_text().splitlines()
    labels = set()
    for line in lines:
        digest, label = line.split("  ", 1)
        target = (output_dir / label).resolve(strict=True)
        assert target.is_relative_to(output_dir.resolve())
        assert digest == _sha256(target)
        labels.add(label)
    assert labels == {
        "profile_report.md", "result.json", "reproducibility/commands.sh",
        "reproducibility/environment.yml", "reproducibility/inputs.json",
    }


def test_demo_records_its_input_and_complete_reproducibility_bundle(tmp_path):
    output_dir = tmp_path / "demo output"
    _run_cli(["--demo", "--output", str(output_dir)], tmp_path)

    repro = output_dir / "reproducibility"
    assert (repro / "commands.sh").is_file()
    assert (repro / "environment.yml").is_file()
    expected_digest = _sha256(SKILL_DIR / "demo_full_profile.json")
    inputs = json.loads((repro / "inputs.json").read_text())
    assert inputs == {
        "input_sha256": expected_digest,
        "input_kind": "demo-file",
        "checksum_kind": "file-bytes",
    }
    result = json.loads((output_dir / "result.json").read_text())
    assert result["input_checksum"] == f"sha256:{expected_digest}"
    commands = (repro / "commands.sh").read_text()
    assert "--demo" in commands and "--profile" not in commands
    assert "$OUTPUT_DIR" in commands and str(output_dir) not in commands
    environment = (repro / "environment.yml").read_text()
    assert f"python={sys.version_info.major}.{sys.version_info.minor}" in environment
    _assert_output_checksums(output_dir)


def test_profile_replays_safely_after_relocating_output(tmp_path):
    profile_path = tmp_path / "patient's profile $(touch REPLAY_MARKER).json"
    shutil.copyfile(SKILL_DIR / "tests" / "fixtures" / "mock_profile.json", profile_path)
    original_bytes = profile_path.read_bytes()
    output_dir = tmp_path / "original output"
    _run_cli(["--profile", str(profile_path), "--output", str(output_dir)], tmp_path)

    inputs = json.loads((output_dir / "reproducibility" / "inputs.json").read_text())
    assert inputs == {
        "input_sha256": _sha256(profile_path),
        "input_kind": "profile-file",
        "checksum_kind": "file-bytes",
    }
    _assert_output_checksums(output_dir)

    relocated = tmp_path / "relocated output $(touch OUTPUT_MARKER)"
    repro = relocated / "reproducibility"
    repro.mkdir(parents=True)
    shutil.copyfile(output_dir / "reproducibility" / "commands.sh", repro / "commands.sh")
    env = dict(os.environ)
    env.pop("CLAWBIO_OTLP_ENDPOINT", None)
    env.update(CLAWBIO_ROOT=str(PROJECT_ROOT), PYTHON=sys.executable,
               CLAWBIO_AUDIT_LOG=str(tmp_path / "replay-audit.jsonl"))
    replay = subprocess.run(
        ["bash", str(repro / "commands.sh")], cwd=tmp_path, env=env,
        capture_output=True, text=True, timeout=30, check=False,
    )
    assert replay.returncode == 0, replay.stderr
    assert not (tmp_path / "REPLAY_MARKER").exists()
    assert not (tmp_path / "OUTPUT_MARKER").exists()
    assert profile_path.read_bytes() == original_bytes
    original = json.loads((output_dir / "result.json").read_text())
    repeated = json.loads((relocated / "result.json").read_text())
    for key in ("summary", "data", "input_checksum"):
        assert repeated[key] == original[key]
    assert repeated["input_checksum"] == f"sha256:{_sha256(profile_path)}"
    assert (relocated / "profile_report.md").read_text() == (output_dir / "profile_report.md").read_text()
    _assert_output_checksums(relocated)


def test_repo_profile_uses_checkout_relative_replay_path(tmp_path):
    profile_path = SKILL_DIR / "tests" / "fixtures" / "mock_profile.json"
    output_dir = tmp_path / "repo profile"
    _run_cli(["--profile", str(profile_path), "--output", str(output_dir)], tmp_path)
    commands = (output_dir / "reproducibility" / "commands.sh").read_text()
    assert '"$CLAWBIO_ROOT"/skills/profile-report/tests/fixtures/mock_profile.json' in commands
    assert str(profile_path) not in commands
    _assert_output_checksums(output_dir)


@pytest.mark.parametrize("with_base_profile", [False, True])
def test_generated_demo_hashes_effective_profile_without_claiming_file_identity(
    tmp_path, monkeypatch, with_base_profile,
):
    sys.path.insert(0, str(SKILL_DIR))
    import profile_report as pr

    class FrozenDatetime(datetime):
        @classmethod
        def now(cls, tz=None):
            return cls(2026, 1, 1, tzinfo=timezone.utc)

    demo_root = tmp_path / "demo checkout"
    demo_root.mkdir()
    if with_base_profile:
        (demo_root / "profiles").mkdir()
        shutil.copyfile(SKILL_DIR / "tests" / "fixtures" / "mock_profile.json",
                        demo_root / "profiles" / "DEMO001.json")
    monkeypatch.setattr(pr, "_SCRIPT_DIR", demo_root / "skills" / "profile-report")
    monkeypatch.setattr(pr, "_PROJECT_ROOT", demo_root)
    monkeypatch.setattr(pr, "datetime", FrozenDatetime)
    profile = pr.build_demo_profile()
    canonical = json.dumps(profile, sort_keys=True, separators=(",", ":"), ensure_ascii=False)
    expected_digest = hashlib.sha256(canonical.encode("utf-8")).hexdigest()
    output_dir = tmp_path / "generated output"
    monkeypatch.setattr(sys, "argv", ["profile_report.py", "--demo", "--output", str(output_dir)])
    pr.main()

    inputs = json.loads((output_dir / "reproducibility" / "inputs.json").read_text())
    assert inputs == {
        "input_sha256": expected_digest,
        "input_kind": "generated-demo",
        "checksum_kind": "canonical-json",
    }
    result = json.loads((output_dir / "result.json").read_text())
    assert result["input_checksum"] == f"sha256:{expected_digest}"
    _assert_output_checksums(output_dir)
