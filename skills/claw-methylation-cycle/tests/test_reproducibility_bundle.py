"""Demo run writes a reproducibility bundle whose checksums verify."""

import hashlib
import subprocess
import sys
from pathlib import Path

SKILL_DIR = Path(__file__).resolve().parent.parent
REPO_ROOT = SKILL_DIR.parent.parent


def test_demo_writes_verifiable_repro_bundle(tmp_path):
    out = tmp_path / "out"
    subprocess.run(
        [sys.executable, str(SKILL_DIR / "methylation_cycle.py"), "--demo", "--output", str(out)],
        check=True, cwd=REPO_ROOT,
    )
    repro = out / "reproducibility"
    for name in ("commands.sh", "environment.yml", "checksums.sha256"):
        assert (repro / name).is_file(), name

    commands = (repro / "commands.sh").read_text()
    assert "python skills/claw-methylation-cycle/methylation_cycle.py --demo" in commands
    assert str(SKILL_DIR) not in commands

    lines = (repro / "checksums.sha256").read_text().splitlines()
    labels = {line.split("  ", 1)[1] for line in lines}
    assert {"report.md", "result.json", "reproducibility/environment.yml"} <= labels
    for line in lines:
        digest, label = line.split("  ", 1)
        assert hashlib.sha256((out / label).read_bytes()).hexdigest() == digest, label
