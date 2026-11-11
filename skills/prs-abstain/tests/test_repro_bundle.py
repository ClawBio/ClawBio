"""The demo run writes a reproducibility bundle that verifies from the output dir."""
from __future__ import annotations

import sys
from pathlib import Path

SKILL_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SKILL_DIR))
sys.path.insert(0, str(SKILL_DIR.parents[1]))

from prs_abstain import main

from clawbio.common.checksums import sha256_file


def test_demo_writes_verifiable_bundle(tmp_path):
    out = tmp_path / "out"
    assert main(["--demo", "--output", str(out)]) == 0
    repro = out / "reproducibility"
    for name in ("commands.sh", "environment.yml", "checksums.sha256"):
        assert (repro / name).is_file(), name

    commands = (repro / "commands.sh").read_text()
    assert "--demo" in commands and str(out) in commands

    lines = (repro / "checksums.sha256").read_text().splitlines()
    labels = {line.split("  ", 1)[1] for line in lines}
    assert {"result.json", "tables/decisions.csv", "reproducibility/commands.sh"} <= labels
    for line in lines:
        digest, label = line.split("  ", 1)
        assert not Path(label).is_absolute(), label
        assert sha256_file(out / label) == digest, label
