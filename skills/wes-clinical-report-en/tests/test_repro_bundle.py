"""The demo writes a reproducibility bundle that verifies from the output dir."""

import subprocess
import sys
from pathlib import Path

import pytest

from clawbio.common.checksums import sha256_file

pytest.importorskip("reportlab")

SCRIPT = Path(__file__).resolve().parents[1] / "wes_clinical_report_en.py"


def test_demo_writes_verifiable_bundle(tmp_path):
    out = tmp_path / "out"
    subprocess.run(
        [sys.executable, str(SCRIPT), "--demo", "--output", str(out)],
        check=True, cwd=tmp_path,
    )
    repro = out / "reproducibility"
    for name in ("commands.sh", "environment.yml", "checksums.sha256"):
        assert (repro / name).is_file(), f"missing {name}"

    commands = (repro / "commands.sh").read_text()
    assert "--demo" in commands
    assert '"$CLAWBIO_ROOT/skills/wes-clinical-report-en/' in commands
    assert '--output \\\n  "$OUTPUT_DIR"' in commands
    assert str(tmp_path) not in commands, "commands.sh must not embed this run's paths"

    lines = (repro / "checksums.sha256").read_text().splitlines()
    labels = {line.split("  ", 1)[1] for line in lines}
    assert {"reproducibility/commands.sh", "reproducibility/environment.yml"} <= labels
    for line in lines:
        digest, label = line.split("  ", 1)
        assert not Path(label).is_absolute()
        assert digest == sha256_file(out / label), label
