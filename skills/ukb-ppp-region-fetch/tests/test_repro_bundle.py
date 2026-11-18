"""Reproducibility bundle: every run writes <output>/reproducibility/ via
clawbio.common.reproducibility, and checksums.sha256 verifies from <output>."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
_PROJECT_ROOT = Path(__file__).resolve().parents[3]
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from ukb_ppp_region_fetch import main

from clawbio.common.checksums import sha256_file


def test_demo_writes_verifiable_bundle(tmp_path):
    out = tmp_path / "out"
    assert main(["--demo", "--output", str(out)]) == 0

    repro = out / "reproducibility"
    for name in ("commands.sh", "environment.yml", "checksums.sha256"):
        assert (repro / name).is_file(), name

    lines = (repro / "checksums.sha256").read_text().splitlines()
    labels = {line.split("  ", 1)[1] for line in lines}
    assert {"variants.tsv", "report.md"} <= labels
    for line in lines:
        digest, label = line.split("  ", 1)
        assert not Path(label).is_absolute()
        assert digest == sha256_file(out / label), label

    commands = (repro / "commands.sh").read_text()
    assert '"$CLAWBIO_ROOT/skills/ukb-ppp-region-fetch/ukb_ppp_region_fetch.py"' in commands
    assert "--demo" in commands and '"$OUTPUT_DIR"' in commands
    assert str(out) not in commands
