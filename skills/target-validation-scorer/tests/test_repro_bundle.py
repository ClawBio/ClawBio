"""Reproducibility bundle tests for target-validation-scorer.

Every run writes <output_dir>/reproducibility/ (commands.sh, environment.yml,
checksums.sha256) via clawbio.common.reproducibility, and the manifest
verifies from the output directory (`sha256sum -c` contract).
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from target_validation_scorer import main

from clawbio.common.checksums import sha256_file


def read_checksums(out):
    lines = (out / "reproducibility" / "checksums.sha256").read_text().splitlines()
    return dict(reversed(line.split("  ", 1)) for line in lines)


def test_demo_writes_bundle_that_verifies(tmp_path):
    out = tmp_path / "out"
    main(["--demo", "--output", str(out)])
    repro = out / "reproducibility"
    for name in ("commands.sh", "environment.yml", "checksums.sha256"):
        assert (repro / name).is_file(), name

    entries = read_checksums(out)
    assert {"report.md", "validation_report.json", "figures/scoring_summary.png"} <= set(entries)
    for label, digest in entries.items():
        assert not Path(label).is_absolute()
        assert digest == sha256_file(out / label), label


def test_commands_sh_records_this_run(tmp_path):
    out = tmp_path / "out"
    main(["--demo", "--output", str(out)])
    content = (out / "reproducibility" / "commands.sh").read_text()
    assert '"$CLAWBIO_ROOT/skills/target-validation-scorer/target_validation_scorer.py"' in content
    assert "--demo" in content
    assert '"$OUTPUT_DIR"' in content
    assert str(out) not in content
