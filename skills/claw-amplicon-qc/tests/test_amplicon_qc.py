"""Demo-mode output-contract test for claw-amplicon-qc.

Runs the skill in --demo mode (pure Python, no external tools) and asserts the
documented output contract plus the automatic-flag behaviour on the synthetic
low-quality sample.
"""
import json
import subprocess
import sys
from pathlib import Path

SKILL_DIR = Path(__file__).resolve().parent.parent
SCRIPT = SKILL_DIR / "amplicon_qc.py"


def test_demo_output_contract(tmp_path):
    out = tmp_path / "demo"
    res = subprocess.run(
        [sys.executable, str(SCRIPT), "--demo", "--output", str(out)],
        capture_output=True, text=True,
    )
    assert res.returncode == 0, res.stderr

    # Output contract (every file the SKILL.md tree promises for --demo)
    for rel in [
        "report.md", "qc_summary.json", "cutadapt.log",
        "stats/stats_raw.tsv", "stats/stats_trimmed.tsv",
        "reproducibility/commands.sh", "reproducibility/checksums.sha256",
    ]:
        assert (out / rel).exists(), f"missing {rel}"

    summary = json.loads((out / "qc_summary.json").read_text())
    assert summary["skill"] == "claw-amplicon-qc"
    assert summary["overall"]["n_samples"] == 3

    # The deliberately-degraded sample must be flagged on both per-sample checks.
    sample_c = next(s for s in summary["samples"] if s["sample_id"] == "sampleC")
    assert any("LOW_RETENTION" in f for f in sample_c["flags"])
    assert any("LOW_READ_COUNT" in f for f in sample_c["flags"])

    # A clean sample must not be flagged.
    sample_a = next(s for s in summary["samples"] if s["sample_id"] == "sampleA")
    assert sample_a["flags"] == []
