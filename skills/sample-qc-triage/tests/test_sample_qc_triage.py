import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import pytest
from inline_snapshot import snapshot

from clawbio.common import demo_runner

SKILL_DIR = Path(__file__).resolve().parents[1]
MODULE_PATH = SKILL_DIR / "sample_qc_triage.py"
DISCLAIMER = (
    "ClawBio is a research and educational tool. It is not a medical device "
    "and does not provide clinical diagnoses. Consult a healthcare professional "
    "before making any medical decisions."
)


def load_module():
    spec = importlib.util.spec_from_file_location("sample_qc_triage", MODULE_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_parse_metrics_rejects_missing_columns(tmp_path):
    module = load_module()
    bad = tmp_path / "bad.csv"
    bad.write_text("sample_id,total_reads\nS1,100\n", encoding="utf-8")
    try:
        module.load_metrics(bad)
    except ValueError as exc:
        assert "missing required columns" in str(exc).lower()
    else:
        raise AssertionError("missing columns should fail")


def test_empty_input_handled(tmp_path):
    module = load_module()
    empty = tmp_path / "empty.csv"
    empty.write_text(
        "sample_id,batch,total_reads,mapped_pct,duplicate_pct,mitochondrial_pct,contamination_pct,complexity_score\n",
        encoding="utf-8",
    )
    try:
        module.load_metrics(empty)
    except ValueError as exc:
        assert "no samples" in str(exc).lower()
    else:
        raise AssertionError("empty input should fail")


def test_cli_rejects_malformed_input_without_traceback(tmp_path):
    bad = tmp_path / "bad.csv"
    bad.write_text("sample_id,total_reads\nS1,100\n", encoding="utf-8")
    completed = subprocess.run(
        [sys.executable, str(MODULE_PATH), "--input", str(bad), "--output", str(tmp_path / "out")],
        text=True,
        capture_output=True,
        check=False,
    )
    assert completed.returncode == 2
    assert "ERROR:" in completed.stderr
    assert "Traceback" not in completed.stderr


def test_demo_analysis_flags_qc_outliers():
    module = load_module()
    records = module.load_metrics(SKILL_DIR / "demo_qc_metrics.csv")
    result = module.analyze_records(records)
    assert result["summary"]["sample_count"] == 5
    assert result["summary"]["flagged_count"] == 3
    flagged = {row["sample_id"]: row for row in result["samples"] if row["status"] == "flagged"}
    assert flagged["CB_QC_003"]["dominant_issue"] == "low_complexity"
    assert flagged["CB_QC_004"]["dominant_issue"] == "contamination"
    assert flagged["CB_QC_005"]["dominant_issue"] == "sex_mismatch"
    assert "identity_mismatch" in flagged["CB_QC_005"]["issues"]


def test_demo_cli_writes_report_json_tables_and_reproducibility(tmp_path):
    out = tmp_path / "qc_out"
    completed = subprocess.run(
        [sys.executable, str(MODULE_PATH), "--demo", "--output", str(out)],
        text=True,
        capture_output=True,
        check=True,
    )
    assert "Sample QC Triage" in completed.stdout
    report = (out / "report.md").read_text(encoding="utf-8")
    assert DISCLAIMER in report
    assert "CB_QC_003" in report
    assert "F->M" in report
    assert "Synthetic demo data" in report
    result = json.loads((out / "result.json").read_text(encoding="utf-8"))
    assert result["skill"] == "sample-qc-triage"
    assert result["summary"]["flagged_count"] == 3
    assert (out / "tables" / "sample_flags.csv").exists()
    assert (out / "reproducibility" / "commands.sh").exists()


def test_demo_reproducibility_bundle_is_complete(tmp_path):
    out = tmp_path / "qc_out"
    subprocess.run(
        [sys.executable, str(MODULE_PATH), "--demo", "--output", str(out)],
        text=True,
        capture_output=True,
        check=True,
    )
    repro = out / "reproducibility"

    commands = repro / "commands.sh"
    assert commands.exists()
    commands_text = commands.read_text(encoding="utf-8")
    assert "CLAWBIO_ROOT" in commands_text
    assert "$OUTPUT_DIR" in commands_text
    assert str(out) not in commands_text

    environment = (repro / "environment.yml").read_text(encoding="utf-8")
    assert "name: clawbio-sample-qc-triage" in environment

    checksum_lines = [
        line
        for line in (repro / "checksums.sha256").read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    assert checksum_lines
    labels = set()
    for line in checksum_lines:
        digest, label = line.split("  ", 1)
        assert len(digest) == 64
        labels.add(label)
        assert (out / label).exists()
    assert {"report.md", "result.json", "tables/sample_flags.csv"} <= labels


def test_demo_outputs_are_byte_stable_across_hash_seeds(tmp_path):
    """Checksums are only meaningful if two identical runs produce identical bytes."""
    import os

    digests = []
    for seed in ("0", "12345"):
        out = tmp_path / f"run_{seed}"
        env = dict(os.environ, PYTHONHASHSEED=seed)
        subprocess.run(
            [sys.executable, str(MODULE_PATH), "--demo", "--output", str(out)],
            text=True,
            capture_output=True,
            check=True,
            env=env,
        )
        digests.append((out / "reproducibility" / "checksums.sha256").read_text(encoding="utf-8"))
    assert digests[0] == digests[1]


def test_bundle_replays_a_non_demo_run_from_a_path_with_spaces(tmp_path):
    """Every other bundle test runs --demo. A real input, a directory with a
    space in it, and an actual replay are what the bundle claims to support."""
    import shutil

    workspace = tmp_path / "my qc runs"
    workspace.mkdir()
    source = workspace / "metrics copy.csv"
    shutil.copy(SKILL_DIR / "demo_qc_metrics.csv", source)
    out = workspace / "first run"

    subprocess.run(
        [sys.executable, str(MODULE_PATH), "--input", str(source), "--output", str(out)],
        text=True, capture_output=True, check=True,
    )
    first = (out / "reproducibility" / "checksums.sha256").read_text(encoding="utf-8")

    replay = workspace / "replayed"
    shutil.copytree(out, replay)
    done = subprocess.run(
        ["bash", str(replay / "reproducibility" / "commands.sh")],
        text=True, capture_output=True,
    )
    assert done.returncode == 0, done.stderr
    assert (replay / "reproducibility" / "checksums.sha256").read_text(encoding="utf-8") == first


@pytest.mark.parametrize("report", ["report.md", "result.json", "tables/sample_flags.csv"])
def test_reports_do_not_record_the_input_directory(tmp_path, report):
    """Reports are shared; an absolute input path leaks the user's name and folder layout."""
    import shutil

    private_dir = tmp_path / "private_user_dir"
    private_dir.mkdir()
    source = private_dir / "metrics.csv"
    shutil.copy(SKILL_DIR / "demo_qc_metrics.csv", source)
    out = tmp_path / "out"
    subprocess.run(
        [sys.executable, str(MODULE_PATH), "--input", str(source), "--output", str(out)],
        text=True, capture_output=True, check=True,
    )
    text = (out / report).read_text(encoding="utf-8")
    assert str(private_dir) not in text and str(private_dir.resolve()) not in text


# ── Demo output snapshot ──────────────────────────────────────────────────────
# Fails when a code change alters any demo number. If the change is intended,
# run `pytest skills/sample-qc-triage/tests --inline-snapshot=fix`; it rewrites the
# snapshot() values below, and the PR diff shows which numbers moved.

@pytest.mark.demo_snapshot
def test_demo_output_matches_snapshot(tmp_path):
    assert demo_runner.demo_outputs("sample-qc-triage", tmp_path) == snapshot({'report.md': """\
# Sample QC Triage Report

**Input**: `demo_qc_metrics.csv`
**Mode**: Synthetic demo data
**Samples**: 5
**Flagged samples**: 3

| Sample | Batch | Status | Dominant issue | Sex check | Fingerprint % | Mapped % | Duplicate % | Contamination % |
|---|---|---|---|---|---:|---:|---:|---:|
| CB_QC_001 | A | pass | none | F->F | 99.8 | 96.2 | 12.0 | 0.3 |
| CB_QC_002 | A | pass | none | M->M | 99.6 | 95.8 | 13.5 | 0.5 |
| CB_QC_003 | A | flagged | low_complexity | F->F | 98.9 | 88.5 | 44.0 | 0.8 |
| CB_QC_004 | B | flagged | contamination | M->M | 98.2 | 91.0 | 16.0 | 8.5 |
| CB_QC_005 | C | flagged | sex_mismatch | F->M | 91.4 | 72.0 | 28.0 | 2.1 |

## Interpretation

Flagged samples are candidates for local review. Identity and sex mismatch flags are operational QC signals, not forensic proof.

ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions.
""", 'result.json': """\
{
  "skill": "sample-qc-triage",
  "summary": {
    "sample_count": 5,
    "flagged_count": 3,
    "median_reads": 48000000.0,
    "issue_counts": {
      "low_complexity": 1,
      "batch_shift": 2,
      "contamination": 1,
      "sex_mismatch": 1,
      "identity_mismatch": 1
    }
  },
  "samples": [
    {
      "sample_id": "CB_QC_001",
      "batch": "A",
      "complexity_score": 0.91,
      "contamination_pct": 0.3,
      "duplicate_pct": 12.0,
      "mapped_pct": 96.2,
      "mitochondrial_pct": 4.1,
      "total_reads": 50000000.0,
      "fingerprint_match_pct": 99.8,
      "expected_sex": "F",
      "observed_sex": "F",
      "read_depth_delta_pct": 4.17,
      "issues": [],
      "dominant_issue": "none",
      "status": "pass"
    },
    {
      "sample_id": "CB_QC_002",
      "batch": "A",
      "complexity_score": 0.89,
      "contamination_pct": 0.5,
      "duplicate_pct": 13.5,
      "mapped_pct": 95.8,
      "mitochondrial_pct": 4.8,
      "total_reads": 48000000.0,
      "fingerprint_match_pct": 99.6,
      "expected_sex": "M",
      "observed_sex": "M",
      "read_depth_delta_pct": 0.0,
      "issues": [],
      "dominant_issue": "none",
      "status": "pass"
    },
    {
      "sample_id": "CB_QC_003",
      "batch": "A",
      "complexity_score": 0.42,
      "contamination_pct": 0.8,
      "duplicate_pct": 44.0,
      "mapped_pct": 88.5,
      "mitochondrial_pct": 13.2,
      "total_reads": 22000000.0,
      "fingerprint_match_pct": 98.9,
      "expected_sex": "F",
      "observed_sex": "F",
      "read_depth_delta_pct": -54.17,
      "issues": [
        "low_complexity",
        "batch_shift"
      ],
      "dominant_issue": "low_complexity",
      "status": "flagged"
    },
    {
      "sample_id": "CB_QC_004",
      "batch": "B",
      "complexity_score": 0.82,
      "contamination_pct": 8.5,
      "duplicate_pct": 16.0,
      "mapped_pct": 91.0,
      "mitochondrial_pct": 7.4,
      "total_reads": 51000000.0,
      "fingerprint_match_pct": 98.2,
      "expected_sex": "M",
      "observed_sex": "M",
      "read_depth_delta_pct": 6.25,
      "issues": [
        "contamination"
      ],
      "dominant_issue": "contamination",
      "status": "flagged"
    },
    {
      "sample_id": "CB_QC_005",
      "batch": "C",
      "complexity_score": 0.65,
      "contamination_pct": 2.1,
      "duplicate_pct": 28.0,
      "mapped_pct": 72.0,
      "mitochondrial_pct": 18.5,
      "total_reads": 15000000.0,
      "fingerprint_match_pct": 91.4,
      "expected_sex": "F",
      "observed_sex": "M",
      "read_depth_delta_pct": -68.75,
      "issues": [
        "sex_mismatch",
        "identity_mismatch",
        "batch_shift"
      ],
      "dominant_issue": "sex_mismatch",
      "status": "flagged"
    }
  ],
  "disclaimer": "ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions."
}
""", 'tables/sample_flags.csv': """\
sample_id,batch,status,dominant_issue,expected_sex,observed_sex,fingerprint_match_pct,total_reads,mapped_pct,duplicate_pct,mitochondrial_pct,contamination_pct,complexity_score,read_depth_delta_pct
CB_QC_001,A,pass,none,F,F,99.8,50000000.0,96.2,12.0,4.1,0.3,0.91,4.17
CB_QC_002,A,pass,none,M,M,99.6,48000000.0,95.8,13.5,4.8,0.5,0.89,0.0
CB_QC_003,A,flagged,low_complexity,F,F,98.9,22000000.0,88.5,44.0,13.2,0.8,0.42,-54.17
CB_QC_004,B,flagged,contamination,M,M,98.2,51000000.0,91.0,16.0,7.4,8.5,0.82,6.25
CB_QC_005,C,flagged,sex_mismatch,F,M,91.4,15000000.0,72.0,28.0,18.5,2.1,0.65,-68.75
"""})
