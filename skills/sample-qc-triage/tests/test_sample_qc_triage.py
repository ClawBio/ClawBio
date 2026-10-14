import importlib.util
import json
import subprocess
import sys
from pathlib import Path


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
