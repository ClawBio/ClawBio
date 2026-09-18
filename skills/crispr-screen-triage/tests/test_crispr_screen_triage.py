import importlib.util
import json
import subprocess
import sys
from pathlib import Path


SKILL_DIR = Path(__file__).resolve().parents[1]
MODULE_PATH = SKILL_DIR / "crispr_screen_triage.py"
DISCLAIMER = (
    "ClawBio is a research and educational tool. It is not a medical device "
    "and does not provide clinical diagnoses. Consult a healthcare professional "
    "before making any medical decisions."
)


def load_module():
    spec = importlib.util.spec_from_file_location("crispr_screen_triage", MODULE_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_load_counts_rejects_non_numeric_counts(tmp_path):
    module = load_module()
    bad = tmp_path / "bad.csv"
    bad.write_text("guide_id,gene,control_count,treatment_count,essentiality,druggability\ng1,TP53,10,nope,0.5,0.5\n", encoding="utf-8")
    try:
        module.load_counts(bad)
    except ValueError as exc:
        assert "numeric" in str(exc).lower()
    else:
        raise AssertionError("non-numeric counts should fail")


def test_empty_input_handled(tmp_path):
    module = load_module()
    empty = tmp_path / "empty.csv"
    empty.write_text("guide_id,gene,control_count,treatment_count,essentiality,druggability\n", encoding="utf-8")
    try:
        module.load_counts(empty)
    except ValueError as exc:
        assert "no guides" in str(exc).lower()
    else:
        raise AssertionError("empty input should fail")


def test_cli_rejects_malformed_input_without_traceback(tmp_path):
    bad = tmp_path / "bad.csv"
    bad.write_text(
        "guide_id,gene,control_count,treatment_count,essentiality,druggability\ng1,TP53,10,nope,0.5,0.5\n",
        encoding="utf-8",
    )
    completed = subprocess.run(
        [sys.executable, str(MODULE_PATH), "--input", str(bad), "--output", str(tmp_path / "out")],
        text=True,
        capture_output=True,
        check=False,
    )
    assert completed.returncode == 2
    assert "ERROR:" in completed.stderr
    assert "Traceback" not in completed.stderr


def test_triage_ranks_depleted_druggable_genes():
    module = load_module()
    rows = module.load_counts(SKILL_DIR / "demo_screen_counts.csv")
    result = module.triage_genes(rows)
    assert result["summary"]["gene_count"] == 6
    assert result["summary"]["guide_count"] == 12
    assert result["top_hits"][0]["gene"] == "BRCA1"
    assert result["top_hits"][0]["priority"] == "high"
    assert result["genes"][0]["median_log2_fold_change"] < -1.0
    assert result["genes"][0]["guide_count"] == 2
    assert {hit["gene"] for hit in result["top_hits"][:3]} == {"BRCA1", "PARP1", "ATR"}


def test_demo_cli_writes_expected_outputs(tmp_path):
    out = tmp_path / "crispr_out"
    completed = subprocess.run(
        [sys.executable, str(MODULE_PATH), "--demo", "--output", str(out)],
        text=True,
        capture_output=True,
        check=True,
    )
    assert "CRISPR Screen Triage" in completed.stdout
    report = (out / "report.md").read_text(encoding="utf-8")
    assert DISCLAIMER in report
    assert "BRCA1" in report
    assert "Synthetic demo data" in report
    result = json.loads((out / "result.json").read_text(encoding="utf-8"))
    assert result["skill"] == "crispr-screen-triage"
    assert result["top_hits"][0]["gene"] == "BRCA1"
    assert result["summary"]["guide_count"] == 12
    assert (out / "tables" / "triaged_genes.csv").exists()
    assert (out / "tables" / "guide_metrics.csv").exists()
    assert (out / "reproducibility" / "commands.sh").exists()


def test_demo_reproducibility_bundle_is_complete(tmp_path):
    out = tmp_path / "repro_out"
    subprocess.run(
        [sys.executable, str(MODULE_PATH), "--demo", "--output", str(out)],
        text=True,
        capture_output=True,
        check=True,
    )
    repro = out / "reproducibility"

    commands_text = (repro / "commands.sh").read_text(encoding="utf-8")
    assert "CLAWBIO_ROOT" in commands_text
    assert "$OUTPUT_DIR" in commands_text
    assert str(out) not in commands_text

    environment = (repro / "environment.yml").read_text(encoding="utf-8")
    assert "name: clawbio-crispr-screen-triage" in environment

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
    assert {"report.md", "result.json", "tables/triaged_genes.csv", "tables/guide_metrics.csv"} <= labels


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
