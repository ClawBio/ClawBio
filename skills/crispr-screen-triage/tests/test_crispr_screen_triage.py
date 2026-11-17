import importlib.util
import json
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
from hypothesis import given
from hypothesis import strategies as st


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


def _guide(gene, control, treated, essentiality=0.5, druggability=0.5, guide_id=None):
    return {"guide_id": guide_id or f"g{gene}_{control}_{treated}", "gene": gene,
            "control_count": float(control), "treatment_count": float(treated),
            "essentiality": essentiality, "druggability": druggability}


@pytest.mark.parametrize("score, median_log2fc, expected", [
    (1.35, -1.0, "high"),     # both "high" thresholds are inclusive
    (1.34, -1.0, "medium"),
    (1.35, -0.99, "medium"),  # a high score alone is not enough: the depletion gate decides
    (2.0, -0.5, "medium"),
    (0.75, 0.0, "medium"),
    (0.74, 0.0, "watch"),
])
def test_priority_boundaries(score, median_log2fc, expected):
    assert load_module()._priority(score, median_log2fc) == expected


def test_median_of_odd_and_even_guide_counts():
    median = load_module()._median
    assert median([5.0]) == 5.0
    assert median([3.0, 1.0, 2.0]) == 2.0
    assert median([4.0, 1.0, 3.0, 2.0]) == 2.5


def test_low_count_flag_is_below_100_total_reads():
    result = load_module().triage_genes([_guide("A", 49, 50), _guide("A", 50, 50)])
    assert [g["low_count"] for g in result["guides"]] == [True, False]
    assert result["genes"][0]["low_count_guides"] == 1


def test_genes_with_equal_scores_are_ordered_by_name():
    result = load_module().triage_genes([_guide("ZZZ", 500, 100), _guide("AAA", 500, 100)])
    assert [g["gene"] for g in result["genes"]] == ["AAA", "ZZZ"]


def test_load_counts_names_each_missing_column(tmp_path):
    module = load_module()
    for column in sorted(module.REQUIRED_COLUMNS):
        kept = sorted(module.REQUIRED_COLUMNS - {column})
        path = tmp_path / f"no_{column}.csv"
        path.write_text(",".join(kept) + "\n" + ",".join("1" for _ in kept) + "\n", encoding="utf-8")
        with pytest.raises(ValueError, match=f"missing required columns: .*'{column}'"):
            module.load_counts(path)


_RANK = {"watch": 0, "medium": 1, "high": 2}


@given(
    control=st.integers(0, 5000),
    treated=st.integers(0, 5000),
    fewer=st.integers(0, 5000),
    essentiality=st.floats(0, 1),
    druggability=st.floats(0, 1),
)
def test_more_depletion_never_lowers_priority(control, treated, fewer, essentiality, druggability):
    module = load_module()
    depleted = min(treated, fewer)

    def gene(t):
        return module.triage_genes([_guide("G", control, t, essentiality, druggability)])["genes"][0]

    base, more = gene(treated), gene(depleted)
    assert more["priority_score"] >= base["priority_score"]
    assert _RANK[more["priority"]] >= _RANK[base["priority"]]


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


@pytest.mark.parametrize(
    "report",
    ["report.md", "result.json", "tables/triaged_genes.csv", "tables/guide_metrics.csv"],
)
def test_reports_do_not_record_the_input_directory(tmp_path, report):
    """Reports are shared; an absolute input path leaks the user's name and folder layout."""
    private_dir = tmp_path / "private_user_dir"
    private_dir.mkdir()
    counts = shutil.copy(SKILL_DIR / "demo_screen_counts.csv", private_dir / "screen.csv")
    out = tmp_path / "out"
    subprocess.run(
        [sys.executable, str(MODULE_PATH), "--input", str(counts), "--output", str(out)],
        text=True,
        capture_output=True,
        check=True,
    )
    text = (out / report).read_text(encoding="utf-8")
    assert str(private_dir) not in text and str(private_dir.resolve()) not in text


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
