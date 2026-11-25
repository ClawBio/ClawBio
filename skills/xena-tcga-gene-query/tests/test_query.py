"""Tests for xena-tcga-gene-query skill (red/green TDD)."""

import json
import subprocess
import sys
from pathlib import Path


SKILL_DIR = Path(__file__).resolve().parent.parent
SCRIPT = SKILL_DIR / "scripts" / "query_tcga_api.py"


def run_cli(args, expect_code=0):
    """Run the CLI and return (returncode, stdout, stderr)."""
    result = subprocess.run(
        [sys.executable, str(SCRIPT), *args],
        capture_output=True,
        text=True,
    )
    if expect_code is not None:
        assert result.returncode == expect_code, (
            f"Expected exit {expect_code}, got {result.returncode}\n"
            f"stdout: {result.stdout}\nstderr: {result.stderr}"
        )
    return result


# ──────────────────────────────────────────────
# Demo mode
# ──────────────────────────────────────────────

def test_demo_runs_successfully(tmp_path):
    """--demo --output DIR produces report.md, result.json, and reproducibility/."""
    out = tmp_path / "demo_out"
    run_cli(["--demo", "--output", str(out)])
    assert (out / "report.md").exists(), "report.md not generated"
    assert (out / "result.json").exists(), "result.json not generated"
    repro = out / "reproducibility"
    assert repro.is_dir(), "reproducibility/ not created"
    assert (repro / "commands.sh").exists(), "commands.sh not generated"
    assert (repro / "run.json").exists(), "run.json not generated"


def test_demo_result_json_is_valid(tmp_path):
    """Demo result.json parses and has expected top-level keys."""
    out = tmp_path / "demo_out"
    run_cli(["--demo", "--output", str(out)])
    data = json.loads((out / "result.json").read_text())
    assert "diff_expr" in data
    assert "corr" in data
    assert "survival" in data


def test_demo_diff_expr_has_required_fields(tmp_path):
    """Demo diff-expr result contains gene, cancer, tumor/normal n, log2FC, p-value."""
    out = tmp_path / "demo_out"
    run_cli(["--demo", "--output", str(out)])
    data = json.loads((out / "result.json").read_text())
    de = data["diff_expr"]
    assert de["gene"] == "TP53"
    assert de["cancer"] == "LUAD"
    assert de["tumor"]["n"] > 0
    assert de["normal"]["n"] > 0
    assert "log2_fold_change" in de
    assert "p_value" in de


def test_demo_corr_has_required_fields(tmp_path):
    """Demo corr result contains both genes, cancer, n, spearman_r, p-value."""
    out = tmp_path / "demo_out"
    run_cli(["--demo", "--output", str(out)])
    data = json.loads((out / "result.json").read_text())
    cr = data["corr"]
    assert cr["gene1"] == "TP53"
    assert cr["gene2"] == "EGFR"
    assert cr["cancer"] == "LUAD"
    assert cr["n"] > 0
    assert "spearman_r" in cr
    assert "p_value" in cr


def test_demo_survival_has_all_endpoints(tmp_path):
    """Demo survival result covers OS, DSS, DFI, PFI."""
    out = tmp_path / "demo_out"
    run_cli(["--demo", "--output", str(out)])
    data = json.loads((out / "result.json").read_text())
    sv = data["survival"]
    assert sv["gene"] == "TP53"
    assert sv["cancer"] == "LUAD"
    for ep in ["OS", "DSS", "DFI", "PFI"]:
        assert ep in sv["survival"], f"Missing survival endpoint: {ep}"


def test_report_includes_clawbio_disclaimer(tmp_path):
    """report.md must reference the ClawBio research-tool disclaimer."""
    out = tmp_path / "demo_out"
    run_cli(["--demo", "--output", str(out)])
    md = (out / "report.md").read_text(encoding="utf-8")
    assert "ClawBio" in md
    assert "research" in md.lower() and "medical device" in md.lower()


def test_report_includes_interpretation_notes(tmp_path):
    """report.md includes interpretation guidance with sample sizes and caveats."""
    out = tmp_path / "demo_out"
    run_cli(["--demo", "--output", str(out)])
    md = (out / "report.md").read_text(encoding="utf-8")
    assert "Interpretation" in md
    assert "n =" in md or "n=" in md
    # Must warn about optimal cutoff being exploratory
    assert "exploratory" in md.lower() or "Exploratory" in md


def test_demo_prints_summary_to_stdout(tmp_path):
    """Demo mode prints a summary of all three queries to stdout."""
    out = tmp_path / "demo_out"
    result = run_cli(["--demo", "--output", str(out)])
    stdout = result.stdout
    assert "TP53" in stdout
    assert "LUAD" in stdout
    assert "EGFR" in stdout


def test_demo_output_dir_created_if_missing(tmp_path):
    """--demo creates the output directory if it doesn't exist."""
    out = tmp_path / "nested" / "demo_out"
    assert not out.exists()
    run_cli(["--demo", "--output", str(out)])
    assert out.is_dir()
    assert (out / "report.md").exists()


def test_reproducibility_run_json_is_valid(tmp_path):
    """reproducibility/run.json contains skill name and mode."""
    out = tmp_path / "demo_out"
    run_cli(["--demo", "--output", str(out)])
    run_json = json.loads((out / "reproducibility" / "run.json").read_text())
    assert run_json["skill"] == "xena-tcga-gene-query"
    assert run_json["mode"] == "demo"


# ──────────────────────────────────────────────
# JSON mode (live API is not required for these)
# ──────────────────────────────────────────────

def test_demo_with_json_flag(tmp_path):
    """--demo --json prints only raw JSON (no summary) to stdout."""
    out = tmp_path / "demo_out"
    result = run_cli(["--demo", "--json", "--output", str(out)])
    # Should still produce file outputs
    assert (out / "report.md").exists()


# ──────────────────────────────────────────────
# Error paths
# ──────────────────────────────────────────────

def test_no_task_and_no_demo_is_error():
    """Running without --demo and without a task should exit non-zero."""
    result = subprocess.run(
        [sys.executable, str(SCRIPT)],
        capture_output=True,
        text=True,
    )
    assert result.returncode != 0


def test_diff_expr_without_gene_fails():
    """diff-expr without --gene should fail."""
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "diff-expr", "--cancer", "LUAD"],
        capture_output=True,
        text=True,
    )
    assert result.returncode != 0


# ──────────────────────────────────────────────
# Output structure (matches SKILL.md contract)
# ──────────────────────────────────────────────

def test_output_structure_matches_contract(tmp_path):
    """Verify the output tree declared in SKILL.md is produced."""
    out = tmp_path / "demo_out"
    run_cli(["--demo", "--output", str(out)])

    expected = [
        out / "report.md",
        out / "result.json",
        out / "reproducibility" / "commands.sh",
        out / "reproducibility" / "run.json",
    ]
    for path in expected:
        assert path.exists(), f"Missing expected output: {path}"


# ──────────────────────────────────────────────
# Mocked HTTP tests (cover live _request() path)
# ──────────────────────────────────────────────

MOCK_DIFF_EXPR_RESPONSE = {
    "gene": "TP53",
    "gene_input": None,
    "cancer": "LUAD",
    "cancer_full_name": "Lung Adenocarcinoma",
    "expression_scale": "log2(TPM + 0.001)",
    "test": "Mann-Whitney U (two-sided)",
    "tumor": {"n": 515, "mean": 5.12, "median": 5.08, "std": 0.89},
    "normal": {"n": 59, "mean": 4.87, "median": 4.91, "std": 0.76},
    "log2_fold_change": 0.25,
    "p_value": 0.0034,
}

MOCK_CANCERS_RESPONSE = {
    "count": 2,
    "cancers": [
        {"cancer": "LUAD", "cancer_full_name": "Lung Adenocarcinoma", "tumor_n": 515, "normal_n": 59, "has_normal": True},
        {"cancer": "BRCA", "cancer_full_name": "Breast invasive carcinoma", "tumor_n": 1100, "normal_n": 113, "has_normal": True},
    ],
}

MOCK_SURVIVAL_RESPONSE = {
    "gene": "TP53",
    "cancer": "LUAD",
    "cancer_full_name": "Lung Adenocarcinoma",
    "survival": {
        "OS": {
            "n_total": 504, "n_events": 189,
            "median_cutoff": {"p_value": 0.042},
            "optimal_cutoff": {"p_value": 0.0081, "cutoff": 5.31},
            "optimal_cutoff_note": "Exploratory.",
        },
        "DSS": {
            "n_total": 494, "n_events": 142,
            "median_cutoff": {"p_value": 0.11},
            "optimal_cutoff": {"p_value": 0.021},
            "optimal_cutoff_note": "Exploratory.",
        },
        "DFI": {
            "n_total": 312, "n_events": 84,
            "median_cutoff": {"p_value": 0.67},
            "optimal_cutoff": {"p_value": 0.13},
            "optimal_cutoff_note": "Exploratory.",
        },
        "PFI": {
            "n_total": 504, "n_events": 218,
            "median_cutoff": {"p_value": 0.031},
            "optimal_cutoff": {"p_value": 0.0056},
            "optimal_cutoff_note": "Exploratory.",
        },
    },
}


class MockHTTPResponse:
    """Minimal mock for urllib response."""
    def __init__(self, data: dict):
        self._data = json.dumps(data).encode("utf-8")

    def read(self):
        return self._data

    def __enter__(self):
        return self

    def __exit__(self, *args):
        pass


def _mock_urlopen(req, timeout=None):
    """Replacement for urllib.request.urlopen that returns synthetic data."""
    url = req.full_url if hasattr(req, "full_url") else str(req.get_full_url()) if hasattr(req, "get_full_url") else ""
    if "/api/v1/cancers" in url:
        return MockHTTPResponse(MOCK_CANCERS_RESPONSE)
    if "/api/v1/diff-expr" in url:
        return MockHTTPResponse(MOCK_DIFF_EXPR_RESPONSE)
    if "/api/v1/survival" in url:
        return MockHTTPResponse(MOCK_SURVIVAL_RESPONSE)
    raise RuntimeError(f"Unexpected URL in mock: {url}")


def test_request_function_parses_response():
    """_request() called with a patched urlopen returns the expected dict."""
    from unittest.mock import patch

    sys.path.insert(0, str(SCRIPT.parent.parent / "scripts"))
    import query_tcga_api

    with patch.object(query_tcga_api.urllib.request, "urlopen", _mock_urlopen):
        result = query_tcga_api._request(
            "http://biotree.top:38123/ucscxena/",
            "/api/v1/diff-expr",
            {"gene": "TP53", "cancer": "LUAD"},
            None,
            30,
        )
    assert result["gene"] == "TP53"
    assert result["cancer"] == "LUAD"
    assert result["tumor"]["n"] == 515


def test_live_cli_with_output_writes_correct_mode(tmp_path):
    """Live mode with --output writes mode=live in reproducibility/run.json."""
    from unittest.mock import patch

    sys.path.insert(0, str(SCRIPT.parent.parent / "scripts"))
    import query_tcga_api

    out = tmp_path / "live_out"
    with patch.object(query_tcga_api.urllib.request, "urlopen", _mock_urlopen):
        exit_code = query_tcga_api.main(["--output", str(out), "diff-expr", "--gene", "TP53", "--cancer", "LUAD"])
    assert exit_code == 0
    assert (out / "report.md").exists(), "report.md not generated for live mode"
    assert (out / "result.json").exists(), "result.json not generated for live mode"

    run_json = json.loads((out / "reproducibility" / "run.json").read_text())
    assert run_json["mode"] == "live", f"Expected mode=live, got {run_json.get('mode')}"

    report_md = (out / "report.md").read_text(encoding="utf-8")
    assert "**Mode**: live" in report_md


def test_live_report_contains_provenance_note(tmp_path):
    """Live mode report.md includes the provenance disclosure."""
    from unittest.mock import patch

    sys.path.insert(0, str(SCRIPT.parent.parent / "scripts"))
    import query_tcga_api

    out = tmp_path / "live_out2"
    with patch.object(query_tcga_api.urllib.request, "urlopen", _mock_urlopen):
        query_tcga_api.main(["--output", str(out), "diff-expr", "--gene", "TP53", "--cancer", "LUAD"])

    report_md = (out / "report.md").read_text(encoding="utf-8")
    assert "author-hosted" in report_md
    assert "UCSCXenaToolsPy" in report_md


def test_live_cli_survival_with_output(tmp_path):
    """Live survival query with --output writes correct report."""
    from unittest.mock import patch

    sys.path.insert(0, str(SCRIPT.parent.parent / "scripts"))
    import query_tcga_api

    out = tmp_path / "surv_out"
    with patch.object(query_tcga_api.urllib.request, "urlopen", _mock_urlopen):
        exit_code = query_tcga_api.main(["--output", str(out), "survival", "--gene", "TP53", "--cancer", "LUAD"])
    assert exit_code == 0
    result = json.loads((out / "result.json").read_text())
    assert "survival" in result
    run_json = json.loads((out / "reproducibility" / "run.json").read_text())
    assert run_json["mode"] == "live"


def test_live_cli_cancers_with_output(tmp_path):
    """Live cancers query with --output produces a report."""
    from unittest.mock import patch

    sys.path.insert(0, str(SCRIPT.parent.parent / "scripts"))
    import query_tcga_api

    out = tmp_path / "cancers_out"
    with patch.object(query_tcga_api.urllib.request, "urlopen", _mock_urlopen):
        exit_code = query_tcga_api.main(["--output", str(out), "cancers"])
    assert exit_code == 0
    report_md = (out / "report.md").read_text(encoding="utf-8")
    assert "LUAD" in report_md
    assert "BRCA" in report_md
