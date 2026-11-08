"""Curated demo panels must not be reported as PGS Catalog scores (2.0 scripts).

The 2.0 scripts under skills/gwas-prs/scripts/ shadow 1.x module names, so they
are driven by subprocess here. Offline: `--pgs-ids all` uses bundled panels only.
"""
from __future__ import annotations

import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import pytest

SKILL_DIR = Path(__file__).resolve().parents[1]
SCRIPTS = SKILL_DIR / "scripts"
DEMO = SKILL_DIR / "demo_patient_prs.txt"


def _run(script: str, *args: str, cwd: Path) -> None:
    subprocess.run(
        [sys.executable, str(SCRIPTS / script), *args],
        check=True, cwd=cwd, capture_output=True, text=True,
    )


@pytest.fixture()
def outputs(tmp_path):
    out = tmp_path / "out"
    out.mkdir()
    _run("resolve_scores.py", "--pgs-ids", "all", "--skill-dir", str(SKILL_DIR),
         "--out", str(out / "scoring_manifest.json"), cwd=tmp_path)
    _run("parse_genotype.py", "--input", str(DEMO),
         "--out", str(out / "genotype.json"), cwd=tmp_path)
    _run("calculate_prs.py", "--genotype-json", str(out / "genotype.json"),
         "--scoring-manifest", str(out / "scoring_manifest.json"),
         "--out-results", str(out / "prs_results.json"),
         "--out-variants", str(out / "prs_variants.csv"), cwd=tmp_path)
    _run("report.py", "--results", str(out / "prs_results.json"),
         "--genotype-json", str(out / "genotype.json"),
         "--out-report", str(out / "report.md"),
         "--out-result-json", str(out / "result.json"), cwd=tmp_path)
    return out


def test_report_labels_curated_panel(outputs):
    report = (outputs / "report.md").read_text()
    assert "curated demo panel" in report
    assert "NOT a PGS Catalog score" in report
    assert "**Scoring files**: PGS Catalog (" not in report
    assert "| **Publication** |" not in report
    assert "| **PGS Catalog ID** |" not in report
    assert "| **Loci reference** |" in report
    assert "| **Historical alias** |" in report
    assert "| **Source** |" in report


def test_result_json_flags_curated(outputs):
    data = json.loads((outputs / "result.json").read_text())["data"]["results"]
    assert data
    for r in data:
        assert r["curated_demo_panel"] is True
        assert r["curated_panel_id"].startswith("CLAWBIO-")
        assert r["legacy_pgs_id"]


def test_non_curated_still_renders_pgs_catalog_wording():
    spec = importlib.util.spec_from_file_location(
        "gwas_prs_scripts_report", SCRIPTS / "report.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    prs = {"raw_score": 1.0, "variants_used": 1, "variants_total": 1,
           "variants_missing": 0, "overlap_fraction": 1.0, "per_variant": []}
    pct = {"percentile": None, "risk_category": None, "z_score": None, "method": None}
    result = {"pgs_id": "PGS000001", "trait": "T", "metadata": {"publication": "X 2020"},
              "prs": prs, "percentile_info": pct}
    text = mod.generate_report([result], {"filepath": "f", "format": "x", "total_snps": 1})
    assert "**Scoring files**: PGS Catalog (" in text
    assert "| **Publication** | X 2020 |" in text
    assert "curated demo panel" not in text
