"""Drive skills/nutrigx/pipeline.yaml end to end through clawbio.workflow.commands.

start -> (open, run rendered command, done) x 4 -> bundle, on the synthetic
patient, which is heterozygous (GA) at rs4988235 (MCM6, dominant_protective).
A het there must leave the lactose domain at Low (regression for 6b427989).
"""
from __future__ import annotations

import importlib.util
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

from clawbio.workflow import bundle, commands
from clawbio.workflow import state as st

SKILL_DIR = Path(__file__).resolve().parents[1]
PATIENT = SKILL_DIR / "tests" / "synthetic_patient.csv"

# emit name -> filename each step actually writes (see pipeline.yaml templates
# and generate_report.py, which names its files nutrigx_*).
STEP_OUTPUTS = {
    "parse_input": {"genotype_table_json": "genotype_table.json"},
    "extract_genotypes": {"snp_calls_json": "snp_calls.json"},
    "score_variants": {"risk_scores_json": "risk_scores.json"},
    "generate_report": {
        "report_md": "nutrigx_report.md",
        "radar_png": "nutrigx_radar.png",
        "heatmap_png": "nutrigx_heatmap.png",
        "result_json": "result.json",
    },
}

# generate_report.py imports seaborn for the heatmap; it is only in the
# affinity-proteomics extra, so CI (--all-extras) runs this and a bare sync skips.
needs_seaborn = pytest.mark.skipif(
    importlib.util.find_spec("seaborn") is None,
    reason="needs seaborn (uv sync --extra affinity-proteomics)",
)


@needs_seaborn
def test_pipeline_end_to_end(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    # Rendered commands call python3: make it the interpreter running pytest.
    monkeypatch.setenv(
        "PATH", str(Path(sys.executable).parent) + os.pathsep + os.environ["PATH"]
    )
    assert "rs4988235\t2\t136616754\tGA" in PATIENT.read_text()

    commands.cmd_start(str(SKILL_DIR), [f"genetic_data={PATIENT}"])
    out_dir = Path(st.run_dir()) / "out"

    for step_id, emits in STEP_OUTPUTS.items():
        commands.cmd_open(step_id)
        rendered = st.read()["steps"][step_id]["command"]
        subprocess.run(rendered, shell=True, check=True, executable="/bin/bash")
        commands.cmd_done(step_id, [f"{k}={out_dir / v}" for k, v in emits.items()])

    state = st.read()
    assert all(s["status"] == st.DONE for s in state["steps"].values())

    report = out_dir / "nutrigx_report.md"
    result_path = out_dir / "result.json"
    assert report.is_file() and result_path.is_file()
    assert "does not constitute medical advice" in report.read_text()

    scores = json.loads((out_dir / "risk_scores.json").read_text())
    assert scores["lactose"]["category"] == "Low"
    result = json.loads(result_path.read_text())
    assert result["data"]["risk_scores"]["lactose"]["category"] == "Low"
    assert result["data"]["snp_calls"]["rs4988235"]["genotype"] in ("GA", "AG")

    bundle.cmd_bundle()
    run = st.run_dir()
    for name in ("commands.sh", "checksums.sha256", "environment.yml"):
        assert (run / name).is_file(), name
    assert st.read()["bundled_at"] is not None
