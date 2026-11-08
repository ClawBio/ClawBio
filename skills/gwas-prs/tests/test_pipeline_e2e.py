"""Drive skills/gwas-prs/pipeline.yaml end to end through clawbio.workflow.commands.

start -> (resolve_scores, parse_genotype) -> calculate_prs -> report -> bundle,
running each rendered command as a real subprocess on the bundled demo genotype
with `pgs_ids=all` (bundled curated panels only).

Offline guarantee: the child processes get a sitecustomize.py on PYTHONPATH that
makes every socket connect raise, so a fetch attempt would fail loudly. The test
also asserts no manifest entry has source "fetched" or "error".
"""
from __future__ import annotations

import json
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

from clawbio.workflow import commands
from clawbio.workflow import state as st

SKILL_DIR = Path(__file__).resolve().parents[1]
DEMO = SKILL_DIR / "demo_patient_prs.txt"

# step id -> {emit name: filename written under <run_dir>/out}
STEP_OUTPUTS = {
    "resolve_scores": {"scoring_manifest_json": "scoring_manifest.json"},
    "parse_genotype": {"genotype_json": "genotype.json"},
    "calculate_prs": {"results_json": "prs_results.json", "variants_csv": "prs_variants.csv"},
    "report": {"report_md": "report.md", "result_json": "result.json"},
}
ORDER = ["resolve_scores", "parse_genotype", "calculate_prs", "report"]


@pytest.fixture
def offline_env(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """Make `python3` resolve to this interpreter and block sockets in children."""
    block = tmp_path / "_netblock"
    block.mkdir()
    (block / "sitecustomize.py").write_text(textwrap.dedent("""
        import socket

        def _no_net(*a, **k):
            raise OSError("network disabled by test")

        socket.socket.connect = _no_net
        socket.socket.connect_ex = _no_net
        socket.create_connection = _no_net
    """))
    monkeypatch.setenv("PATH", f"{Path(sys.executable).parent}:{__import__('os').environ['PATH']}")
    monkeypatch.setenv("PYTHONPATH", str(block))
    return block


def test_pipeline_yaml_end_to_end(tmp_path: Path, monkeypatch, offline_env):
    work = tmp_path / "work"
    work.mkdir()
    monkeypatch.chdir(work)

    # Sanity: the child really cannot open sockets.
    probe = subprocess.run(
        ["python3", "-c", "import socket; socket.socket().connect(('127.0.0.1', 9))"],
        capture_output=True, text=True,
    )
    assert "network disabled by test" in probe.stderr

    commands.cmd_start(str(SKILL_DIR), [f"genotype_file={DEMO}", "pgs_ids=all"])

    for step_id in ORDER:
        commands.cmd_open(step_id)
        rendered = st.read()["steps"][step_id]["command"]
        subprocess.run(rendered, shell=True, check=True, executable="/bin/bash")
        out_dir = st.run_dir() / "out"
        commands.cmd_done(
            step_id,
            [f"{k}={out_dir / fname}" for k, fname in STEP_OUTPUTS[step_id].items()],
        )

    state = st.read()
    assert all(e["status"] == st.DONE for e in state["steps"].values())

    out_dir = st.run_dir() / "out"
    manifest = json.loads((out_dir / "scoring_manifest.json").read_text())
    assert manifest and all(m["source"] == "bundled" for m in manifest)  # nothing fetched

    report_path = out_dir / "report.md"
    result_path = out_dir / "result.json"
    assert report_path.is_file() and result_path.is_file()

    report = report_path.read_text()
    assert "curated demo panel" in report
    assert "**Scoring files**: PGS Catalog (" not in report
    assert "| **PGS Catalog ID** |" not in report
    assert "| **Publication** |" not in report

    results = json.loads(result_path.read_text())["data"]["results"]
    assert results and all(r["curated_demo_panel"] is True for r in results)

    from clawbio.workflow import bundle
    bundle.cmd_bundle()
    run = st.run_dir()
    assert (run / "commands.sh").is_file()
    assert (run / "checksums.sha256").is_file()
    assert (run / "environment.yml").is_file()
    assert st.read()["bundled_at"] is not None
