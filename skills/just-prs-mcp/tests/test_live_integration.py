from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

SKILL_DIR = Path(__file__).resolve().parents[1]
DEMO_VCF = SKILL_DIR / "examples" / "demo_patient.vcf"


@pytest.mark.integration
@pytest.mark.network
@pytest.mark.mcp
@pytest.mark.slow
def test_real_uvx_mcp_scores_public_pgs_model(tmp_path: Path) -> None:
    """Launch the pinned PyPI MCP and perform a real PGS Catalog score request."""
    if os.environ.get("CLAWBIO_RUN_LIVE_JUST_PRS") != "1":
        pytest.skip("Set CLAWBIO_RUN_LIVE_JUST_PRS=1 to run the live PGS Catalog test.")
    if shutil.which("uvx") is None:
        pytest.skip("uvx is not installed; a local stdio MCP process cannot be launched.")

    output_dir = tmp_path / "live"
    completed = subprocess.run(
        [
            sys.executable,
            str(SKILL_DIR / "just_prs_mcp_bridge.py"),
            "--input",
            str(DEMO_VCF),
            "--trait",
            "type 2 diabetes mellitus",
            "--superpopulation",
            "EUR",
            "--limit",
            "1",
            "--top-n",
            "1",
            "--output",
            str(output_dir),
            "--timeout-seconds",
            "600",
        ],
        capture_output=True,
        text=True,
        timeout=720,
        check=False,
    )

    assert completed.returncode == 0, completed.stderr
    result = json.loads((output_dir / "result.json").read_text(encoding="utf-8"))
    assert result["ok"] is True
    assert result["datasets"]["just-prs-mcp"] == "0.2.0"
    assert result["datasets"]["PGS Catalog"] == "live public metadata/scoring files"
    assert result["data"]["trait_id"] == "MONDO_0005148"
    assert result["data"]["n_requested"] == 1
    assert result["data"]["n_scored"] + result["data"]["n_failed"] == 1
    if result["data"]["scores"]:
        assert result["data"]["scores"][0]["pgs_id"].startswith("PGS")
    else:
        assert result["data"]["n_filtered"] == 1
        assert result["data"]["filter_summary"]
    serialized = json.dumps(result)
    assert "#CHROM" not in serialized
    assert "\tGT\t" not in serialized


@pytest.mark.integration
@pytest.mark.network
@pytest.mark.mcp
@pytest.mark.slow
def test_real_uvx_mcp_withholds_unreliable_absolute_risk(tmp_path: Path) -> None:
    """A real zero-coverage score must not be presented as an absolute-risk estimate."""
    if os.environ.get("CLAWBIO_RUN_LIVE_JUST_PRS") != "1":
        pytest.skip("Set CLAWBIO_RUN_LIVE_JUST_PRS=1 to run the live PGS Catalog test.")
    if shutil.which("uvx") is None:
        pytest.skip("uvx is not installed; a local stdio MCP process cannot be launched.")

    output_dir = tmp_path / "single"
    completed = subprocess.run(
        [
            sys.executable,
            str(SKILL_DIR / "just_prs_mcp_bridge.py"),
            "--input",
            str(DEMO_VCF),
            "--pgs-id",
            "PGS000014",
            "--superpopulation",
            "EUR",
            "--output",
            str(output_dir),
            "--timeout-seconds",
            "600",
        ],
        capture_output=True,
        text=True,
        timeout=720,
        check=False,
    )

    assert completed.returncode == 0, completed.stderr
    result = json.loads((output_dir / "result.json").read_text(encoding="utf-8"))
    score = result["data"]["scores"][0]
    assert result["data"]["trait"].startswith("Type 2 diabetes")
    assert score["weight_mass_coverage"] == 0.0
    assert score["absolute_risk"]["status"] == "unavailable"
    assert "coverage" in score["absolute_risk"]["reason"].lower()
