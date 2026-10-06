"""The one test that runs real snakemake. Skipped unless the host can run it.

Everything else in this suite uses a stub, which proves the wrapper's plumbing
and nothing about the workflow. This one asks the pinned upstream checkout to
actually build its conda environments and run QC, assembly and reporting on
synthetic reads — minutes of compute and a network connection, so it is marked
`integration` and `slow`, and it skips rather than fails when snakemake or
conda is absent.

    pytest skills/metagenome-mag-pipeline/tests/ -m "integration and slow"

The synthetic reads have no biological truth: this asserts execution, never
accuracy.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

SKILL_DIR = Path(__file__).resolve().parents[1]
SCRIPT = SKILL_DIR / "metagenome_mag_pipeline.py"
REPO_ROOT = SKILL_DIR.parent.parent

pytestmark = [pytest.mark.integration, pytest.mark.slow]

# Section 3.5 of the plan: the files a real `assembly`-preset demo must produce.
EXPECTED_DEMO_OUTPUTS = (
    "results/final/qc/read_counts.tsv",
    "results/final/assembly/quast.tsv",
    "results/reports/assembly/quast.html",
    "results/final/diagnostic/bottleneck.tsv",
    "results/final/benchmarks/runtime.tsv",
    "results/provenance/software_versions.tsv",
    "results/provenance/pipeline_revision.tsv",
    "results/provenance/reference_db.md5",
    "results/provenance/run_manifest.json",
    "report.md",
    "result.json",
    "reproducibility/checksums.sha256",
)


@pytest.mark.skipif(shutil.which("snakemake") is None, reason="snakemake is not on PATH")
@pytest.mark.skipif(
    shutil.which("conda") is None and shutil.which("mamba") is None,
    reason="neither conda nor mamba is on PATH, so no tool environment can be built",
)
def test_real_demo_runs_the_pinned_workflow(tmp_path: Path) -> None:
    output = tmp_path / "demo"
    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--demo",
            "--output",
            str(output),
            "--timeout-hours",
            "2",
        ],
        capture_output=True,
        text=True,
        env={**os.environ, "PYTHONPATH": str(REPO_ROOT)},
        timeout=3 * 3600,
    )
    assert result.returncode == 0, f"demo failed:\n{result.stdout[-4000:]}\n{result.stderr[-4000:]}"
    data = json.loads((output / "result.json").read_text())
    assert data["status"] == "ok"
    assert data["mode"] == "demo"
    missing = [name for name in EXPECTED_DEMO_OUTPUTS if not (output / name).exists()]
    assert not missing, f"a real demo run did not produce: {missing}"
    assert data["findings"]["qc"]["read_counts"]
    assert data["findings"]["assembly"]["quast"]
