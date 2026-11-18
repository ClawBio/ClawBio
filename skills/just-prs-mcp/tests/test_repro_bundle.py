"""The demo writes the shared reproducibility bundle and its checksums verify."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

pytest.importorskip("typer", reason="needs the just-prs extra (uv sync --extra just-prs)")

SKILL_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SKILL_DIR))

import just_prs_mcp_bridge as bridge  # noqa: E402

from clawbio.common.checksums import sha256_file  # noqa: E402


def test_demo_writes_verifiable_repro_bundle(tmp_path: Path) -> None:
    out = tmp_path / "out"
    bridge.run_demo(out)
    repro = out / "reproducibility"
    for name in ("commands.sh", "environment.yml", "checksums.sha256"):
        assert (repro / name).is_file(), name

    entries = dict(
        reversed(line.split("  ", 1))
        for line in (repro / "checksums.sha256").read_text().splitlines()
    )
    assert {"report.md", "tables/scores.csv", "reproducibility/environment.yml"} <= set(entries)
    for label, digest in entries.items():
        assert not Path(label).is_absolute(), label
        assert digest == sha256_file(out / label), label

    env = (repro / "environment.yml").read_text()
    assert "name: clawbio-just-prs-mcp" in env
    assert "fastmcp" in env and "typer" in env
