"""Reproducibility bundle tests for drug-repurposing-screen.

Every run writes <output_dir>/reproducibility/ via clawbio.common.reproducibility:
commands.sh, environment.yml and checksums.sha256, where
`cd <output_dir> && sha256sum -c reproducibility/checksums.sha256` passes.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

SKILL_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(SKILL_DIR))
sys.path.insert(0, str(SKILL_DIR.parents[1]))

from drug_repurposing_screen import main as cli_main

from clawbio.common.checksums import sha256_file

pytestmark = pytest.mark.skipif(
    importlib.util.find_spec("pyarrow") is None,
    reason="needs pyarrow (uv sync --extra drug-repurposing-screen)",
)


def run_demo(out: Path) -> Path:
    assert cli_main(["--demo", "--output", str(out)]) == 0
    return out


def read_checksums(out: Path) -> dict[str, str]:
    lines = (out / "reproducibility" / "checksums.sha256").read_text().strip().splitlines()
    return {label: digest for digest, label in (line.split("  ", 1) for line in lines)}


def test_demo_writes_bundle(tmp_path):
    repro = run_demo(tmp_path / "out") / "reproducibility"
    for name in ("commands.sh", "environment.yml", "checksums.sha256"):
        assert (repro / name).is_file(), name


def test_checksums_verify_from_output_dir(tmp_path):
    out = run_demo(tmp_path / "out")
    entries = read_checksums(out)
    assert "tables/priority_table.csv" in entries
    for label, digest in entries.items():
        assert not Path(label).is_absolute(), label
        assert digest == sha256_file(out / label), label


def test_commands_sh_records_this_run_portably(tmp_path):
    out = run_demo(tmp_path / "out")
    content = (out / "reproducibility" / "commands.sh").read_text()
    assert '"$CLAWBIO_ROOT/skills/drug-repurposing-screen/drug_repurposing_screen.py"' in content
    assert "--demo" in content
    assert '"$OUTPUT_DIR"' in content
    assert str(out) not in content


def test_environment_yml_is_shared_format(tmp_path):
    content = (run_demo(tmp_path / "out") / "reproducibility" / "environment.yml").read_text()
    assert content.startswith("name: clawbio-drug-repurposing-screen")
    assert "pyarrow>=14.0" in content


def test_checksums_reproduce_on_replay(tmp_path):
    a = run_demo(tmp_path / "a") / "reproducibility" / "checksums.sha256"
    b = run_demo(tmp_path / "b") / "reproducibility" / "checksums.sha256"
    assert a.read_text() == b.read_text()
