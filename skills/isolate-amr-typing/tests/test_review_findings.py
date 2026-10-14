"""Regression tests for review findings on PR #540 (camlloyd).

Each one failed on the pre-fix head and now pins the fix in place:
MLST multi-line ties, commands.sh override quoting, --sample-name escaping,
and the 999_950 bp rounding boundary.
"""

from __future__ import annotations

import importlib.util
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

SKILL_DIR = Path(__file__).resolve().parents[1]
SCRIPT = SKILL_DIR / "isolate_amr_typing.py"
EXAMPLES = SKILL_DIR / "examples"


@pytest.fixture(scope="module")
def mod():
    spec = importlib.util.spec_from_file_location("isolate_amr_typing", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_three_way_mlst_tie_keeps_every_scheme(mod):
    # mlst 2.35.0 bin/mlst prints one wrn() line per scheme tied with the top scorer.
    stderr = (
        "[mlst] WARNING: salmonella(3529)==ecoli_achtman_4(131) score=100 x.fasta\n"
        "[mlst] WARNING: salmonella(3529)==klebsiella(11) score=100 x.fasta\n"
    )
    ties = mod.parse_mlst_ties(stderr)
    assert sorted(ties) == [("ecoli_achtman_4", "131"), ("klebsiella", "11"), ("salmonella", "3529")]
    assert mod.resolve_tie(ties, "Klebsiella_pneumoniae") == "klebsiella"


def test_fmt_bp_never_shows_1000_kb(mod):
    assert mod._fmt_bp(999_950) == "1.00 Mb"


def _replay(output: Path, env: dict) -> subprocess.CompletedProcess:
    return subprocess.run(["bash", str(output / "reproducibility" / "commands.sh")],
                          capture_output=True, text=True, env={**env, "PYTHON": sys.executable})


def test_tool_override_with_quoted_path_replays(tmp_path):
    tools = tmp_path / "my tools"
    tools.mkdir()
    fake = tools / "mlst"
    fake.write_text('#!/bin/bash\n[ "$1" = --version ] && { echo "mlst 2.35.0"; exit; }\n'
                    'echo -e "$1\\tecoli_achtman_4\\t131\\tadk(53)"\n')
    fake.chmod(0o755)
    base = {"PATH": "/usr/bin:/bin", "HOME": os.environ.get("HOME", "/tmp")}
    out = tmp_path / "out"
    first = subprocess.run(
        [sys.executable, str(SCRIPT), "--input", str(EXAMPLES / "demo_isolate.fasta"),
         "--amrfinder", str(EXAMPLES / "demo_amrfinder.tsv"),
         "--plasmidfinder", str(EXAMPLES / "demo_plasmidfinder.tsv"), "--output", str(out)],
        capture_output=True, text=True, env={**base, "CLAWBIO_MLST_CMD": f'"{fake}"'})
    assert first.returncode == 0, first.stderr
    # Replay on a machine without the override set: commands.sh must carry it intact.
    replay = _replay(out, base)
    assert replay.returncode == 0, replay.stderr


def test_sample_name_survives_replay(tmp_path):
    name = 'E. coli "JJ1886" $HOME'
    out = tmp_path / "out"
    first = subprocess.run([sys.executable, str(SCRIPT), "--mlst", str(EXAMPLES / "demo_mlst.tsv"),
                            "--sample-name", name, "--output", str(out)], capture_output=True, text=True)
    assert first.returncode == 0, first.stderr
    replay = _replay(out, dict(os.environ))
    assert replay.returncode == 0, replay.stderr
    assert json.loads((out / "result.json").read_text())["sample"] == name
