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


# ── Review of c3e156a (manuelcorpas): every user string reaches commands.sh inert ──


def _fake_mlst(tmp_path: Path) -> Path:
    tools = tmp_path / "bin"
    tools.mkdir()
    fake = tools / "mlst"
    # Prints a row for the assembly (always the last argument), whatever --scheme says.
    fake.write_text('#!/bin/bash\n[ "$1" = --version ] && { echo "mlst 2.35.0"; exit; }\n'
                    'echo -e "${@: -1}\\tecoli_achtman_4\\t131\\tadk(53)"\n')
    fake.chmod(0o755)
    return fake


def test_organism_survives_replay_without_running_it(tmp_path):
    marker = tmp_path / "pwned"
    organism = f"Escherichia; touch {marker}; :"
    out = tmp_path / "out"
    first = subprocess.run([sys.executable, str(SCRIPT), "--amrfinder", str(EXAMPLES / "demo_amrfinder.tsv"),
                            "--organism", organism, "--output", str(out)], capture_output=True, text=True)
    assert first.returncode == 0, first.stderr
    replay = _replay(out, dict(os.environ))
    assert replay.returncode == 0, replay.stderr
    assert not marker.exists()
    assert json.loads((out / "result.json").read_text())["amr"]["organism"] == organism


def test_mlst_scheme_survives_replay_without_running_it(tmp_path):
    marker = tmp_path / "pwned"
    scheme = f"ecoli_achtman_4; touch {marker}; :"
    fake = _fake_mlst(tmp_path)
    base = {"PATH": "/usr/bin:/bin", "HOME": os.environ.get("HOME", "/tmp"), "CLAWBIO_MLST_CMD": str(fake)}
    out = tmp_path / "out"
    first = subprocess.run(
        [sys.executable, str(SCRIPT), "--input", str(EXAMPLES / "demo_isolate.fasta"),
         "--amrfinder", str(EXAMPLES / "demo_amrfinder.tsv"),
         "--plasmidfinder", str(EXAMPLES / "demo_plasmidfinder.tsv"),
         "--mlst-scheme", scheme, "--output", str(out)],
        capture_output=True, text=True, env=base)
    assert first.returncode == 0, first.stderr
    assert not marker.exists()
    commands = (out / "reproducibility" / "commands.sh").read_text()
    assert "--mlst-scheme" in commands
    replay = _replay(out, base)
    assert replay.returncode == 0, replay.stderr
    assert not marker.exists()


@pytest.mark.parametrize("name", ['odd "quoted" mlst.tsv', "dollar $(touch PWNED) mlst.tsv",
                                  "back`tick` mlst.tsv"])
def test_input_path_with_shell_characters_survives_replay(tmp_path, name):
    src = tmp_path / name
    src.write_text((EXAMPLES / "demo_mlst.tsv").read_text())
    out = tmp_path / "out"
    first = subprocess.run([sys.executable, str(SCRIPT), "--mlst", str(src), "--output", str(out)],
                           capture_output=True, text=True, cwd=tmp_path)
    assert first.returncode == 0, first.stderr
    replay = subprocess.run(["bash", str(out / "reproducibility" / "commands.sh")], capture_output=True,
                            text=True, cwd=tmp_path, env={**os.environ, "PYTHON": sys.executable})
    assert replay.returncode == 0, replay.stderr
    assert not (tmp_path / "PWNED").exists()
    assert not list(tmp_path.glob("tick*"))


def test_tie_rerun_hint_quotes_paths(mod, tmp_path, capsys, monkeypatch):
    # The hint is meant to be pasted into a shell; a path with a space must stay one word.
    import shlex
    asm = tmp_path / "my isolate.fasta"
    asm.write_text((EXAMPLES / "demo_isolate.fasta").read_text())
    fake = tmp_path / "mlst"
    fake.write_text('#!/bin/bash\n[ "$1" = --version ] && { echo "mlst 2.35.0"; exit; }\n'
                    'echo "[mlst] WARNING: salmonella(3529)==ecoli_achtman_4(131) score=100 x" >&2\n'
                    'echo -e "${@: -1}\\tecoli_achtman_4\\t131\\tadk(53)"\n')
    fake.chmod(0o755)
    monkeypatch.setenv("CLAWBIO_MLST_CMD", str(fake))
    out = tmp_path / "out dir"
    rc = mod.main(["--input", str(asm), "--amrfinder", str(EXAMPLES / "demo_amrfinder.tsv"),
                   "--plasmidfinder", str(EXAMPLES / "demo_plasmidfinder.tsv"), "--output", str(out)])
    assert rc == 0
    hints = [l.strip() for l in capsys.readouterr().err.splitlines() if l.strip().startswith("python ")]
    assert hints
    for hint in hints:
        argv = shlex.split(hint)
        assert argv[argv.index("--input") + 1] == str(asm)
        assert argv[argv.index("--output") + 1].startswith(str(out) + "_")


def test_relative_paths_through_clawbio_run_resolve_against_cwd(tmp_path):
    # clawbio.py runs the script from the skill folder; without extra_path_flags a relative
    # --mlst would be looked up there instead of where the user typed it.
    (tmp_path / "mlst.tsv").write_text((EXAMPLES / "demo_mlst.tsv").read_text())
    (tmp_path / "amr.tsv").write_text((EXAMPLES / "demo_amrfinder.tsv").read_text())
    (tmp_path / "plasmid.tsv").write_text((EXAMPLES / "demo_plasmidfinder.tsv").read_text())
    r = subprocess.run([sys.executable, str(SKILL_DIR.parents[1] / "clawbio.py"), "run", "isolate-amr",
                        "--input", str(EXAMPLES / "demo_isolate.fasta"), "--mlst", "mlst.tsv",
                        "--amrfinder", "amr.tsv", "--plasmidfinder", "plasmid.tsv",
                        "--output", str(tmp_path / "out")],
                       capture_output=True, text=True, cwd=tmp_path)
    assert r.returncode == 0, r.stdout + r.stderr
    result = json.loads((tmp_path / "out" / "result.json").read_text())
    assert result["mlst"]["st"] == "131"
    assert result["amr"]["status"] != "not_run"
    assert result["plasmids"]["status"] != "not_run"
