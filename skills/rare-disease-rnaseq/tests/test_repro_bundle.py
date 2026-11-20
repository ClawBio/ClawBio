"""rare-disease-rnaseq writes the shared reproducibility bundle, and it replays."""

from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
from hypothesis import HealthCheck, example, given, settings
from hypothesis import strategies as st

SKILL_DIR = Path(__file__).resolve().parents[1]
SCRIPT = SKILL_DIR / "rare_disease_rnaseq.py"
sys.path.insert(0, str(SKILL_DIR))

from rare_disease_rnaseq import generate_demo  # noqa: E402

from clawbio.common.checksums import sha256_file  # noqa: E402


@pytest.fixture(autouse=True)
def _cwd_in_tmp(tmp_path: Path, monkeypatch) -> None:
    """Run each test, and any subprocess it starts, from its own tmp_path (AGENTS.md)."""
    monkeypatch.chdir(tmp_path)


def _run(*args: str, cwd: Path) -> None:
    subprocess.run([sys.executable, str(SCRIPT), *args], check=True,
                   capture_output=True, text=True, cwd=cwd)


def _verify_checksums(out: Path) -> set[str]:
    lines = (out / "reproducibility" / "checksums.sha256").read_text().splitlines()
    labels = set()
    for digest, label in (line.split("  ", 1) for line in lines):
        assert not Path(label).is_absolute(), label
        assert digest == sha256_file(out / label), label
        labels.add(label)
    return labels


@pytest.fixture(scope="module")
def demo_out(tmp_path_factory) -> Path:
    out = tmp_path_factory.mktemp("demo") / "out"
    _run("--demo", "--output", str(out), cwd=out.parent)
    return out


def test_demo_checksums_cover_every_output(demo_out: Path) -> None:
    labels = _verify_checksums(demo_out)
    assert {
        "report.md", "result.json", "tables/outlier_calls.csv",
        "figures/case_outlier_heatmap.png", "demo_data/demo_counts.csv",
        "reproducibility/commands.sh", "reproducibility/environment.yml",
        "reproducibility/inputs.json",
    } <= labels


def test_demo_commands_sh_is_portable(demo_out: Path) -> None:
    commands = (demo_out / "reproducibility" / "commands.sh").read_text()
    assert '"$CLAWBIO_ROOT/skills/rare-disease-rnaseq/rare_disease_rnaseq.py"' in commands
    assert '--output \\\n  "$OUTPUT_DIR"' in commands
    assert str(demo_out) not in commands and sys.executable not in commands


def test_demo_replays_from_a_copy(demo_out: Path, tmp_path: Path) -> None:
    copy = tmp_path / "copy"
    shutil.copytree(demo_out, copy)
    elsewhere = tmp_path / "elsewhere"  # replay from a directory unrelated to the output
    elsewhere.mkdir()
    subprocess.run(["bash", str(copy / "reproducibility" / "commands.sh")], check=True,
                   capture_output=True, cwd=elsewhere,
                   env={"PATH": "/usr/bin:/bin", "PYTHON": sys.executable})
    _verify_checksums(copy)
    assert not any(elsewhere.iterdir()), "replay wrote into the working directory"


def test_inputs_json_records_the_panel(demo_out: Path) -> None:
    inputs = json.loads((demo_out / "reproducibility" / "inputs.json").read_text())
    panel = SKILL_DIR / "data" / "disease_panel.csv"
    assert inputs == {"panel": {"file": panel.name, "sha256": sha256_file(panel)}}


def _replayed_argv(out: Path) -> list[str]:
    """Run commands.sh with a stand-in python that prints the argv bash built."""
    fake = out.parent / "fake_python"
    fake.write_text('#!/usr/bin/env bash\nprintf "%s\\0" "$@"\n')
    fake.chmod(0o755)
    before = set(out.parent.iterdir())
    run = subprocess.run(["bash", str(out / "reproducibility" / "commands.sh")],
                         env={"PATH": "/usr/bin:/bin", "PYTHON": str(fake)},
                         capture_output=True, check=True, cwd=out.parent)
    assert set(out.parent.iterdir()) == before, "replay created files next to the output"
    return run.stdout.decode("utf-8", "surrogateescape").split("\0")[1:-1]


def test_user_run_records_inputs_and_options(tmp_path: Path) -> None:
    data = tmp_path / "my data"
    counts, cases, controls = generate_demo(seed=3, output_dir=data)
    out = tmp_path / "out"
    _run("--counts", str(counts), "--cases", str(cases), "--controls", str(controls),
         "--z-threshold", "2.5", "--output", str(out), cwd=tmp_path)

    assert _replayed_argv(out) == [
        "--counts", str(counts), "--cases", str(cases), "--controls", str(controls),
        "--z-threshold", "2.5", "--output", str(out),
    ]
    inputs = json.loads((out / "reproducibility" / "inputs.json").read_text())
    assert {k: v["sha256"] for k, v in inputs.items()} == {
        "counts": sha256_file(counts), "cases": sha256_file(cases),
        "controls": sha256_file(controls),
        "panel": sha256_file(SKILL_DIR / "data" / "disease_panel.csv"),
    }
    assert str(tmp_path) not in (out / "reproducibility" / "inputs.json").read_text()
    _verify_checksums(out)
    result = json.loads((out / "result.json").read_text())
    assert result["input_checksum"].removeprefix("sha256:") == sha256_file(counts)


# Any POSIX file name: no NUL or "/". "\r" is left out because the shared LF
# writer behind commands.sh rewrites it to "\n" (clawbio/common/textio.py).
_name = st.text(st.characters(blacklist_categories=("Cs",), blacklist_characters="\x00/\r"),
                min_size=1, max_size=30).filter(lambda s: s not in (".", ".."))


@settings(max_examples=25, deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture])
@given(name=_name, z=st.floats(min_value=0.1, max_value=50, allow_nan=False))
# The old commands.sh joined sys.argv unquoted, so shell syntax in a path ran on replay.
@example(name="my data>x;$(touch y)'\"`", z=2.5)
@example(name="0", z=3.0)  # the default threshold, which commands.sh omits
def test_replay_hands_paths_and_threshold_back_unchanged(tmp_path_factory, name: str, z: float) -> None:
    import rare_disease_rnaseq as rdr

    root = tmp_path_factory.mktemp("prop")
    data = root / name
    counts, cases, controls = generate_demo(seed=1, output_dir=data)
    out = root / "out"
    args = rdr.parse_args(["--counts", str(counts), "--cases", str(cases),
                           "--controls", str(controls), "--z-threshold", repr(z),
                           "--output", str(out)])
    out.mkdir()
    rdr.write_reproducibility(out, args, demo=False)
    argv = _replayed_argv(out)
    assert argv[argv.index("--counts") + 1] == str(counts)
    # The default threshold is left out, and the replay falls back to it.
    replayed_z = float(argv[argv.index("--z-threshold") + 1]) if "--z-threshold" in argv else rdr.DEFAULT_Z_THRESHOLD
    assert replayed_z == z
    assert argv[-2:] == ["--output", str(out)]
