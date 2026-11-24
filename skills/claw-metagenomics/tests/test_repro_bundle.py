"""claw-metagenomics writes the shared reproducibility bundle, and it replays."""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
from hypothesis import HealthCheck, example, given, settings
from hypothesis import strategies as st

SKILL_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SKILL_DIR))

import metagenomics_profiler as mp  # noqa: E402

from clawbio.common.checksums import sha256_file  # noqa: E402

SCRIPT = SKILL_DIR / "metagenomics_profiler.py"


@pytest.fixture(autouse=True)
def _cwd_in_tmp(tmp_path: Path, monkeypatch) -> None:
    """Run each test, and any subprocess it starts, from its own tmp_path (AGENTS.md)."""
    monkeypatch.chdir(tmp_path)


def _checksums(out: Path) -> dict[str, str]:
    lines = (out / "reproducibility" / "checksums.sha256").read_text().splitlines()
    return {label: digest for digest, label in (line.split("  ", 1) for line in lines)}


def _assert_checksums_verify(out: Path) -> dict[str, str]:
    entries = _checksums(out)
    for label, digest in entries.items():
        assert not Path(label).is_absolute(), label
        assert digest == sha256_file(out / label), label
    return entries


def _run_demo(out: Path) -> None:
    subprocess.run([sys.executable, str(SCRIPT), "--demo", "--output", str(out)],
                   check=True, capture_output=True, text=True, cwd=out.parent)


@pytest.fixture(scope="module")
def demo_out(tmp_path_factory) -> Path:
    out = tmp_path_factory.mktemp("demo") / "out"
    _run_demo(out)
    return out


def test_demo_bundle_checksums_every_output(demo_out: Path) -> None:
    entries = _assert_checksums_verify(demo_out)
    assert {"report.md", "reproducibility/commands.sh", "reproducibility/environment.yml"} <= set(entries)
    assert any(label.startswith("figures/") for label in entries)
    assert any(label.startswith("tables/") for label in entries)


def test_demo_commands_sh_is_portable(demo_out: Path) -> None:
    commands = (demo_out / "reproducibility" / "commands.sh").read_text()
    assert '"$CLAWBIO_ROOT/skills/claw-metagenomics/metagenomics_profiler.py"' in commands
    assert '--output \\\n  "$OUTPUT_DIR"' in commands
    assert str(demo_out) not in commands
    # The demo needs no external tools, so the replay does not check for them.
    assert "kraken2 --version" not in commands


def test_demo_replays_from_a_copy(demo_out: Path, tmp_path: Path) -> None:
    copy = tmp_path / "copy"
    shutil.copytree(demo_out, copy)
    elsewhere = tmp_path / "elsewhere"  # replay from a directory unrelated to the output
    elsewhere.mkdir()
    subprocess.run(["bash", str(copy / "reproducibility" / "commands.sh")],
                   check=True, capture_output=True, cwd=elsewhere,
                   env={"PATH": "/usr/bin:/bin", "PYTHON": sys.executable})
    _assert_checksums_verify(copy)
    assert not any(elsewhere.iterdir()), "replay wrote into the working directory"


def test_environment_yml_lists_the_pipeline_tools(demo_out: Path) -> None:
    env = (demo_out / "reproducibility" / "environment.yml").read_text()
    assert "name: clawbio-metagenomics" in env
    assert "bioconda" in env
    for tool in ("kraken2", "bracken", "rgi", "humann"):
        assert f"  - {tool}" in env, tool


def _pipeline_args(**overrides) -> argparse.Namespace:
    base = dict(r1=None, r2=None, input=None, output=None, skip_functional=False,
                kraken2_db=None, humann_db=None, read_length=mp.DEFAULT_READ_LENGTH)
    base.update(overrides)
    return argparse.Namespace(**base)


def _replayed_argv(out: Path) -> list[str]:
    """Run commands.sh with a stand-in python that prints the argv bash built."""
    fake = out.parent / "fake_python"
    fake.write_text('#!/usr/bin/env bash\nprintf "%s\\0" "$@"\n')
    fake.chmod(0o755)
    stub_bin = out.parent / "stub_bin"  # the replay checks tool versions first
    stub_bin.mkdir(exist_ok=True)
    for tool in ("kraken2", "bracken", "rgi", "humann"):
        (stub_bin / tool).write_text("#!/bin/sh\nexit 0\n")
        (stub_bin / tool).chmod(0o755)
    before = set(out.parent.iterdir())
    run = subprocess.run(
        ["bash", str(out / "reproducibility" / "commands.sh")],
        env={"PATH": f"{stub_bin}:/usr/bin:/bin", "PYTHON": str(fake)},
        capture_output=True, check=True, cwd=out.parent,
    )
    assert set(out.parent.iterdir()) == before, "replay created files next to the output"
    return run.stdout.decode("utf-8", "surrogateescape").split("\0")[1:-1]  # drop the script path


def test_pipeline_bundle_records_inputs_and_options(tmp_path: Path) -> None:
    out = tmp_path / "out"
    out.mkdir()
    (out / "report.md").write_text("# report\n")
    r1, r2 = tmp_path / "s R1.fq.gz", tmp_path / "s R2.fq.gz"
    r1.write_bytes(b"r1"), r2.write_bytes(b"r2")
    args = _pipeline_args(r1=str(r1), r2=str(r2), output=str(out), skip_functional=True,
                          kraken2_db="/db/k2", read_length=100)
    mp.write_reproducibility_bundle(out, [r1, r2], args, is_demo=False)

    assert _replayed_argv(out) == [
        "--r1", str(r1), "--r2", str(r2), "--skip-functional",
        "--kraken2-db", "/db/k2", "--read-length", "100", "--output", str(out),
    ]
    inputs = json.loads((out / "reproducibility" / "inputs.json").read_text())
    assert inputs == {
        "r1": {"file": r1.name, "sha256": sha256_file(r1)},
        "r2": {"file": r2.name, "sha256": sha256_file(r2)},
    }
    entries = _assert_checksums_verify(out)
    assert "reproducibility/inputs.json" in entries
    assert str(tmp_path) not in (out / "reproducibility" / "checksums.sha256").read_text()


# Any POSIX file name: no NUL or "/". "\r" is left out because the shared LF
# writer behind commands.sh rewrites it to "\n" (clawbio/common/textio.py).
_name = st.text(st.characters(blacklist_categories=("Cs",), blacklist_characters="\x00/\r"),
                min_size=1, max_size=30).filter(lambda s: s not in (".", ".."))


@settings(max_examples=40, deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture])
@given(r1=_name, db=_name)
# The old commands.sh pasted paths unquoted, so shell syntax in a name ran on replay.
@example(r1="s R1>x;$(touch y)'\"`.fq", db="db dir")
def test_replay_hands_paths_back_unchanged(tmp_path_factory, r1: str, db: str) -> None:
    root = tmp_path_factory.mktemp("prop")
    out = root / "out"
    out.mkdir()
    reads = root / r1
    reads.write_bytes(b"reads")
    args = _pipeline_args(input=str(reads), output=str(out), kraken2_db=str(root / db),
                          humann_db=str(root / db))
    mp.write_reproducibility_bundle(out, [reads], args, is_demo=False)
    argv = _replayed_argv(out)
    assert argv[argv.index("--input") + 1] == str(reads)
    assert argv[argv.index("--kraken2-db") + 1] == str(root / db)
    assert argv[argv.index("--humann-db") + 1] == str(root / db)
    assert argv[-2:] == ["--output", str(out)]
