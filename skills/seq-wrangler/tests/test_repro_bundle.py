"""seq-wrangler writes the shared reproducibility bundle, and it replays."""

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
SCRIPT = SKILL_DIR / "seq_wrangler.py"
sys.path.insert(0, str(SKILL_DIR))

import seq_wrangler as sw  # noqa: E402

from clawbio.common.checksums import sha256_file  # noqa: E402

TOOLS = ("samtools", "bwa", "bowtie2", "minimap2", "fastqc", "fastp", "multiqc")


@pytest.fixture(autouse=True)
def _cwd_in_tmp(tmp_path: Path, monkeypatch) -> None:
    """Run each test, and any subprocess it starts, from its own tmp_path (AGENTS.md)."""
    monkeypatch.chdir(tmp_path)


def _verify_checksums(out: Path) -> set[str]:
    lines = (out / "reproducibility" / "checksums.sha256").read_text().splitlines()
    assert lines, "checksums.sha256 is empty"
    labels = set()
    for digest, label in (line.split("  ", 1) for line in lines):
        assert not Path(label).is_absolute(), label
        assert digest == sha256_file(out / label), label
        labels.add(label)
    return labels


@pytest.fixture(scope="module")
def demo_out(tmp_path_factory) -> Path:
    out = tmp_path_factory.mktemp("demo") / "out"
    subprocess.run([sys.executable, str(SCRIPT), "--demo", "--output", str(out)],
                   check=True, capture_output=True, text=True, cwd=out.parent)
    return out


def test_demo_checksums_cover_every_output(demo_out: Path) -> None:
    labels = _verify_checksums(demo_out)
    assert {"report.md", "summary.json", "reproducibility/commands.sh",
            "reproducibility/environment.yml"} <= labels


def test_demo_commands_sh_is_portable(demo_out: Path) -> None:
    commands = (demo_out / "reproducibility" / "commands.sh").read_text()
    assert '"$CLAWBIO_ROOT/skills/seq-wrangler/seq_wrangler.py"' in commands
    assert '--output \\\n  "$OUTPUT_DIR"' in commands
    assert str(demo_out) not in commands
    # The demo calls no external tools, so its replay does not probe them.
    assert "samtools --version" not in commands


def test_demo_replays_from_a_copy(demo_out: Path, tmp_path: Path) -> None:
    copy = tmp_path / "copy"
    shutil.copytree(demo_out, copy)
    elsewhere = tmp_path / "elsewhere"  # replay from a directory unrelated to the output
    elsewhere.mkdir()
    subprocess.run(["bash", str(copy / "reproducibility" / "commands.sh")], check=True,
                   capture_output=True, cwd=elsewhere,
                   env={"PATH": "/usr/bin:/bin", "PYTHON": sys.executable})
    report = (copy / "report.md").read_text()
    assert "Demo (synthetic data)" in report
    assert not any(elsewhere.iterdir()), "replay wrote into the working directory"


def _replayed_argv(out: Path) -> list[str]:
    """Run commands.sh with stand-ins: tools that exit 0, a python that prints its argv."""
    stub = out.parent / "stub_bin"
    stub.mkdir(exist_ok=True)
    for tool in TOOLS:
        (stub / tool).write_text("#!/bin/sh\nexit 0\n")
        (stub / tool).chmod(0o755)
    fake = stub / "fake_python"
    fake.write_text('#!/usr/bin/env bash\nprintf "%s\\0" "$@"\n')
    fake.chmod(0o755)
    before = set(out.parent.iterdir())
    run = subprocess.run(["bash", str(out / "reproducibility" / "commands.sh")],
                         env={"PATH": f"{stub}:/usr/bin:/bin", "PYTHON": str(fake)},
                         capture_output=True, check=True, cwd=out.parent)
    assert set(out.parent.iterdir()) == before, "replay created files next to the output"
    return run.stdout.decode("utf-8", "surrogateescape").split("\0")[1:-1]


def _args(**overrides):
    args = sw.build_parser().parse_args(["--output", "unused"])
    for key, value in overrides.items():
        setattr(args, key, value)
    return args


def test_pipeline_bundle_records_inputs_and_options(tmp_path: Path) -> None:
    out = tmp_path / "out"
    out.mkdir()
    (out / "report.md").write_text("# report\n")
    r1, r2 = tmp_path / "s R1.fq.gz", tmp_path / "s R2.fq.gz"
    r1.write_bytes(b"r1"), r2.write_bytes(b"r2")
    args = _args(r1=str(r1), r2=str(r2), index="/ref/hg38 idx", aligner="bowtie2",
                 threads=4, mapq=30, trim=True, remove_duplicates=True, output=str(out))
    sw.write_reproducibility_bundle(out, [r1, r2], args, is_demo=False)

    assert _replayed_argv(out) == [
        "--r1", str(r1), "--r2", str(r2), "--index", "/ref/hg38 idx",
        "--aligner", "bowtie2", "--mapq", "30", "--threads", "4",
        "--trim", "--remove-duplicates", "--output", str(out),
    ]
    inputs = json.loads((out / "reproducibility" / "inputs.json").read_text())
    assert inputs == {"files": [
        {"file": r1.name, "sha256": sha256_file(r1)},
        {"file": r2.name, "sha256": sha256_file(r2)},
    ]}
    labels = _verify_checksums(out)
    assert {"report.md", "reproducibility/inputs.json"} <= labels


def test_samplesheet_run_records_the_sheet(tmp_path: Path) -> None:
    out = tmp_path / "out"
    out.mkdir()
    fq = tmp_path / "a.fq.gz"
    fq.write_bytes(b"a")
    sheet = tmp_path / "samples.csv"
    sheet.write_text(f"sample,fastq_1,fastq_2\nA,{fq},\n")
    args = _args(samplesheet=str(sheet), index="/ref/idx", threads=1, output=str(out))
    sw.write_reproducibility_bundle(out, [fq], args, is_demo=False)

    argv = _replayed_argv(out)
    assert argv[:2] == ["--samplesheet", str(sheet)]
    inputs = json.loads((out / "reproducibility" / "inputs.json").read_text())
    assert inputs["samplesheet"] == {"file": "samples.csv", "sha256": sha256_file(sheet)}


# Any POSIX file name: no NUL or "/". "\r" is left out because the shared LF
# writer behind commands.sh rewrites it to "\n" (clawbio/common/textio.py).
_name = st.text(st.characters(blacklist_categories=("Cs",), blacklist_characters="\x00/\r"),
                min_size=1, max_size=30).filter(lambda s: s not in (".", ".."))


@settings(max_examples=40, deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture])
@given(r1=_name, index=_name)
# The old commands.sh pasted paths unquoted, so shell syntax in a name ran on replay.
@example(r1="s R1>x;$(touch y)'\"`.fq", index="--aligner")
def test_replay_hands_paths_back_unchanged(tmp_path_factory, r1: str, index: str) -> None:
    root = tmp_path_factory.mktemp("prop")
    out = root / "out"
    out.mkdir()
    reads = root / r1
    reads.write_bytes(b"reads")
    args = _args(r1=str(reads), index=str(root / index), threads=1, output=str(out))
    sw.write_reproducibility_bundle(out, [reads], args, is_demo=False)
    argv = _replayed_argv(out)
    assert argv[argv.index("--r1") + 1] == str(reads)
    assert argv[argv.index("--index") + 1] == str(root / index)
    assert argv[-2:] == ["--output", str(out)]
