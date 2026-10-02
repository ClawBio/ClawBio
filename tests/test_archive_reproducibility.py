"""An archive skill's reproducibility/commands.sh must actually reproduce the run.

It used to record only `--command <first command> --output <dir>`: no
accession, no flags, and `--command runs` for a `--demo` run. Replaying it
either failed ("needs --accession") or ran something else. These tests replay
the recorded command with bash and compare the output with the original.
"""
from __future__ import annotations

import os
import shlex
import subprocess
import sys
from pathlib import Path

import pytest

from tests.test_archive_command_coverage import SKILLS, _load_entry_point


def _recorded_argv(output_dir: Path) -> list[str]:
    line = (output_dir / "reproducibility" / "commands.sh").read_text().splitlines()[1]
    return shlex.split(line)[2:]  # drop `python <script>`


@pytest.mark.parametrize("skill", SKILLS)
def test_demo_run_is_recorded_as_a_demo_run(skill, tmp_path):
    app = _load_entry_point(skill)
    app.main(["--demo", "--output", str(tmp_path)])
    assert _recorded_argv(tmp_path) == ["--demo", "--output", str(tmp_path)]


def test_every_flag_of_a_command_run_is_recorded(tmp_path):
    app = _load_entry_point("ena-fetch")
    argv = ["--demo", "--command", "samplesheet", "--assay", "bulk",
            "--fastq-dir", "/data/fastq", "--output", str(tmp_path)]
    app.main(argv)
    assert _recorded_argv(tmp_path) == argv


@pytest.mark.parametrize("skill", SKILLS)
def test_replaying_commands_sh_reproduces_the_report(skill, tmp_path):
    app = _load_entry_point(skill)
    app.main(["--demo", "--output", str(tmp_path)])
    original = (tmp_path / "report.md").read_text()
    (tmp_path / "report.md").unlink()

    env = dict(os.environ,
               PATH=f"{Path(sys.executable).parent}{os.pathsep}{os.environ['PATH']}")
    result = subprocess.run(["bash", str(tmp_path / "reproducibility" / "commands.sh")],
                            env=env, capture_output=True, text=True, timeout=120)
    assert result.returncode == 0, result.stderr
    assert (tmp_path / "report.md").read_text() == original
