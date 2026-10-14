"""The launcher must choose a destination accepted by the nf-core wrappers."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

from clawbio import cli

PIPELINES = ("rnaseq-pipeline", "scrnaseq-pipeline", "sarek-pipeline")


@pytest.mark.parametrize("skill", PIPELINES)
def test_default_pipeline_output_is_outside_source_tree(monkeypatch, tmp_path, skill):
    repository = tmp_path / "repository"
    home = tmp_path / "home"
    monkeypatch.setattr(cli, "CLAWBIO_DIR", repository)
    monkeypatch.setattr(cli, "DEFAULT_OUTPUT_ROOT", repository / "output")
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: home))
    captured = []
    monkeypatch.setattr(
        cli.subprocess,
        "run",
        lambda cmd, **kw: (
            captured.append(cmd) or SimpleNamespace(returncode=0, stdout="", stderr="")
        ),
    )
    cli.run_skill(skill, demo=True)
    output = Path(captured[0][captured[0].index("--output") + 1])
    assert not output.is_relative_to(repository)
    assert output.parent == home / "clawbio-runs"


@pytest.mark.parametrize("skill", PIPELINES)
def test_explicit_pipeline_destination_is_preserved(monkeypatch, tmp_path, skill):
    destination = tmp_path / "selected"
    captured = []
    monkeypatch.setattr(
        cli.subprocess,
        "run",
        lambda cmd, **kw: (
            captured.append(cmd) or SimpleNamespace(returncode=0, stdout="", stderr="")
        ),
    )
    cli.run_skill(skill, demo=True, output_dir=str(destination))
    assert captured[0][captured[0].index("--output") + 1] == str(destination)


def test_non_pipeline_output_policy_is_preserved(monkeypatch, tmp_path):
    repository = tmp_path / "repository"
    monkeypatch.setattr(cli, "CLAWBIO_DIR", repository)
    monkeypatch.setattr(cli, "DEFAULT_OUTPUT_ROOT", repository / "output")
    captured = []
    monkeypatch.setattr(
        cli.subprocess,
        "run",
        lambda cmd, **kw: (
            captured.append(cmd) or SimpleNamespace(returncode=0, stdout="", stderr="")
        ),
    )
    cli.run_skill("pharmgx", demo=True)
    output = Path(captured[0][captured[0].index("--output") + 1])
    assert output.parent == repository / "output"
