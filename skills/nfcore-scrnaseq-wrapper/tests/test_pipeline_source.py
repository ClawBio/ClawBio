from __future__ import annotations

import shutil
from pathlib import Path
import subprocess
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from errors import SkillError
from pipeline_source import resolve_pipeline_source

_requires_git = pytest.mark.skipif(shutil.which("git") is None, reason="git not on PATH")


def test_resolves_to_remote_when_no_local_dir(tmp_path):
    absent = tmp_path / "nonexistent_scrnaseq"
    result = resolve_pipeline_source(
        requested_version="4.1.0",
        local_pipeline_dir=absent,
    )
    assert result["source_kind"] == "remote_repo"
    assert result["resolved_version"] == "4.1.0"
    assert result["source_ref"] == "nf-core/scrnaseq"
    assert result["dirty"] is False


def test_resolves_to_local_when_valid_checkout(tmp_path):
    local = tmp_path / "scrnaseq"
    local.mkdir()
    (local / "main.nf").write_text("// main", encoding="utf-8")
    (local / "nextflow.config").write_text("// config", encoding="utf-8")
    assets = local / "assets"
    assets.mkdir()
    (assets / "schema_input.json").write_text("{}", encoding="utf-8")
    result = resolve_pipeline_source(
        requested_version="4.1.0",
        local_pipeline_dir=local,
    )
    assert result["source_kind"] == "local_checkout"
    assert result["source_ref"] == str(local.resolve())
    # The test directory is not a git repo; all _git_stdout calls return "".
    assert result["resolved_version"] == "4.1.0"
    assert result["dirty"] is False
    assert result["branch"] == ""


def test_raises_when_local_dir_missing_required_files(tmp_path):
    local = tmp_path / "scrnaseq"
    local.mkdir()
    (local / "main.nf").write_text("// main", encoding="utf-8")
    with pytest.raises(SkillError) as exc:
        resolve_pipeline_source(
            requested_version="4.1.0",
            local_pipeline_dir=local,
        )
    assert exc.value.error_code == "PIPELINE_SOURCE_INVALID"
    assert "missing_files" in exc.value.details


@_requires_git
def test_local_checkout_dirty_includes_untracked_files(tmp_path):
    local = tmp_path / "scrnaseq"
    local.mkdir()
    (local / "main.nf").write_text("// main", encoding="utf-8")
    (local / "nextflow.config").write_text("// config", encoding="utf-8")
    assets = local / "assets"
    assets.mkdir()
    (assets / "schema_input.json").write_text("{}", encoding="utf-8")
    subprocess.run(["git", "init"], cwd=local, check=True, capture_output=True, text=True)
    result = resolve_pipeline_source(
        requested_version="4.1.0",
        local_pipeline_dir=local,
    )
    assert result["dirty"] is True
