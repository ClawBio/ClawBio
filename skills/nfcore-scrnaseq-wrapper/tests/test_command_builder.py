from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from command_builder import build_nextflow_command


def test_local_checkout_command(tmp_path):
    params_path = tmp_path / "params.yaml"
    params_path.write_text("aligner: star\n", encoding="utf-8")
    source = {
        "source_kind": "local_checkout",
        "source_ref": "/abs/path/to/scrnaseq",
        "resolved_version": "abc1234",
    }
    cmd, cmd_str = build_nextflow_command(
        pipeline_source=source,
        profile="docker",
        params_path=params_path,
        resume=False,
    )
    assert cmd[0] == "nextflow"
    assert cmd[1] == "run"
    assert cmd[2] == "/abs/path/to/scrnaseq"
    assert "-profile" in cmd
    assert "docker" in cmd
    assert "-params-file" in cmd
    assert params_path.as_posix() in cmd  # forward slashes on all platforms
    assert "-resume" not in cmd
    assert "-r" not in cmd


def test_work_dir_is_under_wrapper_output_when_provided(tmp_path):
    params_path = tmp_path / "params.yaml"
    params_path.write_text("aligner: star\n", encoding="utf-8")
    work_dir = tmp_path / "upstream" / "work"
    source = {
        "source_kind": "local_checkout",
        "source_ref": "/abs/path/to/scrnaseq",
        "resolved_version": "abc1234",
    }
    cmd, _cmd_str = build_nextflow_command(
        pipeline_source=source,
        profile="docker",
        params_path=params_path,
        resume=False,
        work_dir=work_dir,
    )
    assert "-work-dir" in cmd
    assert cmd[cmd.index("-work-dir") + 1] == work_dir.as_posix()


def test_remote_repo_command_includes_version(tmp_path):
    params_path = tmp_path / "params.yaml"
    params_path.write_text("aligner: star\n", encoding="utf-8")
    source = {
        "source_kind": "remote_repo",
        "source_ref": "nf-core/scrnaseq",
        "resolved_version": "4.1.0",
    }
    cmd, cmd_str = build_nextflow_command(
        pipeline_source=source,
        profile="singularity",
        params_path=params_path,
        resume=False,
    )
    assert "nf-core/scrnaseq" in cmd
    assert "-r" in cmd
    assert "4.1.0" in cmd
    assert "singularity" in cmd


def test_resume_flag_appended(tmp_path):
    params_path = tmp_path / "params.yaml"
    params_path.write_text("aligner: star\n", encoding="utf-8")
    source = {
        "source_kind": "local_checkout",
        "source_ref": "/abs/path/to/scrnaseq",
        "resolved_version": "abc1234",
    }
    cmd, _cmd_str = build_nextflow_command(
        pipeline_source=source,
        profile="docker",
        params_path=params_path,
        resume=True,
    )
    assert "-resume" in cmd


def test_command_str_is_shell_safe(tmp_path):
    params_path = tmp_path / "params with spaces.yaml"
    params_path.write_text("aligner: star\n", encoding="utf-8")
    source = {
        "source_kind": "local_checkout",
        "source_ref": "/abs/path to/scrnaseq",
        "resolved_version": "abc1234",
    }
    _cmd, cmd_str = build_nextflow_command(
        pipeline_source=source,
        profile="docker",
        params_path=params_path,
        resume=False,
    )
    assert "'" in cmd_str or '"' in cmd_str or "\\ " in cmd_str


def test_extra_configs_added_to_command(tmp_path):
    params_path = tmp_path / "params.yaml"
    params_path.write_text("aligner: star\n", encoding="utf-8")
    cfg1 = tmp_path / "custom.config"
    cfg1.write_text("process { }\n", encoding="utf-8")
    cfg2 = tmp_path / "extra.config"
    cfg2.write_text("docker { enabled = true }\n", encoding="utf-8")
    source = {
        "source_kind": "local_checkout",
        "source_ref": "/abs/path/to/scrnaseq",
        "resolved_version": "abc1234",
    }
    cmd, cmd_str = build_nextflow_command(
        pipeline_source=source,
        profile="docker",
        params_path=params_path,
        resume=False,
        extra_configs=[cfg1, cfg2],
    )
    assert cmd.count("-c") == 2
    assert cfg1.as_posix() in cmd  # forward slashes
    assert cfg2.as_posix() in cmd


def test_no_extra_configs_by_default(tmp_path):
    params_path = tmp_path / "params.yaml"
    params_path.write_text("aligner: star\n", encoding="utf-8")
    source = {
        "source_kind": "local_checkout",
        "source_ref": "/abs/path/to/scrnaseq",
        "resolved_version": "abc1234",
    }
    cmd, _cmd_str = build_nextflow_command(
        pipeline_source=source,
        profile="docker",
        params_path=params_path,
        resume=False,
    )
    assert "-c" not in cmd


def test_resume_comes_after_configs(tmp_path):
    params_path = tmp_path / "params.yaml"
    params_path.write_text("aligner: star\n", encoding="utf-8")
    cfg = tmp_path / "custom.config"
    cfg.write_text("process { }\n", encoding="utf-8")
    source = {
        "source_kind": "local_checkout",
        "source_ref": "/abs/path/to/scrnaseq",
        "resolved_version": "abc1234",
    }
    cmd, _cmd_str = build_nextflow_command(
        pipeline_source=source,
        profile="docker",
        params_path=params_path,
        resume=True,
        extra_configs=[cfg],
    )
    assert cmd.index("-resume") > cmd.index("-c")


def test_command_uses_posix_paths(tmp_path):
    """All file paths in the command list must use forward slashes (no backslashes)."""
    params_path = tmp_path / "params.yaml"
    params_path.write_text("aligner: star\n", encoding="utf-8")
    cfg = tmp_path / "extra.config"
    cfg.write_text("process { }\n", encoding="utf-8")
    source = {
        "source_kind": "local_checkout",
        "source_ref": str(tmp_path / "scrnaseq"),
        "resolved_version": "abc1234",
    }
    cmd, _cmd_str = build_nextflow_command(
        pipeline_source=source,
        profile="docker",
        params_path=params_path,
        resume=False,
        extra_configs=[cfg],
    )
    for part in cmd:
        assert "\\" not in part, f"Backslash found in command part: {part!r}"


def test_local_checkout_source_ref_uses_posix_path(tmp_path):
    """Local checkout source_ref is converted to forward slashes in the command."""
    params_path = tmp_path / "params.yaml"
    params_path.write_text("aligner: star\n", encoding="utf-8")
    pipeline_dir = tmp_path / "scrnaseq"
    source = {
        "source_kind": "local_checkout",
        "source_ref": str(pipeline_dir),
        "resolved_version": "abc",
    }
    cmd, _ = build_nextflow_command(
        pipeline_source=source,
        profile="docker",
        params_path=params_path,
        resume=False,
    )
    assert cmd[2] == pipeline_dir.as_posix()
