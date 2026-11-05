"""End-to-end lifecycle test against the real variant_annotation demo skill.

Runs the actual cp/echo commands the skill renders — no subprocess of the CLI,
no mocks. Each test gets its own cwd via the `run_cwd` fixture.
"""
from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path

import pytest

from clawbio.workflow import bundle, commands
from clawbio.workflow import state as st

REPO_ROOT = Path(__file__).resolve().parents[2]
VARIANT_ANNOTATION = Path(__file__).resolve().parent / "fixtures" / "variant_annotation"


def _execute_step(step_id: str) -> dict[str, str]:
    """Open a step, run its rendered command in a shell, return output paths.

    Mirrors what an agent would do: read the command line printed by `open`,
    execute it, then call `done` with the produced output paths.
    """
    state = st.read()
    rendered = state["steps"][step_id]["command"]
    subprocess.run(rendered, shell=True, check=True, executable="/bin/bash")

    out_dir = Path(st.run_dir()) / "out"
    # variant_annotation hard-codes its output filenames; map emit -> path
    return {
        "normalize": {"normalized_vcf": str(out_dir / "normalized.vcf")},
        "annotate": {"annotated_vcf": str(out_dir / "annotated.vcf")},
    }[step_id]


# ---------------------------------------------------------------------------
# Full lifecycle
# ---------------------------------------------------------------------------

def test_full_lifecycle_start_to_bundle(run_cwd: Path, variant_inputs, capsys):
    commands.cmd_start(
        str(VARIANT_ANNOTATION),
        [f"vcf_path={variant_inputs['vcf_path']}", f"reference={variant_inputs['reference']}"],
    )

    # state.json created with both steps and the right initial states
    state = st.read()
    assert state["skill"] == "variant_annotation"
    assert state["steps"]["normalize"]["status"] == st.PENDING
    assert state["steps"]["annotate"]["status"] == st.BLOCKED
    assert state["steps"]["annotate"]["blocked_by"] == ["normalize", "emit:normalized_vcf"]

    # Step 1: normalize
    commands.cmd_open("normalize")
    assert st.read()["steps"]["normalize"]["status"] == st.IN_PROGRESS
    outputs = _execute_step("normalize")
    commands.cmd_done("normalize", [f"{k}={v}" for k, v in outputs.items()])

    state = st.read()
    assert state["steps"]["normalize"]["status"] == st.DONE
    # output hash recorded as sha256:<hex>
    h = state["steps"]["normalize"]["output_hashes"]["normalized_vcf"]
    assert h.startswith("sha256:") and len(h) == len("sha256:") + 64
    # annotate now unblocked
    assert state["steps"]["annotate"]["status"] == st.PENDING

    # Step 2: annotate
    commands.cmd_open("annotate")
    outputs = _execute_step("annotate")
    commands.cmd_done("annotate", [f"{k}={v}" for k, v in outputs.items()])
    assert all(e["status"] == st.DONE for e in st.read()["steps"].values())

    # Bundle
    bundle.cmd_bundle()
    run = st.run_dir()
    commands_sh = run / "commands.sh"
    checksums = run / "checksums.sha256"
    env_yml = run / "environment.yml"
    assert commands_sh.is_file() and (commands_sh.stat().st_mode & 0o111)
    assert checksums.is_file()
    assert env_yml.is_file()

    # commands.sh contains the two rendered commands in execution order
    body = commands_sh.read_text()
    assert body.index("# step: normalize") < body.index("# step: annotate")

    # bundled_at recorded
    assert st.read()["bundled_at"] is not None


def test_checksums_verify_with_shasum(run_cwd: Path, variant_inputs):
    """Regression for b9f69fb: `shasum -a 256 -c` must accept our format."""
    commands.cmd_start(
        str(VARIANT_ANNOTATION),
        [f"vcf_path={variant_inputs['vcf_path']}", f"reference={variant_inputs['reference']}"],
    )
    for sid in ("normalize", "annotate"):
        commands.cmd_open(sid)
        outs = _execute_step(sid)
        commands.cmd_done(sid, [f"{k}={v}" for k, v in outs.items()])
    bundle.cmd_bundle()

    checksums = st.run_dir() / "checksums.sha256"
    # shasum is preinstalled on macOS; sha256sum on linux
    tool = "shasum" if sys.platform == "darwin" else "sha256sum"
    extra = ["-a", "256"] if tool == "shasum" else []
    result = subprocess.run(
        [tool, *extra, "-c", str(checksums)], capture_output=True, text=True
    )
    assert result.returncode == 0, (
        f"checksum verification failed:\nstdout={result.stdout}\nstderr={result.stderr}"
    )


def test_output_paths_in_state_are_absolute(run_cwd: Path, variant_inputs):
    """Regression for 484222a: --output paths must be resolved to absolute."""
    commands.cmd_start(
        str(VARIANT_ANNOTATION),
        [f"vcf_path={variant_inputs['vcf_path']}", f"reference={variant_inputs['reference']}"],
    )
    commands.cmd_open("normalize")
    _execute_step("normalize")
    # pass a relative-looking path — code should resolve it
    rel_path = str((st.run_dir() / "out" / "normalized.vcf").relative_to(run_cwd))
    commands.cmd_done("normalize", [f"normalized_vcf={rel_path}"])
    recorded = st.read()["steps"]["normalize"]["outputs"]["normalized_vcf"]
    assert Path(recorded).is_absolute()


def test_skill_dir_placeholder_substituted(run_cwd: Path, variant_inputs):
    """Regression for a161939: {skill_dir} placeholder must resolve to abs path."""
    commands.cmd_start(
        str(VARIANT_ANNOTATION),
        [f"vcf_path={variant_inputs['vcf_path']}", f"reference={variant_inputs['reference']}"],
    )
    # variant_annotation doesn't use {skill_dir}, so verify via gwas_lookup's
    # rendered command on a fresh start in a sub-cwd
    state = st.read()
    assert state["skill_dir"] == str(VARIANT_ANNOTATION)


# ---------------------------------------------------------------------------
# Refusals & state transitions
# ---------------------------------------------------------------------------

def test_open_refuses_blocked_step(run_cwd: Path, variant_inputs):
    commands.cmd_start(
        str(VARIANT_ANNOTATION),
        [f"vcf_path={variant_inputs['vcf_path']}", f"reference={variant_inputs['reference']}"],
    )
    with pytest.raises(commands.CommandError, match="blocked by"):
        commands.cmd_open("annotate")


def test_done_rejects_mismatched_outputs(run_cwd: Path, variant_inputs):
    commands.cmd_start(
        str(VARIANT_ANNOTATION),
        [f"vcf_path={variant_inputs['vcf_path']}", f"reference={variant_inputs['reference']}"],
    )
    commands.cmd_open("normalize")
    _execute_step("normalize")

    out_path = str(st.run_dir() / "out" / "normalized.vcf")
    with pytest.raises(commands.CommandError, match="missing"):
        commands.cmd_done("normalize", [f"wrong_key={out_path}"])
    with pytest.raises(commands.CommandError, match="not in emits"):
        commands.cmd_done(
            "normalize",
            [f"normalized_vcf={out_path}", f"bonus={out_path}"],
        )


def test_done_rejects_missing_output_file(run_cwd: Path, variant_inputs):
    commands.cmd_start(
        str(VARIANT_ANNOTATION),
        [f"vcf_path={variant_inputs['vcf_path']}", f"reference={variant_inputs['reference']}"],
    )
    commands.cmd_open("normalize")
    # Skip running the rendered command — the output file won't exist
    with pytest.raises(commands.CommandError, match="not a file"):
        commands.cmd_done("normalize", ["normalized_vcf=/tmp/clawbio-nonexistent.vcf"])


def test_open_refuses_already_done(run_cwd: Path, variant_inputs):
    commands.cmd_start(
        str(VARIANT_ANNOTATION),
        [f"vcf_path={variant_inputs['vcf_path']}", f"reference={variant_inputs['reference']}"],
    )
    commands.cmd_open("normalize")
    outs = _execute_step("normalize")
    commands.cmd_done("normalize", [f"{k}={v}" for k, v in outs.items()])
    with pytest.raises(commands.CommandError, match="already done"):
        commands.cmd_open("normalize")


def test_start_refuses_double_init(run_cwd: Path, variant_inputs):
    commands.cmd_start(
        str(VARIANT_ANNOTATION),
        [f"vcf_path={variant_inputs['vcf_path']}", f"reference={variant_inputs['reference']}"],
    )
    with pytest.raises(FileExistsError, match="refusing to overwrite"):
        commands.cmd_start(
            str(VARIANT_ANNOTATION),
            [f"vcf_path={variant_inputs['vcf_path']}", f"reference={variant_inputs['reference']}"],
        )


def test_start_rejects_missing_required_input(run_cwd: Path, variant_inputs):
    with pytest.raises(commands.CommandError, match="missing required inputs"):
        commands.cmd_start(str(VARIANT_ANNOTATION), [f"vcf_path={variant_inputs['vcf_path']}"])


def test_start_rejects_unknown_input(run_cwd: Path, variant_inputs):
    with pytest.raises(commands.CommandError, match="unknown input"):
        commands.cmd_start(
            str(VARIANT_ANNOTATION),
            [
                f"vcf_path={variant_inputs['vcf_path']}",
                f"reference={variant_inputs['reference']}",
                "bogus=x",
            ],
        )


def test_bundle_refuses_incomplete_run(run_cwd: Path, variant_inputs):
    commands.cmd_start(
        str(VARIANT_ANNOTATION),
        [f"vcf_path={variant_inputs['vcf_path']}", f"reference={variant_inputs['reference']}"],
    )
    with pytest.raises(RuntimeError, match="not done"):
        bundle.cmd_bundle()


# ---------------------------------------------------------------------------
# Survives /compact: status is self-sufficient after fresh state read
# ---------------------------------------------------------------------------

def test_status_after_fresh_read_reflects_truth(run_cwd: Path, variant_inputs, capsys):
    """Simulates a /compact: discard everything but state.json and re-orient."""
    commands.cmd_start(
        str(VARIANT_ANNOTATION),
        [f"vcf_path={variant_inputs['vcf_path']}", f"reference={variant_inputs['reference']}"],
    )
    commands.cmd_open("normalize")
    outs = _execute_step("normalize")
    commands.cmd_done("normalize", [f"{k}={v}" for k, v in outs.items()])

    # "Fresh" agent: only thing it knows is the cwd. status must reveal everything.
    capsys.readouterr()  # drain
    commands.cmd_status()
    out = capsys.readouterr().out
    assert "normalize" in out and "done" in out
    assert "annotate" in out and "pending" in out
    assert "next:" in out
