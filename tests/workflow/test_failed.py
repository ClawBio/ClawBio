"""Tests for `clawbio-flow done --exit-code N` with N != 0.

A nonzero exit code marks the step `failed` (not done). A failed step:
- appends a `fail` event and no `done` event,
- needs no outputs and contributes nothing to emitted outputs,
- keeps downstream steps blocked,
- can be retried with `open`, which moves it back to `in_progress`,
- makes `bundle` refuse the run, and
- makes the Stop hook (`status --check-bundle`) print a notice without blocking.

Runs the real commands against the variant_annotation fixture skill.
"""
from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from clawbio.workflow import bundle, commands
from clawbio.workflow import state as st

VARIANT_ANNOTATION = Path(__file__).resolve().parent / "fixtures" / "variant_annotation"


def _start(variant_inputs: dict[str, str]) -> None:
    commands.cmd_start(
        str(VARIANT_ANNOTATION),
        [f"vcf_path={variant_inputs['vcf_path']}", f"reference={variant_inputs['reference']}"],
    )


def _run_rendered(step_id: str) -> None:
    rendered = st.read()["steps"][step_id]["command"]
    subprocess.run(rendered, shell=True, check=True, executable="/bin/bash")


def _out(name: str) -> str:
    return str(Path(st.run_dir()) / "out" / name)


def _fail_normalize() -> None:
    commands.cmd_open("normalize")
    commands.cmd_done("normalize", [], exit_code=1)


# ---------------------------------------------------------------------------
# (a) nonzero exit code marks the step failed
# ---------------------------------------------------------------------------

def test_nonzero_exit_code_marks_step_failed(run_cwd: Path, variant_inputs, capsys):
    _start(variant_inputs)
    commands.cmd_open("normalize")
    capsys.readouterr()

    # No --output at all: a failing step often produces nothing.
    commands.cmd_done("normalize", [], exit_code=1)

    out = capsys.readouterr().out
    assert "step 'normalize' marked failed" in out
    assert "marked done" not in out

    state = st.read()
    entry = state["steps"]["normalize"]
    assert entry["status"] == st.FAILED
    assert entry["exit_code"] == 1
    assert entry["failed_at"] is not None
    assert entry["outputs"] == {}
    assert entry["output_hashes"] == {}

    kinds = [e["kind"] for e in state["events"] if e.get("step") == "normalize"]
    assert kinds == ["open", "fail"]
    fail_event = [e for e in state["events"] if e["kind"] == "fail"][0]
    assert fail_event["exit_code"] == 1

    assert st.emitted_outputs(state) == {}
    assert state["steps"]["annotate"]["status"] == st.BLOCKED
    assert "normalize" in state["steps"]["annotate"]["blocked_by"]


def test_nonzero_exit_code_accepts_nonexistent_outputs(run_cwd: Path, variant_inputs):
    _start(variant_inputs)
    commands.cmd_open("normalize")

    missing = _out("normalized.vcf")
    assert not Path(missing).exists()
    commands.cmd_done("normalize", [f"normalized_vcf={missing}"], exit_code=3)

    state = st.read()
    entry = state["steps"]["normalize"]
    assert entry["status"] == st.FAILED
    assert entry["exit_code"] == 3
    assert st.emitted_outputs(state) == {}
    assert "done" not in [e["kind"] for e in state["events"]]
    assert state["steps"]["annotate"]["status"] == st.BLOCKED


def test_nonzero_exit_code_partial_outputs_never_emitted(run_cwd: Path, variant_inputs):
    """An output file that does exist is still not made available downstream."""
    _start(variant_inputs)
    commands.cmd_open("normalize")
    _run_rendered("normalize")  # produces normalized.vcf
    produced = _out("normalized.vcf")
    assert Path(produced).is_file()

    commands.cmd_done("normalize", [f"normalized_vcf={produced}"], exit_code=1)

    state = st.read()
    assert state["steps"]["normalize"]["status"] == st.FAILED
    assert st.emitted_outputs(state) == {}
    assert state["steps"]["annotate"]["status"] == st.BLOCKED
    # recorded on the fail event for diagnosis only
    fail_event = [e for e in state["events"] if e["kind"] == "fail"][0]
    assert fail_event["outputs"] == {"normalized_vcf": str(Path(produced).resolve())}

    # the downstream step cannot be opened
    with pytest.raises(commands.CommandError, match="blocked by"):
        commands.cmd_open("annotate")


# ---------------------------------------------------------------------------
# (b) retry: open a failed step again
# ---------------------------------------------------------------------------

def test_open_retries_failed_step_and_bundle_lists_it_once(run_cwd: Path, variant_inputs, capsys):
    _start(variant_inputs)
    _fail_normalize()
    capsys.readouterr()

    commands.cmd_open("normalize")
    out = capsys.readouterr().out
    state = st.read()
    entry = state["steps"]["normalize"]
    assert entry["status"] == st.IN_PROGRESS
    assert entry["command"] in out  # command re-rendered and printed
    assert entry["exit_code"] is None
    assert entry["failed_at"] is None

    _run_rendered("normalize")
    commands.cmd_done("normalize", [f"normalized_vcf={_out('normalized.vcf')}"])
    assert st.read()["steps"]["normalize"]["status"] == st.DONE
    assert st.read()["steps"]["annotate"]["status"] == st.PENDING

    commands.cmd_open("annotate")
    _run_rendered("annotate")
    commands.cmd_done("annotate", [f"annotated_vcf={_out('annotated.vcf')}"])

    bundle.cmd_bundle()
    commands_sh = (st.run_dir() / "commands.sh").read_text()
    assert commands_sh.count("# step: normalize") == 1
    assert commands_sh.count("# step: annotate") == 1

    kinds = [e["kind"] for e in st.read()["events"] if e.get("step") == "normalize"]
    assert kinds == ["open", "fail", "open", "done"]


def test_open_still_rejects_done_step(run_cwd: Path, variant_inputs):
    _start(variant_inputs)
    commands.cmd_open("normalize")
    _run_rendered("normalize")
    commands.cmd_done("normalize", [f"normalized_vcf={_out('normalized.vcf')}"])
    with pytest.raises(commands.CommandError, match="already done"):
        commands.cmd_open("normalize")


# ---------------------------------------------------------------------------
# (c) bundle refuses a run with failed steps
# ---------------------------------------------------------------------------

def test_bundle_refuses_failed_steps(run_cwd: Path, variant_inputs):
    _start(variant_inputs)
    _fail_normalize()

    with pytest.raises(RuntimeError, match=r"failed steps: \['normalize'\]"):
        bundle.cmd_bundle()

    run = st.run_dir()
    assert not (run / "checksums.sha256").exists()
    assert not (run / "commands.sh").exists()
    assert st.read()["bundled_at"] is None


# ---------------------------------------------------------------------------
# (d) Stop hook is non-blocking on failed steps
# ---------------------------------------------------------------------------

def test_check_bundle_with_failed_step_is_nonblocking_notice(run_cwd: Path, variant_inputs, capsys):
    _start(variant_inputs)
    _fail_normalize()
    capsys.readouterr()

    # Must return normally: no SystemExit (exit 0 for the Stop hook).
    commands.cmd_status(check_bundle=True)

    out = capsys.readouterr().out
    assert "failed step(s): normalize" in out


def test_cli_check_bundle_with_failed_step_exits_zero(run_cwd: Path, variant_inputs, capsys):
    from clawbio.workflow.__main__ import main

    _start(variant_inputs)
    _fail_normalize()
    capsys.readouterr()

    assert main(["status", "--check-bundle"]) == 0
    assert "normalize" in capsys.readouterr().out


def test_status_shows_failed_marker(run_cwd: Path, variant_inputs, capsys):
    _start(variant_inputs)
    _fail_normalize()
    capsys.readouterr()

    commands.cmd_status()
    out = capsys.readouterr().out
    assert "[!] normalize" in out


def test_next_names_failed_step_when_nothing_ready(run_cwd: Path, variant_inputs):
    _start(variant_inputs)
    _fail_normalize()
    with pytest.raises(commands.CommandError, match=r"failed steps: \['normalize'\]"):
        commands.cmd_next()
