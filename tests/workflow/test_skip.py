"""Tests for `optional: true` steps and `clawbio-flow skip`.

Covers the five cases from optional-steps-proposal.md: skipping an optional
step, refusing to skip a required step, refusing to skip a blocked step,
bundle accepting a mixed done/skipped run, and a downstream step that gates
on a skipped emit ending up permanently blocked.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from clawbio.workflow import bundle, commands
from clawbio.workflow import state as st


def _finish_step(step_id: str, emits: dict[str, str]) -> None:
    """Open a step, fabricate its declared output files, mark it done."""
    commands.cmd_open(step_id)
    out_dir = Path(st.run_dir()) / "out"
    out_dir.mkdir(parents=True, exist_ok=True)
    pairs = []
    for key, filename in emits.items():
        p = out_dir / filename
        p.write_text("x")
        pairs.append(f"{key}={p}")
    commands.cmd_done(step_id, pairs)


def _ordering_only_pipeline(run_cwd: Path) -> Path:
    """a -> b (optional) -> c, where c only cares about DAG order, not b's emit."""
    skill_dir = run_cwd / "src_ordering"
    skill_dir.mkdir()
    (skill_dir / "pipeline.yaml").write_text("""
skill: skip_ordering_demo
version: 1
steps:
  - {id: a, summary: "produce a", command_template: "true", emits: [a_out]}
  - id: b
    summary: "optional step"
    optional: true
    prereqs: [a]
    requires: {a_out: present}
    command_template: "true"
    emits: [b_out]
  - id: c
    summary: "ordering-only downstream of b"
    prereqs: [a, b]
    command_template: "true"
""")
    return skill_dir


def _emit_gated_pipeline(run_cwd: Path) -> Path:
    """a -> b (optional) -> d, where d gates on b's emit."""
    skill_dir = run_cwd / "src_emit_gated"
    skill_dir.mkdir()
    (skill_dir / "pipeline.yaml").write_text("""
skill: skip_emit_gated_demo
version: 1
steps:
  - {id: a, summary: "produce a", command_template: "true", emits: [a_out]}
  - id: b
    summary: "optional step"
    optional: true
    prereqs: [a]
    requires: {a_out: present}
    command_template: "true"
    emits: [b_out]
  - id: d
    summary: "gates on b's emit"
    prereqs: [a, b]
    requires: {b_out: present}
    command_template: "true"
""")
    return skill_dir


def test_skip_marks_skipped_and_unblocks_ordering_only_downstream(run_cwd: Path):
    skill_dir = _ordering_only_pipeline(run_cwd)
    commands.cmd_start(str(skill_dir), [])
    _finish_step("a", {"a_out": "a.txt"})

    assert st.read()["steps"]["b"]["status"] == st.PENDING
    commands.cmd_skip("b", "single batch, nothing to correct for")

    state = st.read()
    assert state["steps"]["b"]["status"] == st.SKIPPED
    assert state["steps"]["b"]["skip_reason"] == "single batch, nothing to correct for"
    assert state["steps"]["b"]["skipped_at"] is not None
    # c only orders after b — it should unblock even though b never emitted
    assert state["steps"]["c"]["status"] == st.PENDING


def test_downstream_gated_on_skipped_emit_blocks_forever(run_cwd: Path):
    skill_dir = _emit_gated_pipeline(run_cwd)
    commands.cmd_start(str(skill_dir), [])
    _finish_step("a", {"a_out": "a.txt"})
    commands.cmd_skip("b", "no doublets expected")

    state = st.read()
    assert state["steps"]["d"]["status"] == st.BLOCKED
    assert "emit:b_out" in state["steps"]["d"]["blocked_by"]


def test_skip_refuses_non_optional_step(run_cwd: Path):
    skill_dir = _ordering_only_pipeline(run_cwd)
    commands.cmd_start(str(skill_dir), [])
    with pytest.raises(commands.CommandError, match="not optional"):
        commands.cmd_skip("a", "trying to skip a required step")


def test_skip_refuses_blocked_step(run_cwd: Path):
    skill_dir = _ordering_only_pipeline(run_cwd)
    commands.cmd_start(str(skill_dir), [])
    # b is blocked: its prereq `a` hasn't run yet
    with pytest.raises(commands.CommandError, match="blocked by"):
        commands.cmd_skip("b", "too early")


def test_skip_refuses_already_skipped_step(run_cwd: Path):
    skill_dir = _ordering_only_pipeline(run_cwd)
    commands.cmd_start(str(skill_dir), [])
    _finish_step("a", {"a_out": "a.txt"})
    commands.cmd_skip("b", "first skip")
    with pytest.raises(commands.CommandError, match="not pending"):
        commands.cmd_skip("b", "second skip")


def test_bundle_accepts_mixed_done_and_skipped(run_cwd: Path):
    skill_dir = _ordering_only_pipeline(run_cwd)
    commands.cmd_start(str(skill_dir), [])
    _finish_step("a", {"a_out": "a.txt"})
    commands.cmd_skip("b", "single batch")
    _finish_step("c", {})

    bundle.cmd_bundle()
    assert st.read()["bundled_at"] is not None

    commands_sh = (st.run_dir() / "commands.sh").read_text()
    assert "# step: a" in commands_sh
    assert "# step: c" in commands_sh
    assert "# step: b" not in commands_sh  # skipped steps aren't invoked
    assert "# skipped: b (reason: single batch)" in commands_sh
