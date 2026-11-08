"""Command implementations: start, status, next, open, done."""
from __future__ import annotations

import hashlib
import sys
from pathlib import Path

from clawbio.workflow import pipeline as pl
from clawbio.workflow import state as st


class CommandError(RuntimeError):
    """Raised when a command should exit non-zero with a clean message."""


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return "sha256:" + h.hexdigest()


def _out_dir(cwd: Path | None = None) -> str:
    return str(st.run_dir(cwd) / "out")


# ---------------------------------------------------------------------------
# start
# ---------------------------------------------------------------------------

def cmd_start(skill_path: str, input_pairs: list[str]) -> None:
    pipeline = pl.load(Path(skill_path))
    inputs = _parse_kv(input_pairs)

    declared = {i.name: i for i in pipeline.inputs}
    for name in inputs:
        if name not in declared:
            raise CommandError(f"unknown input {name!r}; expected one of {sorted(declared)}")
    missing = [i.name for i in pipeline.inputs if i.required and i.name not in inputs]
    if missing:
        raise CommandError(f"missing required inputs: {missing}")

    step_defs = [
        {"id": s.id, "prereqs": s.prereqs, "requires": s.requires, "emits": s.emits}
        for s in pipeline.steps
    ]

    st.init(
        skill=pipeline.skill,
        skill_dir=str(pipeline.skill_dir),
        inputs=inputs,
        steps=step_defs,
    )
    Path(_out_dir()).mkdir(parents=True, exist_ok=True)
    print(f"started clawbio run for skill {pipeline.skill!r} with {len(pipeline.steps)} step(s)")
    cmd_next(quiet_if_none=True)


# ---------------------------------------------------------------------------
# status
# ---------------------------------------------------------------------------

def cmd_status(check_bundle: bool = False) -> None:
    if not st.exists():
        if check_bundle:
            return
        print("no active clawbio run in this directory")
        return

    state = st.read()
    if check_bundle:
        if state.get("bundled_at"):
            return
        failed = _failed_ids(state)
        if failed:
            # Non-blocking: a failed step needs a retry or a human, not a bundle.
            print(
                f"clawbio-flow: run has failed step(s): {', '.join(failed)} — "
                "retry with `clawbio-flow open <step>`; the run cannot be bundled until they succeed"
            )
            return
        all_done = all(e["status"] in (st.DONE, st.SKIPPED) for e in state["steps"].values())
        if all_done:
            print("active clawbio run is complete but `clawbio-flow bundle` has not been called")
            sys.exit(2)
        return

    print(f"skill:       {state['skill']}")
    print(f"started:     {state['started_at']}")
    print(f"bundled:     {state.get('bundled_at') or '—'}")
    print("steps:")
    for sid, entry in state["steps"].items():
        marker = {
            st.DONE: "[x]",
            st.IN_PROGRESS: "[~]",
            st.BLOCKED: "[ ]",
            st.PENDING: "[ ]",
            st.FAILED: "[!]",
            st.SKIPPED: "[-]",
        }.get(entry["status"], "[?]")
        extra = ""
        if entry["status"] == st.BLOCKED and entry.get("blocked_by"):
            extra = f"  (blocked by: {', '.join(entry['blocked_by'])})"
        elif entry["status"] == st.SKIPPED and entry.get("skip_reason"):
            extra = f"  (reason: {entry['skip_reason']})"
        print(f"  {marker} {sid:<24} {entry['status']}{extra}")

    nxt = _next_ready_id(state)
    print(f"next:        {nxt or '— (run `clawbio-flow bundle` if all steps are done)'}")


# ---------------------------------------------------------------------------
# next
# ---------------------------------------------------------------------------

def cmd_next(quiet_if_none: bool = False) -> None:
    state = st.read()
    pipeline = pl.load(Path(state["skill_dir"]))
    st.recompute_blocks(state, pipeline)
    st.write(state)

    sid = _next_ready_id(state)
    if sid is None:
        if quiet_if_none:
            return
        all_done = all(e["status"] in (st.DONE, st.SKIPPED) for e in state["steps"].values())
        if all_done:
            print("all steps done — run `clawbio-flow bundle`")
            return
        failed = _failed_ids(state)
        if failed:
            raise CommandError(
                f"no ready steps; failed steps: {failed} — retry with `clawbio-flow open <step>`"
            )
        raise CommandError("no ready steps (everything is blocked or in_progress)")

    summary = pipeline.step(sid).summary
    print(f"{sid}\t{summary}")


def _failed_ids(state: dict) -> list[str]:
    return [sid for sid, e in state["steps"].items() if e["status"] == st.FAILED]


def _next_ready_id(state: dict) -> str | None:
    for sid, entry in state["steps"].items():
        if entry["status"] == st.PENDING:
            return sid
    return None


# ---------------------------------------------------------------------------
# open
# ---------------------------------------------------------------------------

def cmd_open(step_id: str) -> None:
    with st.state_lock():
        state = st.read()
        pipeline = pl.load(Path(state["skill_dir"]))
        st.recompute_blocks(state, pipeline)

        if step_id not in state["steps"]:
            raise CommandError(f"unknown step {step_id!r}")

        entry = state["steps"][step_id]
        status = entry["status"]
        if status == st.DONE:
            raise CommandError(f"step {step_id!r} is already done")
        if status == st.IN_PROGRESS:
            raise CommandError(f"step {step_id!r} is already in_progress")
        if status == st.BLOCKED:
            raise CommandError(
                f"step {step_id!r} is blocked by: {', '.join(entry['blocked_by'])} "
                "— complete prerequisites first"
            )
        if status not in (st.PENDING, st.FAILED):
            raise CommandError(f"step {step_id!r} has unexpected status {status!r}")

        step = pipeline.step(step_id)
        inputs = state["inputs"]
        emitted = st.emitted_outputs(state)
        out_dir = _out_dir()
        Path(out_dir).mkdir(parents=True, exist_ok=True)
        rendered = pl.render_command(step, inputs, emitted, out_dir, state["skill_dir"])

        # Re-opening a failed step is a retry: clear the failed attempt's
        # fields. The attempt itself stays in the event log as a `fail` event.
        entry["status"] = st.IN_PROGRESS
        entry["started_at"] = st.now_iso()
        entry["command"] = rendered
        entry["exit_code"] = None
        entry["failed_at"] = None
        st.append_event(state, "open", step=step_id, command=rendered)
        st.write(state)

    print(f"# step:    {step.id}")
    print(f"# summary: {step.summary}")
    print(f"# emits:   {', '.join(step.emits) or '(none)'}")
    print(rendered)


# ---------------------------------------------------------------------------
# done
# ---------------------------------------------------------------------------

def cmd_done(step_id: str, output_pairs: list[str], exit_code: int = 0) -> None:
    outputs = _parse_kv(output_pairs)

    with st.state_lock():
        state = st.read()
        pipeline = pl.load(Path(state["skill_dir"]))

        if step_id not in state["steps"]:
            raise CommandError(f"unknown step {step_id!r}")

        entry = state["steps"][step_id]
        if entry["status"] != st.IN_PROGRESS:
            raise CommandError(
                f"step {step_id!r} is not in_progress (status={entry['status']!r}); "
                "call `clawbio-flow open` first"
            )

        step = pipeline.step(step_id)

        if exit_code != 0:
            extra = [k for k in outputs if k not in step.emits]
            if extra:
                raise CommandError(
                    f"step {step_id!r} produced --output keys not in emits: {extra} (declared: {step.emits})"
                )

            # A failing step often produces nothing: skip the emits-present and
            # is-a-file checks. Outputs that do exist are recorded on the `fail`
            # event for diagnosis only; they never become emitted outputs, so
            # downstream steps stay blocked.
            hashes, resolved = {}, {}
            for key, path_str in outputs.items():
                p = Path(path_str)
                if p.is_file():
                    hashes[key] = _sha256(p)
                    resolved[key] = str(p.resolve())

            entry["status"] = st.FAILED
            entry["failed_at"] = st.now_iso()
            entry["completed_at"] = None
            entry["outputs"] = {}
            entry["output_hashes"] = {}
            entry["exit_code"] = exit_code

            st.append_event(
                state, "fail", step=step_id, outputs=resolved, hashes=hashes, exit_code=exit_code
            )
            st.recompute_blocks(state, pipeline)
            st.write(state)

            print(f"step {step_id!r} marked failed (exit code {exit_code})")
            print(f"retry with `clawbio-flow open {step_id}`; downstream steps stay blocked")
            return

        missing_emits = [e for e in step.emits if e not in outputs]
        if missing_emits:
            raise CommandError(
                f"step {step_id!r} declares emits {step.emits} but --output missing: {missing_emits}"
            )
        extra = [k for k in outputs if k not in step.emits]
        if extra:
            raise CommandError(
                f"step {step_id!r} produced --output keys not in emits: {extra} (declared: {step.emits})"
            )

        hashes: dict[str, str] = {}
        resolved: dict[str, str] = {}
        for key, path_str in outputs.items():
            p = Path(path_str)
            if not p.is_file():
                raise CommandError(f"output {key}={path_str!r} is not a file")
            hashes[key] = _sha256(p)
            resolved[key] = str(p.resolve())

        entry["status"] = st.DONE
        entry["completed_at"] = st.now_iso()
        entry["outputs"] = resolved
        entry["output_hashes"] = hashes
        entry["exit_code"] = exit_code

        st.append_event(
            state, "done", step=step_id, outputs=resolved, hashes=hashes, exit_code=exit_code
        )
        st.recompute_blocks(state, pipeline)
        st.write(state)

    print(f"step {step_id!r} marked done")
    nxt = _next_ready_id(state)
    if nxt:
        print(f"next: {nxt}\t{pipeline.step(nxt).summary}")
    elif all(e["status"] in (st.DONE, st.SKIPPED) for e in state["steps"].values()):
        print("all steps done — run `clawbio-flow bundle`")


# ---------------------------------------------------------------------------
# skip
# ---------------------------------------------------------------------------

def cmd_skip(step_id: str, reason: str) -> None:
    with st.state_lock():
        state = st.read()
        pipeline = pl.load(Path(state["skill_dir"]))
        st.recompute_blocks(state, pipeline)

        if step_id not in state["steps"]:
            raise CommandError(f"unknown step {step_id!r}")

        step = pipeline.step(step_id)
        if not step.optional:
            raise CommandError(f"step {step_id!r} is not optional; cannot skip")

        entry = state["steps"][step_id]
        status = entry["status"]
        if status == st.BLOCKED:
            raise CommandError(
                f"step {step_id!r} is blocked by: {', '.join(entry['blocked_by'])} "
                "— resolve upstream before skipping"
            )
        if status != st.PENDING:
            raise CommandError(f"step {step_id!r} is not pending (status={status!r}); cannot skip")

        entry["status"] = st.SKIPPED
        entry["skipped_at"] = st.now_iso()
        entry["skip_reason"] = reason
        st.append_event(state, "skip", step=step_id, reason=reason)
        st.recompute_blocks(state, pipeline)
        st.write(state)

    print(f"step {step_id!r} marked skipped (reason: {reason})")
    nxt = _next_ready_id(state)
    if nxt:
        print(f"next: {nxt}\t{pipeline.step(nxt).summary}")
    elif all(e["status"] in (st.DONE, st.SKIPPED) for e in state["steps"].values()):
        print("all steps done — run `clawbio-flow bundle`")


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def _parse_kv(pairs: list[str]) -> dict[str, str]:
    out: dict[str, str] = {}
    for pair in pairs:
        if "=" not in pair:
            raise CommandError(f"expected key=value, got {pair!r}")
        k, v = pair.split("=", 1)
        if not k:
            raise CommandError(f"empty key in {pair!r}")
        out[k] = v
    return out
