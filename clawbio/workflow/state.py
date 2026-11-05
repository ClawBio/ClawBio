"""state.json read/write helpers. Atomic writes; append-only events log."""
from __future__ import annotations

import contextlib
import datetime as dt
import fcntl
import json
import os
import tempfile
from pathlib import Path
from typing import Any

SCHEMA_VERSION = 1
RUN_DIR_NAME = ".clawbio-run"
STATE_FILE = "state.json"

# Step statuses
PENDING = "pending"
BLOCKED = "blocked"
IN_PROGRESS = "in_progress"
DONE = "done"
FAILED = "failed"
SKIPPED = "skipped"


def now_iso() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z")


def run_dir(cwd: Path | None = None) -> Path:
    return (Path(cwd) if cwd else Path.cwd()) / RUN_DIR_NAME


def state_path(cwd: Path | None = None) -> Path:
    return run_dir(cwd) / STATE_FILE


def exists(cwd: Path | None = None) -> bool:
    return state_path(cwd).is_file()


def read(cwd: Path | None = None) -> dict[str, Any]:
    p = state_path(cwd)
    if not p.is_file():
        raise FileNotFoundError(f"no active clawbio run: {p} not found (run `clawbio-flow start` first)")
    with p.open() as f:
        return json.load(f)


@contextlib.contextmanager
def state_lock(cwd: Path | None = None):
    """Exclusive flock on a sibling .clawbio-run/state.lock file.

    Wrap the entire read-modify-write window of cmd_open/cmd_done so two
    parallel processes cannot lose updates. Lockfile is separate from
    state.json because os.replace swaps inodes; locking state.json itself
    would race the atomic write.
    """
    lock_path = run_dir(cwd) / "state.lock"
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    with lock_path.open("a+") as f:
        fcntl.flock(f.fileno(), fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(f.fileno(), fcntl.LOCK_UN)


def write(state: dict[str, Any], cwd: Path | None = None) -> None:
    """Atomic write: tmpfile in same dir + os.replace."""
    p = state_path(cwd)
    p.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w", dir=p.parent, prefix=".state-", suffix=".json.tmp", delete=False
    ) as tf:
        json.dump(state, tf, indent=2, sort_keys=False)
        tf.flush()
        os.fsync(tf.fileno())
        tmp_name = tf.name
    os.replace(tmp_name, p)


def init(
    skill: str,
    skill_dir: str,
    inputs: dict[str, str],
    steps: list[dict[str, Any]],
    cwd: Path | None = None,
) -> dict[str, Any]:
    """Create a fresh state.json. Step entries already prepared by caller."""
    if exists(cwd):
        raise FileExistsError(
            f"{state_path(cwd)} already exists — refusing to overwrite "
            "(remove .clawbio-run/ to start a new run)"
        )
    state = {
        "schema_version": SCHEMA_VERSION,
        "skill": skill,
        "skill_dir": skill_dir,
        "started_at": now_iso(),
        "inputs": inputs,
        "steps": {s["id"]: _new_step_entry(s) for s in steps},
        "bundled_at": None,
        "events": [{"ts": now_iso(), "kind": "start", "skill": skill, "inputs": inputs}],
    }
    _recompute_blocks(state, steps)
    write(state, cwd)
    return state


def _new_step_entry(step_def: dict[str, Any]) -> dict[str, Any]:
    return {
        "status": PENDING,
        "prereqs": list(step_def.get("prereqs", [])),
        "requires": dict(step_def.get("requires", {})),
        "emits": list(step_def.get("emits", [])),
        "optional": bool(step_def.get("optional", False)),
        "blocked_by": [],
        "started_at": None,
        "completed_at": None,
        "skipped_at": None,
        "skip_reason": None,
        "command": None,
        "outputs": {},
        "output_hashes": {},
        "exit_code": None,
    }


def _recompute_blocks(state: dict[str, Any], step_defs: list[dict[str, Any]]) -> None:
    """Mark steps blocked if any prereq is not done/skipped or any required emit is unsatisfied.

    A skipped step satisfies downstream `prereqs` (DAG ordering — it ran its
    course, just chose not to execute) but never satisfies `requires` emit
    gates (it produced no output). This is the "(a) strict by default" model
    from optional-steps-proposal.md: a downstream step that gates on a
    skipped step's emit blocks forever — the agent must skip it too.
    """
    all_emitted: dict[str, str] = {}
    for sid, entry in state["steps"].items():
        if entry["status"] == DONE:
            all_emitted.update(entry["outputs"])

    for sdef in step_defs:
        sid = sdef["id"]
        entry = state["steps"][sid]
        if entry["status"] in (DONE, IN_PROGRESS, FAILED, SKIPPED):
            continue
        missing_prereqs = [
            p for p in entry["prereqs"] if state["steps"][p]["status"] not in (DONE, SKIPPED)
        ]
        missing_requires = [k for k in entry["requires"] if k not in all_emitted]
        if missing_prereqs or missing_requires:
            entry["status"] = BLOCKED
            entry["blocked_by"] = missing_prereqs + [f"emit:{k}" for k in missing_requires]
        else:
            entry["status"] = PENDING
            entry["blocked_by"] = []


def recompute_blocks(state: dict[str, Any], pipeline) -> None:
    """Public wrapper accepting a Pipeline object."""
    step_defs = [
        {"id": s.id, "prereqs": s.prereqs, "requires": s.requires, "emits": s.emits}
        for s in pipeline.steps
    ]
    _recompute_blocks(state, step_defs)


def append_event(state: dict[str, Any], kind: str, **payload: Any) -> None:
    event = {"ts": now_iso(), "kind": kind, **payload}
    state.setdefault("events", []).append(event)


def emitted_outputs(state: dict[str, Any]) -> dict[str, str]:
    out: dict[str, str] = {}
    for entry in state["steps"].values():
        if entry["status"] == DONE:
            out.update(entry["outputs"])
    return out
