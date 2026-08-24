#!/usr/bin/env python3
"""Approve queued CI runs for fork PRs that touch only `skills/`.

Issue #360. This repository requires approval for first-time contributors'
workflow runs. Those runs park in `action_required` and surface through
`gh pr checks` as "no checks reported", which reads exactly like a pass. Five
PRs sat in that state in a single week, and #329 went five remediation rounds
with CI having never run against it at all, so its 17 assertions executed
nowhere but the contributor's laptop.

GitHub offers no path-scoped approval policy: `approval_policy` is one
repository-wide value with three options, none path-aware. So the narrowing
lives here instead of in a setting.

WHAT THIS BOUNDARY IS WORTH, stated plainly so nobody mistakes it for more.
`skills/` is precisely where contributor code lives, and CI already executes
it: pytest imports a fork's test files, and the skill-harness job runs a
fork's shell. Auto-approving skills-only PRs therefore does NOT mean untrusted
code stops running unreviewed. What it does mean is that a PR which runs
untrusted code cannot in the same breath rewrite the pipeline that runs it
(`.github/`), the core package (`clawbio/`), the tooling (`scripts/`), or
dependency resolution (`pyproject.toml`, `uv.lock`). Anything touching those
still waits for a human.

The residual exposure is runner compute on a token with `contents: read` and
no secrets in any workflow. That is the trade this script encodes.

Usage:
    python scripts/approve_skills_only_runs.py            # act
    python scripts/approve_skills_only_runs.py --dry-run  # report only
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys

REPO = "ClawBio/ClawBio"

# The only prefix admitted. Trailing slash is load-bearing: a bare
# `startswith("skills")` would admit `skills-evil/payload.py`.
ALLOWED_PREFIX = "skills/"


def is_skills_only(files: list[str]) -> bool:
    """True only if every path is a plain file under `skills/`.

    Fails closed on an empty list: no files means the diff could not be read,
    which is not the same as a harmless PR.
    """
    if not files:
        return False

    for raw in files:
        path = (raw or "").strip()
        if not path:
            return False
        # Absolute paths, traversal, and backslash separators are all outside
        # what a well-formed git pathname looks like here. Refuse rather than
        # normalise: normalising is where the bypasses live.
        if path.startswith("/") or "\\" in path:
            return False
        if ".." in path.split("/"):
            return False
        if "." in path.split("/"):
            return False
        # Case-sensitive on purpose. A case-insensitive checkout could resolve
        # `SKILLS/` to the same directory, but guessing which filesystem is in
        # play is worse than refusing an unusual spelling.
        if not path.startswith(ALLOWED_PREFIX):
            return False
        if len(path) <= len(ALLOWED_PREFIX):
            return False
    return True


def pair_runs_to_prs(runs: list[dict], prs: list[dict]) -> list[tuple[int, int]]:
    """Pair each queued run with the open PR whose CURRENT head it is.

    GitHub leaves `pull_requests` empty on cross-repository runs, so the join
    is by head SHA. Matching the current head also drops superseded runs: all
    four pending on this repo when this was written were stale commits on one
    branch, and approving those spends runner minutes on code nobody is
    reviewing.
    """
    by_sha = {pr["headRefOid"]: pr["number"] for pr in prs if pr.get("headRefOid")}
    pairs = {
        (int(run["id"]), by_sha[run["head_sha"]])
        for run in runs
        if run.get("head_sha") in by_sha
    }
    return sorted(pairs)


def _gh_json(args: list[str]) -> object:
    result = subprocess.run(
        ["gh", *args], capture_output=True, text=True, check=False
    )
    if result.returncode != 0:
        raise RuntimeError(f"gh {' '.join(args)} failed: {result.stderr.strip()}")
    return json.loads(result.stdout or "null")


def queued_runs() -> list[dict]:
    data = _gh_json([
        "api", f"repos/{REPO}/actions/runs?status=action_required&per_page=50"
    ])
    return (data or {}).get("workflow_runs", [])


def open_prs() -> list[dict]:
    return _gh_json([
        "pr", "list", "--repo", REPO, "--state", "open",
        "--limit", "100", "--json", "number,headRefOid,author",
    ]) or []


def pr_files(number: int) -> list[str]:
    data = _gh_json(["api", f"repos/{REPO}/pulls/{number}/files?per_page=300"])
    return [f["filename"] for f in (data or [])]


def approve(run_id: int) -> None:
    subprocess.run(
        ["gh", "api", "-X", "POST",
         f"repos/{REPO}/actions/runs/{run_id}/approve"],
        check=True, capture_output=True, text=True,
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--dry-run", action="store_true",
        help="report what would be approved without approving it",
    )
    args = parser.parse_args()

    runs = queued_runs()
    if not runs:
        print("No runs awaiting approval.")
        return 0

    pairs = pair_runs_to_prs(runs, open_prs())
    stale = len(runs) - len({r for r, _ in pairs})
    if stale:
        print(f"Ignoring {stale} run(s) not at any open PR's current head.")

    approved = 0
    for run_id, pr_number in pairs:
        try:
            files = pr_files(pr_number)
        except RuntimeError as exc:
            print(f"  PR #{pr_number}: could not read files, leaving queued ({exc})")
            continue

        if not is_skills_only(files):
            outside = sorted({
                f for f in files if not f.startswith(ALLOWED_PREFIX)
            })[:5]
            print(f"  PR #{pr_number}: run {run_id} left for a human "
                  f"(touches {', '.join(outside) or 'unreadable paths'})")
            continue

        if args.dry_run:
            print(f"  PR #{pr_number}: would approve run {run_id} "
                  f"({len(files)} files, all under skills/)")
            continue

        try:
            approve(run_id)
        except subprocess.CalledProcessError as exc:
            print(f"  PR #{pr_number}: approve failed for run {run_id}: "
                  f"{exc.stderr.strip()}")
            continue
        approved += 1
        print(f"  PR #{pr_number}: approved run {run_id} "
              f"({len(files)} files, all under skills/)")

    print(f"{'Would approve' if args.dry_run else 'Approved'} {approved} run(s).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
