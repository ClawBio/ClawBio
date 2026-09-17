#!/usr/bin/env python3
"""Diff MAINTAINERS.md against the live GitHub state.

MAINTAINERS.md is a claim; the API is the fact. GOVERNANCE.md says the file is
read from the API rather than written from memory, and this is the script that
makes that true on demand. Run before each release.

Exit 0 when they agree, 1 when an account with access is not declared, 2 when
the API could not be read (INCOMPLETE, never a pass).

Usage: python scripts/check_maintainers.py
"""
from __future__ import annotations

import json
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ORG = "ClawBio"
REPO = "ClawBio/ClawBio"

# Handles appear as markdown profile links: [@handle](https://github.com/handle)
_HANDLE = re.compile(r"\[@([A-Za-z0-9-]+)\]\(https://github\.com/[A-Za-z0-9-]+\)")


def _run(args: list[str]) -> str | None:
    try:
        out = subprocess.run(args, capture_output=True, text=True, timeout=60)
    except (OSError, subprocess.TimeoutExpired):
        return None
    return out.stdout if out.returncode == 0 else None


def _gh(path: str) -> list | dict | None:
    """Read every page. `gh api` stops at 30 items without --paginate, and a
    truncated read here would print OK while missing the accounts that matter.
    --slurp wraps the pages into one array so this parses as a single document.
    """
    raw = _run(["gh", "api", "--paginate", "--slurp", path])
    if raw is None:
        return None
    try:
        return json.loads(raw)
    except json.JSONDecodeError:
        return None


def _logins(payload: list | dict) -> set[str]:
    """Logins from a --slurp response (list of pages) or a plain list."""
    if isinstance(payload, dict):
        payload = [payload]
    logins: set[str] = set()
    for item in payload:
        if isinstance(item, list):
            logins |= {entry["login"] for entry in item}
        else:
            logins.add(item["login"])
    return logins


def declared_handles(text: str) -> set[str]:
    return set(_HANDLE.findall(text))


def main() -> int:
    text = (ROOT / "MAINTAINERS.md").read_text()
    declared = declared_handles(text)

    members = _gh(f"orgs/{ORG}/members")
    collabs = _gh(f"repos/{REPO}/collaborators?affiliation=direct")
    if members is None or collabs is None:
        print("INCOMPLETE: could not read the GitHub API (gh missing, offline, or unauthenticated)")
        return 2
    if not _logins(members):
        print(
            f"INCOMPLETE: orgs/{ORG}/members came back empty. An organisation always has "
            "an owner, so this token cannot see its members (it needs read:org and "
            "membership of the org). Treating this as unread rather than as agreement."
        )
        return 2

    live = _logins(members) | _logins(collabs)
    admins = sorted(
        c["login"]
        for page in (collabs if collabs and isinstance(collabs[0], list) else [collabs])
        for c in page
        if c.get("permissions", {}).get("admin")
    )

    undeclared = sorted(live - declared)
    # No account is exempt. The lead maintainer's row going stale is the one
    # worth hearing about, and this is a warning, not a failure.
    stale = sorted(declared - live)
    print(f"live accounts : {sorted(live)}")
    print(f"direct admins : {admins}")
    print(f"declared      : {sorted(declared)}")
    if undeclared:
        print(f"FAIL: accounts with access not listed in MAINTAINERS.md: {undeclared}")
    if stale:
        print(f"WARN: handles listed without live access (stale row?): {stale}")
    if not undeclared:
        print("OK: every account with access is declared")
    return 1 if undeclared else 0


if __name__ == "__main__":
    sys.exit(main())
