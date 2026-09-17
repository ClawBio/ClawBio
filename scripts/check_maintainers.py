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


def _gh(path: str) -> list | dict | None:
    try:
        out = subprocess.run(["gh", "api", path], capture_output=True, text=True, timeout=30)
    except (OSError, subprocess.TimeoutExpired):
        return None
    if out.returncode != 0:
        return None
    try:
        return json.loads(out.stdout)
    except json.JSONDecodeError:
        return None


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

    live = {m["login"] for m in members} | {c["login"] for c in collabs}
    admins = sorted(c["login"] for c in collabs if c.get("permissions", {}).get("admin"))

    undeclared = sorted(live - declared)
    stale = sorted(declared - live - {"manuelcorpas"})
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
