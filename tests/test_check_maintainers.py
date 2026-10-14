"""The script that proves MAINTAINERS.md matches the API must not pass on a
truncated read: `gh api` stops at 30 items without --paginate, and a silent
`OK: every account with access is declared` is the one failure mode that
defeats the script's whole purpose."""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import check_maintainers as cm  # noqa: E402


def test_paginated_pages_are_flattened():
    """--slurp returns a list of pages; every account across them must survive."""
    pages = [[{"login": f"user{i}"} for i in range(30)], [{"login": "user30"}]]
    assert cm._logins(pages) == {f"user{i}" for i in range(31)}


def test_single_unpaginated_page_still_works():
    assert cm._logins([{"login": "solo"}]) == {"solo"}


def test_api_calls_request_every_page():
    """Without --paginate the read stops at 30 and the script passes on a lie."""
    calls = []
    cm._run = lambda args: calls.append(args) or "[]"
    cm._gh("orgs/ClawBio/members")
    assert "--paginate" in calls[0] and "--slurp" in calls[0]


def test_empty_org_members_is_incomplete_not_a_pass(capsys, monkeypatch, tmp_path):
    """An org always has at least an owner. An empty read means the token cannot
    see members, which must not be reported as agreement."""
    monkeypatch.setattr(cm, "ROOT", tmp_path)
    (tmp_path / "MAINTAINERS.md").write_text("[@someone](https://github.com/someone)\n")
    monkeypatch.setattr(cm, "_gh", lambda path: [[]] if "members" in path else [[{"login": "someone"}]])

    assert cm.main() == 2
    assert "INCOMPLETE" in capsys.readouterr().out
