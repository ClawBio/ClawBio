"""MAINTAINERS.md is a claim about who controls the project. These tests keep
it structurally honest: the documented rights tables exist, every handle named
in CODEOWNERS is declared there, and the two-person rule on the controls that
matter is asserted as an xfail that turns into a visible XPASS the day it is
met, so the marker is removed rather than the gap being forgotten.
"""
from __future__ import annotations

import re
from pathlib import Path
from urllib.parse import urlparse

import pytest

ROOT = Path(__file__).resolve().parents[1]
MAINTAINERS = ROOT / "MAINTAINERS.md"
CODEOWNERS = ROOT / ".github" / "CODEOWNERS"

REQUIRED_SECTIONS = [
    "## Lead maintainer",
    "## Repository permissions",
    "## Organisation",
    "## Other controls",
]

# Controls that GOVERNANCE.md's thirty-day rule and the release path depend on.
SINGLE_PERSON_CONTROLS = ["PyPI project", "Domain `clawbio.ai`"]

_HANDLE = re.compile(r"\[@([A-Za-z0-9-]+)\]\(https://github\.com/[A-Za-z0-9-]+\)")


def _text() -> str:
    return MAINTAINERS.read_text()


def test_maintainers_file_names_a_lead_with_orcid():
    t = _text()
    assert "## Lead maintainer" in t
    assert "[@manuelcorpas](https://github.com/manuelcorpas)" in t
    urls = re.findall(r"\[[^\]]+\]\((https?://[^)\s]+)\)", t)
    assert any(
        (host := urlparse(url).hostname) and (host == "orcid.org" or host.endswith(".orcid.org"))
        for url in urls
    )


def test_every_documented_section_is_present():
    t = _text()
    missing = [s for s in REQUIRED_SECTIONS if s not in t]
    assert missing == [], f"sections missing from MAINTAINERS.md: {missing}"


def test_codeowners_handles_are_declared_in_maintainers():
    owners = set(re.findall(r"@([A-Za-z0-9-]+)", CODEOWNERS.read_text()))
    declared = set(_HANDLE.findall(_text()))
    undeclared = sorted(owners - declared)
    assert undeclared == [], f"CODEOWNERS handles absent from MAINTAINERS.md: {undeclared}"


def test_every_handle_link_points_at_its_own_profile():
    """A link whose text and target disagree is a copy-paste error waiting to
    misdirect a security contact."""
    bad = [(h, url) for h, url in re.findall(
        r"\[@([A-Za-z0-9-]+)\]\(https://github\.com/([A-Za-z0-9-]+)\)", _text()) if h != url]
    assert bad == [], f"handle/link mismatches: {bad}"


def test_critical_control_rows_exist():
    """Outside the xfail below on purpose: under strict xfail a missing row is
    indistinguishable from a genuine single holder, so renaming a control would
    silently retire the check instead of failing."""
    missing = [c for c in SINGLE_PERSON_CONTROLS
               if not any(c in ln for ln in _text().splitlines())]
    assert missing == [], f"controls named in the test but absent from MAINTAINERS.md: {missing}"


@pytest.mark.xfail(strict=True, reason="PyPI and DNS are held by one person; remove this marker when a second holder is listed")
def test_critical_controls_have_two_holders():
    t = _text()
    single = []
    for control in SINGLE_PERSON_CONTROLS:
        row = next((ln for ln in t.splitlines() if control in ln), None)
        assert row is not None, (
            f"no row for {control!r} in MAINTAINERS.md; this test cannot report on a "
            "control it cannot find, and an unfound row must not read as a single holder")
        holders = [c.strip() for c in row.split("|")[2:3]]
        names = re.split(r",| and ", holders[0]) if holders else []
        if len([n for n in names if n.strip()]) < 2:
            single.append(control)
    assert single == [], f"controls with a single holder: {single}"
