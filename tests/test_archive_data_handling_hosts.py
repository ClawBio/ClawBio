"""Every host an archive skill contacts must be named in its data-handling row.

docs/data-handling.md tells network administrators which hosts to allowlist
for each skill. A host missing from a row is worse than no row: someone who
allowlists exactly the listed hosts sees some commands fail, which is the
confusion the page's allowlisting section exists to prevent. That happened to
geo-fetch: `samplesheet`, `runtable` and `metadata-table` call
www.ncbi.nlm.nih.gov, and the row did not list it.
"""
from __future__ import annotations

import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
DOC = ROOT / "docs" / "data-handling.md"
SKILLS = ["arrayexpress-fetch", "biostudies-fetch", "ena-fetch",
          "geo-fetch", "pride-fetch"]
HOST = re.compile(r"https?://([a-z0-9.-]+)")


def hosts_in_code(skill: str) -> set[str]:
    """URL hosts written in the skill's source, ignoring comment lines."""
    hosts = set()
    for path in (ROOT / "skills" / skill).glob("*.py"):
        for line in path.read_text().splitlines():
            if not line.lstrip().startswith("#"):
                hosts.update(HOST.findall(line))
    return hosts


def doc_row(skill: str) -> str:
    rows = [ln for ln in DOC.read_text().splitlines() if ln.startswith(f"| `{skill}` |")]
    assert len(rows) == 1, f"expected one data-handling row for {skill}, found {len(rows)}"
    return rows[0]


@pytest.mark.parametrize("skill", SKILLS)
def test_every_contacted_host_is_in_the_skills_row(skill):
    hosts = hosts_in_code(skill)
    assert hosts, f"found no hosts in {skill}; the scan is broken"
    row = doc_row(skill)
    missing = sorted(h for h in hosts if h not in row)
    assert not missing, f"docs/data-handling.md row for {skill} omits {missing}"


def test_every_contacted_host_is_in_the_allowlisting_table():
    text = DOC.read_text().split("### Allowlisting for the public-archive skills", 1)[1]
    table = {m for m in re.findall(r"^\| `([a-z0-9.-]+)` \|", text, re.M)}
    wanted = set().union(*(hosts_in_code(s) for s in SKILLS))
    assert not sorted(wanted - table), f"allowlisting table omits {sorted(wanted - table)}"
