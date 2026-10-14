"""The homepage must agree with itself and with skills/catalog.json.

The automated site-stats routine has repeatedly opened PRs whose index.html
refreshed the hero numbers and left other count sites stale (ClawBio#576,
#589, #609: demo 92 against a catalog 99, CLI 91 against 98, contributors 50
in one place and 70 in another). A reader sees both numbers on the same page.

Every count site is found by a fixed anchor. An anchor that matches nothing
FAILS ("anchor missing"): a redesign must update this test, never silently
skip a check. No network: stars and contributors cannot be checked against
GitHub here, so they are checked for agreement across every site instead.

Sources of truth:
  skills    skills/catalog.json  skill_count
  demo      count of skills with has_demo
  CLI       count of skills with has_script (the page says "with a
            deterministic CLI entry point", which is a runnable script; the
            narrower cli_alias field counts registered `clawbio.py run`
            aliases, a different and smaller number)
  Galaxy    skills/catalog.json  galaxy_tool_count (only if the page states one)
"""

from __future__ import annotations

import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
INDEX = ROOT / "index.html"
CATALOG = ROOT / "skills" / "catalog.json"

EXPECTED_SKILL_CARDS = 5
N = r"(\d[\d,]*)"

HERO = re.compile(
    r'class="hero-badge">[^<]*?' + N + r" stars · " + N + r" skills · "
    + N + r" contributors · " + N + r" tests"
)
LIST_ITEM = re.compile(N + r" skills</strong>, " + N + r" with a deterministic CLI entry point")
SECTION = re.compile(
    r'class="section-title">' + N + r" skills · " + N + r" with runnable demo data · " + N + r" tests"
)
MORE = re.compile(r"\+ " + N + r" more skills")
COMPARISON_ROW = re.compile(
    r'href="https://github\.com/ClawBio/ClawBio"[^>]*>ClawBio</a></td>'
    r"<td[^>]*>[^<]*</td>"            # OpenClaw column
    r"<td[^>]*>" + N + r"</td>"       # stars
    r"<td[^>]*>" + N + r"</td>"       # forks
    r"<td[^>]*>" + N + r"</td>"       # skills
)
COMMUNITY = re.compile(r"More than " + N + r" contributors")
GALAXY = re.compile(N + r" Galaxy tools")  # optional: no anchor required
SKILL_CARD = 'class="skill-card"'


def _int(s: str) -> int:
    return int(s.replace(",", ""))


def homepage_problems(html: str, catalog: dict) -> list[str]:
    """Every inconsistency between the page and the catalog, as readable strings."""
    skills_cat = catalog["skill_count"]
    demo_cat = sum(1 for s in catalog["skills"] if s.get("has_demo"))
    cli_cat = sum(1 for s in catalog["skills"] if s.get("has_script"))
    problems: list[str] = []

    def one(name: str, pattern: re.Pattern) -> re.Match | None:
        hits = list(pattern.finditer(html))
        if not hits:
            problems.append(f"anchor missing: {name} (pattern {pattern.pattern!r}); update this test with the redesign")
            return None
        if len(hits) > 1:
            problems.append(f"anchor ambiguous: {name} matched {len(hits)} times")
        return hits[0]

    skills: dict[str, int] = {}
    stars: dict[str, int] = {}
    contributors: dict[str, int] = {}
    tests: dict[str, int] = {}

    if m := one("hero badge", HERO):
        stars["hero"], skills["hero"] = _int(m[1]), _int(m[2])
        contributors["hero"], tests["hero"] = _int(m[3]), _int(m[4])
    if m := one("skills/CLI list item", LIST_ITEM):
        skills["list item"] = _int(m[1])
        if _int(m[2]) != cli_cat:
            problems.append(f"CLI figure {_int(m[2])} != catalog has_script count {cli_cat}")
    if m := one("skills section title", SECTION):
        skills["section title"], tests["section title"] = _int(m[1]), _int(m[3])
        if _int(m[2]) != demo_cat:
            problems.append(f"demo figure {_int(m[2])} != catalog has_demo count {demo_cat}")
    if m := one("comparison table ClawBio row", COMPARISON_ROW):
        stars["comparison row"], skills["comparison row"] = _int(m[1]), _int(m[3])
    if m := one("community contributors", COMMUNITY):
        contributors["community"] = _int(m[1])

    for site, value in skills.items():
        if value != skills_cat:
            problems.append(f"skills at {site} = {value} != catalog skill_count {skills_cat}")

    cards = html.count(SKILL_CARD)
    if cards != EXPECTED_SKILL_CARDS:
        problems.append(f"skill cards = {cards}, expected exactly {EXPECTED_SKILL_CARDS}")
    if m := one("+ N more skills", MORE):
        if _int(m[1]) != skills_cat - cards:
            problems.append(f"'+ {_int(m[1])} more skills' != skill_count {skills_cat} - {cards} cards = {skills_cat - cards}")

    for label, sites in (("stars", stars), ("contributors", contributors), ("tests", tests)):
        if len(set(sites.values())) > 1:
            problems.append(f"{label} disagree across the page: {sites}")

    for m in GALAXY.finditer(html):
        if _int(m[1]) != catalog.get("galaxy_tool_count"):
            problems.append(f"Galaxy tools {_int(m[1])} != catalog galaxy_tool_count {catalog.get('galaxy_tool_count')}")

    return problems


def test_homepage_counts_are_consistent():
    problems = homepage_problems(INDEX.read_text(), json.loads(CATALOG.read_text()))
    assert not problems, "index.html is inconsistent:\n  " + "\n  ".join(problems)


# Guards: prove each check still fires, so a loosened pattern cannot pass silently.

def _fixture_catalog() -> dict:
    skills = [{"has_demo": i < 9, "has_script": i < 8} for i in range(10)]
    return {"skill_count": 10, "galaxy_tool_count": 7, "skills": skills}


GOOD = (
    '<div class="hero-badge">x · 1,000 stars · 10 skills · 70 contributors · 5,000 tests</div>'
    "<li><strong>10 skills</strong>, 8 with a deterministic CLI entry point</li>"
    '<div class="section-title">10 skills · 9 with runnable demo data · 5,000 tests · 1 / 2 passing</div>'
    + '<div class="skill-card"></div>' * 5
    + "<p>+ 5 more skills on GitHub</p>"
    '<tr><td><a href="https://github.com/ClawBio/ClawBio" target="_blank">ClawBio</a></td>'
    "<td>&#10003;</td><td>1,000</td><td>200</td><td>10</td></tr>"
    "<p>More than 70 contributors from labs</p>"
)


def test_guard_good_fixture_passes():
    assert homepage_problems(GOOD, _fixture_catalog()) == []


def test_guard_fires_on_each_stale_site():
    cases = {
        "demo": GOOD.replace("9 with runnable", "8 with runnable"),
        "CLI": GOOD.replace("8 with a deterministic", "7 with a deterministic"),
        "community": GOOD.replace("More than 70", "More than 50"),
        "stars": GOOD.replace("<td>1,000</td>", "<td>999</td>"),
        "skills": GOOD.replace("<td>10</td></tr>", "<td>9</td></tr>"),
        "more": GOOD.replace("+ 5 more", "+ 4 more"),
        "cards": GOOD.replace('<div class="skill-card"></div>', "", 1),
        "tests": GOOD.replace("9 with runnable demo data · 5,000", "9 with runnable demo data · 4,000"),
        "galaxy": GOOD + "<p>6 Galaxy tools</p>",
    }
    for name, html in cases.items():
        assert homepage_problems(html, _fixture_catalog()), f"guard disarmed: {name}"


def test_guard_missing_anchor_fails_not_passes():
    problems = homepage_problems(GOOD.replace("More than 70 contributors", "Many contributors"), _fixture_catalog())
    assert any(p.startswith("anchor missing") for p in problems)
