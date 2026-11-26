import importlib.util
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("build_site", ROOT / "scripts" / "build_site.py")
build_site = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(build_site)

PAGE = (
    '<div class="skill-card"></div><div class="skill-card"></div>'
    '<p><span data-stat="stars">1</span> stars · <span data-stat="skills">2</span> skills'
    ' · + <span data-stat="more_skills">0</span> more · <span data-stat="forks">7</span></p>'
)


def test_render_fills_formats_and_derives_more_skills():
    out = build_site.render(PAGE, {"stars": 1234, "skills": 9})
    assert '<span data-stat="stars">1,234</span>' in out
    assert '<span data-stat="skills">9</span>' in out
    assert '<span data-stat="more_skills">7</span>' in out  # 9 skills - 2 cards


def test_render_keeps_page_value_for_missing_fact():
    out = build_site.render(PAGE, {"skills": 9})
    assert '<span data-stat="forks">7</span>' in out


def test_render_rejects_unknown_key():
    with pytest.raises(ValueError, match="starz"):
        build_site.render('<span data-stat="starz">1</span>', {"skills": 1})


def test_github_failure_keeps_local_facts(capsys):
    def broken():
        raise subprocess.CalledProcessError(1, ["gh"])

    facts = build_site.site_facts(github=broken)
    assert {"skills", "demo", "cli", "tests"} <= facts.keys()
    assert "stars" not in facts
    assert "GitHub lookup failed" in capsys.readouterr().err


def test_homepage_marks_every_fact_and_only_known_ones():
    page = (ROOT / "index.html").read_text(encoding="utf-8")
    keys = {m[2] for m in build_site.STAT.finditer(page)}
    assert keys == build_site.KNOWN
    # A full render with every fact present leaves no value unfilled.
    facts = {k: 1 for k in build_site.KNOWN - {"more_skills"}} | {"skills": 105}
    out = build_site.render(page, facts)
    assert '<span data-stat="more_skills">100</span>' in out
