"""Public claims must match the repository they describe.

Every number a reader can quote (skill count, demo count, CLI count, release
version) is derived here from its source of truth and compared against the
files where it is repeated by hand. The plan that produced this gate found the
README, CITATION.cff, .zenodo.json and llms.txt disagreeing with each other by
up to 55 skills and two minor releases, and a "first" claim that no evidence
supports. A stale public claim is a defect, so this test fails CI on drift.

Sources of truth:
  skill counts   skills/catalog.json (written by scripts/generate_catalog.py)
  version        clawbio/__init__.py

How to fix a failure: regenerate the catalogue (`python scripts/generate_catalog.py`)
and copy the numbers it reports into the file named in the assertion; never
hand-type a count. To bump the version, change clawbio/__init__.py first.
"""

from __future__ import annotations

import json
import re
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
CATALOG = ROOT / "skills" / "catalog.json"

# Files that repeat the skill count by hand. Each `<number> skills` (or
# `<number> Agent Skills`) in them must equal catalog.json's skill_count.
COUNT_FILES = [
    "README.md",
    "llms.txt",
    "CITATION.cff",
    ".zenodo.json",
    ".claude-plugin/plugin.json",
    ".claude-plugin/marketplace.json",
]

# Files that state the current release. Each must carry the package version.
VERSION_FILES = {
    "README.md": r"The current release is \*\*v([0-9][0-9a-z.\-]*)\*\*",
    "llms.txt": r"Current public release: `v([0-9][0-9a-z.\-]*)`",
    "CITATION.cff": r"^version: ([0-9][0-9a-z.\-]*)$",
    ".zenodo.json": r'"version": "([0-9][0-9a-z.\-]*)"',
    ".claude-plugin/plugin.json": r'"version": "([0-9][0-9a-z.\-]*)"',
    ".claude-plugin/marketplace.json": r'"version": "([0-9][0-9a-z.\-]*)"',
}

# A primacy claim the project cannot evidence. Scanned across every tracked
# text file; the allowlist holds dated records of past talks and submissions,
# which are history and are not rewritten.
PRIMACY = re.compile(
    r"\bfirst\s+bioinformatics[\s-]native\b|\bfirst\s+(?:AI\s+)?agent\s+skill\s+library\b",
    re.IGNORECASE,
)
PRIMACY_ALLOWLIST_PREFIXES = (
    "slides/",  # London Bioinformatics Meetup deck, February 2026, as delivered
    "docs/dorahacks-buidl-submission.md",  # submission text as filed, 5 March 2026
    "tests/test_public_claims.py",
)
TEXT_SUFFIXES = {".md", ".txt", ".toml", ".json", ".cff", ".html", ".py", ".yml", ".yaml", ".cfg", ".ini"}

COUNT_PATTERN = re.compile(r"(?<![\d.,])(\d[\d,]*)\s+(?:Agent\s+)?skills\b", re.IGNORECASE)
DEMO_PATTERN = re.compile(r"(\d+) with runnable demo data")
CLI_PATTERN = re.compile(r"(\d+) (?:with a deterministic CLI entry point|run deterministically from the CLI)")


def _catalog() -> dict:
    return json.loads(CATALOG.read_text())


def _package_version() -> str:
    text = (ROOT / "clawbio" / "__init__.py").read_text()
    match = re.search(r'^__version__\s*=\s*"([^"]+)"', text, re.MULTILINE)
    assert match, "clawbio/__init__.py has no __version__"
    return match.group(1)


def _tracked_text_files() -> list[Path]:
    out = subprocess.run(
        ["git", "ls-files", "-z"], cwd=ROOT, capture_output=True, check=True
    ).stdout
    paths = [ROOT / p for p in out.decode().split("\0") if p]
    return [p for p in paths if p.suffix in TEXT_SUFFIXES and p.is_file()]


def test_catalog_skill_count_matches_skill_entries():
    cat = _catalog()
    assert cat["skill_count"] == len(cat["skills"]), "catalog.json is internally inconsistent; regenerate it"


@pytest.mark.parametrize("relpath", COUNT_FILES)
def test_skill_counts_match_catalog(relpath: str):
    expected = _catalog()["skill_count"]
    text = (ROOT / relpath).read_text()
    found = [(m.group(1), m.group(0)) for m in COUNT_PATTERN.finditer(text)]
    assert found, f"{relpath} states no skill count; add one from catalog.json or drop it from COUNT_FILES"
    wrong = [phrase for number, phrase in found if int(number.replace(",", "")) != expected]
    assert not wrong, f"{relpath} disagrees with catalog.json skill_count={expected}: {wrong}"


def test_demo_and_cli_counts_match_catalog():
    skills = _catalog()["skills"]
    demo = sum(1 for s in skills if s.get("has_demo"))
    cli = sum(1 for s in skills if s.get("cli_alias"))
    problems = []
    for relpath in COUNT_FILES:
        text = (ROOT / relpath).read_text()
        for m in DEMO_PATTERN.finditer(text):
            if int(m.group(1)) != demo:
                problems.append(f"{relpath}: '{m.group(0)}' but catalog has_demo={demo}")
        for m in CLI_PATTERN.finditer(text):
            if int(m.group(1)) != cli:
                problems.append(f"{relpath}: '{m.group(0)}' but catalog cli_alias={cli}")
    assert not problems, problems


@pytest.mark.parametrize("relpath,pattern", sorted(VERSION_FILES.items()))
def test_stated_release_matches_package_version(relpath: str, pattern: str):
    expected = _package_version()
    text = (ROOT / relpath).read_text()
    found = re.findall(pattern, text, re.MULTILINE)
    assert found, f"{relpath} states no release version (pattern {pattern!r})"
    wrong = sorted(set(v for v in found if v != expected))
    assert not wrong, f"{relpath} states release {wrong}, package is {expected}"


def test_no_primacy_claim_anywhere():
    hits = []
    for path in _tracked_text_files():
        rel = path.relative_to(ROOT).as_posix()
        if rel.startswith(PRIMACY_ALLOWLIST_PREFIXES):
            continue
        for lineno, line in enumerate(path.read_text(errors="ignore").splitlines(), 1):
            if PRIMACY.search(line):
                hits.append(f"{rel}:{lineno}: {line.strip()[:100]}")
    assert not hits, "primacy claim found (the project does not evidence it):\n" + "\n".join(hits)
