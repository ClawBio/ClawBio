"""Fill index.html's live numbers at deploy (.github/workflows/pages.yml).

Every live number on the page sits in `<span data-stat="KEY">N</span>`; this rewrites N
from its source. The committed N is the fallback: if the GitHub lookup fails, the page
keeps it. Benchmark figures are not filled here: bench_runs/latest/ is not refreshed by
the weekly job (b52ac21a), so they stay hand-typed with their date and scope.

    python scripts/build_site.py          # rewrite index.html in place
    python scripts/build_site.py --check  # print the facts, write nothing
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from collections.abc import Callable
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REPO = "ClawBio/ClawBio"

sys.path.insert(0, str(ROOT / "scripts"))
from generate_catalog import build_catalog

KNOWN = {"stars", "forks", "contributors", "skills", "demo", "cli", "tests", "more_skills"}
STAT = re.compile(r'(<span data-stat="([a-z_]+)">)[^<]*(</span>)')
TEST_DEF = re.compile(r"^\s*def test_", re.MULTILINE)


def _gh(*args: str) -> str:
    return subprocess.run(["gh", *args], capture_output=True, text=True,
                          cwd=ROOT, check=True).stdout


def github_facts() -> dict[str, int]:
    repo = json.loads(_gh("api", f"repos/{REPO}"))
    # The contributors API misses anyone without an attributed commit, so union it
    # with merged-PR authors (the rule settled in b52ac21a).
    api = set(_gh("api", f"repos/{REPO}/contributors", "--paginate",
                  "-q", ".[].login").split())
    prs = set(_gh("pr", "list", "--repo", REPO, "--state", "merged", "--limit", "10000",
                  "--json", "author", "-q", ".[].author.login").split())
    return {"stars": repo["stargazers_count"], "forks": repo["forks_count"],
            "contributors": len(api | prs)}


def count_tests(root: Path = ROOT) -> int:
    files = subprocess.run(["git", "ls-files", "*.py"], capture_output=True, text=True,
                           cwd=root, check=True).stdout.splitlines()
    return sum(len(TEST_DEF.findall((root / f).read_text(encoding="utf-8", errors="replace")))
               for f in files if (root / f).is_file())


def site_facts(github: Callable[[], dict[str, int]] = github_facts) -> dict[str, int]:
    """Every live number the site shows, each from its one source."""
    catalog = build_catalog()
    facts = {
        "skills": len(catalog),
        "demo": sum(1 for s in catalog if s["has_demo"]),
        "cli": sum(1 for s in catalog if s["has_script"]),
        "tests": count_tests(),
    }
    try:
        facts |= github()
    except (subprocess.CalledProcessError, OSError, ValueError, KeyError) as exc:
        print(f"warning: GitHub lookup failed, keeping the page's stars, forks and "
              f"contributors: {exc}", file=sys.stderr)
    return facts


def render(page: str, facts: dict[str, int]) -> str:
    """Fill each data-stat span from facts; a fact that is missing keeps the page's value."""
    facts = {**facts, "more_skills": facts["skills"] - page.count('class="skill-card"')}

    def fill(m: re.Match) -> str:
        if m[2] not in KNOWN:
            raise ValueError(f"index.html: unknown data-stat {m[2]!r}")
        return m[1] + f"{facts[m[2]]:,}" + m[3] if m[2] in facts else m[0]

    return STAT.sub(fill, page)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--check", action="store_true", help="print the facts, write nothing")
    args = ap.parse_args()
    facts = site_facts()
    if args.check:
        print(json.dumps(facts, indent=2))
        return
    page_path = ROOT / "index.html"
    page_path.write_text(render(page_path.read_text(encoding="utf-8"), facts), encoding="utf-8")
    print(f"index.html filled: {facts}")


if __name__ == "__main__":
    main()
