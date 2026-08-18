"""Check whether the ACMG secondary-findings list we screen against is current.

Why this is not decoration
--------------------------
The list this skill screens against is `ACMG_SF_V32_GENES` — **v3.2**, 81 genes —
because that is what the library ships. Screening against a superseded list is
exactly the class of error the ledger exists to name, so leaving the version
unexamined would be the project failing its own standard.

But the honest way to establish "the current version is v3.3" is not to assert it
from memory. It is to fetch evidence and record what came back. This script asks
Tavily, prints the URLs it got, and states plainly which of three outcomes
occurred:

    CONFIRMED      evidence names a newer version and the added genes
    UP_TO_DATE     evidence indicates v3.2 is still current
    UNVERIFIED     no key, fetch failed, or the evidence did not settle it

`UNVERIFIED` is a normal outcome, not an error. It is reported as
`SF_LIST_VERSION_UNVERIFIED` rather than silently falling back to a remembered
answer — which is the whole argument of this skill applied to this skill.

Reads TAVILY_API_KEY from the environment or from a local .env. The key is never
printed, and neither is the URL that carries it.
"""

from __future__ import annotations

import argparse
import json
import os
import pathlib
import re
import sys
import urllib.error
import urllib.request

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

TAVILY_SEARCH = "https://api.tavily.com/search"
QUERY = (
    "ACMG SF v3.3 secondary findings gene list 2025 number of genes "
    "added ABCD1 CYP27A1 PLN policy statement"
)
BUNDLED_VERSION = "v3.2"


def load_key() -> str | None:
    key = os.environ.get("TAVILY_API_KEY")
    if key:
        return key
    for candidate in (
        pathlib.Path(__file__).resolve().parents[2] / ".env",
        pathlib.Path.cwd() / ".env",
    ):
        if candidate.exists():
            for line in candidate.read_text().splitlines():
                line = line.strip()
                if line.startswith("TAVILY_API_KEY="):
                    return line.split("=", 1)[1].strip() or None
    return None


def search(key: str, query: str, *, timeout: int = 60) -> dict:
    req = urllib.request.Request(
        TAVILY_SEARCH,
        data=json.dumps(
            {"query": query, "max_results": 6, "search_depth": "advanced",
             "include_answer": True}
        ).encode(),
        headers={
            "Content-Type": "application/json",
            "Authorization": f"Bearer {key}",
        },
        method="POST",
    )
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return json.load(resp)


def assess(blob: dict) -> dict:
    """Read a verdict out of the evidence, without reaching past it."""
    corpus = " ".join(
        [blob.get("answer") or ""]
        + [r.get("content", "") for r in blob.get("results", [])]
        + [r.get("title", "") for r in blob.get("results", [])]
    )

    versions = sorted(set(re.findall(r"\bv?3\.(\d)\b", corpus)))
    counts = sorted({int(c) for c in re.findall(r"\b(8[0-9])\s+genes\b", corpus, re.I)})
    added = sorted({g for g in ("ABCD1", "CYP27A1", "PLN") if re.search(rf"\b{g}\b", corpus)})

    newer = [v for v in versions if int(v) > 2]
    if newer and counts:
        verdict = "CONFIRMED"
        note = (
            f"evidence names ACMG SF v3.{max(newer)} with {max(counts)} genes; "
            f"the bundled list is {BUNDLED_VERSION} with 81. "
            f"Genes named as added: {', '.join(added) if added else 'none identified in the snippets'}."
        )
    elif versions and not newer:
        verdict = "UP_TO_DATE"
        note = f"evidence mentions only v3.{max(versions)}; the bundled list appears current."
    else:
        verdict = "UNVERIFIED"
        note = (
            "the retrieved snippets did not settle the version. Reported as "
            "SF_LIST_VERSION_UNVERIFIED rather than filled in from memory."
        )

    return {
        "verdict": verdict,
        "note": note,
        "versions_seen": [f"v3.{v}" for v in versions],
        "gene_counts_seen": counts,
        "added_genes_named": added,
        "sources": [
            {"title": r.get("title", ""), "url": r.get("url", "")}
            for r in blob.get("results", [])
        ],
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--output", type=pathlib.Path, required=True)
    args = ap.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    from gates import ACMG_SF_GENES, ACMG_SF_SOURCE  # noqa: PLC0415

    payload: dict = {
        "bundled_version": BUNDLED_VERSION,
        "bundled_gene_count": len(ACMG_SF_GENES) if ACMG_SF_GENES else 0,
        "bundled_source": ACMG_SF_SOURCE,
        "query": QUERY,
    }

    key = load_key()
    if not key:
        payload |= {
            "verdict": "UNVERIFIED",
            "note": "no TAVILY_API_KEY available; the version was not checked.",
            "sources": [],
        }
    else:
        try:
            payload |= assess(search(key, QUERY))
        except (urllib.error.URLError, TimeoutError, OSError, json.JSONDecodeError) as exc:
            payload |= {
                "verdict": "UNVERIFIED",
                "note": f"fetch failed ({type(exc).__name__}); the version was not checked.",
                "sources": [],
            }

    (args.output / "sf_version_check.json").write_text(json.dumps(payload, indent=1))

    L = [
        "# Is the secondary-findings list we screen against current?",
        "",
        f"**Verdict: `{payload['verdict']}`**",
        "",
        payload["note"],
        "",
        f"- Bundled list: ACMG SF **{payload['bundled_version']}**, "
        f"{payload['bundled_gene_count']} genes",
        f"- Provenance: {payload['bundled_source']}",
        "",
    ]
    if payload.get("sources"):
        L += ["## Evidence retrieved", ""]
        L += [f"- [{s['title']}]({s['url']})" for s in payload["sources"] if s.get("url")]
        L += [""]
    L += [
        "## Why this matters here rather than in general",
        "",
        "This skill withholds records for screening against stale evidence. Screening "
        "them against a stale gene list would be the same error one level up. The "
        "version therefore travels with every result, and where it could not be "
        "established the report says so instead of assuming.",
        "",
    ]
    (args.output / "sf_version_check.md").write_text("\n".join(L))

    print(f"{payload['verdict']}: {payload['note']}")
    print(f"-> {args.output}/sf_version_check.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
