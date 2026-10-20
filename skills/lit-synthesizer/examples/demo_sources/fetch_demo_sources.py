#!/usr/bin/env python3
"""Re-fetch the raw source records behind the lit-synthesizer demo.

The demo papers are real. This script saves the exact responses they are
parsed from, so `tests/test_demo_provenance.py` can check every demo field
against the source offline, then rebuilds examples/demo_papers.json from them with the
skill's own parsers. Refresh the sources (network), or rebuild only:

    python skills/lit-synthesizer/examples/demo_sources/fetch_demo_sources.py
    python skills/lit-synthesizer/examples/demo_sources/fetch_demo_sources.py --offline
"""
import importlib.util
import json
import sys
import urllib.parse
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
PMIDS = ["27096365", "29160308", "31634902"]
BIORXIV_DOI = "10.64898/2026.07.06.736660"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def _get(url: str) -> bytes:
    with urllib.request.urlopen(url, timeout=30) as resp:
        return resp.read()


def fetch() -> None:
    ids = ",".join(PMIDS)
    efetch = f"{EUTILS}/efetch.fcgi?" + urllib.parse.urlencode(
        {"db": "pubmed", "id": ids, "retmode": "xml", "rettype": "abstract"}
    )
    (HERE / "pubmed_efetch.xml").write_bytes(_get(efetch))
    # One linkset per PMID: `id` repeated, not comma-joined.
    elink = f"{EUTILS}/elink.fcgi?" + urllib.parse.urlencode(
        [("dbfrom", "pubmed"), ("linkname", "pubmed_pubmed_refs"), ("retmode", "json")]
        + [("id", p) for p in PMIDS]
    )
    refs = json.loads(_get(elink))
    (HERE / "pubmed_elink_refs.json").write_text(json.dumps(refs, indent=1) + "\n")
    details = json.loads(_get(f"https://api.biorxiv.org/details/biorxiv/{BIORXIV_DOI}"))
    (HERE / "biorxiv_details.json").write_text(json.dumps(details, indent=1) + "\n")


def build() -> list[dict]:
    """Parse the saved sources into demo papers, citations within the set only."""
    spec = importlib.util.spec_from_file_location("lit_synthesizer", HERE.parents[1] / "lit_synthesizer.py")
    ls = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(ls)
    by_pmid = {p["pmid"]: p for p in ls.parse_pubmed_xml((HERE / "pubmed_efetch.xml").read_bytes())}
    refs = {}
    for linkset in json.loads((HERE / "pubmed_elink_refs.json").read_text())["linksets"]:
        refs[linkset["ids"][0]] = {x for db in linkset.get("linksetdbs", []) for x in db["links"]}
    papers = [by_pmid[p] for p in PMIDS]
    for paper in papers:
        paper["citations"] = [p for p in PMIDS if p in refs.get(paper["pmid"], set())]
    items = json.loads((HERE / "biorxiv_details.json").read_text())["collection"]
    papers.append(ls.biorxiv_item_to_paper(items[-1]))
    out = HERE.parent / "demo_papers.json"
    out.write_text(json.dumps(papers, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return papers


if __name__ == "__main__":
    if "--offline" not in sys.argv[1:]:
        fetch()
    build()
