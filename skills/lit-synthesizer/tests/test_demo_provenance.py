"""Every demo paper must be a real record, field for field.

The original demo (Feb 2026) shipped three invented papers: PMIDs that
resolve to unrelated papers (37001234 is a d-lyxose synthesis), DOIs that
404, and real authors attached to titles they never wrote. These tests check
each demo field against the saved PubMed and bioRxiv responses in
examples/demo_sources/, parsed by the skill's own parsers, so a made-up or
hand-edited record cannot pass. Refresh the sources with
examples/demo_sources/fetch_demo_sources.py.
"""
import json
import re
import sys
from pathlib import Path

SKILL = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(SKILL))
import lit_synthesizer as ls  # noqa: E402

SOURCES = SKILL / "examples" / "demo_sources"
FIELDS = ("pmid", "title", "authors", "journal", "year", "doi", "abstract", "url")


def _demo():
    return ls.load_demo_papers()


def test_demo_has_pubmed_and_biorxiv_papers():
    sources = [p["source"] for p in _demo()]
    assert sources.count("PubMed") >= 2 and sources.count("bioRxiv") >= 1


def test_every_pubmed_demo_paper_matches_the_saved_pubmed_record():
    parsed = {p["pmid"]: p for p in ls.parse_pubmed_xml((SOURCES / "pubmed_efetch.xml").read_bytes())}
    for paper in (p for p in _demo() if p["source"] == "PubMed"):
        record = parsed[paper["pmid"]]
        for field in FIELDS:
            assert paper[field] == record[field], (paper["pmid"], field)


def test_every_biorxiv_demo_paper_matches_the_saved_biorxiv_record():
    items = json.loads((SOURCES / "biorxiv_details.json").read_text())["collection"]
    parsed = {ls.biorxiv_item_to_paper(i)["doi"]: ls.biorxiv_item_to_paper(i) for i in items}
    for paper in (p for p in _demo() if p["source"] == "bioRxiv"):
        record = parsed[paper["doi"]]
        for field in FIELDS:
            assert paper[field] == record[field], (paper["doi"], field)


def test_every_citation_edge_is_in_pubmeds_reference_list():
    refs = {}
    for linkset in json.loads((SOURCES / "pubmed_elink_refs.json").read_text())["linksets"]:
        links = {x for db in linkset.get("linksetdbs", []) for x in db["links"]}
        refs[linkset["ids"][0]] = links
    edges = 0
    for paper in _demo():
        for cited in paper["citations"]:
            assert cited in refs.get(paper["pmid"], set()), (paper["pmid"], cited)
            edges += 1
    assert edges >= 1


def test_abstract_keeps_text_after_inline_markup():
    # PubMed wraps reference numbers in <sup>; ElementTree's .text stops at
    # the first child, which cut Anzalone 2019 to "Most genetic variants that
    # contribute to disease".
    xml = (
        b"<PubmedArticleSet><PubmedArticle><MedlineCitation><PMID>1</PMID><Article>"
        b"<ArticleTitle>T</ArticleTitle><Abstract>"
        b"<AbstractText>Most variants<sup>1</sup> are hard to fix.</AbstractText>"
        b"<AbstractText>Second part.</AbstractText>"
        b"</Abstract></Article></MedlineCitation></PubmedArticle></PubmedArticleSet>"
    )
    assert ls.parse_pubmed_xml(xml)[0]["abstract"] == "Most variants1 are hard to fix. Second part."


def test_skill_md_example_titles_are_real_demo_titles():
    titles = {p["title"] for p in _demo()}
    text = (SKILL / "SKILL.md").read_text()
    example = text.split("## Example Output", 1)[1].split("## Output Structure", 1)[0]
    shown = re.findall(r"^### \d+\. (.+)$", example, flags=re.M)
    assert shown and set(shown) <= titles


def test_committed_demo_report_lists_the_demo_titles():
    report = (SKILL / "examples" / "demo_output" / "report.md").read_text()
    for paper in _demo():
        assert paper["title"] in report
