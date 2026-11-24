"""
Tests for ClawBio Lit Synthesizer skill.
Run with: pytest skills/lit-synthesizer/tests/test_lit_synthesizer.py -v
"""

import json
import sys
import tempfile
from pathlib import Path

import pytest
from inline_snapshot import snapshot
from clawbio.common import demo_runner

# Add skill directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))
from lit_synthesizer import (
    load_demo_papers,
    build_citation_graph,
    generate_report,
)


# ── Citation graph tests ───────────────────────────────────────────────────────

class TestBuildCitationGraph:
    def test_returns_nodes_and_edges_keys(self):
        graph = build_citation_graph(load_demo_papers())
        assert "nodes" in graph
        assert "edges" in graph

    def test_node_count_matches_papers(self):
        graph = build_citation_graph(load_demo_papers())
        assert len(graph["nodes"]) == len(load_demo_papers())

    def test_nodes_have_required_fields(self):
        graph = build_citation_graph(load_demo_papers())
        for node in graph["nodes"]:
            assert "id" in node
            assert "title" in node
            assert "source" in node
            assert "year" in node

    def test_edges_reference_known_nodes(self):
        graph = build_citation_graph(load_demo_papers())
        node_ids = {n["id"] for n in graph["nodes"]}
        for edge in graph["edges"]:
            assert "from" in edge
            assert "to" in edge
            # 'to' should be a PMID from our demo set or an external reference
            assert edge["from"] in node_ids

    def test_empty_papers_returns_empty_graph(self):
        graph = build_citation_graph([])
        assert graph["nodes"] == []
        assert graph["edges"] == []

    def test_paper_without_pmid_uses_doi_as_id(self):
        paper = {
            "source": "bioRxiv",
            "pmid": None,
            "title": "Test Paper",
            "year": "2025",
            "doi": "10.1101/test.001",
            "citations": [],
        }
        graph = build_citation_graph([paper])
        assert graph["nodes"][0]["id"] == "10.1101/test.001"


# ── Report generation tests ────────────────────────────────────────────────────

class TestGenerateReport:
    def test_creates_report_md(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            out = Path(tmpdir)
            graph = build_citation_graph(load_demo_papers())
            generate_report("CRISPR", load_demo_papers(), graph, out)
            assert (out / "report.md").exists()

    def test_creates_results_json(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            out = Path(tmpdir)
            graph = build_citation_graph(load_demo_papers())
            generate_report("CRISPR", load_demo_papers(), graph, out)
            assert (out / "results.json").exists()

    def test_creates_citation_graph_json(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            out = Path(tmpdir)
            graph = build_citation_graph(load_demo_papers())
            generate_report("CRISPR", load_demo_papers(), graph, out)
            assert (out / "citation_graph.json").exists()

    def test_creates_csv_table(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            out = Path(tmpdir)
            graph = build_citation_graph(load_demo_papers())
            generate_report("CRISPR", load_demo_papers(), graph, out)
            assert (out / "tables" / "papers.csv").exists()

    def test_creates_reproducibility_bundle(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            out = Path(tmpdir)
            graph = build_citation_graph(load_demo_papers())
            generate_report("CRISPR", load_demo_papers(), graph, out)
            assert (out / "reproducibility" / "commands.sh").exists()
            assert (out / "reproducibility" / "environment.yml").exists()
            assert (out / "reproducibility" / "checksums.sha256").exists()

    def test_results_json_is_valid(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            out = Path(tmpdir)
            graph = build_citation_graph(load_demo_papers())
            generate_report("CRISPR", load_demo_papers(), graph, out)
            data = json.loads((out / "results.json").read_text())
            assert isinstance(data, list)
            assert len(data) == len(load_demo_papers())

    def test_report_contains_query(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            out = Path(tmpdir)
            graph = build_citation_graph(load_demo_papers())
            generate_report("CRISPR genome editing", load_demo_papers(), graph, out)
            report_text = (out / "report.md").read_text()
            assert "CRISPR genome editing" in report_text

    def test_report_contains_disclaimer(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            out = Path(tmpdir)
            graph = build_citation_graph(load_demo_papers())
            generate_report("CRISPR", load_demo_papers(), graph, out)
            report_text = (out / "report.md").read_text()
            assert "research tool" in report_text.lower()

    def test_report_contains_paper_titles(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            out = Path(tmpdir)
            graph = build_citation_graph(load_demo_papers())
            generate_report("CRISPR", load_demo_papers(), graph, out)
            report_text = (out / "report.md").read_text()
            for paper in load_demo_papers():
                assert paper["title"][:30] in report_text

    def test_checksums_cover_all_files(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            out = Path(tmpdir)
            graph = build_citation_graph(load_demo_papers())
            generate_report("CRISPR", load_demo_papers(), graph, out)
            checksums_text = (out / "reproducibility" / "checksums.sha256").read_text()
            # report.md should be checksummed
            assert "report.md" in checksums_text

    def test_empty_papers_does_not_crash(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            out = Path(tmpdir)
            graph = build_citation_graph([])
            generate_report("empty query", [], graph, out)
            assert (out / "report.md").exists()


# ── Demo data integrity tests ──────────────────────────────────────────────────

class TestDemoData:
    def test_demo_papers_have_required_fields(self):
        required = ["source", "title", "authors", "journal", "year", "abstract", "doi", "url", "citations"]
        for paper in load_demo_papers():
            for field in required:
                assert field in paper, f"Demo paper missing field: {field}"

    def test_demo_papers_sources_are_valid(self):
        valid_sources = {"PubMed", "bioRxiv"}
        for paper in load_demo_papers():
            assert paper["source"] in valid_sources

    def test_demo_papers_authors_are_list(self):
        for paper in load_demo_papers():
            assert isinstance(paper["authors"], list)

    def test_demo_papers_citations_are_list(self):
        for paper in load_demo_papers():
            assert isinstance(paper["citations"], list)


# ── Demo output snapshot ──────────────────────────────────────────────────────
# Fails when a code change alters any demo number. If the change is intended,
# run `pytest skills/lit-synthesizer/tests --inline-snapshot=fix`; it rewrites the
# snapshot() values below, and the PR diff shows which numbers moved.

@pytest.mark.demo_snapshot
def test_demo_output_matches_snapshot(tmp_path):
    assert demo_runner.demo_outputs("lit-synthesizer", tmp_path) == snapshot({'citation_graph.json': """\
{
  "nodes": [
    {
      "id": "27096365",
      "title": "Programmable editing of a target base in genomic DNA without double-stranded DNA cleavage.",
      "source": "PubMed",
      "year": "2016"
    },
    {
      "id": "29160308",
      "title": "Programmable base editing of A\\u2022T to G\\u2022C in genomic DNA without DNA cleavage.",
      "source": "PubMed",
      "year": "2017"
    },
    {
      "id": "31634902",
      "title": "Search-and-replace genome editing without double-strand breaks or donor DNA.",
      "source": "PubMed",
      "year": "2019"
    },
    {
      "id": "10.64898/2026.07.06.736660",
      "title": "Saturated mutagenesis screen of M-MLV reverse transcriptase identifies variants enhancing prime editing efficiency",
      "source": "bioRxiv",
      "year": "2026"
    }
  ],
  "edges": [
    {
      "from": "29160308",
      "to": "27096365"
    },
    {
      "from": "31634902",
      "to": "27096365"
    },
    {
      "from": "31634902",
      "to": "29160308"
    }
  ]
}\
""", 'report.md': """\
# 🦖 ClawBio Lit Synthesizer Report

**Query**: `CRISPR genome editing`  \n\
**Date**: <TIMESTAMP> UTC  \n\
**Sources**: PubMed (3 results) · bioRxiv (1 results)  \n\
**Total papers**: 4

---

## Summary

Across 4 retrieved papers, recurring themes include: **crispr**, **genome editing**, **off-target**, **base editing**, **prime editing**, **cas9**. The literature spans 2016 to 2026.

---

## Papers

### 1. Programmable editing of a target base in genomic DNA without double-stranded DNA cleavage.

| Field | Value |
|-------|-------|
| **Source** | PubMed |
| **Authors** | Komor AC, Kim YB, Packer MS et al. |
| **Journal** | Nature |
| **Year** | 2016 |
| **DOI** | 10.1038/nature17946 |
| **URL** | [https://pubmed.ncbi.nlm.nih.gov/27096365](https://pubmed.ncbi.nlm.nih.gov/27096365) |

**Abstract**: Current genome-editing technologies introduce double-stranded (ds) DNA breaks at a target locus as the first step to gene correction. Although most genetic diseases arise from point mutations, current approaches to point mutation correction are inefficient and typically induce an abundance of random insertions and deletions (indels) at the target locus resulting from the cellular response to dsDNA...

### 2. Programmable base editing of A•T to G•C in genomic DNA without DNA cleavage.

| Field | Value |
|-------|-------|
| **Source** | PubMed |
| **Authors** | Gaudelli NM, Komor AC, Rees HA et al. |
| **Journal** | Nature |
| **Year** | 2017 |
| **DOI** | 10.1038/nature24644 |
| **URL** | [https://pubmed.ncbi.nlm.nih.gov/29160308](https://pubmed.ncbi.nlm.nih.gov/29160308) |

**Abstract**: The spontaneous deamination of cytosine is a major source of transitions from C•G to T•A base pairs, which account for half of known pathogenic point mutations in humans. The ability to efficiently convert targeted A•T base pairs to G•C could therefore advance the study and treatment of genetic diseases. The deamination of adenine yields inosine, which is treated as guanine by polymerases, but no ...

### 3. Search-and-replace genome editing without double-strand breaks or donor DNA.

| Field | Value |
|-------|-------|
| **Source** | PubMed |
| **Authors** | Anzalone AV, Randolph PB, Davis JR et al. |
| **Journal** | Nature |
| **Year** | 2019 |
| **DOI** | 10.1038/s41586-019-1711-4 |
| **URL** | [https://pubmed.ncbi.nlm.nih.gov/31634902](https://pubmed.ncbi.nlm.nih.gov/31634902) |

**Abstract**: Most genetic variants that contribute to disease1 are challenging to correct efficiently and without excess byproducts2-5. Here we describe prime editing, a versatile and precise genome editing method that directly writes new genetic information into a specified DNA site using a catalytically impaired Cas9 endonuclease fused to an engineered reverse transcriptase, programmed with a prime editing g...

### 4. Saturated mutagenesis screen of M-MLV reverse transcriptase identifies variants enhancing prime editing efficiency

| Field | Value |
|-------|-------|
| **Source** | bioRxiv |
| **Authors** | Li, H., Wang, Y., Zhang, C. et al. |
| **Journal** | bioRxiv (preprint) |
| **Year** | 2026 |
| **DOI** | 10.64898/2026.07.06.736660 |
| **URL** | [https://biorxiv.org/content/10.64898/2026.07.06.736660](https://biorxiv.org/content/10.64898/2026.07.06.736660) |

**Abstract**: Prime editing enables the precise modification of genomes, thereby holding great potential for the treatment of genetic diseases. Despite substantial advancements in prime editing technology and the initiation of the first clinical trial for treating chronic granulomatous disease, further enhancement of editing efficiency across edit types is still urgently needed. Here, we developed a compact pri...

---

## Citation Graph

**Nodes**: 4  \n\
**Edges** (internal citations): 3

```json
{
  "nodes": [
    {
      "id": "27096365",
      "title": "Programmable editing of a target base in genomic DNA without double-stranded DNA cleavage.",
      "source": "PubMed",
      "year": "2016"
    },
    {
      "id": "29160308",
      "title": "Programmable base editing of A\\u2022T to G\\u2022C in genomic DNA without DNA cleavage.",
      "source": "PubMed",
      "year": "2017"
    },
    {
      "id": "31634902",
      "title": "Search-and-replace genome editing without double-strand breaks or donor DNA.",
      "source": "PubMed",
      "year": "2019"
    },
    {
      "id": "10.64898/2026.07.06.736660",
      "title": "Saturated mutagenesis screen of M-MLV reverse transcriptase identifies variants enhancing prime editing efficiency",
      "source": "bioRxiv",
      "year": "2
... (see citation_graph.json for full graph)
```

---

---
*ClawBio Lit Synthesizer is a research tool. Always verify citations independently before use in publications.*
""", 'results.json': """\
[
  {
    "source": "PubMed",
    "pmid": "27096365",
    "title": "Programmable editing of a target base in genomic DNA without double-stranded DNA cleavage.",
    "authors": [
      "Komor AC",
      "Kim YB",
      "Packer MS",
      "Zuris JA",
      "Liu DR"
    ],
    "journal": "Nature",
    "year": "2016",
    "abstract": "Current genome-editing technologies introduce double-stranded (ds) DNA breaks at a target locus as the first step to gene correction. Although most genetic diseases arise from point mutations, current approaches to point mutation correction are inefficient and typically induce an abundance of random insertions and deletions (indels) at the target locus resulting from the cellular response to dsDNA breaks. Here we report the development of 'base editing', a new approach to genome editing that enables the direct, irreversible conversion of one target DNA base into another in a programmable manner, without requiring dsDNA backbone cleavage or a donor template. We engineered fusions of CRISPR/Cas9 and a cytidine deaminase enzyme that retain the ability to be programmed with a guide RNA, do not induce dsDNA breaks, and mediate the direct conversion of cytidine to uridine, thereby effecting a C\\u2192T (or G\\u2192A) substitution. The resulting 'base editors' convert cytidines within a window of approximately five nucleotides, and can efficiently correct a variety of point mutations relevant to human disease. In four transformed human and murine cell lines, second- and third-generation base editors that fuse uracil glycosylase inhibitor, and that use a Cas9 nickase targeting the non-edited strand, manipulate the cellular DNA repair response to favour desired base-editing outcomes, resulting in permanent correction of ~15-75% of total cellular DNA with minimal (typically \\u22641%) indel formation. Base editing expands the scope and efficiency of genome editing of point mutations.",
    "doi": "10.1038/nature17946",
    "url": "https://pubmed.ncbi.nlm.nih.gov/27096365",
    "citations": []
  },
  {
    "source": "PubMed",
    "pmid": "29160308",
    "title": "Programmable base editing of A\\u2022T to G\\u2022C in genomic DNA without DNA cleavage.",
    "authors": [
      "Gaudelli NM",
      "Komor AC",
      "Rees HA",
      "Packer MS",
      "Badran AH",
      "Bryson DI",
      "Liu DR"
    ],
    "journal": "Nature",
    "year": "2017",
    "abstract": "The spontaneous deamination of cytosine is a major source of transitions from C\\u2022G to T\\u2022A base pairs, which account for half of known pathogenic point mutations in humans. The ability to efficiently convert targeted A\\u2022T base pairs to G\\u2022C could therefore advance the study and treatment of genetic diseases. The deamination of adenine yields inosine, which is treated as guanine by polymerases, but no enzymes are known to deaminate adenine in DNA. Here we describe adenine base editors (ABEs) that mediate the conversion of A\\u2022T to G\\u2022C in genomic DNA. We evolved a transfer RNA adenosine deaminase to operate on DNA when fused to a catalytically impaired CRISPR-Cas9 mutant. Extensive directed evolution and protein engineering resulted in seventh-generation ABEs that convert targeted A\\u2022T base pairs efficiently to G\\u2022C (approximately 50% efficiency in human cells) with high product purity (typically at least 99.9%) and low rates of indels (typically no more than 0.1%). ABEs introduce point mutations more efficiently and cleanly, and with less off-target genome modification, than a current Cas9 nuclease-based method, and can install disease-correcting or disease-suppressing mutations in human cells. Together with previous base editors, ABEs enable the direct, programmable introduction of all four transition mutations without double-stranded DNA cleavage.",
    "doi": "10.1038/nature24644",
    "url": "https://pubmed.ncbi.nlm.nih.gov/29160308",
    "citations": [
      "27096365"
    ]
  },
  {
    "source": "PubMed",
    "pmid": "31634902",
    "title": "Search-and-replace genome editing without double-strand breaks or donor DNA.",
    "authors": [
      "Anzalone AV",
      "Randolph PB",
      "Davis JR",
      "Sousa AA",
      "Koblan LW",
      "Levy JM",
      "Chen PJ",
      "Wilson C",
      "Newby GA",
      "Raguram A",
      "Liu DR"
    ],
    "journal": "Nature",
    "year": "2019",
    "abstract": "Most genetic variants that contribute to disease1 are challenging to correct efficiently and without excess byproducts2-5. Here we describe prime editing, a versatile and precise genome editing method that directly writes new genetic information into a specified DNA site using a catalytically impaired Cas9 endonuclease fused to an engineered reverse transcriptase, programmed with a prime editing guide RNA (pegRNA) that both specifies the target site and encodes the desired edit. We performed more than 175 edits in human cells, including targeted insertions, deletions, and all 12 types of point mutation, without requiring double-strand breaks or donor DNA templates. We used prime editing in human cells to correct, efficiently and with few byproducts, the primary genetic causes of sickle cell disease (requiring a transversion in HBB) and Tay-Sachs disease (requiring a deletion in HEXA); to install a protective transversion in PRNP; and to insert various tags and epitopes precisely into target loci. Four human cell lines and primary post-mitotic mouse cortical neurons support prime editing with varying efficiencies. Prime editing shows higher or similar efficiency and fewer byproducts than homology-directed repair, has complementary strengths and weaknesses compared to base editing, and induces much lower off-target editing than Cas9 nuclease at known Cas9 off-target sites. Prime editing substantially expands the scope and capabilities of genome editing, and in principle could correct up to 89% of known genetic variants associated with human diseases.",
    "doi": "10.1038/s41586-019-1711-4",
    "url": "https://pubmed.ncbi.nlm.nih.gov/31634902",
    "citations": [
      "27096365",
      "29160308"
    ]
  },
  {
    "source": "bioRxiv",
    "pmid": null,
    "title": "Saturated mutagenesis screen of M-MLV reverse transcriptase identifies variants enhancing prime editing efficiency",
    "authors": [
      "Li, H.",
      "Wang, Y.",
      "Zhang, C.",
      "Tun, T. T.",
      "Yu, S.",
      "Hu, C.",
      "Yu, H."
    ],
    "journal": "bioRxiv (preprint)",
    "year": "2026",
    "abstract": "Prime editing enables the precise modification of genomes, thereby holding great potential for the treatment of genetic diseases. Despite substantial advancements in prime editing technology and the initiation of the first clinical trial for treating chronic granulomatous disease, further enhancement of editing efficiency across edit types is still urgently needed. Here, we developed a compact prime editor, PE2{Delta}R, by deleting the RNase H domain of the MMLV reverse transcriptase (MMLV-RT). We then conducted a saturated mutagenesis screen targeting two DNA interacting regions within the PE2{Delta}R-RT Fingers domain. By integrating three highly effective mutations (I61R, V101R, S67W) into PEmax lacking RNase H domain (termed PEmax{Delta}RM3), we achieved up to a 90% increase in editing efficiency across editing types compared to PEmax. Structural modelling using AlphaFold 3 suggests that these mutations enhance primer-template stabilization and guide the RNA/DNA hybrid into a catalytically favourable trajectory, providing a mechanistic explanation for the enhanced activity. Taken together, our study demonstrates proof-of-concept for the application of unbiased mutagenesis screen to identify novel mutations that enhance prime editor performance. Furthermore, we discovered that RT variants (I61R, V101R, S67W) synergize with PEmax and epegRNA to improve prime editing efficiency across edit types, with the strongest improvement observed in introducing small deletions.",
    "doi": "10.64898/2026.07.06.736660",
    "url": "https://biorxiv.org/content/10.64898/2026.07.06.736660",
    "citations": []
  }
]\
""", 'tables/papers.csv': """\
source,pmid,title,authors,journal,year,doi,url
PubMed,27096365,Programmable editing of a target base in genomic DNA without double-stranded DNA cleavage.,Komor AC; Kim YB; Packer MS; Zuris JA; Liu DR,Nature,2016,10.1038/nature17946,https://pubmed.ncbi.nlm.nih.gov/27096365
PubMed,29160308,Programmable base editing of A•T to G•C in genomic DNA without DNA cleavage.,Gaudelli NM; Komor AC; Rees HA; Packer MS; Badran AH; Bryson DI; Liu DR,Nature,2017,10.1038/nature24644,https://pubmed.ncbi.nlm.nih.gov/29160308
PubMed,31634902,Search-and-replace genome editing without double-strand breaks or donor DNA.,Anzalone AV; Randolph PB; Davis JR; Sousa AA; Koblan LW; Levy JM; Chen PJ; Wilson C; Newby GA; Raguram A; Liu DR,Nature,2019,10.1038/s41586-019-1711-4,https://pubmed.ncbi.nlm.nih.gov/31634902
bioRxiv,,Saturated mutagenesis screen of M-MLV reverse transcriptase identifies variants enhancing prime editing efficiency,"Li, H.; Wang, Y.; Zhang, C.; Tun, T. T.; Yu, S.; Hu, C.; Yu, H.",bioRxiv (preprint),2026,10.64898/2026.07.06.736660,https://biorxiv.org/content/10.64898/2026.07.06.736660
"""})
