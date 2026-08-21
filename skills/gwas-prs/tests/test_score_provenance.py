"""Provenance tests for the bundled curated demo panels.

Context (issue #356). The six files under ``skills/gwas-prs/data/`` are
ClawBio-curated illustrative panels of well-established trait-associated loci.
They are *not* PGS Catalog scoring files, even though they were named with the
PGS Catalog harmonised convention and keyed by real PGS accessions.

Two classes of defect followed from that, and these tests pin both shut:

1. Citation drift. ``CURATED_SCORES`` in ``gwas_prs.py`` and the standalone
   ``curated_scores.json`` carried different publications for the same panel,
   and between them named three papers that do not exist as cited (a
   ventilator-pneumonia nursing paper, the ExAC paper, and a body-fat
   distribution meta-analysis).

2. Silent score substitution. The non-demo fetch path consulted
   ``DATA_DIR/{pgs_id}_hmPOS_{build}.txt`` before contacting the API, so
   ``--pgs-id PGS000013`` returned the bundled 8-variant type 2 diabetes panel
   instead of the real PGS000013 (Khera 2018, coronary artery disease,
   6,630,150 variants) and never touched the network.
"""

from __future__ import annotations

import json
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import pytest

SKILL_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = SKILL_DIR / "data"


def _load_engine():
    spec = spec_from_file_location("gwas_prs", SKILL_DIR / "gwas_prs.py")
    mod = module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


ENGINE = _load_engine()
CURATED = ENGINE.CURATED_SCORES

# Every field here was verified against PubMed E-utilities on 2026-08-21.
# Do not edit without re-resolving the PMID; a citation that names the wrong
# paper is the defect this file exists to prevent.
VERIFIED = {
    "PGS000013": ("Type 2 diabetes", 8, "24520119", "Vassy", "Diabetes"),
    "PGS000011": ("Atrial fibrillation", 12, "25123217", "Tada", "Stroke"),
    "PGS000004": ("Coronary artery disease", 46, "27655226", "Abraham", "Eur Heart J"),
    "PGS000001": ("Breast cancer", 77, "25855707", "Mavaddat", "J Natl Cancer Inst"),
    "PGS000057": ("Prostate cancer", 147, "29892016", "Schumacher", "Nat Genet"),
    "PGS000039": ("BMI", 97, "25673413", "Locke", "Nature"),
}

PANEL_IDS = sorted(VERIFIED)


def _panel_path(pgs_id: str) -> Path:
    return DATA_DIR / f"{pgs_id}_hmPOS_GRCh37.txt"


def _headers(path: Path) -> dict[str, str]:
    out = {}
    for line in path.read_text().splitlines():
        if not line.startswith("#"):
            break
        if "=" in line:
            key, _, value = line[1:].partition("=")
            out[key.strip()] = value.strip()
    return out


class TestCitationsAreReal:
    @pytest.mark.parametrize("pgs_id", PANEL_IDS)
    def test_pmid_matches_the_verified_paper(self, pgs_id):
        expected_pmid = VERIFIED[pgs_id][2]
        assert CURATED[pgs_id]["pmid"] == expected_pmid

    @pytest.mark.parametrize("pgs_id", PANEL_IDS)
    def test_publication_names_the_right_author_and_journal(self, pgs_id):
        _, _, _, author, journal = VERIFIED[pgs_id]
        publication = CURATED[pgs_id]["publication"]
        assert author in publication
        assert journal in publication

    def test_no_retired_pmid_survives_anywhere_in_the_skill(self):
        """The three PMIDs that named unrelated papers must be gone.

        25087147 is a ventilator-associated-pneumonia nursing paper,
        27535533 is the ExAC paper, 30239722 is a body-fat meta-analysis.
        A banned-string check is the only guard that fires for a stray copy
        nobody edited.
        """
        retired = ("25087147", "27535533", "30239722")
        offenders = []
        for path in SKILL_DIR.rglob("*"):
            if not path.is_file() or ".git" in path.parts:
                continue
            if path.name == Path(__file__).name:
                continue
            try:
                text = path.read_text()
            except (UnicodeDecodeError, OSError):
                continue
            for pmid in retired:
                if pmid in text:
                    offenders.append(f"{path.relative_to(SKILL_DIR)}:{pmid}")
        assert offenders == []


class TestTheTwoSourcesAgree:
    """`curated_scores.json` is not read by any code, so it drifted."""

    def test_json_and_engine_cover_the_same_panels(self):
        data = json.loads((SKILL_DIR / "curated_scores.json").read_text())
        assert sorted(data) == sorted(CURATED)

    @pytest.mark.parametrize("pgs_id", PANEL_IDS)
    def test_json_and_engine_agree_on_every_shared_field(self, pgs_id):
        data = json.loads((SKILL_DIR / "curated_scores.json").read_text())
        j, c = data[pgs_id], CURATED[pgs_id]
        assert j["trait"] == c["trait"]
        assert j["variants_count"] == c["variants_count"]
        assert j["pmid"] == c["pmid"]
        assert j["reference_distribution"] == c["reference_distribution"]


class TestPanelFilesDeclareWhatTheyAre:
    @pytest.mark.parametrize("pgs_id", PANEL_IDS)
    def test_file_is_marked_as_a_curated_panel_not_a_catalog_download(self, pgs_id):
        headers = _headers(_panel_path(pgs_id))
        assert headers.get("clawbio_panel") == "curated_demo"

    @pytest.mark.parametrize("pgs_id", PANEL_IDS)
    def test_file_header_agrees_with_the_engine(self, pgs_id):
        trait, count, pmid, _, _ = VERIFIED[pgs_id]
        headers = _headers(_panel_path(pgs_id))
        assert headers["trait_reported"] == trait
        assert int(headers["variants_number"]) == count
        assert headers["clawbio_pmid"] == pmid

    @pytest.mark.parametrize("pgs_id", PANEL_IDS)
    def test_declared_variant_count_matches_the_actual_rows(self, pgs_id):
        path = _panel_path(pgs_id)
        rows = [
            ln for ln in path.read_text().splitlines()
            if ln and not ln.startswith("#")
        ]
        # first non-comment line is the column header
        assert len(rows) - 1 == VERIFIED[pgs_id][1]

    @pytest.mark.parametrize("pgs_id", PANEL_IDS)
    def test_file_does_not_claim_to_be_a_pgs_catalog_score(self, pgs_id):
        text = _panel_path(pgs_id).read_text()
        assert "#pgs_id=" not in text, (
            "a curated panel must not present a PGS Catalog accession as its own id"
        )


class TestCuratedPanelsNeverShadowARealScore:
    """The substitution bug from issue #356.

    A curated panel sitting in the download cache path must never be served
    in answer to a real request for that accession.
    """

    @pytest.mark.parametrize("pgs_id", PANEL_IDS)
    def test_bundled_panels_are_recognised_as_curated(self, pgs_id):
        assert ENGINE.is_curated_demo_panel(_panel_path(pgs_id)) is True

    def test_a_genuine_catalog_file_is_not_flagged(self, tmp_path):
        genuine = tmp_path / "PGS000013_hmPOS_GRCh37.txt"
        genuine.write_text(
            "#pgs_id=PGS000013\n"
            "#trait_reported=Coronary artery disease\n"
            "rsID\tchr_name\tchr_position\teffect_allele\tother_allele\teffect_weight\n"
            "rs1234\t1\t1000\tA\tG\t0.1\n"
        )
        assert ENGINE.is_curated_demo_panel(genuine) is False

    def test_missing_file_is_not_flagged(self, tmp_path):
        assert ENGINE.is_curated_demo_panel(tmp_path / "absent.txt") is False

    def test_guard_would_fail_if_deleted(self):
        """If `is_curated_demo_panel` always returned False the bundled panels
        would be served as catalog scores again, so assert it discriminates."""
        genuine_like = DATA_DIR / "PGS000013_hmPOS_GRCh37.txt"
        assert ENGINE.is_curated_demo_panel(genuine_like) != ENGINE.is_curated_demo_panel(
            SKILL_DIR / "demo_patient_prs.txt"
        )
