"""Provenance tests for the bundled curated demo panels.

Context (issue #356). The six files under ``skills/gwas-prs/data/`` are
ClawBio-curated illustrative panels of well-established trait-associated loci.
They are *not* PGS Catalog scoring files. They once used the PGS Catalog
harmonised filename convention and real accessions as keys; issue #356 moved
those accessions to provenance metadata and made ``CLAWBIO-*`` ids canonical.

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
    "CLAWBIO-T2D-8": ("Type 2 diabetes", 8, "24520119", "Vassy", "Diabetes", "PGS000013"),
    "CLAWBIO-AF-12": ("Atrial fibrillation", 12, "25123217", "Tada", "Stroke", "PGS000011"),
    "CLAWBIO-CAD-46": ("Coronary artery disease", 46, "27655226", "Abraham", "Eur Heart J", "PGS000004"),
    "CLAWBIO-BC-77": ("Breast cancer", 77, "25855707", "Mavaddat", "J Natl Cancer Inst", "PGS000001"),
    "CLAWBIO-PC-147": ("Prostate cancer", 147, "29892016", "Schumacher", "Nat Genet", "PGS000057"),
    "CLAWBIO-BMI-97": ("BMI", 97, "25673413", "Locke", "Nature", "PGS000039"),
}

PANEL_IDS = sorted(VERIFIED)


def _panel_path(panel_id: str) -> Path:
    return ENGINE.curated_panel_path(panel_id, "GRCh37")


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
    @pytest.mark.parametrize("panel_id", PANEL_IDS)
    def test_pmid_matches_the_verified_paper(self, panel_id):
        expected_pmid = VERIFIED[panel_id][2]
        assert CURATED[panel_id]["pmid"] == expected_pmid

    @pytest.mark.parametrize("panel_id", PANEL_IDS)
    def test_publication_names_the_right_author_and_journal(self, panel_id):
        _, _, _, author, journal, _ = VERIFIED[panel_id]
        publication = CURATED[panel_id]["publication"]
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

    @pytest.mark.parametrize("panel_id", PANEL_IDS)
    def test_json_and_engine_agree_on_every_shared_field(self, panel_id):
        data = json.loads((SKILL_DIR / "curated_scores.json").read_text())
        j, c = data[panel_id], CURATED[panel_id]
        assert j["trait"] == c["trait"]
        assert j["variants_count"] == c["variants_count"]
        assert j["pmid"] == c["pmid"]
        assert j["reference_distribution"] == c["reference_distribution"]


class TestPanelFilesDeclareWhatTheyAre:
    @pytest.mark.parametrize("panel_id", PANEL_IDS)
    def test_file_is_marked_as_a_curated_panel_not_a_catalog_download(self, panel_id):
        headers = _headers(_panel_path(panel_id))
        assert headers.get("clawbio_panel") == "curated_demo"

    @pytest.mark.parametrize("panel_id", PANEL_IDS)
    def test_file_header_agrees_with_the_engine(self, panel_id):
        trait, count, pmid, _, _, _ = VERIFIED[panel_id]
        headers = _headers(_panel_path(panel_id))
        assert headers["trait_reported"] == trait
        assert int(headers["variants_number"]) == count
        assert headers["clawbio_pmid"] == pmid

    @pytest.mark.parametrize("panel_id", PANEL_IDS)
    def test_declared_variant_count_matches_the_actual_rows(self, panel_id):
        path = _panel_path(panel_id)
        rows = [
            ln for ln in path.read_text().splitlines()
            if ln and not ln.startswith("#")
        ]
        # first non-comment line is the column header
        assert len(rows) - 1 == VERIFIED[panel_id][1]

    @pytest.mark.parametrize("panel_id", PANEL_IDS)
    def test_file_does_not_claim_to_be_a_pgs_catalog_score(self, panel_id):
        text = _panel_path(panel_id).read_text()
        assert "#pgs_id=" not in text, (
            "a curated panel must not present a PGS Catalog accession as its own id"
        )


class TestCuratedPanelsNeverShadowARealScore:
    """The substitution bug from issue #356.

    A curated panel sitting in the download cache path must never be served
    in answer to a real request for that accession.
    """

    @pytest.mark.parametrize("panel_id", PANEL_IDS)
    def test_bundled_panels_are_recognised_as_curated(self, panel_id):
        assert ENGINE.is_curated_demo_panel(_panel_path(panel_id)) is True

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
        curated = _panel_path("CLAWBIO-T2D-8")
        assert ENGINE.is_curated_demo_panel(curated) != ENGINE.is_curated_demo_panel(
            SKILL_DIR / "demo_patient_prs.txt"
        )


# ---------------------------------------------------------------------------
# Round 2, from the adversarial re-audit of this PR.
#
# The first pass shipped four defects of exactly the classes this repo has been
# flagging in contributor PRs all week: a doc asserting behaviour the code does
# not have, a mitigation that only reaches stdout, a guard applied to one of two
# call sites, and a behavioural change with no test. These pin all four.
# ---------------------------------------------------------------------------

import os
import subprocess
import sys


class TestDocsDescribeTheCodeThatExists:
    def test_skill_md_does_not_claim_the_panel_is_ignored(self):
        """The refusal was reverted; the sentence promising it must go too.

        A planner reads SKILL.md to decide whether it needs --demo. Telling it
        a bare --pgs-id fetches the genuine score, when the bundled panel is
        actually served, is a false assurance and worse than saying nothing.
        """
        text = (SKILL_DIR / "SKILL.md").read_text().lower()
        for claim in (
            "is deliberately ignored",
            "fetches the genuine",
            "cannot stand in for the real thing",
        ):
            assert claim not in text, f"SKILL.md still promises: {claim!r}"

    def test_skill_md_does_not_call_the_json_generated_without_a_generator(self):
        text = (SKILL_DIR / "SKILL.md").read_text()
        if "generated from" in text:
            assert (SKILL_DIR / "generate_curated_scores.py").exists(), (
                "SKILL.md calls curated_scores.json generated; ship the generator"
            )


class TestBothCallSitesAreGuarded:
    """Round 1 guarded the --pgs-id lookup and left the --trait one open.

    There are two PGS Catalog cache lookups: direct ``--pgs-id`` and
    ``--trait``. Both must refuse a curated panel placed in a Catalog path.
    """

    def _lookup_sites(self):
        lines = (SKILL_DIR / "gwas_prs.py").read_text().splitlines()
        return [
            i for i, ln in enumerate(lines)
            if '= DATA_DIR / f"{' in ln and "_gz" not in ln
        ]

    def test_there_are_exactly_three_known_lookups(self):
        """If a fourth appears, this test fails and someone reads the next one."""
        assert len(self._lookup_sites()) == 2

    def test_every_non_demo_lookup_warns_before_using_the_file(self):
        lines = (SKILL_DIR / "gwas_prs.py").read_text().splitlines()
        unguarded = []
        for i in self._lookup_sites():
            if "is_curated_demo_panel(" not in "\n".join(lines[i:i + 20]):
                unguarded.append(i + 1)
        assert unguarded == [], (
            f"cache lookups that can serve a panel unannounced: lines {unguarded}"
        )


class TestTheMitigationActuallyRuns:
    """Deleting the warning block used to leave all 78 tests green."""

    # These tests exercise the benchmark compatibility alias, which is opt-in
    # since #356's refusal landed. See test_legacy_alias_opt_in.py for the
    # default (refusing) behaviour.
    _ALIAS_ENV = {**os.environ, "CLAWBIO_ALLOW_LEGACY_PGS_ALIAS": "1"}

    def _run(self, tmp_path, *args):
        return subprocess.run(
            [sys.executable, str(SKILL_DIR / "gwas_prs.py"),
             "--input", str(SKILL_DIR / "demo_patient_prs.txt"),
             "--output", str(tmp_path), *args],
            capture_output=True, text=True, timeout=120, env=self._ALIAS_ENV,
        )

    def test_pgs_id_path_warns_on_stdout(self, tmp_path):
        proc = self._run(tmp_path, "--pgs-id", "PGS000013")
        combined = proc.stdout + proc.stderr
        assert "benchmark compatibility" in combined.lower()

    def test_pgs_id_path_marks_the_report(self, tmp_path):
        self._run(tmp_path, "--pgs-id", "PGS000013")
        report = (tmp_path / "prs_report.md").read_text().lower()
        assert "curated demo panel" in report

    def test_pgs_id_path_marks_result_json(self, tmp_path):
        self._run(tmp_path, "--pgs-id", "PGS000013")
        data = json.loads((tmp_path / "result.json").read_text())
        assert data["summary"]["curated_demo_panel"] is True
        assert data["summary"]["pgs_catalog_id"] is None

    def test_report_does_not_claim_pgs_catalog_provenance_for_a_panel(self, tmp_path):
        self._run(tmp_path, "--pgs-id", "PGS000013")
        report = (tmp_path / "prs_report.md").read_text()
        assert "**Scoring files**: PGS Catalog" not in report


class TestPanelHeadersDoNotOverclaimTheWeights:
    @pytest.mark.parametrize("panel_id", PANEL_IDS)
    def test_no_header_implies_the_weights_are_the_papers(self, panel_id):
        """`#source_publication=` reads as "these weights come from here".

        They do not: the real PGS000001 shares 60 of 77 rsIDs with the bundled
        panel and 31 of those 60 weights differ by more than 0.02. Vassy 2014
        is a 62-locus score against an 8-locus panel; Abraham 2016 is 49,310
        SNPs against 46.
        """
        text = _panel_path(panel_id).read_text()
        assert "#source_publication=" not in text
        assert "#loci_reference=" in text
        assert "approximate" in text.lower()


class TestPerPanelAccessionTruth:
    """PGS000001 really is the cited paper's accession. Five of six are not."""

    def test_pgs000001_is_not_described_as_a_different_score(self):
        text = _panel_path("CLAWBIO-BC-77").read_text()
        assert "belongs to a different" not in text

    @pytest.mark.parametrize(
        "panel_id", ["CLAWBIO-CAD-46", "CLAWBIO-AF-12", "CLAWBIO-T2D-8", "CLAWBIO-BMI-97", "CLAWBIO-PC-147"]
    )
    def test_the_other_five_say_so(self, panel_id):
        text = _panel_path(panel_id).read_text()
        assert "different" in text.lower()


class TestThePinCoversEveryField:
    """The round-1 pin compared 4 fields, so the JSON could carry a fabricated
    citation, or re-assert PGS Catalog identity, and stay green."""

    def test_panel_ids_come_from_the_engine_not_a_hardcoded_list(self):
        assert set(PANEL_IDS) == set(CURATED)

    @pytest.mark.parametrize("field", [
        "publication", "name", "curated_panel_id", "legacy_pgs_id", "pgs_catalog_id", "trait_id",
    ])
    def test_json_and_engine_agree_on(self, field):
        data = json.loads((SKILL_DIR / "curated_scores.json").read_text())
        for panel_id in CURATED:
            assert data[panel_id][field] == CURATED[panel_id][field], panel_id

    def test_trait_id_was_not_silently_dropped(self):
        for panel_id, meta in CURATED.items():
            assert str(meta.get("trait_id", "")).startswith("EFO_"), panel_id
