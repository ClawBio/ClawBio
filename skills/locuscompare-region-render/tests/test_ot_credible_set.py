"""Tests for ot_credible_set.py: choosing the QTL and GWAS from an Open Targets credible set.

Offline tests replay Open Targets API responses and GWAS Catalog FTP listings
recorded on 2026-10-13 (tests/fixtures/open_targets/). The one live test is
skipped unless RUN_LIVE_TESTS=1.
"""
from __future__ import annotations

import ast
import copy
import importlib.util
import json
import os
import sys
import urllib.error
from pathlib import Path

import pytest
import yaml

SKILL_DIR = Path(__file__).resolve().parent.parent
FIX = Path(__file__).resolve().parent / "fixtures" / "open_targets"


def _load():
    name = "ot_credible_set_under_test"
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, SKILL_DIR / "ot_credible_set.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


ocs = _load()

IBD = "2f6755df0a7737586c9ff36fa2643c61"
CHOL = "e723d7b2d797c4c3d9e3889ef3e79997"
PCSK9_LUNG_GE = "227c96a4630d2dbb1cb3cf1e637ebc49"


def _fixture(kind: str, slid: str) -> dict:
    return json.loads((FIX / f"{kind}__{slid}.json").read_text())


def _listing(files: list[str]) -> str:
    return "".join(f'<tr><td><a href="{f}">{f}</a></td></tr>' for f in files)


class FakeOT:
    """Replays recorded responses; `overrides` replaces whole responses by key."""

    def __init__(self, probe_overrides: dict | None = None, coloc_overrides: dict | None = None):
        self.probes = json.loads((FIX / "ftp_probe_2026-10-13.json").read_text())
        self.probe_overrides = probe_overrides or {}
        self.coloc_overrides = coloc_overrides or {}
        self.get_calls: list[str] = []

    def post(self, query, variables):
        slid = variables["id"]
        types = variables.get("types")
        if types is None:
            return _fixture("credibleSet", slid)
        kind = "coloc_gwas" if types == ["gwas"] else "coloc_molqtl"
        resp = self.coloc_overrides.get((kind, slid)) or _fixture(kind, slid)
        if variables.get("index", 0) > 0:
            resp = copy.deepcopy(resp)
            resp["data"]["credibleSet"]["colocalisation"]["rows"] = []
        return resp

    def get(self, url, timeout):
        self.get_calls.append(url)
        acc = url.rstrip("/").split("/")[-2]
        if acc in self.probe_overrides:
            o = self.probe_overrides[acc]
            if isinstance(o, Exception):
                raise o
            return o
        if acc in self.probes:
            return 200, _listing(self.probes[acc]["files"])
        raise TimeoutError("not recorded")


@pytest.fixture
def fake(monkeypatch):
    f = FakeOT()
    monkeypatch.setattr(ocs, "_post", f.post)
    monkeypatch.setattr(ocs, "_get", f.get)
    return f


def _synthetic_qtl(i: int, **kw) -> dict:
    study = {
        "id": kw.get("study_id", f"gtex-v10_ge_lung_gene{i}"), "studyType": "eqtl",
        "projectId": kw.get("project", "GTEx-v10"), "traitFromSource": f"ENSG{i:011d}",
        "publicationFirstAuthor": None, "publicationDate": None, "nSamples": 652,
        "hasSumstats": True,
        "summarystatsLocation": kw.get("loc", "gs://bucket/QTS000015/QTD000271/QTD000271.cc.parquet"),
        "condition": "naive",
        "target": {"id": kw.get("target", f"ENSG{i:011d}"), "approvedSymbol": kw.get("symbol", f"GENE{i}")},
        "biosample": {"biosampleId": kw.get("bios", "UBERON_0008952"), "biosampleName": "lung"},
        "diseases": [],
    }
    return {"h4": kw.get("h4", 0.9 - i * 0.001), "betaRatioSignAverage": 1,
            "colocalisationMethod": "COLOC_PIP_ECAVIAR",
            "otherStudyLocus": {"studyLocusId": kw.get("slid", f"{i:032x}"), "study": study}}


def _molqtl_response(rows: list[dict], count: int | None = None) -> dict:
    return {"data": {"credibleSet": {"colocalisation": {
        "count": len(rows) if count is None else count, "rows": rows}}}}


# --------------------------------------------------------------------------- 2f6755df: one QTL
def test_single_qtl_page_gives_the_example_08_config(fake, tmp_path, capsys):
    assert ocs.main(["options", f"https://platform.opentargets.org/credible-set/{IBD}"]) == 0
    text = capsys.readouterr().out
    assert "Page GWAS file check (GCST90225550): yes" in text
    assert "could stand in" not in text  # page study readable: no substitutes offered

    out = tmp_path / "config.yaml"
    assert ocs.main(["config", IBD, "--qtl", "1", "--out", str(out)]) == 0
    got = yaml.safe_load(out.read_text())
    want = yaml.safe_load((SKILL_DIR / "examples" / "08_open_targets_credible_set" / "config.yaml").read_text())
    # Example 08 predates lead candidates and pins its lead; the helper now also writes
    # the credible set's 15 members, PIP high to low, for the render-time fallback.
    cands = got["lead"].pop("candidates")
    assert got["lead"] == want["lead"]
    assert len(cands) == 15 and cands[0]["variant_id"] == "1_67240275_G_A"
    assert [c["pip"] for c in cands] == sorted((c["pip"] for c in cands), reverse=True)
    assert got["exposure"]["fetch"] == want["exposure"]["fetch"]
    assert got["outcome"]["fetch"] == want["outcome"]["fetch"]
    assert got["provenance"]["substitutions"] == []
    assert "caveats" not in got


# --------------------------------------------------------------------------- e723d7b2: many QTLs
def test_many_qtls_table_keeps_kinds_apart_and_recommends_pcsk9_lung(fake):
    a = ocs.analyse(CHOL)
    assert a.qtl_count == 10
    # Gene-level and exon-level rows of the same gene and tissue are different
    # measurements, so nothing folds here.
    assert a.qtl_tally == {"shown": 10, "hidden_by_cap": 0, "excluded": 0}
    rec = a.qtl_rec.row
    assert rec.study_locus_id == PCSK9_LUNG_GE
    assert rec.kind == "eqtl" and rec.qtd == "QTD000271"
    assert rec.older_release
    assert any("GTEx v8" in r and "older release" in r for r in a.qtl_rec.reasons)
    assert a.selected is None and a.page_probe is None  # nothing chosen yet, nothing probed
    text = ocs.render_text(a)
    assert "exon eQTL" in text and "eQTL (gene)" in text
    assert "older release: GTEx v8, n=510" in text


def test_fold_only_rows_of_the_same_kind():
    index = ocs.load_eqtl_index()
    ge_a = _synthetic_qtl(1, study_id="gtex-v10_ge_lung_x", target="ENSG1", h4=0.99, slid="a" * 32)
    ge_b = _synthetic_qtl(2, study_id="gtex-v10_ge_lung_y", target="ENSG1", h4=0.98, slid="b" * 32)
    exon = _synthetic_qtl(3, study_id="gtex-v10_exon_lung_x", target="ENSG1", h4=0.97, slid="c" * 32)
    rows = [ocs.qtl_row(r, index) for r in (ge_a, ge_b, exon)]
    kept = ocs.build_qtl_table(rows)
    assert [r.study_id for r in kept] == ["gtex-v10_ge_lung_x", "gtex-v10_exon_lung_x"]
    assert rows[1].status == "excluded:collapsed"
    assert ocs.check_partition(rows, 3, "t") == {"shown": 2, "hidden_by_cap": 0, "excluded": 1}


def test_named_gene_is_listed_first_and_preferred(fake):
    a = ocs.analyse(CHOL, gene="BSND")
    shown = [r for r in a.qtl_kept if r.status == "shown"]
    assert [r.symbol for r in shown[:2]] == ["BSND", "BSND"]
    assert a.qtl_rec.row.symbol == "BSND"
    assert any("the gene you named" in r for r in a.qtl_rec.reasons)


def test_display_cap_states_the_hidden_count(fake, monkeypatch):
    rows = [_synthetic_qtl(i) for i in range(20)]
    fake.coloc_overrides[("coloc_molqtl", CHOL)] = _molqtl_response(rows)
    a = ocs.analyse(CHOL)
    assert a.qtl_tally == {"shown": 15, "hidden_by_cap": 5, "excluded": 0}
    assert "5 more not shown" in ocs.render_text(a)


def test_partition_mismatch_is_raised():
    index = ocs.load_eqtl_index()
    rows = [ocs.qtl_row(_synthetic_qtl(i), index) for i in range(3)]
    ocs.build_qtl_table(rows)
    with pytest.raises(ocs.PartitionError):
        ocs.check_partition(rows, 4, "t")  # Open Targets said 4, three were read
    rows[0].status = ""
    with pytest.raises(ocs.PartitionError):
        ocs.check_partition(rows, 3, "t")  # a row with no state


def test_pagination_reads_every_page(fake):
    rows = [_synthetic_qtl(i) for i in range(4)]
    pages = {0: rows[:2], 1: rows[2:]}

    def post(query, variables):
        if variables.get("types") is None:
            return _fixture("credibleSet", CHOL)
        return _molqtl_response(pages.get(variables["index"], []), count=4)

    ocs_post = ocs._post
    try:
        ocs._post = post
        count, got = ocs.fetch_colocalisations(CHOL, ocs.MOLQTL_STUDY_TYPES)
    finally:
        ocs._post = ocs_post
    assert count == 4 and len(got) == 4


@pytest.mark.parametrize("loc", [
    "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/QTD000271",
    "gs://eqtl_catalogue_processed_otar_internal_data/otar2077/GTEx_v10_majiq/sumstats/QTS000015/QTD000271/QTD000271.cc.parquet",
])
def test_dataset_id_read_from_both_location_shapes(loc):
    r = ocs.qtl_row(_synthetic_qtl(1, loc=loc), ocs.load_eqtl_index())
    assert r.qtd == "QTD000271" and r.exposure == "eqtl_catalogue"


# --------------------------------------------------------------------------- GWAS side
def test_unreadable_page_study_lists_substitutes_and_recommends_same_disease(fake):
    a = ocs.analyse(CHOL, qtl=PCSK9_LUNG_GE)
    assert a.page_probe.status == "no_index"
    assert a.gwas_count == 123
    page_rows = [r for r in a.gwas_rows if r.study_id == "GCST90013932"]
    assert [r.status for r in page_rows] == ["excluded:is_page_study"]
    assert a.gwas_rec.row.study_id == "GCST90475721"
    # GCST90013882 has a higher H4 and the same disease but no indexed file
    firth = next(r for r in a.gwas_kept if r.study_id == "GCST90013882")
    assert firth.status == "shown" and firth.probe.status == "no_index"
    assert sum(a.gwas_tally.values()) == 123
    text = ocs.render_text(a)
    assert "the plot would show GCST90475721, not the page study GCST90013932" in text


def test_unknown_probe_is_never_recommended(fake):
    fake.probe_overrides["GCST90475721"] = urllib.error.URLError("timed out")
    a = ocs.analyse(CHOL, qtl=PCSK9_LUNG_GE)
    r = next(r for r in a.gwas_kept if r.study_id == "GCST90475721")
    assert r.probe.status == "unknown"
    assert a.gwas_rec.row.study_id == "GCST90435755"


def test_failed_page_probe_keeps_the_page_study(fake, tmp_path, capsys):
    fake.probe_overrides["GCST90013932"] = TimeoutError("timed out")
    a = ocs.analyse(CHOL, qtl=PCSK9_LUNG_GE)
    assert a.page_probe.status == "unknown"
    assert a.gwas_rec is None  # no substitutes on a failed check
    assert "TimeoutError" in ocs.render_text(a)

    out = tmp_path / "c.yaml"
    assert ocs.main(["config", CHOL, "--qtl", PCSK9_LUNG_GE, "--out", str(out)]) == 0
    got = yaml.safe_load(out.read_text())
    assert got["outcome"]["fetch"]["accession"] == "GCST90013932"
    assert got["provenance"]["page_study_file_check"] == "unknown"
    assert any("is unknown" in c and "GCST90013932" in c for c in got["caveats"])

    b = ocs.analyse(CHOL, qtl=PCSK9_LUNG_GE, show_substitutes=True)
    assert b.gwas_rec is not None and b.gwas_rec.row.study_id == "GCST90475721"


@pytest.mark.parametrize("response,expected", [
    ((200, _listing(["GCST1.h.tsv.gz", "GCST1.h.tsv.gz.tbi"])), "indexed"),
    ((200, _listing(["33959723-GCST1-EFO_0003767.h.tsv.gz"])), "no_index"),
    ((200, _listing(["GCST1.h.tsv.gz"])), "no_index"),
    ((404, ""), "no_harmonised_dir"),
    ((503, ""), "unknown"),
    (urllib.error.URLError("down"), "unknown"),
    (TimeoutError("slow"), "unknown"),
])
def test_probe_states(monkeypatch, response, expected):
    def get(url, timeout):
        assert url.endswith("/GCST1/harmonised/")
        if isinstance(response, Exception):
            raise response
        return response

    monkeypatch.setattr(ocs, "_get", get)
    assert ocs.probe_gwas("GCST1").status == expected


def test_probe_of_non_gwas_catalog_study():
    assert ocs.probe_gwas("FINNGEN_R12_RX_STATIN").status == "not_gcst"


def test_harmonised_dir_url():
    assert ocs.harmonised_dir_url("GCST90475721").endswith(
        "/GCST90475001-GCST90476000/GCST90475721/harmonised/")
    assert ocs.harmonised_dir_url("GCST005194").endswith("/GCST005001-GCST006000/GCST005194/harmonised/")


# --------------------------------------------------------------------------- config
def test_config_records_both_substitutions_and_matches_example_09(fake, tmp_path, capsys):
    out = tmp_path / "c.yaml"
    assert ocs.main(["config", CHOL, "--qtl", PCSK9_LUNG_GE, "--gwas", "GCST90475721",
                     "--out", str(out)]) == 0
    got = yaml.safe_load(out.read_text())
    assert got["exposure"]["trait_label"] == "PCSK9 expression - lung (GTEx v8 via eQTL Catalogue r7)"
    subs = {s["side"]: s for s in got["provenance"]["substitutions"]}
    assert subs["outcome"]["requested"] == "GCST90013932" and subs["outcome"]["used"] == "GCST90475721"
    assert "QTD000271" in subs["exposure"]["used"]
    # The configured lead stays the page's lead; whether it is usable is decided at
    # render time, from the credible-set members written as candidates (PIP high to low).
    assert "lead" not in subs
    assert got["lead"]["variant_id"] == "1_55025188_T_C" and got["lead"]["rs_id"] == "rs4500361"
    assert [c["variant_id"] for c in got["lead"]["candidates"]] == [
        "1_55025188_T_C", "1_55021673_C_G", "1_55023869_C_T", "1_55026242_C_T", "1_55027911_C_A"]
    assert got["provenance"]["outcome_study_file_check"] == "indexed"
    assert len(got["caveats"]) == 2
    want = yaml.safe_load(
        (SKILL_DIR / "examples" / "09_open_targets_credible_set_many_qtls" / "config.yaml").read_text())
    for key in ("lead", "exposure", "outcome"):
        assert got[key] == want[key]
    assert "caveats:" in capsys.readouterr().out


def test_config_refusals(fake, tmp_path, capsys):
    out = str(tmp_path / "c.yaml")
    assert ocs.main(["config", CHOL, "--out", out]) == 2  # no --qtl
    assert ocs.main(["config", IBD, "--out", out]) == 2  # no --qtl, even when only one QTL is listed
    assert ocs.main(["config", CHOL, "--qtl", PCSK9_LUNG_GE, "--out", out]) == 2  # page unreadable, no --gwas
    assert ocs.main(["config", CHOL, "--qtl", PCSK9_LUNG_GE, "--gwas", "GCST000001",
                     "--out", out]) == 2  # not in the QTL's colocalisations
    assert ocs.main(["config", CHOL, "--qtl", PCSK9_LUNG_GE, "--gwas", "GCST90013882",
                     "--out", out]) == 2  # colocalises, but no indexed file
    assert not Path(out).exists()
    assert ocs.main(["options", "not-an-id"]) == 2


def test_json_lists_every_row_with_its_state(fake, tmp_path):
    out = tmp_path / "rows.json"
    assert ocs.main(["options", CHOL, "--qtl", PCSK9_LUNG_GE, "--json", str(out)]) == 0
    d = json.loads(out.read_text())
    assert len(d["qtl"]["rows"]) == 10 and len(d["gwas"]["rows"]) == 123
    assert {r["status"] for r in d["gwas"]["rows"]} == {"shown", "hidden_by_cap", "excluded:is_page_study"}


# --------------------------------------------------------------------------- one rule, one map
def _dict_literal(path: Path, name: str) -> dict:
    tree = ast.parse(path.read_text())
    for node in tree.body:
        target = getattr(node, "target", None) or (node.targets[0] if isinstance(node, ast.Assign) else None)
        if isinstance(target, ast.Name) and target.id == name:
            return ast.literal_eval(node.value)
    raise AssertionError(f"{name} not found in {path.name}")


def test_qtl_kind_map_matches_the_renderer():
    """The helper carries a copy of the renderer's quantification -> QTL-kind map
    (it must stay standard-library only); the copies may not drift."""
    renderer = _dict_literal(SKILL_DIR / "locuscompare_region_render.py", "_QUANT_METHOD_MODALITY")
    assert ocs._QUANT_METHOD_MODALITY == renderer
    tags = _dict_literal(SKILL_DIR / "locuscompare_region_render.py", "_EXPOSURE_MODALITY_TAGS")
    assert set(ocs.KIND_LABELS) == set(tags)


# --------------------------------------------------------------------------- live
@pytest.mark.live
@pytest.mark.skipif(os.environ.get("RUN_LIVE_TESTS") != "1", reason="set RUN_LIVE_TESTS=1")
def test_live_open_targets_and_gwas_catalog(capsys):
    """Open Targets GraphQL + GWAS Catalog FTP: the IL23R x IBD page still resolves."""
    a = ocs.analyse(IBD)
    assert a.selected is not None and a.selected.qtd == "QTD000584"
    assert a.selected.target_id == "ENSG00000162594"
    assert a.page_probe.status == "indexed"
    assert a.credible_set["variant"]["id"] == "1_67240275_G_A"


# --------------------------------------------------------------------------- review rework
@pytest.mark.parametrize("body", ["", "<html><body>Too many requests</body></html>"])
def test_a_200_without_a_listing_is_unknown_not_no(monkeypatch, body):
    monkeypatch.setattr(ocs, "_get", lambda url, timeout: (200, body))
    assert ocs.probe_gwas("GCST1").status == "unknown"


def test_named_gene_reason_counts_only_that_gene(fake):
    a = ocs.analyse(CHOL, gene="BSND", probe=False)
    assert a.qtl_rec.row.symbol == "BSND"
    reason = [r for r in a.qtl_rec.reasons if r.startswith("highest")][0]
    assert "among the 2 QTLs" in reason and reason.endswith("for BSND")


def test_substitute_with_an_unchecked_file_says_so(monkeypatch, tmp_path):
    f = FakeOT(probe_overrides={"GCST90475721": TimeoutError("slow")})
    monkeypatch.setattr(ocs, "_post", f.post)
    monkeypatch.setattr(ocs, "_get", f.get)
    out = tmp_path / "c.yaml"
    assert ocs.main(["config", CHOL, "--qtl", PCSK9_LUNG_GE, "--gwas", "GCST90475721",
                     "--out", str(out)]) == 0
    got = yaml.safe_load(out.read_text())
    assert got["provenance"]["outcome_study_file_check"] == "unknown"
    assert any("Whether GCST90475721 can be read by region is unknown" in c for c in got["caveats"])


def test_exon_qtl_config_states_the_gene_pooling(fake, tmp_path):
    out = tmp_path / "c.yaml"
    assert ocs.main(["config", CHOL, "--qtl", "gtex-v10_exon_lung_ensg00000169174_11_1_55063529_55064797",
                     "--gwas", "GCST90475721", "--out", str(out)]) == 0
    got = yaml.safe_load(out.read_text())
    assert any("pools every exon eQTL feature of PCSK9" in c for c in got["caveats"])


def _mutated_credible_set(slid: str, **changes) -> dict:
    resp = copy.deepcopy(_fixture("credibleSet", slid))
    cs = resp["data"]["credibleSet"]
    for k, v in changes.items():
        if k == "study_id":
            cs["study"]["id"] = v
        else:
            cs[k] = v
    return resp


def test_a_qtl_credible_set_is_refused(monkeypatch):
    f = FakeOT()
    resp = _mutated_credible_set(IBD, studyType="pqtl")
    monkeypatch.setattr(ocs, "_post", lambda q, v: resp if v.get("types") is None else f.post(q, v))
    assert ocs.main(["options", IBD]) == 2


def test_a_non_gwas_catalog_page_study_is_refused(monkeypatch, tmp_path):
    f = FakeOT()
    resp = _mutated_credible_set(IBD, study_id="FINNGEN_R12_K11_IBD_STRICT")
    monkeypatch.setattr(ocs, "_post", lambda q, v: resp if v.get("types") is None else f.post(q, v))
    monkeypatch.setattr(ocs, "_get", f.get)
    out = tmp_path / "c.yaml"
    # not_gcst is a definite no: substitutes are offered, and keeping the page study is refused.
    assert ocs.main(["config", IBD, "--qtl", "1", "--out", str(out)]) == 2
    assert not out.exists()


def test_an_unreadable_qtl_is_refused(monkeypatch, tmp_path):
    f = FakeOT()
    resp = copy.deepcopy(_fixture("coloc_molqtl", IBD))
    resp["data"]["credibleSet"]["colocalisation"]["rows"][0]["otherStudyLocus"]["study"][
        "summarystatsLocation"] = "https://example.org/not-a-catalogue-path"
    f.coloc_overrides[("coloc_molqtl", IBD)] = resp
    monkeypatch.setattr(ocs, "_post", f.post)
    monkeypatch.setattr(ocs, "_get", f.get)
    out = tmp_path / "c.yaml"
    assert ocs.main(["config", IBD, "--qtl", "1", "--out", str(out)]) == 2
    a = ocs.analyse(IBD, qtl="1")
    assert a.page_probe is None  # nothing on the GWAS side is checked for an unreadable QTL
    assert ocs.main(["options", IBD, "--qtl", "1"]) == 0  # prints, does not crash


def test_index_override_comes_first(monkeypatch, tmp_path):
    alt = tmp_path / "idx.tsv"
    alt.write_text("dataset_id\tstudy_label\tquant_method\tsample_size\nQTD999999\tX\tge\t1\n")
    monkeypatch.setenv("EQTL_CATALOGUE_INDEX", str(alt))
    assert set(ocs.load_eqtl_index()) == {"QTD999999"}


def test_a_description_selects_or_lists(fake, tmp_path):
    a = ocs.analyse(CHOL, qtl="PCSK9 lung gene", probe=False)
    assert a.selected is not None and a.selected.study_id == "gtex-v10_ge_lung_ensg00000169174"
    a = ocs.analyse(CHOL, qtl="PCSK9 lung", probe=False)
    assert a.selected is None
    assert {r.study_id for r in a.matches} == {
        "gtex-v10_ge_lung_ensg00000169174", "gtex-v10_exon_lung_ensg00000169174_11_1_55063529_55064797"}
    assert a.match_rec.row.study_id == "gtex-v10_ge_lung_ensg00000169174"
    out = str(tmp_path / "c.yaml")
    assert ocs.main(["config", CHOL, "--qtl", "PCSK9 lung", "--gwas", "GCST90475721", "--out", out]) == 2
    assert ocs.main(["options", CHOL, "--qtl", "liver"]) == 2  # nothing matches: said, not guessed


def test_lead_override_must_be_in_the_page_credible_set(fake, tmp_path):
    out = tmp_path / "c.yaml"
    base = ["config", CHOL, "--qtl", PCSK9_LUNG_GE, "--gwas", "GCST90475721", "--out", str(out)]
    assert ocs.main(base + ["--lead", "1_55000000_A_G"]) == 2
    assert ocs.main(base + ["--lead", "1_55023869_C_T"]) == 0
    got = yaml.safe_load(out.read_text())
    assert got["lead"]["variant_id"] == "1_55023869_C_T"
    assert "rs_id" not in got["lead"]  # the page's rs id names the page's lead
    assert {s["side"] for s in got["provenance"]["substitutions"]} >= {"lead"}
    assert len(got["lead"]["candidates"]) == 5  # --lead sets the configured lead; fallback stays


def test_a_truncated_credible_set_says_so(fake, monkeypatch, tmp_path):
    resp = copy.deepcopy(_fixture("credibleSet", CHOL))
    resp["data"]["credibleSet"]["locus"]["count"] = 700
    real_post = fake.post
    monkeypatch.setattr(ocs, "_post", lambda q, v: resp if v.get("types") is None else real_post(q, v))
    out = tmp_path / "c.yaml"
    assert ocs.main(["config", CHOL, "--qtl", PCSK9_LUNG_GE, "--gwas", "GCST90475721",
                     "--out", str(out)]) == 0
    got = yaml.safe_load(out.read_text())
    assert any("lists 700 members" in c and "5 it returned" in c for c in got["caveats"])


def test_an_empty_open_targets_answer_is_not_an_absence(monkeypatch, capsys):
    monkeypatch.setattr(ocs, "_post", lambda q, v: {"data": None})
    assert ocs.main(["options", CHOL]) == 3
    assert "no credible set" not in capsys.readouterr().err


def test_candidates_are_written_in_pip_order_whatever_order_open_targets_returns(
        fake, monkeypatch, tmp_path):
    resp = copy.deepcopy(_fixture("credibleSet", CHOL))
    resp["data"]["credibleSet"]["locus"]["rows"].reverse()  # lowest PIP first
    real_post = fake.post
    monkeypatch.setattr(ocs, "_post", lambda q, v: resp if v.get("types") is None else real_post(q, v))
    out = tmp_path / "c.yaml"
    assert ocs.main(["config", CHOL, "--qtl", PCSK9_LUNG_GE, "--gwas", "GCST90475721",
                     "--out", str(out)]) == 0
    pips = [c["pip"] for c in yaml.safe_load(out.read_text())["lead"]["candidates"]]
    assert pips == sorted(pips, reverse=True) and len(pips) == 5


def test_lead_override_is_worded_as_requested_not_final(fake, tmp_path):
    out = tmp_path / "c.yaml"
    assert ocs.main(["config", CHOL, "--qtl", PCSK9_LUNG_GE, "--gwas", "GCST90475721",
                     "--lead", "1_55023869_C_T", "--out", str(out)]) == 0
    caveats = yaml.safe_load(out.read_text())["caveats"]
    assert any(c.startswith("The lead requested is 1_55023869_C_T") for c in caveats)
    assert not any("centred on" in c for c in caveats)


def test_an_unwritable_output_is_not_called_an_open_targets_failure(fake, tmp_path, capsys):
    out = tmp_path / "missing_dir" / "c.yaml"
    assert ocs.main(["config", CHOL, "--qtl", PCSK9_LUNG_GE, "--gwas", "GCST90475721",
                     "--out", str(out)]) == 4
    assert "could not read Open Targets" not in capsys.readouterr().err


def test_a_small_float_is_written_as_a_yaml_float():
    assert yaml.safe_load(f"x: {ocs._y(1e-05)}")["x"] == pytest.approx(1e-05)
    assert yaml.safe_load(f"x: {ocs._y(0.25)}")["x"] == 0.25


def test_a_study_colocalising_twice_records_its_best_h4(fake, tmp_path):
    resp = copy.deepcopy(_fixture("coloc_gwas", PCSK9_LUNG_GE))
    rows = resp["data"]["credibleSet"]["colocalisation"]["rows"]
    best = next(r for r in rows if r["otherStudyLocus"]["study"]["id"] == "GCST90475721")
    worse = copy.deepcopy(best)
    worse["h4"] = 0.5
    worse["otherStudyLocus"]["studyLocusId"] = "0" * 32
    rows.append(worse)  # listed after the better row
    resp["data"]["credibleSet"]["colocalisation"]["count"] += 1
    fake.coloc_overrides[("coloc_gwas", PCSK9_LUNG_GE)] = resp
    out = tmp_path / "c.yaml"
    assert ocs.main(["config", CHOL, "--qtl", PCSK9_LUNG_GE, "--gwas", "GCST90475721",
                     "--out", str(out)]) == 0
    got = yaml.safe_load(out.read_text())
    assert got["provenance"]["coloc_qtl_vs_outcome_h4"] == pytest.approx(best["h4"])
