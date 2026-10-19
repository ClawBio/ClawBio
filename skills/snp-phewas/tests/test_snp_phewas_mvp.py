"""Million Veteran Program (GWAS Catalog) offline tests: the catalog build and its write-time
assertions, the reader's three format facts (CRLF, odds ratio -> beta / SE, p underflow), the
harmonizer's not-orientated records, the autosome rule, unreadable studies, the crosswalk's
MVP rows, and the Manhattan category of an MVP outcome. Demo +/-10 kb slices of GCST90475929
(ischemic heart disease, phecode 411, odds ratios) and GCST90475416 (LDL cholesterol mean,
beta) around rs12740374."""

from __future__ import annotations

import copy
import gzip
import hashlib
import math
from statistics import NormalDist
import os

import pytest

import snp_phewas_mvp_catalog as M
import snp_phewas_schema as S
from snp_phewas_catalog import OutcomeSpec, get_dataset, load_outcome_catalog
from snp_phewas_tabix import read_header
from snp_phewas_lookup import lookup_variant
from snp_phewas_sources import (
    MVP_COLUMNS_BINARY,
    MVP_COLUMNS_QUANT,
    P_UNDERFLOW_DETAIL,
    Z_95,
    MvpReader,
    neglog10_p_from_z,
    record_from_fields,
)
from snp_phewas_demo import MVP_DIR, MVP_OUTCOMES, demo_map as fake_map, demo_rsid as fake_rsid
from snp_phewas_variant import Placement

IHD, LDL = MVP_OUTCOMES
IHD_FILE = os.path.join(MVP_DIR, "GCST90475001-GCST90476000/GCST90475929/harmonised/GCST90475929.h.tsv.gz")
LDL_FILE = os.path.join(MVP_DIR, "GCST90475001-GCST90476000/GCST90475416/harmonised/GCST90475416.h.tsv.gz")
TARGET = Placement("GRCh38", "1", 109274968, "G", "T")


def _reader():
    return MvpReader(get_dataset("mvp_eur"), data_base=MVP_DIR, workers=2)


def _spec(oid):
    return {o.outcome_id: o for o in load_outcome_catalog("mvp_eur")}[oid]


def _records(path):
    with gzip.open(path, "rt", newline="") as fh:
        lines = fh.read().split("\r\n")
    header = lines[0].split("\t")
    return header, [dict(zip(header, ln.split("\t"))) for ln in lines[1:] if ln]


# --- the catalog ---------------------------------------------------------------------------

def test_the_committed_catalog_rebuilds_byte_for_byte_from_the_pinned_snapshots(tmp_path):
    from snp_phewas_build_catalogs import build_mvp

    rows = build_mvp()  # checks each snapshot against its md5 pin
    out = tmp_path / "mvp.tsv.gz"
    assert M.write_catalog(rows, str(out)) == 7537
    committed = os.path.join(os.path.dirname(get_dataset("mvp_eur").outcome_catalog_path),
                             "mvp_outcomes.tsv.gz")
    assert hashlib.md5(out.read_bytes()).hexdigest() == hashlib.md5(
        open(committed, "rb").read()).hexdigest()


def _snapshot():
    from snp_phewas_build_catalogs import _SKILL_DIR, _load_yaml

    pin = _load_yaml()["families"]["mvp"]["release_pin"]

    def rd(k):
        return open(os.path.join(_SKILL_DIR, pin[k]), "rb").read()

    return (M.read_studies(rd("study_list_snapshot")),
            M.read_harmonised_paths(rd("harmonised_subset_snapshot")),
            M.read_meta(rd("harmonised_meta_snapshot")))


def test_each_write_time_assertion_fires_on_a_constructed_failure():
    studies, harm, meta = _snapshot()
    rows = M.build_rows(studies, harm, meta)
    ids = {r["outcome_id"] for r in rows}
    with pytest.raises(M.MvpCatalogError, match="repeats or omits"):
        M.build_rows(studies + studies[:1], harm, meta)
    with pytest.raises(M.MvpCatalogError, match="not in the list"):
        M.build_rows(studies[1:], harm, meta) if studies[0]["accession_id"] in harm else \
            M.build_rows([s for s in studies if s["accession_id"] != next(iter(harm))], harm, meta)
    acc = next(iter(harm))
    bad_meta = copy.deepcopy(meta)
    bad_meta[acc]["data_file_md5sum"] = ""
    with pytest.raises(M.MvpCatalogError, match="no metadata md5"):
        M.build_rows(studies, harm, bad_meta)
    bad_meta = copy.deepcopy(meta)
    bad_meta[acc]["genome_assembly"] = "GRCh37"
    with pytest.raises(M.MvpCatalogError, match="not GRCh38"):
        M.build_rows(studies, harm, bad_meta)
    # a filename on a study that is not harmonized, and a lost row
    raw = next(r for r in rows if r["access_status"] != M.ACCESS_HARMONISED)
    with pytest.raises(M.MvpCatalogError, match="set exactly when harmonised"):
        M.assert_catalog([{**r, "filename": "x"} if r is raw else r for r in rows], ids)
    with pytest.raises(M.MvpCatalogError, match="do not partition"):
        M.assert_catalog(rows[1:], ids)
    odd = copy.deepcopy(studies)
    odd[0]["discovery_ancestry"] = ["100 Martian (U.S.)"]
    with pytest.raises(M.MvpCatalogError, match="unknown MVP discovery ancestry"):
        M.build_rows(odd, harm, meta)


def test_trait_population_and_sample_size_parsing():
    t = M.parse_trait("Ischemic heart disease (PheCode 411)")
    assert (t["trait_kind"], t["phecode"], t["label"]) == ("phecode", "411", "Ischemic heart disease")
    assert M.parse_trait("Bacterial enteritis (PheCode 008.5)")["phecode"] == "8.5"
    q = M.parse_trait("low density lipoprotein cholesterol (LDLC, mean, inv-norm transformed)")
    assert (q["trait_kind"], q["quant_base"], q["quant_qualifier"], q["quant_stat"]) == (
        "quantitative", "low density lipoprotein cholesterol", "LDLC", "mean")
    assert M.parse_trait("high density lipoprotein cholesterol (HDLC, minimm, inv-norm "
                         "transformed)")["quant_stat"] == "minimum"
    assert M.parse_trait("Age at menopause (inv-normal transformed)")["quant_stat"] == "single"
    assert M.parse_trait("Biological Father: Asthma")["trait_kind"] == ""  # decided by sample size
    assert M.parse_population(["417274 European (U.S.)"]) == ("EUR", "")
    assert M.parse_population(["1 East Asian (U.S.)", "2 European (U.S.)"]) == ("META", "EUR+EAS")
    assert M.parse_sample_size("11,954 Hispanic or Latin American cases, 43,516 Hispanic or "
                               "Latin American controls") == (11954, 43516, 55470)
    assert M.parse_sample_size("404,741 European ancestry individuals") == (None, None, 404741)
    with pytest.raises(M.MvpCatalogError):
        M.parse_sample_size("about 400k people")


def test_catalog_counts_on_the_snapshot():
    rows = list(load_outcome_catalog("mvp_eur"))
    kinds = {}
    for o in rows:
        kinds[(o.trait_type, o.access_status == "harmonised")] = kinds.get(
            (o.trait_type, o.access_status == "harmonised"), 0) + 1
    # counted on the 2026-10-19 GWAS Catalog snapshot
    assert kinds == {("phecode", True): 930, ("phecode", False): 686,
                     ("quantitative", True): 107, ("quantitative", False): 107,
                     ("survey_binary", True): 184, ("survey_binary", False): 56}
    assert all(o.outcome_id.startswith("GCST") for o in rows)
    labels = [o.outcome_label for o in rows]
    assert len(set(labels)) == len(labels)


# --- the reader ----------------------------------------------------------------------------

def test_crlf_lines_are_stripped_from_the_header_and_the_records():
    with gzip.open(IHD_FILE, "rb") as fh:
        assert fh.readline().endswith(b"variant_id\r\n")  # the slice keeps the source's CRLF
    assert read_header(IHD_FILE)[-1] == "variant_id"
    rec = _reader().read_outcome(_spec(IHD), TARGET)
    assert rec.status == S.ST_OK and rec.variant_id_native == "1_109274968_G_T"
    assert "\r" not in rec.variant_id_native


def test_odds_ratio_becomes_beta_and_an_se_that_reproduces_the_files_p():
    header, recs = _records(IHD_FILE)
    assert "odds_ratio" in header and "beta" not in header
    n = n_interval_off = 0
    for r in recs:
        tgt = Placement("GRCh38", "1", int(r["base_pair_location"]), r["other_allele"],
                        r["effect_allele"])
        raw = record_from_fields([r], MVP_COLUMNS_BINARY, tgt, "tabix")
        assert raw.status == S.ST_OK, r
        p_file = float(r["p_value"])
        p_from_se = math.erfc(abs(raw.beta / raw.se) / math.sqrt(2))
        assert abs(p_from_se - p_file) < 1e-3 and \
            abs(math.log10(p_from_se) - math.log10(p_file)) < 1e-6, (r["rsid"], p_from_se, p_file)
        assert raw.neglog10_p == pytest.approx(-math.log10(p_file))
        assert raw.beta == pytest.approx(math.log(float(r["odds_ratio"])))
        se_ci = (math.log(float(r["ci_upper"])) - math.log(float(r["ci_lower"]))) / (2 * Z_95)
        p_ci = math.erfc(abs(raw.beta / se_ci) / math.sqrt(2))
        n_interval_off += (abs(p_ci - p_file) > 1e-3
                           or abs(math.log10(max(p_ci, 1e-300)) - math.log10(p_file)) > 0.05)
        n += 1
    # every record of the slice; the rounded interval alone would have missed the file's p
    # on 85 of them (why the SE is taken from p where p can give it)
    assert (n, n_interval_off) == (150, 85)
    # p cannot give it (p underflowed to 0): the interval's SE, and -log10 p from it
    r = dict(recs[0], p_value="0.0", odds_ratio="2.0", ci_lower="1.8", ci_upper="2.2222")
    tgt = Placement("GRCh38", "1", int(r["base_pair_location"]), r["other_allele"], r["effect_allele"])
    raw = record_from_fields([r], MVP_COLUMNS_BINARY, tgt, "tabix")
    assert raw.se == pytest.approx((math.log(2.2222) - math.log(1.8)) / (2 * Z_95))
    assert raw.detail == P_UNDERFLOW_DETAIL
    # an odds ratio and a usable p with NO interval is still a value: the interval is only
    # the fallback (requiring it made such a record "not tested")
    r = dict(recs[0], p_value="0.01", odds_ratio="1.2", ci_lower="NA", ci_upper="NA")
    raw = record_from_fields([r], MVP_COLUMNS_BINARY, tgt, "tabix")
    assert raw.status == S.ST_OK
    assert raw.se == pytest.approx(math.log(1.2) / -NormalDist().inv_cdf(0.005))
    # and with neither a usable p nor an interval there is no SE: unparsable, never a guess
    r = dict(recs[0], p_value="0.0", odds_ratio="1.2", ci_lower="NA", ci_upper="NA")
    assert record_from_fields([r], MVP_COLUMNS_BINARY, tgt, "tabix").status == \
        S.ST_RECORD_UNPARSABLE


def test_p_underflow_is_derived_from_beta_over_se_and_said():
    raw = _reader().read_outcome(_spec(LDL), TARGET)
    assert raw.status == S.ST_OK and (raw.beta, raw.se) == (-0.1561, 0.003009)
    assert raw.detail == P_UNDERFLOW_DETAIL
    assert raw.neglog10_p == pytest.approx(neglog10_p_from_z(-0.1561 / 0.003009))
    assert raw.neglog10_p > 300  # p underflows a double; -log10 p does not
    scipy = pytest.importorskip("scipy.stats")
    ref = -(scipy.norm.logsf(0.1561 / 0.003009) + math.log(2)) / math.log(10)
    assert raw.neglog10_p == pytest.approx(ref, abs=1e-9)
    # no usable SE with p = 0 is unparsable, never a made-up p
    r = {"chromosome": "1", "base_pair_location": "5", "other_allele": "A", "effect_allele": "G",
         "beta": "0.1", "standard_error": "0", "p_value": "0.0", "effect_allele_frequency": "0.1",
         "hm_code": "10"}
    assert record_from_fields([r], MVP_COLUMNS_QUANT, Placement("GRCh38", "1", 5, "A", "G"),
                              "tabix").status == S.ST_RECORD_UNPARSABLE


def test_neglog10_p_from_z_matches_scipy_across_the_range():
    scipy = pytest.importorskip("scipy.stats")
    for z in (0.1, 1, 3, 10, 35, 35.4, 36, 40, 60, 100, 300):
        ref = -(scipy.norm.logsf(z) + math.log(2)) / math.log(10)
        assert abs(neglog10_p_from_z(z) - ref) < 2e-11, z
        assert neglog10_p_from_z(-z) == neglog10_p_from_z(z)


@pytest.mark.parametrize("code,withheld", [("9", True), ("14", True), ("18", True),
                                           ("10", False), ("11", False), ("6", False)])
def test_a_record_the_harmonizer_could_not_orient_keeps_its_p_and_withholds_its_direction(code, withheld):
    r = {"chromosome": "1", "base_pair_location": "5", "other_allele": "A", "effect_allele": "G",
         "beta": "0.2", "standard_error": "0.01", "p_value": "1e-20",
         "effect_allele_frequency": "0.3", "hm_code": code}
    raw = record_from_fields([r], MVP_COLUMNS_QUANT, Placement("GRCh38", "1", 5, "A", "G"), "tabix")
    assert raw.status == S.ST_OK and raw.neglog10_p == pytest.approx(20)
    assert raw.direction_withheld == (S.WITHHELD_SOURCE_NOT_ORIENTATED if withheld else None)


def test_lookup_on_the_slices_orients_withholds_and_accounts_for_every_mvp_study(monkeypatch):
    import snp_phewas_sources as src

    r = lookup_variant("rs12740374", datasets="mvp_eur", readers={"mvp_eur": _reader()},
                       outcomes={"mvp_eur": MVP_OUTCOMES}, rsid_resolver=fake_rsid,
                       map_fn=fake_map)
    rows = r.dataset_rows("mvp_eur")
    assert len(rows) == 2070  # the partition: every study of the population is one row
    by = {x.outcome_id: x for x in rows}
    assert by[IHD].effect_allele == "T" and by[IHD].beta == pytest.approx(math.log(0.9107468123861566))
    assert by[IHD].harmonized_outcome_id == "PHECODE:411"
    assert by[LDL].status_detail == P_UNDERFLOW_DETAIL
    counts = r.status_counts("mvp_eur")
    assert counts[S.ST_NOT_READABLE] == 849 and counts[S.ST_OK] == 2
    assert counts[S.ST_NOT_QUERIED] == 2070 - 849 - 2
    # not readable is outside the significance denominator, like not looked up
    assert S.phenome_wide_n(rows) == 2
    assert any("1,221 of 2,070 outcomes are readable" in n and "448 not harmonized" in n
               and "401 deposited without summary statistics" in n for n in r.notes)
    # a hm_code 9 record at the target: same p, no direction, no beta
    real = src.record_from_fields

    def not_orientated(records, cols, tgt, route):
        rec = real(records, cols, tgt, route)
        rec.direction_withheld = S.WITHHELD_SOURCE_NOT_ORIENTATED
        return rec

    monkeypatch.setattr(src, "record_from_fields", not_orientated)
    w = lookup_variant("rs12740374", datasets="mvp_eur", readers={"mvp_eur": _reader()},
                       outcomes={"mvp_eur": [IHD]}, rsid_resolver=fake_rsid, map_fn=fake_map)
    ihd = {x.outcome_id: x for x in w.rows}[IHD]
    assert ihd.beta is None and ihd.eaf is None and ihd.neglog10_p == by[IHD].neglog10_p
    assert ihd.direction_withheld_reason == S.WITHHELD_SOURCE_NOT_ORIENTATED
    assert any("could not orient" in n for n in w.notes)


def test_the_palindromic_variant_in_the_mvp_slices_keeps_its_direction_in_order():
    raw = _reader().read_outcome(_spec(LDL), Placement("GRCh38", "1", 109270888, "A", "T"))
    assert raw.status == S.ST_OK and (raw.src_ref, raw.src_alt) == ("A", "T")
    r = lookup_variant("1_109270888_A_T", datasets="mvp_eur", readers={"mvp_eur": _reader()},
                       outcomes={"mvp_eur": [LDL]}, map_fn=fake_map)
    row = {x.outcome_id: x for x in r.rows}[LDL]
    assert row.palindromic and row.direction_withheld_reason is None and row.beta == -0.05772


def test_an_x_variant_is_absent_by_rule_with_no_file_read(monkeypatch):
    import snp_phewas_sources as src

    def boom(*a, **k):
        raise AssertionError("no file may be read for an X variant in an autosomes-only dataset")

    monkeypatch.setattr(src, "read_outcome_file", boom)
    r = lookup_variant("X_1000000_A_G", datasets="mvp_eur", readers={"mvp_eur": _reader()},
                       outcomes={"mvp_eur": MVP_OUTCOMES}, map_fn=fake_map)
    by = {x.outcome_id: x for x in r.rows}
    assert by[IHD].status == S.ST_ABSENT and by[IHD].status_detail == "released autosomal results only"
    assert r.provenance["mvp_eur"]["autosomes_only"]
    assert any("autosomal results only" in n for n in r.notes)


def test_an_unreadable_study_is_answered_from_the_catalog_with_no_network(monkeypatch):
    import snp_phewas_sources as src

    monkeypatch.setattr(src, "read_outcome_file",
                        lambda *a, **k: (_ for _ in ()).throw(AssertionError("read attempted")))
    o = next(o for o in load_outcome_catalog("mvp_eur") if o.access_status == "raw_only_unindexed")
    rec = _reader().read_outcome(o, TARGET)
    assert rec.status == S.ST_NOT_READABLE and rec.detail.startswith("raw_only_unindexed:")
    none = OutcomeSpec("mvp_eur", "GCSTX", "x", "", "binary", None, 1, 1, None,
                       access_status="no_summary_statistics")
    assert _reader().read_outcome(none, TARGET).detail.startswith("no_summary_statistics:")


def test_an_mvp_phecode_takes_the_pan_ukbb_category_of_the_same_code():
    from snp_phewas_plots import outcome_category, phecode_category

    def row(native, kind="binary"):
        return S.OutcomeRow(
            dataset_id="mvp_eur", dataset_release="r", outcome_id="GCST1", outcome_label="x",
            outcome_category_native=native, outcome_kind=kind, n_cases=1, n_controls=1,
            n_total=2, variant_id_grch38="1_1_A_G", variant_id_native=None, rs_id=None,
            effect_allele="G", other_allele="A", status=S.ST_OK)

    assert outcome_category(row("PheCode 411")) == "circulatory"  # Pan-UKBB: circulatory system
    assert outcome_category(row("PheCode 8.5")) == "infectious"  # Pan-UKBB writes 008.5
    assert outcome_category(row("PheCode 250.2")) == "endocrine & metabolic"
    # a code Pan-UKBB lacks: its numeric neighbors' shared category (both circulatory here)
    assert phecode_category("411.9") == "circulatory"
    assert outcome_category(row(M.TRAIT_KIND_LABELS[M.KIND_QUANT], "continuous")) == "measurements"
    assert outcome_category(row(M.TRAIT_KIND_LABELS[M.KIND_SURVEY])) == "other"


# --- the crosswalk's MVP rows ---------------------------------------------------------------

def test_the_crosswalk_carries_every_mvp_study_by_the_phecode_and_measurement_rules():
    """One row per study (the 7,537 of the catalog), phecode studies keyed on the same code
    as Pan-UKBB's both-sexes phecode outcome, and every MVP measurement classified exactly as
    the vendored reviewed table says (mean / single compared, maximum / minimum listed only)."""
    import csv
    import gzip as gz

    from snp_phewas_crosswalk import (
        DEFAULT_CROSSWALK_PATH,
        SOURCES,
        load_crosswalk,
        normalize_phecode,
        quant_key,
    )

    cw = load_crosswalk()
    skill = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    with gz.open(os.path.join(skill, "catalogs", "mvp_outcomes.tsv.gz"), "rt", newline="") as fh:
        mvp = list(csv.DictReader(fh, delimiter="\t"))
    assert {o for f, o in cw if f == "mvp"} == {r["outcome_id"] for r in mvp}
    assert sum(1 for f, _ in cw if f == "pan_ukbb") == 7223
    assert sum(1 for f, _ in cw if f == "finngen") == 2469
    for name in ("mvp_outcome_catalog", "mvp_measurement_matches"):
        path, md5 = SOURCES[name]
        assert hashlib.md5(open(os.path.join(skill, path), "rb").read()).hexdigest() == md5, name
    lines = [ln for ln in open(os.path.join(skill, "catalogs", "mvp_measurement_matches.tsv"),
                               encoding="utf-8").read().splitlines() if ln and not ln.startswith("#")]
    table = {r["mvp_measurement"]: r for r in csv.DictReader(lines, delimiter="\t")}
    assert len(table) == 76
    pan_phecode = {normalize_phecode(o.split("-")[1]): o for (f, o) in cw
                   if f == "pan_ukbb" and o.startswith("phecode-") and o.endswith("-both_sexes")}
    for r in mvp:
        x = cw[("mvp", r["outcome_id"])]
        if r["trait_kind"] == "phecode":
            if r["phecode"] in pan_phecode:
                assert (x.harmonized_outcome_id, x.harmonization_basis) == (
                    f"PHECODE:{r['phecode']}", "phecode_same_code"), r["outcome_id"]
            else:
                assert x.unmapped_reason == "no_counterpart_code", r["outcome_id"]
        elif r["trait_kind"] == "quantitative":
            t = table[quant_key(r["quant_base"], r["quant_qualifier"])]
            want = {"exact": ("summary_statistic_variant" if r["quant_stat"] in ("maximum", "minimum")
                              else "curated_name_match"),
                    "near_equivalent": "measurement_near_equivalent", "none": "unmapped"}[t["match"]]
            assert x.harmonization_basis == want, r["outcome_id"]
            if t["match"] == "exact":
                assert x.anchor_outcome_id == t["pan_ukbb_outcome_ids"], r["outcome_id"]
        else:
            assert x.unmapped_reason == "survey_item_no_counterpart"
    assert DEFAULT_CROSSWALK_PATH.endswith("catalogs/biobank_outcome_crosswalk.tsv.gz")


def test_a_maximum_or_minimum_keeps_its_id_and_is_never_the_representative(offline_kwargs):
    """MVP's LDL maximum and minimum share the mean's harmonized id; only the mean may stand
    for MVP in the replication view."""
    from snp_phewas_crosswalk import load_crosswalk

    cw = load_crosswalk()
    ldl_id = cw[("mvp", LDL)].harmonized_outcome_id
    variants = [o for o in load_outcome_catalog("mvp_eur")
                if cw[("mvp", o.outcome_id)].harmonized_outcome_id == ldl_id
                and o.outcome_id != LDL]
    assert variants and all(cw[("mvp", o.outcome_id)].harmonization_basis
                            == "summary_statistic_variant" for o in variants)
    r = lookup_variant("rs12740374", **offline_kwargs)
    rows = {x.outcome_id: x for x in r.dataset_rows("mvp_eur")}
    assert rows[LDL].harmonized_representative is True
    assert all(rows[o.outcome_id].harmonized_representative is False for o in variants)


def test_a_tested_maximum_does_not_stand_in_for_an_unread_mean():
    """The maximum / minimum keep the mean's harmonized id but are never compared: when the
    mean could not be read, the maximum must not become MVP's representative."""
    from types import SimpleNamespace

    from snp_phewas_crosswalk import load_crosswalk
    from snp_phewas_harmonize import attach_harmonization

    cw = load_crosswalk()
    hid = cw[("mvp", LDL)].harmonized_outcome_id
    mx = next(o for o in load_outcome_catalog("mvp_eur")
              if o.outcome_id != LDL and cw[("mvp", o.outcome_id)].harmonized_outcome_id == hid
              and "maximum" in o.outcome_label)

    def row(oid, status, n):
        return S.OutcomeRow(
            dataset_id="mvp_eur", dataset_release="r", outcome_id=oid, outcome_label=oid,
            outcome_category_native="", outcome_kind="continuous", n_cases=n, n_controls=None,
            n_total=n, variant_id_grch38="1_1_A_G", variant_id_native=None, rs_id=None,
            effect_allele="G", other_allele="A", status=status)

    mean, maximum = row(LDL, S.ST_UNREACHABLE, 1), row(mx.outcome_id, S.ST_OK, 10**6)
    attach_harmonization(SimpleNamespace(dataset_ids=["mvp_eur"], rows=[mean, maximum],
                                         uncatalogued=[]))
    assert maximum.harmonized_outcome_id == hid and maximum.harmonized_representative is False
    assert mean.harmonized_representative is True



class _Fixed:
    """A reader that answers every outcome with one status (ok carries a weak value)."""

    def __init__(self, status):
        self.status = status

    def read(self, outcomes, tgt):
        from snp_phewas_sources import RawRecord

        rec = (RawRecord(S.ST_OK, src_ref=tgt.ref, src_alt=tgt.alt, beta=0.01, se=0.01,
                         neglog10_p=0.5, eaf=0.2, variant_id_native=tgt.variant_id)
               if self.status == S.ST_OK else RawRecord(S.ST_UNREACHABLE, detail="HTTP 503"))
        return {o.outcome_id: rec for o in outcomes}, {"route": "tabix",
                                                       "n_files_read": len(outcomes),
                                                       "elapsed_s": 0.0}


def test_an_outage_never_changes_the_representative_or_the_threshold():
    """MVP AFR holds two studies for phecode 274.1: GCST90476378 (readable, 6,538 cases) and
    GCST90475734 (no readable file, 14,437 cases). With the readable one unreachable, the
    case-count tie-break picked the unreadable one, which the threshold never counts, so an
    outage cut the shared count (213 -> 212 in this construction) and loosened 0.05 / n. A
    study known unreadable from the catalog is never the representative over a readable one."""
    from snp_phewas_replication import build_replication

    pair = ("pan_ukbb_afr", "mvp_afr")
    shared, reps = {}, {}
    for status in (S.ST_OK, S.ST_UNREACHABLE):
        res = lookup_variant("1_109274968_G_T", datasets=list(pair), map_fn=fake_map,
                             readers={"pan_ukbb_afr": _Fixed(S.ST_OK), "mvp_afr": _Fixed(status)},
                             scope="harmonized")
        shared[status] = build_replication(res).n_shared[pair]
        reps[status] = {r.outcome_id for r in res.rows if r.dataset_id == "mvp_afr"
                        and r.harmonized_outcome_id == "PHECODE:274.1"
                        and r.harmonized_representative}
    assert reps[S.ST_OK] == reps[S.ST_UNREACHABLE] == {"GCST90476378"}
    assert shared[S.ST_OK] == shared[S.ST_UNREACHABLE]


def test_the_meta_fallback_sentence_names_the_group_without_nested_brackets():
    from snp_phewas_catalog import datasets_for_ancestry

    ids, said = datasets_for_ancestry("META")
    assert "mvp_meta" in ids and "finngen_r12" in ids
    finngen = [s for s in said if "FinnGen" in s]
    assert finngen == ["This lookup serves no multi-population meta-analysis for FinnGen, so "
                       "its default group is shown instead: FinnGen R12 (Finnish)."]
