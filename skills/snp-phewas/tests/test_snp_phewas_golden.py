"""Golden parity on rs12740374 (offline replay of vendored source slices).

Asserts the KEY SET of tested rows (not only their values: a golden that checks values
of the keys it expects cannot see an extra or missing row), every value, the status
partition, and API-equals-tabix parity for the FinnGen endpoints.
"""

from __future__ import annotations

import os

import pytest
import yaml

from snp_phewas_lookup import lookup_variant
from snp_phewas_replication import NO_COMPARISON, build_replication
from snp_phewas_demo import demo_map, demo_rsid, demo_readers

GOLD = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fixtures", "golden", "rs12740374")


def _load(name):
    with open(os.path.join(GOLD, name)) as fh:
        return yaml.safe_load(fh)


def _run(route="hybrid"):
    q = _load("query.yaml")
    return lookup_variant(q["variant"], datasets=q["datasets"], outcomes=q["outcomes"],
                          readers=demo_readers(finngen_route=route),
                          rsid_resolver=demo_rsid, map_fn=demo_map)


def test_golden_rs12740374():
    exp = _load("expected.yaml")
    r = _run()
    assert r.variant.variant_id_grch38 == _load("query.yaml")["variant_id_grch38"]
    got = {(x.dataset_id, x.outcome_id): x for x in r.tested}
    want = {(ds, oid) for ds, d in exp["tested"].items() for oid in d}
    assert set(got) == want
    for ds, d in exp["tested"].items():
        for oid, v in d.items():
            row = got[(ds, oid)]
            for k in ("beta", "se", "neglog10_p", "eaf"):
                assert getattr(row, k) == pytest.approx(v[k], rel=1e-9), (ds, oid, k)
            assert row.variant_id_native == v["variant_id_native"]
            assert row.harmonized_outcome_id == v["harmonized_outcome_id"]
            assert row.effect_allele == "T" and not row.allele_flip_applied
    for ds, counts in exp["status_counts"].items():
        assert dict(r.status_counts(ds)) == counts
    assert [u.outcome_id for u in r.uncatalogued] == exp["uncatalogued"]
    t = build_replication(r)
    pair = ("pan_ukbb_eur", "finngen_r12")
    assert {row.harmonized_outcome_id: row.verdicts[pair] for row in t.rows
            if row.verdicts[pair] not in NO_COMPARISON} == exp["replication"]


def test_finngen_api_equals_tabix_on_the_golden():
    api = {x.outcome_id: x for x in _run("api").tested if x.dataset_id == "finngen_r12"}
    tab = {x.outcome_id: x for x in _run("tabix").tested if x.dataset_id == "finngen_r12"}
    assert set(api) == set(tab) == {"I9_CHD", "E4_LIPOPROT", "I9_IHD"}
    for oid in api:
        assert api[oid].source_route == "api" and tab[oid].source_route == "tabix"
        for k in ("beta", "se", "neglog10_p", "eaf"):
            assert getattr(api[oid], k) == pytest.approx(getattr(tab[oid], k), rel=1e-6), (oid, k)


GOLD3 = os.path.join(os.path.dirname(GOLD), "rs12740374_three_biobanks")


def test_golden_rs12740374_three_biobanks():
    """Ischemic heart disease is one harmonized outcome in all three biobanks (Pan-UKBB
    phecode 411, FinnGen I9_IHD by FinnGen's published match, MVP phecode 411 by the same
    code) and replicates in each pair and in all three. LDL cholesterol is a Pan-UKBB x MVP
    pair only. The MVP LDL p underflows to 0 in the file; -log10 p is derived from beta / SE."""
    def load(name):
        with open(os.path.join(GOLD3, name)) as fh:
            return yaml.safe_load(fh)

    q, exp = load("query.yaml"), load("expected.yaml")
    r = lookup_variant(q["variant"], datasets=q["datasets"], outcomes=q["outcomes"],
                       readers=demo_readers(), rsid_resolver=demo_rsid, map_fn=demo_map)
    got = {x.outcome_id: x for x in r.tested if x.dataset_id == "mvp_eur"}
    assert set(got) == set(exp["mvp_tested"])  # the key set, not only the values
    for oid, v in exp["mvp_tested"].items():
        row = got[oid]
        for k in ("beta", "se", "neglog10_p", "eaf"):
            assert getattr(row, k) == pytest.approx(v[k], rel=1e-9), (oid, k)
        for k in ("variant_id_native", "harmonized_outcome_id", "status_detail"):
            assert getattr(row, k) == v[k], (oid, k)
        assert row.effect_allele == "T" and not row.allele_flip_applied
    for ds, counts in exp["status_counts"].items():
        assert dict(r.status_counts(ds)) == counts
    t = build_replication(r)
    assert {f"{a}|{b}": n for (a, b), n in t.n_shared.items()} == exp["n_shared"]
    assert t.n_shared_3 == exp["n_shared_3"] and t.alpha_3 == pytest.approx(0.05 / exp["n_shared_3"])
    comparisons = {}
    for row in t.rows:
        v = {f"{a}|{b}": x for (a, b), x in row.verdicts.items()
             if x not in NO_COMPARISON and x != "no_counterpart"}
        if v:
            comparisons[row.harmonized_outcome_id] = v
    assert comparisons == exp["replication"]
    assert {row.harmonized_outcome_id: row.three_way for row in t.rows
            if row.three_way and row.three_way not in NO_COMPARISON} == exp["three_way"]
