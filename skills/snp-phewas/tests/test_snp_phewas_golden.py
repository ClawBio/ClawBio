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
