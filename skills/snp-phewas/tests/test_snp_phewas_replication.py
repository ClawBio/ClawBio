"""Unit tests: one per verdict class, the n_shared threshold, not-tested vs tested-null,
harmonization attachment, and the two cross-dataset figures."""

from __future__ import annotations

import math

import pytest

import snp_phewas_schema as S
from snp_phewas_lookup import lookup_variant
from snp_phewas_replication import (
    V_BOTH_NULL,
    V_CONCORDANT_ONE,
    V_DIRECTION_UNKNOWN,
    V_DISCORDANT,
    V_LOW_CONFIDENCE,
    V_NOT_READ,
    V_NOT_TESTED,
    V_REPLICATES,
    V_SINGLE,
    VERDICT_ORDER,
    ReplicationRow,
    ReplicationTable,
    build_replication,
    cell_signal,
    verdict,
)

ALPHA = 0.001


def _cell(p=None, beta=None, status=S.ST_OK):
    return S.OutcomeRow(
        dataset_id="x", dataset_release="r", outcome_id="o", outcome_label="o",
        outcome_category_native="", outcome_kind="binary", n_cases=1, n_controls=1, n_total=2,
        variant_id_grch38="1_1_A_G", variant_id_native="1_1_A_G", rs_id=None,
        effect_allele="G", other_allele="A", status=status, beta=beta, se=0.01,
        neglog10_p=None if p is None else -math.log10(p))


def _v(a, b):
    return verdict(cell_signal(a, ALPHA), cell_signal(b, ALPHA))


def test_each_verdict_class():
    sig_neg, sig_pos = _cell(1e-6, -0.1), _cell(1e-6, 0.1)
    null_neg, null_pos = _cell(0.3, -0.01), _cell(0.3, 0.01)
    assert _v(sig_neg, sig_neg) == V_REPLICATES
    assert _v(sig_neg, sig_pos) == V_DISCORDANT
    assert _v(sig_neg, null_neg) == V_CONCORDANT_ONE
    assert _v(null_pos, sig_pos) == V_CONCORDANT_ONE
    assert _v(sig_neg, null_pos) == V_SINGLE
    assert _v(null_neg, null_pos) == V_BOTH_NULL
    withheld = _cell(1e-6, None)
    assert _v(sig_neg, withheld) == V_DIRECTION_UNKNOWN
    assert _v(withheld, null_neg) == V_SINGLE
    for missing in (None, _cell(status=S.ST_ABSENT), _cell(status=S.ST_ALLELE_MISMATCH)):
        assert _v(sig_neg, missing) == V_NOT_TESTED
    # A read that failed or did not happen is never "not tested", and never a null.
    for unread in (S.ST_UNREACHABLE, S.ST_NOT_QUERIED, S.ST_FILE_MISSING, S.ST_RECORD_UNPARSABLE):
        assert _v(sig_neg, _cell(status=unread)) == V_NOT_READ, unread
        assert _v(_cell(status=S.ST_ABSENT), _cell(status=unread)) == V_NOT_READ, unread
    assert _v(sig_neg, _cell(1e-6, -0.1, status=S.ST_LOW_CONFIDENCE)) == V_LOW_CONFIDENCE
    censored = _cell(status=S.ST_API_CENSORED)  # FinnGen API: only "p >= 0.05" is known
    assert _v(sig_neg, censored) == V_SINGLE
    assert _v(null_neg, censored) == V_BOTH_NULL


def test_threshold_is_point05_over_outcomes_tested_in_both():
    def table(cells_a, cells_b):
        class R:
            dataset_ids = ["a", "b"]
            rows = []

        for i, (ca, cb) in enumerate(zip(cells_a, cells_b)):
            for ds, c in (("a", ca), ("b", cb)):
                if c is None:
                    continue
                c = S.OutcomeRow(**{**c.__dict__, "dataset_id": ds, "outcome_id": f"o{i}",
                                    "harmonized_outcome_id": f"H{i}", "harmonized_label": f"h{i}",
                                    "harmonized_representative": True})
                R.rows.append(c)
        return build_replication(R)

    # 4 harmonized outcomes shared by both catalogs: alpha = 0.05 / 4, whatever was read.
    t = table([_cell(0.01, -0.1), _cell(0.5, 0.1), _cell(0.5, 0.1), _cell(0.01, -0.1)],
              [_cell(0.01, -0.1), _cell(0.5, 0.1), _cell(0.5, 0.1), _cell(status=S.ST_ABSENT)])
    pair = ("a", "b")
    assert t.n_shared[pair] == 4 and t.alpha[pair] == pytest.approx(0.05 / 4)
    by = {r.harmonized_outcome_id: r.verdicts[pair] for r in t.rows}
    assert by == {"H0": V_REPLICATES, "H1": V_BOTH_NULL, "H2": V_BOTH_NULL, "H3": V_NOT_TESTED}
    assert set(t.counts(pair)) == set(VERDICT_ORDER)
    # An outage must not loosen the threshold: with every b cell unreachable but one, the
    # denominator stays the shared count (it used to fall to 1, i.e. alpha 0.05).
    out = table([_cell(0.02, -0.1)] + [_cell(0.5, 0.1)] * 3,
                [_cell(0.02, -0.1)] + [_cell(status=S.ST_UNREACHABLE)] * 3)
    assert out.n_shared[pair] == 4 and out.alpha[pair] == pytest.approx(0.05 / 4)
    assert {r.harmonized_outcome_id: r.verdicts[pair] for r in out.rows}["H0"] == V_BOTH_NULL


def _offline_result(**kw):
    from snp_phewas_demo import (
        FG_OUTCOMES, PAN_OUTCOMES, demo_map, demo_rsid, demo_readers)

    return lookup_variant("rs12740374", readers=demo_readers(), rsid_resolver=demo_rsid,
                          map_fn=demo_map,
                          outcomes={"pan_ukbb_eur": PAN_OUTCOMES, "finngen_r12": FG_OUTCOMES}, **kw)


def test_harmonization_attaches_and_unmapped_stays_a_row():
    r = _offline_result()
    by = {(x.dataset_id, x.outcome_id): x for x in r.rows}
    e78p, e78f = by[("pan_ukbb_eur", "icd10-E78-both_sexes")], by[("finngen_r12", "E4_LIPOPROT")]
    assert e78p.harmonized_outcome_id == e78f.harmonized_outcome_id == "ICD10:E78"
    assert e78f.harmonization_basis == "finngen_published_match"
    assert by[("finngen_r12", "I9_IHD")].harmonized_outcome_id == "PANUKBB:phecode-411-both_sexes"
    chd = by[("finngen_r12", "I9_CHD")]
    assert chd.harmonized_outcome_id is None and chd.harmonization_basis == "icd10_broad_only"
    ldl = by[("pan_ukbb_eur", "biomarkers-30780-both_sexes-irnt")]
    assert ldl.harmonization_basis == "unmapped" and ldl.status == "ok"
    # every row carries a basis; none was dropped for lacking a counterpart
    assert all(x.harmonization_basis for x in r.rows)
    reps = [x for x in r.rows if x.harmonized_outcome_id == "ICD10:E78" and x.dataset_id == "finngen_r12"]
    assert sum(1 for x in reps if x.harmonized_representative) == 1


def test_replication_on_the_fixture_variant_and_not_tested_is_not_null():
    r = _offline_result()
    t = build_replication(r)
    pair = ("pan_ukbb_eur", "finngen_r12")
    by = {row.harmonized_outcome_id: row for row in t.rows}
    assert by["ICD10:E78"].verdicts[pair] == V_REPLICATES
    assert by["PANUKBB:phecode-411-both_sexes"].verdicts[pair] == V_REPLICATES
    # a harmonized outcome the lookup did not read in Pan-UKBB is "not read", never null
    others = [row for row in t.rows if row.cells["pan_ukbb_eur"].status == S.ST_NOT_QUERIED]
    assert others and {row.verdicts[pair] for row in others} == {V_NOT_READ}
    # the threshold counts every outcome both catalogs share, not the 2 this lookup read
    assert t.n_shared[pair] == len(t.rows) > 2
    recs = t.to_records()
    assert {"harmonized_outcome_id", "pan_ukbb_eur__beta", "verdict__pan_ukbb_eur__vs__finngen_r12"} <= set(recs[0])


def test_cross_dataset_figures_render(tmp_path):
    import snp_phewas_plots as plots

    r = _offline_result()
    t = build_replication(r)
    for fn, name in ((plots.replication_matrix, "m.png"), (plots.paired_effects, "p.png")):
        out = fn(t, str(tmp_path / name))
        assert (tmp_path / name).stat().st_size > 10_000, out
    empty = ReplicationTable(["a"], [], {}, {})
    with pytest.raises(ValueError):
        plots.paired_effects(empty, str(tmp_path / "x.png"))
    assert isinstance(t.rows[0], ReplicationRow)


def test_the_representative_is_a_tested_endpoint_before_basis_and_cases():
    """An absent or unreachable representative must not hide a tested sibling's value."""
    from types import SimpleNamespace

    from snp_phewas_harmonize import attach_harmonization

    def row(oid, status, n_cases):
        return S.OutcomeRow(
            dataset_id="finngen_r12", dataset_release="r", outcome_id=oid, outcome_label=oid,
            outcome_category_native="", outcome_kind="binary", n_cases=n_cases,
            n_controls=None, n_total=None, variant_id_grch38="1_1_A_G", variant_id_native=None,
            rs_id=None, effect_allele="G", other_allele="A", status=status)

    big = row("C3_BRAIN_EXALLC", S.ST_UNREACHABLE, 5000)
    small = row("C3_ASTROCYTOMA_EXALLC", S.ST_OK, 100)
    res = SimpleNamespace(dataset_ids=["finngen_r12"], rows=[big, small], uncatalogued=[])
    attach_harmonization(res)
    assert big.harmonized_outcome_id == small.harmonized_outcome_id == "ICD10:C71"
    assert small.harmonized_representative and not big.harmonized_representative
    big.status = S.ST_OK  # both tested: the most cases wins, as before
    attach_harmonization(res)
    assert big.harmonized_representative and not small.harmonized_representative
