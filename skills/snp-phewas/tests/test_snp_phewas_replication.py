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
        FG_OUTCOMES, MVP_OUTCOMES, PAN_OUTCOMES, demo_map, demo_rsid, demo_readers)

    kw.setdefault("datasets", ["pan_ukbb_eur", "finngen_r12"])
    return lookup_variant("rs12740374", readers=demo_readers(), rsid_resolver=demo_rsid,
                          map_fn=demo_map,
                          outcomes={"pan_ukbb_eur": PAN_OUTCOMES, "finngen_r12": FG_OUTCOMES,
                                    "mvp_eur": MVP_OUTCOMES}, **kw)


def test_harmonization_attaches_and_unmapped_stays_a_row():
    r = _offline_result()
    by = {(x.dataset_id, x.outcome_id): x for x in r.rows}
    e78p, e78f = by[("pan_ukbb_eur", "icd10-E78-both_sexes")], by[("finngen_r12", "E4_LIPOPROT")]
    assert e78p.harmonized_outcome_id == e78f.harmonized_outcome_id == "ICD10:E78"
    assert e78f.harmonization_basis == "finngen_published_match"
    assert by[("finngen_r12", "I9_IHD")].harmonized_outcome_id == "PHECODE:411"
    chd = by[("finngen_r12", "I9_CHD")]
    assert chd.harmonized_outcome_id is None and chd.harmonization_basis == "icd10_broad_only"
    # LDL has an id (MVP's LDL mean matches it) but no partner in a Pan-UKBB + FinnGen lookup
    ldl = by[("pan_ukbb_eur", "biomarkers-30780-both_sexes-irnt")]
    assert ldl.harmonization_basis == "curated_name_match" and ldl.status == "ok"
    assert any(x.harmonization_basis == "unmapped" for x in r.rows)  # unmapped rows stay rows
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
    assert by["PHECODE:411"].verdicts[pair] == V_REPLICATES
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


# --- three datasets ---------------------------------------------------------------------

from snp_phewas_replication import (  # noqa: E402
    T_ALL_THREE,
    T_DIRECTION_UNKNOWN,
    T_DISCORDANT,
    T_NULL,
    T_ONE,
    T_TWO,
    V_NO_COUNTERPART,
    three_way,
)


def _t(a, b, c):
    return three_way(*(cell_signal(x, ALPHA) for x in (a, b, c)))


def test_each_three_way_class():
    sig_neg, sig_pos = _cell(1e-6, -0.1), _cell(1e-6, 0.1)
    null = _cell(0.3, -0.01)
    assert _t(sig_neg, sig_neg, sig_neg) == (T_ALL_THREE, None)
    assert _t(sig_neg, null, sig_neg) == (T_TWO, (0, 2))  # the replicating pair is named
    assert _t(sig_neg, sig_pos, null) == (T_DISCORDANT, None)
    assert _t(sig_neg, sig_neg, sig_pos) == (T_DISCORDANT, None)  # any two opposite
    assert _t(sig_neg, _cell(1e-6, None), null) == (T_DIRECTION_UNKNOWN, None)
    assert _t(null, sig_pos, null) == (T_ONE, None)
    assert _t(null, null, _cell(status=S.ST_API_CENSORED)) == (T_NULL, None)
    # a missing cell decides, worst first: not read, then low confidence, then not tested
    assert _t(sig_neg, _cell(status=S.ST_ABSENT), _cell(status=S.ST_NOT_READABLE))[0] == V_NOT_READ
    assert _t(sig_neg, _cell(1e-6, -0.1, status=S.ST_LOW_CONFIDENCE),
              _cell(status=S.ST_ABSENT))[0] == V_LOW_CONFIDENCE
    assert _t(sig_neg, sig_neg, _cell(status=S.ST_ABSENT))[0] == V_NOT_TESTED


def _table3(cells):
    """cells: one tuple (a, b, c) per harmonized outcome; None = no such outcome there."""
    class R:
        dataset_ids = ["a", "b", "c"]
        rows = []

    for i, trio in enumerate(cells):
        for ds, c in zip("abc", trio):
            if c is None:
                continue
            R.rows.append(S.OutcomeRow(**{**c.__dict__, "dataset_id": ds, "outcome_id": f"o{i}",
                                         "harmonized_outcome_id": f"H{i}",
                                         "harmonized_label": f"h{i}",
                                         "harmonized_representative": True}))
    return build_replication(R)


def test_three_way_threshold_no_counterpart_and_an_outage_cannot_loosen_it():
    sig, null = _cell(1e-5, -0.1), _cell(0.5, 0.1)
    t = _table3([(sig, sig, sig), (sig, sig, None), (null, None, null), (sig, sig, sig)])
    assert t.n_shared == {("a", "b"): 3, ("a", "c"): 3, ("b", "c"): 2}
    assert t.n_shared_3 == 2 and t.alpha_3 == pytest.approx(0.05 / 2)
    by = {r.harmonized_outcome_id: r for r in t.rows}
    # H1 has no outcome in c: its a/c and b/c verdicts say so, and it has no three-way class
    assert by["H1"].verdicts[("a", "c")] == V_NO_COUNTERPART and by["H1"].three_way is None
    assert by["H2"].verdicts[("a", "b")] == V_NO_COUNTERPART
    # no_counterpart is not counted in a pair: a/c counts H0, H2, H3 only
    assert sum(t.counts(("a", "c")).values()) == 3 and V_NO_COUNTERPART not in t.counts(("a", "c"))
    assert by["H0"].three_way == T_ALL_THREE
    recs = {r["harmonized_outcome_id"]: r for r in t.to_records()}
    assert recs["H1"]["c__status"] == "no_counterpart" and recs["H1"]["three_way"] == V_NO_COUNTERPART
    assert recs["H0"]["three_way_n_shared"] == 2
    # every c cell unreachable: the three-way denominator stays the shared count
    out = _table3([(sig, sig, _cell(status=S.ST_UNREACHABLE))] * 1 +
                  [(null, null, _cell(status=S.ST_UNREACHABLE))] * 3)
    assert out.n_shared_3 == 4 and out.alpha_3 == pytest.approx(0.05 / 4)
    assert t.counts_three_way()[T_ALL_THREE] == 2
    # with two datasets there is no three-way column
    two = build_replication(_offline_result())
    assert two.n_shared_3 is None and "three_way" not in two.to_records()[0]


def test_matrix_pair_labels_are_one_line_and_name_both_biobanks():
    """A pair label split after the "×" would leave "Pan-UKBB EUR ×" alone on a line, read as
    the title of the MVP column under it."""
    import snp_phewas_plots as plots

    t = build_replication(_offline_result(datasets=["pan_ukbb_eur", "finngen_r12", "mvp_eur"]))
    assert [name for name, _ in plots.pair_strips(t)] == [
        "Pan-UKBB × FinnGen", "Pan-UKBB × MVP", "FinnGen × MVP", "all three"]
    two = build_replication(_offline_result())
    assert [name for name, _ in plots.pair_strips(two)] == ["Pan-UKBB × FinnGen"]


def test_three_biobank_figures_render_and_the_matrix_draws_signals_first(tmp_path):
    import snp_phewas_plots as plots

    r = _offline_result(datasets=["pan_ukbb_eur", "finngen_r12", "mvp_eur"])
    t = build_replication(r)
    assert len(t.pairs) == 3 and t.n_shared_3
    drawn, n_sig, n_null_hidden = plots.select_matrix_rows(t, max_rows=40)
    sig_ids = [x.harmonized_outcome_id for x in drawn[:n_sig]]
    assert {"PHECODE:411", "ICD10:E78", "PANUKBB:biomarkers-30780-both_sexes-irnt"} <= set(sig_ids)
    assert all(x.significant_anywhere(t.alpha) for x in drawn[:n_sig])
    assert len(drawn) == 40
    for fn, name in ((plots.replication_matrix, "m3.png"), (plots.paired_effects, "p3.png")):
        fn(t, str(tmp_path / name))
        assert (tmp_path / name).stat().st_size > 10_000
    # null rows: everything tested and significant nowhere is counted as not drawn
    sig, null = _cell(1e-9, -0.1), _cell(0.5, 0.1)
    many = _table3([(sig, sig, sig)] + [(null, null, null)] * 60)
    drawn, n_sig, n_null_hidden = plots.select_matrix_rows(many, max_rows=40)
    assert (len(drawn), n_sig, n_null_hidden) == (40, 1, 21)
    plots.replication_matrix(many, str(tmp_path / "many.png"))


def test_paired_correlation_over_all_and_over_the_significant():
    from snp_phewas_plots import paired_correlations

    sig = [(_cell(1e-9, b), _cell(1e-9, b * 1.1)) for b in (-0.3, -0.2, 0.25)]
    nulls = [(_cell(0.5, b), _cell(0.5, -b)) for b in (0.01, -0.02, 0.03, -0.01, 0.02)]

    class R:
        dataset_ids = ["a", "b"]
        rows = []

    for i, (ca, cb) in enumerate(sig + nulls):
        for ds, c in (("a", ca), ("b", cb)):
            R.rows.append(S.OutcomeRow(**{**c.__dict__, "dataset_id": ds, "outcome_id": f"o{i}",
                                         "harmonized_outcome_id": f"H{i}", "harmonized_label": "h",
                                         "harmonized_representative": True}))
    t = build_replication(R)
    c = paired_correlations(t.rows, "a", "b", "binary")
    assert (c["n_all"], c["n_sig"]) == (8, 3)
    assert c["r_sig"] == pytest.approx(1.0) and c["r_all"] < c["r_sig"]


def _table_for(ids, cells):
    """`_table3` with real dataset ids, so labels and families come from the catalog."""
    class R:
        dataset_ids = list(ids)
        rows = []

    for i, cs in enumerate(cells):
        for ds, c in zip(ids, cs):
            R.rows.append(S.OutcomeRow(**{**c.__dict__, "dataset_id": ds, "outcome_id": f"o{i}",
                                         "harmonized_outcome_id": f"H{i}",
                                         "harmonized_label": f"h{i}",
                                         "harmonized_representative": True}))
    return build_replication(R)


@pytest.mark.parametrize("ids,n_sig,n_null", [
    (["pan_ukbb_eur", "finngen_r12", "mvp_eur"], 3, 0),
    (["pan_ukbb_eur", "finngen_r12", "mvp_eur"], 40, 0),
    (["pan_ukbb_eur", "pan_ukbb_afr"], 3, 0),  # one biobank twice: the longest labels
    (["pan_ukbb_eur", "pan_ukbb_afr", "mvp_eur"], 3, 0),
    # long labels, three pairs, and outcomes not drawn (a three-line title)
    (["pan_ukbb_eur", "pan_ukbb_afr", "mvp_eur"], 3, 60),
])
def test_the_pair_heading_sits_above_every_pair_label_and_below_the_title(
        tmp_path, monkeypatch, ids, n_sig, n_null):
    """A heading placed from a measurement taken before the final layout printed on the
    title (two datasets of one biobank, long labels, outcomes not drawn)."""
    import snp_phewas_plots as plots

    sig, null = _cell(1e-9, -0.1), _cell(0.5, 0.1)
    t = _table_for(ids, [tuple(sig for _ in ids)] * n_sig + [tuple(null for _ in ids)] * n_null)
    figs = []
    monkeypatch.setattr(plots.plt, "close", figs.append)
    plots.replication_matrix(t, str(tmp_path / "m.png"))
    fig = figs[0]
    ax = fig.axes[0]
    r = plots.FigureCanvasAgg(fig).get_renderer()
    texts = [x for x in ax.texts if x.get_text()]
    heading = next(x for x in texts if x.get_text() == "replication by pair")
    pair_labels = [x for x in texts if x.get_rotation() == 90]
    assert pair_labels
    # the title is left-aligned, so it is the axes' left title, not `ax.title` (empty)
    title = next(x for x in (ax.title, ax._left_title, ax._right_title) if x.get_text())
    assert fig.dpi == plots.MATRIX_DPI  # measured at the resolution the PNG is saved at
    hb, tb = heading.get_window_extent(r), title.get_window_extent(r)
    label_top = max(x.get_window_extent(r).y1 for x in pair_labels)
    (line,) = [ln for ln in ax.lines if len(set(ln.get_ydata())) == 1]
    line_y = ax.transData.transform((0, line.get_ydata()[0]))[1]
    assert label_top < line_y < hb.y0  # the bracket clears the labels; the heading sits on it
    assert hb.y1 < tb.y0
    assert hb.y1 <= ax.get_window_extent(r).y1 + 1  # inside the axes, so the save keeps it
    plots.plt.close(fig)


def test_pair_labels_keep_the_ancestry_only_when_two_datasets_share_a_biobank():
    import snp_phewas_plots as plots

    sig = _cell(1e-9, -0.1)
    two = _table_for(["pan_ukbb_eur", "pan_ukbb_afr"], [(sig, sig)])
    assert [n for n, _ in plots.pair_strips(two)] == ["Pan-UKBB EUR × Pan-UKBB AFR"]
