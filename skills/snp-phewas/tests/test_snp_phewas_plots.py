"""The Manhattan's category map and its label rule (every outcome above the line is labeled
or counted as left off)."""

from __future__ import annotations

import math
from types import SimpleNamespace

import pytest

import snp_phewas_plots as plots
import snp_phewas_schema as S


def _row(ds, oid, *, cat="", kind="binary", p=None, beta=-0.1, status=S.ST_OK, label=None):
    return S.OutcomeRow(
        dataset_id=ds, dataset_release="r", outcome_id=oid, outcome_label=label or oid,
        outcome_category_native=cat, outcome_kind=kind, n_cases=None, n_controls=None,
        n_total=None, variant_id_grch38="1_1_A_G", variant_id_native=None, rs_id="rs1",
        effect_allele="G", other_allele="A", status=status, beta=beta, se=0.01,
        neglog10_p=p)


@pytest.mark.parametrize("ds,oid,cat,kind,want", [
    ("pan_ukbb_eur", "phecode-411-both_sexes", "circulatory system", "binary", "circulatory"),
    ("pan_ukbb_eur", "phecode-250-both_sexes", "endocrine/metabolic", "binary",
     "endocrine & metabolic"),
    ("pan_ukbb_eur", "icd10-E78-both_sexes", "Chapter IV Endocrine, nutritional and metabolic "
     "diseases | E70-E90 | E78 Disorders", "binary", "endocrine & metabolic"),
    ("pan_ukbb_eur", "icd10-H40-both_sexes", "Chapter VII Diseases of the eye and adnexa",
     "binary", "eye & ear"),
    ("pan_ukbb_eur", "biomarkers-30780-both_sexes-irnt", "Biological samples > Assay results",
     "continuous", "measurements"),
    ("pan_ukbb_eur", "prescriptions-atorvastatin-both_sexes", "HMG CoA reductase inhibitor",
     "binary", "medications"),
    ("pan_ukbb_eur", "categorical-20003-both_sexes-1140861958",
     "UK Biobank Assessment Centre > Verbal interview > Medications", "binary", "medications"),
    ("pan_ukbb_eur", "categorical-1100-both_sexes-2", "UK Biobank Assessment Centre > Touchscreen",
     "binary", "other"),
    ("finngen_r12", "I9_IHD", "IX Diseases of the circulatory system (I9_)", "binary",
     "circulatory"),
    ("finngen_r12", "C3_BREAST", "II Neoplasms, from cancer register (ICD-O-3)", "binary",
     "neoplasms"),
    # FinnGen files a statin purchase under the circulatory chapter; it is a medication
    ("finngen_r12", "RX_STATIN", "IX Diseases of the circulatory system (I9_)", "binary",
     "medications"),
    ("finngen_r12", "M13_RHEUMA", "Rheuma endpoints", "binary", "musculoskeletal"),
    ("finngen_r12", "J10_COPD_X", "Comorbidities of COPD", "binary", "respiratory"),
    ("finngen_r12", "BMI_IRN", "Quantitative endpoints", "continuous", "measurements"),
    ("finngen_r12", "MISC_X", "Miscellaneous, not yet classified endpoints", "binary", "other"),
])
def test_outcome_category_follows_the_sources_own_grouping(ds, oid, cat, kind, want):
    assert plots.outcome_category(_row(ds, oid, cat=cat, kind=kind)) == want
    assert want in plots.CATEGORIES


def test_every_icd10_chapter_and_phecode_group_maps_to_a_category():
    assert set(plots._ICD_CHAPTER.values()) <= set(plots.CATEGORIES)
    assert set(plots._PHECODE.values()) <= set(plots.CATEGORIES)
    assert len(plots._ICD_CHAPTER) == 22


def _result(rows):
    ids = list(dict.fromkeys(r.dataset_id for r in rows))
    return SimpleNamespace(
        dataset_ids=ids, dataset_rows=lambda d: [r for r in rows if r.dataset_id == d],
        variant=SimpleNamespace(rs_id="rs1", variant_id_grch38="1_1_A_G",
                                grch38=SimpleNamespace(alt="G")))


def _nulls(ds, n, cat):
    return [_row(ds, f"N{i}", cat=cat, p=0.3 + (i % 7) / 10) for i in range(n)]


def test_sparse_hits_are_all_labeled(tmp_path):
    cat = "IX Diseases of the circulatory system (I9_)"
    rows = _nulls("finngen_r12", 300, cat) + [
        _row("finngen_r12", "HIT_A", cat=cat, p=30, label="Hit A"),
        _row("finngen_r12", "HIT_B", cat="IV Endocrine (E4_)", p=12, beta=0.2, label="Hit B")]
    stats = {}
    plots.manhattan(_result(rows), str(tmp_path / "m.png"), stats=stats)
    assert stats["finngen_r12"] == {"n_above_line": 2, "n_labeled": 2, "n_unlabeled": 0}


def test_a_crowd_of_hits_is_labeled_or_counted_never_silently_dropped(tmp_path):
    # 120 hits stacked in one narrow category: they cannot all fit, so some are left off
    # and the count is reported (the panel subtitle prints the same number)
    cat = "IX Diseases of the circulatory system (I9_)"
    hits = [_row("finngen_r12", f"HIT_{i}", cat=cat, p=6 + i * 0.05,
                 label=f"A long outcome label number {i}") for i in range(120)]
    rows = _nulls("finngen_r12", 200, "XI Diseases of the digestive system (K11_)") + hits
    stats = {}
    plots.manhattan(_result(rows), str(tmp_path / "m.png"), stats=stats)
    s = stats["finngen_r12"]
    assert s["n_above_line"] == 120
    assert s["n_labeled"] + s["n_unlabeled"] == s["n_above_line"]
    assert s["n_unlabeled"] > 0 and s["n_labeled"] > 0


def test_line_is_the_phenome_wide_bonferroni_and_low_confidence_is_never_a_hit(tmp_path):
    cat = "IX Diseases of the circulatory system (I9_)"
    rows = _nulls("finngen_r12", 99, cat) + [
        _row("finngen_r12", "LOW", cat=cat, p=40, status=S.ST_LOW_CONFIDENCE)]
    at_line = _row("finngen_r12", "AT_LINE", cat=cat)
    rows.append(at_line)
    at_line.neglog10_p = -math.log10(0.05 / S.phenome_wide_n(rows)) + 1e-9
    stats = {}
    plots.manhattan(_result(rows), str(tmp_path / "m.png"), stats=stats)
    assert stats["finngen_r12"]["n_above_line"] == 1


def test_short_label_cuts_at_a_word_boundary():
    assert plots._short("Disorders of lipoprotein metabolism and other lipidaemias", 44) == \
        "Disorders of lipoprotein … other lipidaemias"
    assert plots._short("Angina pectoris") == "Angina pectoris"


def test_display_label_lowers_a_shouted_word_and_keeps_acronyms():
    from snp_phewas_catalog import load_outcome_catalog

    assert S.display_label("STROKE") == "Stroke"
    for acronym in ("COPD", "GERD", "SHBG", "IGF-1", "Angina pectoris"):
        assert S.display_label(acronym) == acronym
    changed = {o.outcome_label for d in ("pan_ukbb_eur", "finngen_r12")
               for o in load_outcome_catalog(d) if S.display_label(o.outcome_label) != o.outcome_label}
    assert changed == {"STROKE"}  # the measured universe: nothing else moves


def test_long_labels_keep_their_distinguishing_end_and_never_collide():
    """The full UK Biobank panel for rs12740374 showed three identical
    "Medication for cholesterol, blood pressure..." labels: an end-cut hid the answer."""
    labels = [
        "Medication for cholesterol, blood pressure, diabetes, or take exogenous hormones: "
        "Cholesterol lowering medication",
        "Medication for cholesterol, blood pressure or diabetes: Cholesterol lowering medication",
        "Medication for cholesterol, blood pressure or diabetes: Blood pressure medication",
        "LDL direct, adjusted by medication (medication-adjusted, rank-normalized)",
        "LDL direct, adjusted by medication (medication-adjusted, raw units)",
        "Angina pectoris",
    ]
    out = plots._distinct_short(labels)
    assert len(set(out)) == len(out)
    assert out[2].endswith(": Blood pressure medication")
    assert out[5] == "Angina pectoris"
    assert all(s.count("(") >= s.count(")") for s in out)


def test_every_catalog_panel_shortens_to_distinct_labels():
    from snp_phewas_catalog import load_outcome_catalog

    for ds in ("pan_ukbb_eur", "finngen_r12"):
        labels = [S.display_label(o.outcome_label) for o in load_outcome_catalog(ds)]
        out = plots._distinct_short(labels)
        assert len(set(out)) == len(set(labels)) == len(labels), ds


def _drawn(monkeypatch, result, tmp_path, **kw):
    """Render a Manhattan and keep its figure open for inspection."""
    figs = []
    monkeypatch.setattr(plots.plt, "close", figs.append)
    plots.manhattan(result, str(tmp_path / "m.png"), **kw)
    return figs[0]


def test_labels_sit_on_a_background_patch(tmp_path, monkeypatch):
    """A halo covers the glyphs and leaves the gaps between words open, so a line under a
    label showed through ("I21.Acute.myocardial"); the label carries a filled patch."""
    cat = "IX Diseases of the circulatory system (I9_)"
    rows = _nulls("finngen_r12", 100, cat) + [
        _row("finngen_r12", "HIT", cat=cat, p=30, label="Acute myocardial infarction")]
    fig = _drawn(monkeypatch, _result(rows), tmp_path)
    (label,) = [t for t in fig.axes[0].texts if t.get_text() == "Acute myocardial infarction"]
    patch = label.get_bbox_patch()
    assert patch is not None, "the label has no background patch"
    assert patch.get_boxstyle().pad == pytest.approx(0.12)
    assert type(patch.get_boxstyle()).__name__ == "Square"
    assert patch.get_alpha() == pytest.approx(0.85)
    from matplotlib.colors import to_hex

    assert to_hex(patch.get_facecolor(), keep_alpha=False) == plots.SURFACE.lower()
    plots.plt.close(fig)


def _subtitle(fig) -> tuple[str, float, float]:
    from matplotlib.backends.backend_agg import FigureCanvasAgg

    r = FigureCanvasAgg(fig).get_renderer()
    (t,) = [t for t in fig.axes[0].texts if "outcomes" in t.get_text() and "·" in t.get_text()]
    return t.get_text(), t.get_window_extent(r).x1, fig.get_window_extent(r).x1


def test_the_subtitle_names_the_matched_read_only_when_that_was_the_reason(tmp_path, monkeypatch):
    """A read narrowed to the matched outcomes says so; a requested subset (the demo's) is
    not "matched to another biobank"."""
    cat = "IX Diseases of the circulatory system (I9_)"
    rows = _nulls("finngen_r12", 50, cat) + [
        _row("finngen_r12", f"Q{i}", cat=cat, status=S.ST_NOT_QUERIED) for i in range(30)]
    res = _result(rows)
    for scope, says in (("harmonized", True), ("requested", False)):
        res.provenance = {"finngen_r12": {"scope": scope}}
        fig = _drawn(monkeypatch, res, tmp_path)
        text = _subtitle(fig)[0]
        assert text.startswith("50 of 80 outcomes looked up")
        assert ("matched to another biobank" in text) is says, text
        plots.plt.close(fig)


def test_the_mvp_subtitle_fits_inside_the_figure(tmp_path, monkeypatch):
    """The MVP panel says how many outcomes were read, how many have no readable file, how
    many are past the line and that the axis is compressed; it once ran past the edge."""
    from snp_phewas_catalog import load_outcome_catalog
    from snp_phewas_lookup import matched_outcome_ids

    matched = matched_outcome_ids("mvp_eur", ["pan_ukbb_eur", "finngen_r12", "mvp_eur"])
    rows = []
    for i, o in enumerate(load_outcome_catalog("mvp_eur")):
        if o.access_status != "harmonised":
            status, p = S.ST_NOT_READABLE, None
        elif o.outcome_id in matched:
            status, p = S.ST_OK, (300.0 if i % 97 == 0 else 0.2 + (i % 5) / 10)
        else:
            status, p = S.ST_NOT_QUERIED, None
        rows.append(_row("mvp_eur", o.outcome_id, cat=o.outcome_category_native,
                         kind=o.outcome_kind, p=p, status=status, label=o.outcome_label))
    res = _result(rows)
    res.provenance = {"mvp_eur": {"scope": "harmonized"}}
    fig = _drawn(monkeypatch, res, tmp_path)
    text, right, edge = _subtitle(fig)
    assert text.startswith("727 of 1,221 readable outcomes looked up (matched to another biobank)")
    assert "849 with no readable file" in text and "compressed" in text
    assert right <= edge, (text, right, edge)
    plots.plt.close(fig)
