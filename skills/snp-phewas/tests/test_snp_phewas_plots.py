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
