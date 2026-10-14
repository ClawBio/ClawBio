"""The effect-size panel for an outcome that publishes odds ratios, or nothing.

Every test that claims "the panel draws" or "the panel states the reason" goes through
the real render path (`render_locuscompare_for_lead` -> `_render_for_spec` ->
`render_full_locuscompare`) and observes the axes the renderer actually drew, by
wrapping `_render_effect_size_scatter` rather than calling it directly.

The defect being fixed: GCST90476128 (case-control chronic kidney disease) publishes
odds ratios and no beta, the figure never converted them, and the effect-size panel
drew bare axes with no note anywhere.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import pytest

SKILL_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(SKILL_DIR))

import locuscompare_region_render  # noqa: E402  (puts the sibling skills on sys.path)
import regional_plot  # noqa: E402
from _prefetched import PrefetchedEQTLClient, PrefetchedGWASClient  # noqa: E402
from _region_harmonise import (  # noqa: E402
    LOG_ODDS_DERIVED_LABEL,
    LOG_ODDS_FINNGEN_LABEL,
    OR_NO_SE_REASON,
)
from eqtl_catalogue_region_fetch import RegionVariant as EQTLRegionVariant  # noqa: E402
from gwas_catalog_region_fetch import RegionVariant as GWASRegionVariant  # noqa: E402
from locuscompare_region_render import (  # noqa: E402
    LocusCompareSpec,
    render_locuscompare_for_lead,
)

LEAD = "15_45389140_C_T"
CHROM = "15"
LEAD_POS = 45_389_140
GCST = "GCST90476128"

# (variant_id, position, ref, alt, exposure beta, exposure p, outcome OR, CI lo, CI hi, p)
_ROWS = [
    (LEAD, LEAD_POS, "C", "T", 0.40, 1e-20, 1.0516, 1.03, 1.07, 1e-6),
    ("15_45389500_A_G", 45_389_500, "A", "G", 0.30, 1e-12, 1.0400, 1.02, 1.06, 1e-4),
    ("15_45390000_G_A", 45_390_000, "G", "A", -0.10, 1e-3, 0.9800, 0.96, 1.00, 0.04),
    ("15_45391000_C_A", 45_391_000, "C", "A", 0.05, 0.2, 1.0050, 0.99, 1.02, 0.5),
]


def _exposure_rows():
    return [
        EQTLRegionVariant(
            variant_id=vid, chromosome=CHROM, position=pos, ref=ref, alt=alt,
            beta=b, se=0.05, p_value=p, maf=0.3, effect_allele_frequency=0.3,
        )
        for vid, pos, ref, alt, b, p, *_ in _ROWS
    ]


def _or_outcome_rows(*, with_ci: bool = True, with_p: bool = True, flip_lead: bool = False):
    out = []
    for vid, pos, ref, alt, _b, _p, or_, lo, hi, p in _ROWS:
        if flip_lead and vid == LEAD:
            # Same join id, but the outcome's effect allele is the other one: its OR is
            # for C, not T, so the harmoniser must flip the derived beta.
            ref, alt = alt, ref
        out.append(GWASRegionVariant(
            variant_id=vid, chromosome=CHROM, position=pos, ref=ref, alt=alt,
            beta=None, se=None, p_value=p if with_p else None, odds_ratio=or_,
            effect_allele_frequency=0.3,
            ci_lower=lo if with_ci else None, ci_upper=hi if with_ci else None,
        ))
    return out


def _p_only_outcome_rows():
    return [
        GWASRegionVariant(
            variant_id=vid, chromosome=CHROM, position=pos, ref=ref, alt=alt,
            beta=None, se=None, p_value=p, odds_ratio=None, effect_allele_frequency=0.3,
        )
        for vid, pos, ref, alt, _b, _p, _or, _lo, _hi, p in _ROWS
    ]


def _spec(**over) -> LocusCompareSpec:
    base = dict(
        lead_variant_id=LEAD, chromosome=CHROM, lead_position_bp=LEAD_POS,
        window_bp=20_000, eqtl_dataset_id="QTD000341", molecular_trait_id=None,
        gwas_accession=GCST, exposure_gene_symbol="GATM",
        outcome_trait_label="chronic kidney disease", prefetched_gene_track=[],
    )
    base.update(over)
    return LocusCompareSpec(**base)


@pytest.fixture
def es_panel(monkeypatch):
    """Record what the effect-size panel actually drew, on the real render path."""
    seen: dict = {}
    original = regional_plot._render_effect_size_scatter

    def wrapper(ax, pairs, lead_variant_id, **kw):
        n = original(ax, pairs, lead_variant_id, **kw)
        seen.update(
            n_drawn=n,
            pairs=pairs,
            ylabel=ax.get_ylabel(),
            # whitespace-normalised: the renderer wraps the reason onto several lines
            texts=[" ".join(t.get_text().split()) for t in ax.texts],
            line_labels=[ln.get_label() for ln in ax.get_lines()],
        )
        return n

    monkeypatch.setattr(regional_plot, "_render_effect_size_scatter", wrapper)
    return seen


def _render(tmp_path, outcome_rows, *, spec=None, exposure_rows=None):
    return render_locuscompare_for_lead(
        spec or _spec(),
        eqtl_client=PrefetchedEQTLClient(variants=exposure_rows or _exposure_rows()),
        gwas_client=PrefetchedGWASClient(variants=outcome_rows, accession=GCST),
        ld_client=None,
        out_path=tmp_path / "lc.png",
    )


def test_or_only_outcome_draws_points_slope_and_the_derived_label(tmp_path, es_panel):
    result = _render(tmp_path, _or_outcome_rows())

    assert es_panel["n_drawn"] == len(_ROWS)
    assert LOG_ODDS_DERIVED_LABEL in es_panel["ylabel"]
    assert any(lbl.startswith("WR slope (lead)") for lbl in es_panel["line_labels"])
    assert not any("No effect sizes" in t for t in es_panel["texts"])
    lead = next(p for p in es_panel["pairs"] if p.variant_id == LEAD)
    assert lead.beta_outcome == pytest.approx(math.log(1.0516))

    block = result.manifest_block
    assert block["outcome_beta_source"] == "or_derived"
    assert block["outcome_effect_scale_label"] == LOG_ODDS_DERIVED_LABEL
    assert block["effect_size_panel_unavailable_reason"] is None
    assert any("derived from the reported odds ratios" in n for n in result.notes)


def test_conversion_happens_before_the_allele_flip(tmp_path, es_panel):
    """A lead stored with swapped alleles must carry -ln(OR) after harmonisation: the
    derived beta is flipped exactly like a reported one."""
    _render(tmp_path, _or_outcome_rows(flip_lead=True))
    lead = next(p for p in es_panel["pairs"] if p.variant_id == LEAD)
    assert lead.flip_outcome_beta is True
    assert lead.beta_outcome == pytest.approx(-math.log(1.0516))


def test_no_beta_and_no_or_states_the_side_and_reason_everywhere(tmp_path, es_panel):
    result = _render(tmp_path, _p_only_outcome_rows())

    assert es_panel["n_drawn"] == 0
    panel = " ".join(es_panel["texts"])
    assert "No effect sizes to plot" in panel
    assert f"outcome study ({GCST})" in panel and "p-values only" in panel

    reason = result.manifest_block["effect_size_panel_unavailable_reason"]
    assert reason and GCST in reason and "p-values only" in reason
    assert any(n.startswith("effect-size panel not drawn:") and GCST in n for n in result.notes)
    assert result.manifest_block["outcome_beta_source"] is None
    # the p-value panels still had the data: every outcome row joined
    assert result.n_pairs == len(_ROWS)


def test_odds_ratio_without_ci_or_p_states_the_conversion_residual(tmp_path, es_panel):
    result = _render(tmp_path, _or_outcome_rows(with_ci=False, with_p=False))
    assert es_panel["n_drawn"] == 0
    reason = result.manifest_block["effect_size_panel_unavailable_reason"]
    assert "confidence interval" in reason and GCST in reason


def test_exposure_without_effect_sizes_names_the_exposure_side(tmp_path, es_panel):
    exp = [
        EQTLRegionVariant(variant_id=v.variant_id, chromosome=v.chromosome,
                          position=v.position, ref=v.ref, alt=v.alt, beta=None,
                          se=None, p_value=v.p_value, maf=v.maf,
                          effect_allele_frequency=v.effect_allele_frequency)
        for v in _exposure_rows()
    ]
    result = _render(tmp_path, _or_outcome_rows(), exposure_rows=exp)
    reason = result.manifest_block["effect_size_panel_unavailable_reason"]
    assert "exposure study (QTD000341)" in reason
    assert "outcome" not in reason  # the outcome side converted fine


def test_native_beta_outcome_is_unlabelled_and_unchanged(tmp_path, es_panel):
    native = [
        GWASRegionVariant(
            variant_id=vid, chromosome=CHROM, position=pos, ref=ref, alt=alt,
            beta=0.02, se=0.01, p_value=p, odds_ratio=None, effect_allele_frequency=0.3,
        )
        for vid, pos, ref, alt, _b, _p, _or, _lo, _hi, p in _ROWS
    ]
    result = _render(tmp_path, native)
    assert es_panel["ylabel"] == "β (GWAS / outcome)"
    assert result.manifest_block["outcome_beta_source"] == "native"
    assert result.manifest_block["outcome_effect_scale_label"] is None


def test_finngen_outcome_is_labelled_log_odds(tmp_path, es_panel):
    native = [
        GWASRegionVariant(
            variant_id=vid, chromosome=CHROM, position=pos, ref=ref, alt=alt,
            beta=0.02, se=0.01, p_value=p, odds_ratio=None, effect_allele_frequency=0.3,
        )
        for vid, pos, ref, alt, _b, _p, _or, _lo, _hi, p in _ROWS
    ]
    result = _render(
        tmp_path, native,
        spec=_spec(outcome_source_study_id="FINNGEN_R12_N14_CHRONKIDNEYDIS"),
    )
    assert LOG_ODDS_FINNGEN_LABEL in es_panel["ylabel"]
    assert result.manifest_block["outcome_effect_scale_label"] == LOG_ODDS_FINNGEN_LABEL


def test_renderer_never_draws_bare_axes_even_without_a_caller_reason():
    """The pure renderer's own guard, for any caller that passes no reason."""
    import matplotlib.pyplot as plt

    pair = regional_plot.HarmonisedRegionPair(
        variant_id=LEAD, chromosome=CHROM, position=LEAD_POS, ref="C", alt="T",
        beta_exposure=0.4, se_exposure=0.05, p_exposure=1e-20,
        beta_outcome=None, se_outcome=None, p_outcome=1e-6, r2_with_lead=1.0,
    )
    fig, ax = plt.subplots()
    n = regional_plot._render_effect_size_scatter(ax, [pair], LEAD)
    assert n == 0
    texts = [" ".join(t.get_text().split()) for t in ax.texts]
    assert any(regional_plot.EFFECT_SIZE_EMPTY_FALLBACK_REASON in t for t in texts)
    plt.close(fig)


# ---- guards that had no failing test, and the OR = 1 residual ----


def _or_rows_with(ors, *, p=0.01):
    """OR-only outcome rows on the _ROWS positions, no CI, one shared p-value."""
    return [
        GWASRegionVariant(
            variant_id=vid, chromosome=CHROM, position=pos, ref=ref, alt=alt,
            beta=None, se=None, p_value=p, odds_ratio=or_, effect_allele_frequency=0.3,
        )
        for (vid, pos, ref, alt, *_), or_ in zip(_ROWS, ors)
    ]


def test_odds_ratio_of_one_with_a_good_p_is_not_called_a_missing_p_value(tmp_path, es_panel):
    """OR 1.0 and p = 0.5, no CI: the SE cannot be derived because ln(1) = 0, not because
    the p-value is unusable. The reason must say so."""
    result = _render(tmp_path, _or_rows_with([1.0] * len(_ROWS), p=0.5))
    assert es_panel["n_drawn"] == 0
    reason = result.manifest_block["effect_size_panel_unavailable_reason"]
    assert OR_NO_SE_REASON in reason and GCST in reason
    assert "without a confidence interval or usable p-value" not in reason


def test_partly_unconverted_odds_ratios_are_counted_in_a_note(tmp_path, es_panel):
    """Two rows convert (OR 1.2, 1.1) and two cannot (OR 1.0, rounded): the note counts
    the two and gives the true cause."""
    result = _render(tmp_path, _or_rows_with([1.2, 1.0, 1.1, 1.0]))
    assert es_panel["n_drawn"] == 2
    want = (
        f"2 outcome variants report an odds ratio from which {OR_NO_SE_REASON}, so no "
        f"effect size was derived for them; they appear on the p-value panels only"
    )
    assert want in result.notes


def test_a_window_of_only_palindromic_joins_says_so(tmp_path, es_panel):
    """Every joined variant is A/T or C/G: the panel names strand ambiguity, not the
    generic 'no variant carries an effect size' sentence."""
    pal = [
        ("15_45389140_A_T", 45_389_140, "A", "T"),
        ("15_45389500_C_G", 45_389_500, "C", "G"),
        ("15_45390000_T_A", 45_390_000, "T", "A"),
    ]
    exp = [
        EQTLRegionVariant(variant_id=vid, chromosome=CHROM, position=pos, ref=r, alt=a,
                          beta=0.3, se=0.05, p_value=1e-8, maf=0.3,
                          effect_allele_frequency=0.3)
        for vid, pos, r, a in pal
    ]
    out = [
        GWASRegionVariant(variant_id=vid, chromosome=CHROM, position=pos, ref=r, alt=a,
                          beta=0.02, se=0.01, p_value=1e-4, odds_ratio=None,
                          effect_allele_frequency=0.3)
        for vid, pos, r, a in pal
    ]
    result = _render(
        tmp_path, out, exposure_rows=exp,
        spec=_spec(lead_variant_id="15_45389140_A_T"),
    )
    assert result.n_pairs == len(pal)
    reason = result.manifest_block["effect_size_panel_unavailable_reason"]
    assert reason == (
        "every variant joined across the two studies is strand-ambiguous "
        "(palindromic) and is left off this panel"
    )


def test_a_window_mixing_reported_and_derived_effects_is_labelled_mixed(tmp_path, es_panel):
    """One row reports a beta, three only an odds ratio: the derived-log-odds label would
    be false for the reported point, so the axis and the manifest say mixed."""
    rows = _or_outcome_rows()
    first = rows[0]
    rows[0] = GWASRegionVariant(
        variant_id=first.variant_id, chromosome=CHROM, position=first.position,
        ref=first.ref, alt=first.alt, beta=0.9, se=0.1, p_value=first.p_value,
        odds_ratio=None, effect_allele_frequency=0.3,
    )
    result = _render(tmp_path, rows)
    block = result.manifest_block
    assert block["outcome_beta_source"] == "mixed"
    assert LOG_ODDS_DERIVED_LABEL not in es_panel["ylabel"]
    assert "mixed scales" in es_panel["ylabel"]
    assert any(n.startswith("the outcome window mixes effect-size scales: 1 variants")
               for n in result.notes)


def test_a_prefetched_study_is_named_as_a_file_not_by_the_placeholder(tmp_path, es_panel):
    result = _render(
        tmp_path, _p_only_outcome_rows(), spec=_spec(gwas_accession="prefetched"),
    )
    reason = result.manifest_block["effect_size_panel_unavailable_reason"]
    assert "outcome study (supplied as a file)" in reason
    assert "(prefetched)" not in reason


def test_an_ot_finngen_outcome_row_is_labelled_log_odds(tmp_path, es_panel, monkeypatch):
    """Through the Open Targets entry point, the outcome's upstream study id (here a
    FinnGen endpoint mapped to a GWAS Catalog accession) reaches the scale label."""
    monkeypatch.setattr(
        locuscompare_region_render, "fetch_region_genes_remote",
        lambda **kw: ([], {}, []),
    )
    native = [
        GWASRegionVariant(
            variant_id=vid, chromosome=CHROM, position=pos, ref=ref, alt=alt,
            beta=0.02, se=0.01, p_value=p, odds_ratio=None, effect_allele_frequency=0.3,
        )
        for vid, pos, ref, alt, _b, _p, _or, _lo, _hi, p in _ROWS
    ]
    mapping = locuscompare_region_render.StudyIdMapping(
        ot_left_study_id="gtex_ge_kidney_cortex_ensg00000137860",
        ot_right_study_id="FINNGEN_R12_N14_CHRONKIDNEYDIS",
        gwas_catalog_accession=GCST,
        eqtl_catalogue_dataset_id="QTD000341",
    )
    result = locuscompare_region_render.render_tier2_for_lead(
        lead_variant_id=LEAD, chromosome=CHROM, lead_position_bp=LEAD_POS,
        window_bp=20_000, study_mapping=mapping,
        eqtl_client=PrefetchedEQTLClient(variants=_exposure_rows()),
        gwas_client=PrefetchedGWASClient(variants=native, accession=GCST),
        ld_client=None, out_path=tmp_path / "lc.png", ot_release="26.03",
    )
    assert LOG_ODDS_FINNGEN_LABEL in es_panel["ylabel"]
    assert result.manifest_block["outcome_effect_scale_label"] == LOG_ODDS_FINNGEN_LABEL
