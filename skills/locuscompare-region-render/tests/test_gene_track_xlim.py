"""Regression: the gene track must render even when the lead is absent from the
harmonised exposure-intersect-outcome pairs.

Root cause that this guards: `xlim_bp` was derived only from a lead found IN
`inp.pairs`; when the shared lead is missing from one side's summary statistics, the
lead was not in `pairs`, `xlim_bp` came out None, and the gene-track render was gated
out, leaving a blank panel. The window must be resolvable from the lead variant id
(and the data extent) independent of `pairs`.
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

SKILL_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(SKILL_DIR))

import regional_plot  # noqa: E402
from regional_plot import (  # noqa: E402
    GeneTrackEntry,
    HarmonisedRegionPair,
    RegionalLocusCompareInput,
    _compute_xlim_bp,
    _position_from_variant_id,
    render_full_locuscompare,
)


def test_position_from_variant_id():
    assert _position_from_variant_id("7_100482234_A_T") == 100482234
    assert _position_from_variant_id("X_12345_A_C") == 12345
    assert _position_from_variant_id(None) is None
    assert _position_from_variant_id("garbage") is None
    assert _position_from_variant_id("7_notanumber_A_T") is None


def test_xlim_from_lead_id_when_lead_absent_from_pairs():
    # Pairs contain a different variant; the lead is NOT among them (it dropped out
    # of the harmonised outcome). xlim must still center on the lead via its id.
    pair = HarmonisedRegionPair(
        variant_id="7_100100000_C_G", chromosome="7", position=100100000,
        ref="C", alt="G", beta_exposure=0.5, se_exposure=0.04, p_exposure=1e-20,
        beta_outcome=0.1, se_outcome=0.05, p_outcome=1e-3, r2_with_lead=0.4,
        flip_outcome_beta=False, palindromic_excluded=False,
    )
    xlim = _compute_xlim_bp([pair], "7_100482234_A_T", 1_000_000)
    assert xlim == (100482234 - 500_000, 100482234 + 500_000)


def test_xlim_falls_back_to_data_extent_without_lead_or_id():
    pair = HarmonisedRegionPair(
        variant_id="7_100100000_C_G", chromosome="7", position=100100000,
        ref="C", alt="G", beta_exposure=0.5, se_exposure=0.04, p_exposure=1e-20,
        beta_outcome=0.1, se_outcome=0.05, p_outcome=1e-3, r2_with_lead=0.4,
        flip_outcome_beta=False, palindromic_excluded=False,
    )
    # Unparseable lead id + no window -> span the data.
    assert _compute_xlim_bp([pair], "no_position_here", None) == (100100000, 100100000)
    assert _compute_xlim_bp([], None, None) is None


def test_gene_track_renders_when_lead_absent_from_pairs(tmp_path: Path, monkeypatch):
    # The end-to-end guard: lead absent from pairs, but a gene track is supplied.
    # Before the fix the gene-track axis was blank; now the gene track is drawn,
    # centred on the lead.
    drawn: list[tuple] = []
    real_render_gene_track = regional_plot.render_gene_track

    def spy(ax, genes, **kwargs):
        drawn.append((kwargs.get("xlim_bp"), kwargs.get("lead_position")))
        return real_render_gene_track(ax, genes, **kwargs)

    monkeypatch.setattr(regional_plot, "render_gene_track", spy)
    pair = HarmonisedRegionPair(
        variant_id="7_100100000_C_G", chromosome="7", position=100100000,
        ref="C", alt="G", beta_exposure=0.5, se_exposure=0.04, p_exposure=1e-20,
        beta_outcome=0.1, se_outcome=0.05, p_outcome=1e-3, r2_with_lead=0.4,
        flip_outcome_beta=False, palindromic_excluded=False,
    )
    genes = [
        GeneTrackEntry(gene_symbol="NYAP1", start=100482000, end=100506000,
                       strand="+", exons=[(100482000, 100483000)]),
        GeneTrackEntry(gene_symbol="MCM7", start=100092000, end=100101000,
                       strand="-", exons=[(100092000, 100093000)]),
    ]
    inp = RegionalLocusCompareInput(
        pairs=[pair],
        lead_variant_id="7_100482234_A_T",  # NOT in pairs
        chromosome="7",
        window_bp=1_000_000,
        ld_panel_label="(no LD panel)",
        window_label="+/-500 kb",
        exposure_label="NYAP1 (synthetic)",
        outcome_label="Alzheimer (synthetic)",
        provenance_label="synthetic render",
        gene_track=genes,
        focal_gene_symbol="NYAP1",
    )
    out_path = tmp_path / "plot.png"
    render_full_locuscompare(inp, out_path)
    assert out_path.is_file() and out_path.stat().st_size > 0
    assert drawn == [((100482234 - 500_000, 100482234 + 500_000), 100482234)]
    plt.close("all")
