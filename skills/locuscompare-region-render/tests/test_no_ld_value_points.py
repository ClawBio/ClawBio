"""A variant with no r² to the lead is drawn hollow under "no LD value", not as the
grey "r² 0.0 - 0.2" bin. The grey bin claims the variant was measured and found out
of LD; a variant missing from the LD panel, or a render with no LD at all, was not
measured, and the legend must not say otherwise."""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import pytest  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402

SKILL_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(SKILL_DIR))

import regional_plot  # noqa: E402

LEAD = "1_67240275_G_A"
NO_LD_LABEL = "no LD value"  # the legend text is the contract


def _pair(vid: str, pos: int, r2: float | None) -> regional_plot.HarmonisedRegionPair:
    return regional_plot.HarmonisedRegionPair(
        variant_id=vid, chromosome="1", position=pos, ref="G", alt="A",
        beta_exposure=0.2, se_exposure=0.05, p_exposure=1e-6,
        beta_outcome=0.1, se_outcome=0.05, p_outcome=1e-3, r2_with_lead=r2,
        flip_outcome_beta=False, palindromic_excluded=False,
    )


@pytest.mark.parametrize("r2", [None, float("nan")])
def test_a_missing_r2_is_its_own_bin(r2):
    assert regional_plot._r2_color(r2)[1] == NO_LD_LABEL
    assert regional_plot._r2_color(0.1)[1] == "0.0 - 0.2"


def _legend_entries(tmp_path: Path, pairs, monkeypatch) -> list[list[tuple[str, object]]]:
    """(label, collection) pairs per legend, read off the figure at save time."""
    seen: list[list[tuple[str, object]]] = []
    real_savefig = Figure.savefig

    def spy(self, *args, **kwargs):
        for ax in self.axes:
            leg = ax.get_legend()
            if leg is None:
                continue
            handles, labels = ax.get_legend_handles_labels()
            if "lead" in labels:  # the two regional (Manhattan) panels
                seen.append(list(zip(labels, handles)))
        return real_savefig(self, *args, **kwargs)

    monkeypatch.setattr(Figure, "savefig", spy)
    inp = regional_plot.RegionalLocusCompareInput(
        pairs=[_pair(LEAD, 67240275, 1.0)] + pairs,
        lead_variant_id=LEAD, chromosome="1", window_bp=1_000_000,
        ld_panel_label="test", window_label="+/-500 kb",
        exposure_label="x", outcome_label="y", provenance_label="test",
    )
    regional_plot.render_full_locuscompare(inp, tmp_path / "plot.png")
    plt.close("all")
    return seen


def test_unmeasured_points_are_hollow_and_named(tmp_path, monkeypatch):
    legends = _legend_entries(
        tmp_path,
        [_pair("1_67241000_C_T", 67241000, 0.1), _pair("1_67242000_C_T", 67242000, None)],
        monkeypatch,
    )
    assert len(legends) == 2, legends
    for entries in legends:
        by_label = dict(entries)
        assert NO_LD_LABEL in by_label and "r² 0.0 - 0.2" in by_label
        hollow = by_label[NO_LD_LABEL].get_facecolors()
        filled = by_label["r² 0.0 - 0.2"].get_facecolors()
        assert len(hollow) == 0 or all(fc[3] == 0 for fc in hollow), hollow
        assert len(filled) and filled[0][3] > 0


def test_a_render_without_ld_claims_no_low_ld_bin(tmp_path, monkeypatch):
    legends = _legend_entries(
        tmp_path,
        [_pair("1_67241000_C_T", 67241000, None), _pair("1_67242000_C_T", 67242000, None)],
        monkeypatch,
    )
    assert legends
    for entries in legends:
        labels = [label for label, _ in entries]
        assert NO_LD_LABEL in labels
        assert "r² 0.0 - 0.2" not in labels


def test_long_caveats_are_wrapped_not_one_wide_line():
    inp = regional_plot.RegionalLocusCompareInput(
        pairs=[_pair(LEAD, 67240275, 1.0)], lead_variant_id=LEAD, chromosome="1",
        window_bp=1_000_000, ld_panel_label="test", window_label="+/-500 kb",
        exposure_label="x", outcome_label="y", provenance_label="test",
        caveats=["word " * 80, "second caveat"],
    )
    caption = regional_plot._build_caption(inp, 0)
    assert max(len(line) for line in caption.splitlines()) <= regional_plot._CAPTION_WRAP
    assert "\n- second caveat" in caption
