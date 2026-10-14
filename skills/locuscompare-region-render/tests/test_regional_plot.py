"""Smoke + harmoniser tests for `regional_plot.py`.

Audit gap (PR #272): `regional_plot.py` is the largest module (~1.3 kLoC) and
was flagged as untested. The existing test suite already touches the renderer
indirectly via the synthetic-demo end-to-end test and the focal-gene-highlight
test; here we add direct-call smoke coverage for:

  - `harmonise_regions_for_locuscompare` (allele flip / palindromic / drop)
  - `render_full_locuscompare` end-to-end with and without a gene track
  - PNG file exists + has nonzero on-disk size + nonzero pixel dimensions

We do NOT chase 100% coverage in one pass — `regional_plot.py` is mostly
matplotlib choreography. Deeper-branch coverage (no-data panel, no-LD panel,
trans-region downsampling) is left for follow-up PRs per the maintainer
audit's [coverage/L] guidance.
"""
from __future__ import annotations

import struct
import sys
from pathlib import Path

import pytest


SKILL_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(SKILL_DIR))

from regional_plot import (  # noqa: E402
    GeneTrackEntry,
    HarmonisedRegionPair,
    RegionalLocusCompareInput,
    harmonise_regions_for_locuscompare,
    render_full_locuscompare,
)
from _wald_ratio_types import LocusVariant  # noqa: E402


# ----- helpers


def _lv(vid, ref, alt, beta, se=0.05, p=1e-8, pip=0.5):
    return LocusVariant(
        variant_id=vid,
        chromosome="1",
        position=int(vid.split("_")[1]),
        ref=ref,
        alt=alt,
        pip=pip,
        beta=beta,
        se=se,
        p_value=p,
        is95=True,
        is99=True,
        r2_lead=1.0,
    )


def _png_pixel_dims(path: Path) -> tuple[int, int]:
    """Pull the (width, height) from a PNG header without depending on Pillow.

    PNG IHDR chunk lives at bytes 16..23 as big-endian uint32 width + height.
    """
    raw = path.read_bytes()
    assert raw.startswith(b"\x89PNG\r\n\x1a\n"), "not a PNG file"
    width = struct.unpack(">I", raw[16:20])[0]
    height = struct.unpack(">I", raw[20:24])[0]
    return width, height


def _stub_input(*, lead_id="1_500_A_G", with_gene_track=True):
    """Build a minimal but renderable RegionalLocusCompareInput.

    Three joined pairs (one is the lead), plus per-side track variants so the
    Manhattan tracks render via the locus-list code path.
    """
    pairs = [
        HarmonisedRegionPair(
            variant_id=lead_id, chromosome="1", position=500, ref="A", alt="G",
            beta_exposure=0.5, se_exposure=0.04, p_exposure=1e-12,
            beta_outcome=0.4, se_outcome=0.05, p_outcome=1e-10,
            r2_with_lead=1.0,
        ),
        HarmonisedRegionPair(
            variant_id="1_400_C_T", chromosome="1", position=400, ref="C", alt="T",
            beta_exposure=0.30, se_exposure=0.05, p_exposure=1e-6,
            beta_outcome=0.20, se_outcome=0.06, p_outcome=1e-5,
            r2_with_lead=0.65,
        ),
        HarmonisedRegionPair(
            variant_id="1_600_C_T", chromosome="1", position=600, ref="C", alt="T",
            beta_exposure=-0.10, se_exposure=0.05, p_exposure=1e-2,
            beta_outcome=-0.05, se_outcome=0.06, p_outcome=1e-1,
            r2_with_lead=0.25,
        ),
    ]
    exposure_track = [_lv(p.variant_id, p.ref, p.alt, p.beta_exposure) for p in pairs]
    outcome_track = [_lv(p.variant_id, p.ref, p.alt, p.beta_outcome) for p in pairs]
    r2_by_variant = {p.variant_id: p.r2_with_lead or 0.0 for p in pairs}

    gene_track = []
    if with_gene_track:
        gene_track = [
            GeneTrackEntry(gene_symbol="DEMOGENE_A", start=100, end=300, strand="+"),
            GeneTrackEntry(gene_symbol="SORT1",     start=450, end=550, strand="-"),
            GeneTrackEntry(gene_symbol="DEMOGENE_C", start=700, end=900, strand="+"),
        ]

    return RegionalLocusCompareInput(
        pairs=pairs,
        lead_variant_id=lead_id,
        chromosome="1",
        window_bp=1000,
        ld_panel_label="synthetic LD (test)",
        window_label=f"+/-500 bp of lead {lead_id}",
        exposure_label="synthetic exposure",
        outcome_label="synthetic outcome",
        provenance_label="test provenance",
        caveats=["test caveat"],
        title="regional_plot smoke",
        exposure_track_variants=exposure_track,
        outcome_track_variants=outcome_track,
        r2_by_variant=r2_by_variant,
        exposure_short_label="synthetic eQTL",
        outcome_short_label="synthetic GWAS",
        gene_track=gene_track,
        focal_gene_symbol="SORT1" if with_gene_track else None,
    )


# ----- harmonise_regions_for_locuscompare


def test_harmonise_regions_keeps_matching_alleles_unchanged():
    expo = [_lv("1_100_A_G", "A", "G", 0.5)]
    outc = [_lv("1_100_A_G", "A", "G", 0.3)]
    pairs = harmonise_regions_for_locuscompare(expo, outc, {"1_100_A_G": 0.9}, "1_100_A_G")
    assert len(pairs) == 1
    p = pairs[0]
    assert p.flip_outcome_beta is False
    assert p.palindromic_excluded is False
    assert p.beta_outcome == pytest.approx(0.3)
    assert p.r2_with_lead == 1.0  # lead always seeded to 1.0


def test_harmonise_regions_flips_swapped_outcome_beta():
    expo = [_lv("1_100_A_G", "A", "G", 0.5)]
    outc = [_lv("1_100_A_G", "G", "A", 0.3)]  # swapped
    pairs = harmonise_regions_for_locuscompare(expo, outc, {"1_100_A_G": 0.9}, "1_999_T_C")
    assert len(pairs) == 1
    assert pairs[0].flip_outcome_beta is True
    assert pairs[0].beta_outcome == pytest.approx(-0.3)


def test_harmonise_regions_flags_palindromic_without_dropping():
    expo = [_lv("1_100_A_T", "A", "T", 0.5)]
    outc = [_lv("1_100_A_T", "A", "T", 0.3)]
    pairs = harmonise_regions_for_locuscompare(expo, outc, {}, "1_999_T_C")
    assert len(pairs) == 1
    assert pairs[0].palindromic_excluded is True


def test_harmonise_regions_drops_irreconcilable_alleles():
    expo = [_lv("1_100_A_G", "A", "G", 0.5)]
    outc = [_lv("1_100_A_G", "C", "T", 0.3)]  # different SNV at same vid
    pairs = harmonise_regions_for_locuscompare(expo, outc, {}, "1_999_T_C")
    assert pairs == []


def test_harmonise_regions_drops_variants_only_on_one_side():
    expo = [_lv("1_100_A_G", "A", "G", 0.5), _lv("1_200_C_T", "C", "T", 0.1)]
    outc = [_lv("1_200_C_T", "C", "T", 0.05)]
    pairs = harmonise_regions_for_locuscompare(expo, outc, {}, "1_999_T_C")
    assert [p.variant_id for p in pairs] == ["1_200_C_T"]


def test_harmonise_regions_attaches_r2_for_non_lead_variants():
    expo = [_lv("1_400_C_T", "C", "T", 0.2)]
    outc = [_lv("1_400_C_T", "C", "T", 0.1)]
    pairs = harmonise_regions_for_locuscompare(
        expo, outc, r2_by_variant={"1_400_C_T": 0.42}, lead_variant_id="1_500_A_G",
    )
    assert pairs[0].r2_with_lead == pytest.approx(0.42)


# ----- render_full_locuscompare smoke


def test_render_full_locuscompare_with_gene_track(tmp_path):
    inp = _stub_input(with_gene_track=True)
    out_path = tmp_path / "smoke_with_gene_track.png"
    returned = render_full_locuscompare(inp, out_path)

    assert returned == out_path
    assert out_path.is_file()
    assert out_path.stat().st_size > 10_000, "rendered PNG should be > 10 KB"
    width, height = _png_pixel_dims(out_path)
    assert width > 0 and height > 0, "PNG must have nonzero pixel dimensions"


def test_render_full_locuscompare_without_gene_track(tmp_path):
    """When gene_track is empty the renderer must still emit a valid PNG —
    the gene-track row collapses but the four other panels render."""
    inp = _stub_input(with_gene_track=False)
    out_path = tmp_path / "smoke_no_gene_track.png"
    render_full_locuscompare(inp, out_path)
    assert out_path.is_file()
    assert out_path.stat().st_size > 10_000
    width, height = _png_pixel_dims(out_path)
    assert width > 0 and height > 0


def test_render_full_locuscompare_creates_parent_directory(tmp_path):
    """Nested output directory must be created on the fly."""
    inp = _stub_input()
    out_path = tmp_path / "nested" / "sub" / "smoke.png"
    render_full_locuscompare(inp, out_path)
    assert out_path.is_file()


def test_render_full_locuscompare_falls_back_to_joined_pairs_for_manhattan(tmp_path):
    """If the orchestrator passes no per-side track lists, the renderer must
    still produce a valid PNG by rendering Manhattans from the joined pairs."""
    inp = _stub_input(with_gene_track=True)
    inp.exposure_track_variants = []
    inp.outcome_track_variants = []
    out_path = tmp_path / "smoke_joined_only.png"
    render_full_locuscompare(inp, out_path)
    assert out_path.is_file()
    assert out_path.stat().st_size > 10_000


def test_render_default_title_reports_the_half_window(tmp_path, monkeypatch):
    """window_bp is the full width of the fetched region. With no title override
    the suptitle reports half of it, the distance covered on each side of the
    lead: 1,000,000 bp renders as "±500 kb", never "±1000 kb". Captured from the
    Figure at savefig time, since the renderer closes the figure before returning."""
    from matplotlib.figure import Figure

    captured: list[str] = []
    real_savefig = Figure.savefig

    def spy(self, *args, **kwargs):
        captured.append(self.get_suptitle())
        return real_savefig(self, *args, **kwargs)

    monkeypatch.setattr(Figure, "savefig", spy)
    inp = _stub_input(lead_id="1_500000_A_G", with_gene_track=False)
    inp.title = None
    inp.window_bp = 1_000_000
    render_full_locuscompare(inp, tmp_path / "half_window.png")
    assert captured == ["Regional LocusCompare: 3 variants joined (±500 kb of 1_500000_A_G)"]
