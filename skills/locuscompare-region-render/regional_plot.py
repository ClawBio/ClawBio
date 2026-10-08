"""Diagnostic scatter for the lead-variant Wald-ratio.

Produces a per-coloc-pair scatter:
  x = β_exposure (left credible set / QTL)
  y = β_outcome  (right credible set / GWAS or other)
  point area ∝ PIP_left × PIP_right
  lead variant highlighted; Wald-ratio slope drawn through (β_exp_lead, β_out_lead)
  optional 95% CI shading on the slope using SE(WR)

Per pair, this is the LocusCompare-style sanity check that a `betaRatioSignAverage`
of −1 reflects a clean inverse relationship across credible-set members, not a
single outlier.
"""

from __future__ import annotations

import math
import textwrap
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # headless rendering — server / CI safe
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.ticker import FuncFormatter  # noqa: E402

# The shared data classes (HarmonisedPair / LocusVariant / PALINDROMIC_PAIRS /
# WaldRatioResult) live in the bundled `_wald_ratio_types.py` module next to
# this file. Imported here as a sibling.
from _wald_ratio_types import (  # noqa: E402
    HarmonisedPair,
    LocusVariant,
    PALINDROMIC_PAIRS,
    WaldRatioResult,
)

from dataclasses import dataclass, field  # noqa: E402

# LocusZoom r² palette (the conventional 4-bin discrete coloring).
# (lower_bound, upper_bound, color, label).
LD_R2_BINS: list[tuple[float, float, str, str]] = [
    (0.8, 1.000001, "#FF7F0E", "0.8 - 1.0"),  # orange (very strong LD)
    (0.6, 0.8,     "#2CA02C", "0.6 - 0.8"),  # green (strong)
    (0.4, 0.6,     "#87CEEB", "0.4 - 0.6"),  # light blue (moderate)
    (0.2, 0.4,     "#1F3A93", "0.2 - 0.4"),  # navy (weak)
    (0.0, 0.2,     "#999999", "0.0 - 0.2"),  # grey (not in LD)
]
LEAD_COLOR = "#9B30FF"  # purple-square per LocusZoom convention


def _r2_color(r2: float | None) -> tuple[str, str]:
    """Return (hex_color, bin_label) for a given r² value. None -> grey."""
    if r2 is None:
        return LD_R2_BINS[-1][2], LD_R2_BINS[-1][3]
    for lo, hi, color, label in LD_R2_BINS:
        if lo <= r2 < hi:
            return color, label
    return LD_R2_BINS[-1][2], LD_R2_BINS[-1][3]


@dataclass
class HarmonisedRegionPair:
    """One variant joined across exposure + outcome regions, with LD to the lead.

    Variant ids match OT's GRCh38 chr_pos_ref_alt convention. `flip_outcome_beta`
    is True when the outcome row's allele was swapped relative to OT;
    `palindromic_excluded` is True for A/T or G/C SNPs (Hemani 2018 convention
    extended to visualisation here).
    """
    variant_id: str
    chromosome: str
    position: int
    ref: str
    alt: str
    beta_exposure: float | None
    se_exposure: float | None
    p_exposure: float | None
    beta_outcome: float | None
    se_outcome: float | None
    p_outcome: float | None
    r2_with_lead: float | None
    flip_outcome_beta: bool = False
    palindromic_excluded: bool = False


@dataclass
class RegionalLocusCompareInput:
    """Bundled input to `render_full_locuscompare`.

    The renderer is intentionally pure: takes a fully-resolved input with all
    data already joined, harmonised, LD-annotated, and caption-fielded.
    Fetching, harmonising, and LD-computation are upstream (execution skills +
    a thin orchestrator step). This separation keeps the renderer testable
    without external dependencies.
    """
    pairs: list[HarmonisedRegionPair]
    lead_variant_id: str
    chromosome: str
    window_bp: int
    # Caption fields surfaced on the rendered plot:
    ld_panel_label: str          # "1000G Phase 3 v5b GRCh38 (EUR); plink 1.90b6.27"
    window_label: str            # "+/-500 kb of lead 2_36910110_C_T (intersected; PIP_L=0.42 x PIP_R=0.18)"
    exposure_label: str          # "eQTL Catalogue v7.0; study quach_2016_ge_monocyte_iav_ensg00000115808"
    outcome_label: str           # "GWAS Catalog harmonised (snapshot 2026-04-15); study FINNGEN_R12_I9_HEARTFAIL"
    provenance_label: str        # "OT release: 26.03 | Rendered: 2026-05-04T14:32Z"
    caveats: list[str] = field(default_factory=list)
    # Optional title override.
    title: str | None = None
    # Optional full per-side variant lists. When provided, the manhattan tracks
    # render each side's complete region (LocusZoom-style "bottom of plot fills
    # the window") instead of just the joined intersection. r² coloring still
    # comes from the LD lookup against the lead; variants not in the LD result
    # render as grey (r²=None bin) which is correct for "not in LD with lead".
    exposure_track_variants: list[LocusVariant] = field(default_factory=list)
    outcome_track_variants: list[LocusVariant] = field(default_factory=list)
    r2_by_variant: dict[str, float] = field(default_factory=dict)
    # Front-and-center one-line panel titles. Populated by the orchestrator
    # so the plot identifies the outcome and exposure without forcing the
    # reader to parse the small caption. If empty, the renderer falls back
    # to a generic "GWAS outcome" / "eQTL exposure".
    exposure_short_label: str = ""
    outcome_short_label: str = ""
    # Optional gene track: list of GeneTrackEntry. When non-empty, the
    # renderer draws a LocusZoom-convention exon/intron/strand track between
    # the manhattans and the bottom-row scatters. When empty, the gene track
    # row collapses (no panel rendered).
    gene_track: list["GeneTrackEntry"] = field(default_factory=list)
    # Optional focal-gene symbol. When set, the matching entry in `gene_track`
    # renders bold + red on the gene-track panel so the QTL gene is anchored
    # visually in 10-30-gene windows. Threaded from the orchestrator's
    # `spec.exposure_gene_symbol`.
    focal_gene_symbol: str | None = None
    # Scale of the outcome effect sizes on the effect-size panel, when it can be stated
    # (e.g. "log odds ratio, derived from the reported odds ratio"); None leaves the
    # axis as a plain beta.
    outcome_effect_scale_label: str | None = None
    # Why the effect-size panel has nothing to plot, naming the side, e.g. "the outcome
    # study (GCST...) publishes p-values only". The panel states it instead of drawing
    # bare axes. When None and the panel is still empty, a generic reason is drawn.
    effect_size_unavailable_reason: str | None = None


def _build_caption(inp: RegionalLocusCompareInput, n_excluded_palindromic: int) -> str:
    caveats = list(inp.caveats)
    if n_excluded_palindromic:
        caveats.append(f"{n_excluded_palindromic} palindromic-ambiguous variants excluded from LD panel")
    caveats_str = "; ".join(caveats) if caveats else "none"
    return (
        f"LD: {inp.ld_panel_label}\n"
        f"Window: {inp.window_label}\n"
        f"Exposure: {inp.exposure_label}\n"
        f"Outcome: {inp.outcome_label}\n"
        f"{inp.provenance_label}\n"
        f"Caveats: {caveats_str}"
    )


def harmonise_regions_for_locuscompare(
    exposure_variants: list[LocusVariant],
    outcome_variants: list[LocusVariant],
    r2_by_variant: dict[str, float],
    lead_variant_id: str,
) -> list[HarmonisedRegionPair]:
    """Join exposure + outcome rows by variant_id, harmonise alleles, attach r².

    Variants present on only one side are dropped (they cannot appear on the
    LocusCompare scatter or effect-size scatter). For the manhattan tracks,
    the caller can pass the unjoined per-side lists separately if it wants
    side-only variants visible.

    Effect-allele harmonisation rules:
    - Identical (ref, alt) -> no flip needed.
    - Swapped (ref<->alt) -> flip the outcome beta so both are ALT-effect.
    - Palindromic (A/T or G/C) -> mark `palindromic_excluded=True` (caller
      drops from the LD-coloured panel; we keep them in the dataclass for
      auditability).
    - Otherwise (irreconcilable alleles) -> drop.
    """
    by_id_outcome = {v.variant_id: v for v in outcome_variants if v.variant_id}
    out: list[HarmonisedRegionPair] = []
    for ev in exposure_variants:
        if not ev.variant_id:
            continue
        ov = by_id_outcome.get(ev.variant_id)
        if ov is None:
            continue
        if not (ev.ref and ev.alt and ov.ref and ov.alt):
            continue
        ev_ref, ev_alt = ev.ref.upper(), ev.alt.upper()
        ov_ref, ov_alt = ov.ref.upper(), ov.alt.upper()
        is_palindromic = (ev_ref, ev_alt) in PALINDROMIC_PAIRS
        flip = False
        beta_outcome = ov.beta
        if (ev_ref, ev_alt) == (ov_ref, ov_alt):
            pass
        elif (ev_ref, ev_alt) == (ov_alt, ov_ref):
            flip = True
            if beta_outcome is not None:
                beta_outcome = -beta_outcome
        else:
            continue  # irreconcilable; drop
        out.append(HarmonisedRegionPair(
            variant_id=ev.variant_id,
            chromosome=ev.chromosome or "",
            position=ev.position or 0,
            ref=ev_ref,
            alt=ev_alt,
            beta_exposure=ev.beta,
            se_exposure=ev.se,
            p_exposure=ev.p_value,
            beta_outcome=beta_outcome,
            se_outcome=ov.se,
            p_outcome=ov.p_value,
            r2_with_lead=r2_by_variant.get(ev.variant_id) if ev.variant_id != lead_variant_id else 1.0,
            flip_outcome_beta=flip,
            palindromic_excluded=is_palindromic,
        ))
    return out


def _stratified_downsample(
    pairs: list[HarmonisedRegionPair],
    target_max: int = 5000,
) -> list[HarmonisedRegionPair]:
    """Stratified-by-r²-bin sampling.

    Keeps the lead and the highest-r² bin in full; samples from other bins
    proportionally to bin count. Capped at `target_max` total points.
    """
    if len(pairs) <= target_max:
        return pairs
    # Bin pairs by r² range.
    bins: dict[str, list[HarmonisedRegionPair]] = {label: [] for *_, label in LD_R2_BINS}
    for p in pairs:
        _, label = _r2_color(p.r2_with_lead)
        bins[label].append(p)
    # Always keep the highest-r² bin (LD ≥ 0.8).
    highest_label = LD_R2_BINS[0][3]
    keep = list(bins.pop(highest_label, []))
    # Allocate remaining budget across the other bins proportionally.
    remaining = max(0, target_max - len(keep))
    other_total = sum(len(v) for v in bins.values()) or 1
    import random
    rng = random.Random(0)  # deterministic
    for label, members in bins.items():
        share = max(1, int(remaining * len(members) / other_total))
        if len(members) <= share:
            keep.extend(members)
        else:
            keep.extend(rng.sample(members, share))
    return keep


def _render_manhattan_panel(
    ax,
    pairs: list[HarmonisedRegionPair],
    *,
    use_p: str,  # "exposure" or "outcome"
    track_label: str,
    lead_variant_id: str,
    show_xticklabels: bool = True,
    xlim_bp: tuple[int, int] | None = None,
):
    """Render one manhattan track. `xlim_bp` forces the x-axis range to the
    orchestrator's requested window so the lead is visually centered even
    when the data tails off before the window edge.
    """
    bins_data: dict[str, list[HarmonisedRegionPair]] = {label: [] for *_, label in LD_R2_BINS}
    lead_pair: HarmonisedRegionPair | None = None
    for p in pairs:
        if p.variant_id == lead_variant_id:
            lead_pair = p
            continue
        _, label = _r2_color(p.r2_with_lead)
        bins_data[label].append(p)
    # Plot in order: weakest LD first (drawn behind), strongest in front.
    for lo, hi, color, label in reversed(LD_R2_BINS):
        members = bins_data.get(label, [])
        if not members:
            continue
        xs, ys = [], []
        for m in members:
            p_value = m.p_exposure if use_p == "exposure" else m.p_outcome
            nlp = _safe_neg_log10(p_value)
            if nlp is None or m.position is None:
                continue
            xs.append(m.position)
            ys.append(nlp)
        if xs:
            ax.scatter(
                xs, ys, s=18.0, c=color, edgecolors="0.4",
                linewidths=0.3, alpha=0.85,
                label=f"r² {label}", zorder=2,
            )
    if lead_pair is not None:
        p_value = lead_pair.p_exposure if use_p == "exposure" else lead_pair.p_outcome
        nlp = _safe_neg_log10(p_value)
        if nlp is not None:
            ax.scatter(
                [lead_pair.position], [nlp],
                s=140.0, marker="D",
                facecolors=LEAD_COLOR, edgecolors="black", linewidths=1.4,
                label="lead", zorder=4,
            )
            ax.annotate(
                lead_pair.variant_id,
                xy=(lead_pair.position, nlp),
                xytext=(8, 8), textcoords="offset points",
                fontsize=8, color="black",
            )

    ax.set_ylabel(track_label, fontsize=9)
    ax.grid(True, linestyle=":", alpha=0.4)
    # markerscale shrinks the legend lead diamond so it doesn't look like a
    # second lead in the plot. Other entries (r² bins) are already at s=18.
    ax.legend(loc="upper right", fontsize=7, frameon=True, ncol=1, markerscale=0.4)
    ax.xaxis.set_major_formatter(FuncFormatter(_format_position_mb))
    if xlim_bp is not None:
        ax.set_xlim(*xlim_bp)
    if not show_xticklabels:
        ax.tick_params(axis="x", labelbottom=False)


def _shade_not_tested(
    ax,
    data_min_bp: int | None,
    data_max_bp: int | None,
    xlim_bp: tuple[int, int] | None,
) -> None:
    """Shade the regions of the x-axis where the source did not test variants
    (because of cis-window boundaries or sparse GWAS coverage). Done with a
    light grey fill + a small top-edge label so the reader understands why the
    bottom of the plot is empty there.
    """
    if xlim_bp is None or data_min_bp is None or data_max_bp is None:
        return
    lo, hi = xlim_bp
    annotated = False
    if data_min_bp > lo:
        ax.axvspan(lo, data_min_bp, color="0.93", alpha=0.7, zorder=0)
        annotated = True
    if data_max_bp < hi:
        ax.axvspan(data_max_bp, hi, color="0.93", alpha=0.7, zorder=0)
        annotated = True
    if annotated:
        ax.text(
            0.99, 0.97, "shaded = no source data tested",
            transform=ax.transAxes, ha="right", va="top",
            fontsize=6, color="0.45",
        )


def _render_manhattan_panel_from_locus(
    ax,
    locus: list[LocusVariant],
    *,
    track_label: str,
    panel_title: str,
    lead_variant_id: str,
    r2_by_variant: dict[str, float],
    show_xticklabels: bool = True,
    xlim_bp: tuple[int, int] | None = None,
):
    """Draw a manhattan track from a per-side LocusVariant list.

    Variants not in `r2_by_variant` render as grey (r²=None bin). The lead is
    located by `lead_variant_id`; if absent from the list it is not annotated
    (caller is responsible for ensuring the lead variant is present in the
    region fetch).

    `panel_title` is the front-and-center label above the track so the reader
    can identify the source / study without parsing the caption. Region of
    the x-axis where the source did not test variants is shaded.
    """
    bins_data: dict[str, list[tuple[int, float]]] = {label: [] for *_, label in LD_R2_BINS}
    lead_xy: tuple[int, float] | None = None
    data_positions: list[int] = []
    for v in locus:
        if v.position is None:
            continue
        data_positions.append(v.position)
        nlp = _safe_neg_log10(v.p_value)
        if nlp is None:
            continue
        if v.variant_id == lead_variant_id:
            lead_xy = (v.position, nlp)
            continue
        r2 = r2_by_variant.get(v.variant_id) if v.variant_id else None
        _, label = _r2_color(r2)
        bins_data[label].append((v.position, nlp))

    # Mark cis-window bounds with thin dotted verticals so readers see the
    # difference between "no data" and "outside the source's tested window".
    # No shading (per user feedback the panel should look uniformly white);
    # subtle marker plus a small corner annotation.
    if data_positions and xlim_bp is not None:
        data_min_bp = min(data_positions)
        data_max_bp = max(data_positions)
        win_lo, win_hi = xlim_bp
        cis_marker_drawn = False
        if data_min_bp > win_lo:
            ax.axvline(data_min_bp, color="0.55", linewidth=0.8,
                       linestyle=":", alpha=0.7, zorder=1)
            cis_marker_drawn = True
        if data_max_bp < win_hi:
            ax.axvline(data_max_bp, color="0.55", linewidth=0.8,
                       linestyle=":", alpha=0.7, zorder=1)
            cis_marker_drawn = True
        if cis_marker_drawn:
            cis_lo_mb = data_min_bp / 1_000_000
            cis_hi_mb = data_max_bp / 1_000_000
            ax.text(
                0.99, 0.97,
                f"source cis-window: {cis_lo_mb:.2f}–{cis_hi_mb:.2f} Mb",
                transform=ax.transAxes, ha="right", va="top",
                fontsize=6, color="0.45",
            )

    for _, _, color, label in reversed(LD_R2_BINS):
        members = bins_data.get(label, [])
        if not members:
            continue
        ax.scatter(
            [m[0] for m in members],
            [m[1] for m in members],
            s=18.0, c=color, edgecolors="0.4",
            linewidths=0.3, alpha=0.85,
            label=f"r² {label}", zorder=2,
        )
    if lead_xy is not None:
        ax.scatter(
            [lead_xy[0]], [lead_xy[1]],
            s=140.0, marker="D",
            facecolors=LEAD_COLOR, edgecolors="black", linewidths=1.4,
            label="lead", zorder=4,
        )
        ax.annotate(
            lead_variant_id, xy=lead_xy,
            xytext=(8, 8), textcoords="offset points",
            fontsize=8, color="black",
        )
    if panel_title:
        ax.set_title(panel_title, fontsize=10, fontweight="bold", loc="left")
    ax.set_ylabel(track_label, fontsize=9)
    ax.grid(True, linestyle=":", alpha=0.4)
    ax.legend(loc="upper right", fontsize=7, frameon=True, ncol=1, markerscale=0.4)
    ax.xaxis.set_major_formatter(FuncFormatter(_format_position_mb))
    if xlim_bp is not None:
        ax.set_xlim(*xlim_bp)
    if not show_xticklabels:
        ax.tick_params(axis="x", labelbottom=False)


@dataclass
class GeneTrackEntry:
    """Renderer-side representation of one gene. Keeps the renderer decoupled
    from the GENCODE execution skill so callers can pass in any source
    (GENCODE, Ensembl, custom) as long as they map to this shape.
    """
    gene_symbol: str
    start: int                         # 1-based inclusive
    end: int
    strand: str                        # "+" or "-"
    exons: list[tuple[int, int]] = field(default_factory=list)  # list of (start, end)
    biotype: str = "protein_coding"


def _stack_genes_into_rows(
    genes: list[GeneTrackEntry],
    min_gap_bp: int = 20_000,
) -> list[int]:
    """Greedy row assignment to avoid overlap. Returns a row index per gene
    in input order. Two genes share a row only when their spans (with a
    `min_gap_bp` cushion to leave room for labels) do not overlap.
    """
    sorted_idxs = sorted(range(len(genes)), key=lambda i: (genes[i].start, genes[i].end))
    rows: list[int] = [0] * len(genes)
    row_ends: list[int] = []  # last gene-end placed on each row
    for idx in sorted_idxs:
        g = genes[idx]
        placed = False
        for r, rend in enumerate(row_ends):
            if g.start > rend + min_gap_bp:
                rows[idx] = r
                row_ends[r] = g.end
                placed = True
                break
        if not placed:
            rows[idx] = len(row_ends)
            row_ends.append(g.end)
    return rows


FOCAL_GENE_COLOR = "#d62728"  # matplotlib qualitative red — anchors the QTL gene visually


def render_gene_track(
    ax,
    genes: list[GeneTrackEntry],
    *,
    xlim_bp: tuple[int, int],
    lead_position: int | None = None,
    track_label: str = "Genes (GENCODE v39)",
    show_xticklabels: bool = True,
    focal_gene_symbol: str | None = None,
):
    """Render a LocusZoom-style gene track on `ax`.

    Each gene is drawn as: a thin horizontal line spanning [start, end] with
    short triangular arrows along the line indicating transcription direction
    (right for `+`, left for `-`); filled rectangles for each exon stacked on
    the line; gene symbol labeled below the gene with a strand-direction
    indicator (`SYMBOL →` for `+`, `← SYMBOL` for `-`).

    Genes that overlap horizontally are stacked vertically into rows via a
    greedy first-fit algorithm. Y-bounds are snug to the rendered content so
    the panel does not appear taller than its drawing.
    """
    ax.set_xlim(*xlim_bp)
    ax.xaxis.set_major_formatter(FuncFormatter(_format_position_mb))

    if not genes:
        ax.text(
            0.5, 0.5, "no genes in window (after biotype filter)",
            transform=ax.transAxes, ha="center", va="center",
            fontsize=8, color="0.5",
        )
        ax.set_yticks([])
        ax.set_ylabel(track_label, fontsize=9)
        if not show_xticklabels:
            ax.tick_params(axis="x", labelbottom=False)
        return

    rows = _stack_genes_into_rows(genes)
    n_rows = max(rows) + 1
    row_height = 1.0
    exon_height = 0.32
    line_y_offset = 0.55  # gene line slightly above row midline; labels sit below
    label_y_offset = 0.32  # below the gene line, inside [0, n_rows] band

    intron_color = "#4060A0"
    exon_color = "#274A8A"

    for g, r in zip(genes, rows):
        is_focal = bool(focal_gene_symbol) and g.gene_symbol == focal_gene_symbol
        gene_intron_color = FOCAL_GENE_COLOR if is_focal else intron_color
        gene_exon_color = FOCAL_GENE_COLOR if is_focal else exon_color
        gene_label_color = FOCAL_GENE_COLOR if is_focal else "0.15"
        gene_label_weight = "bold" if is_focal else "normal"

        # Top row index (visually) is the LOWEST r; flip so r=0 -> top of axes.
        y_center = (n_rows - 1 - r) * row_height + line_y_offset
        # Intron line.
        ax.hlines(
            y_center, g.start, g.end,
            colors=gene_intron_color, linewidth=(1.4 if is_focal else 1.0), zorder=2,
        )
        # Strand arrows along the intron line. Spaced every ~4% of the window
        # so even short genes get at least one arrow.
        win_lo, win_hi = xlim_bp
        arrow_step = max(int((win_hi - win_lo) * 0.04), 5_000)
        gene_arrow_xs = list(range(g.start + arrow_step // 2, g.end, arrow_step))
        if gene_arrow_xs:
            arrow_marker = ">" if g.strand == "+" else "<"
            ax.scatter(
                gene_arrow_xs, [y_center] * len(gene_arrow_xs),
                marker=arrow_marker, s=14.0,
                facecolors=gene_intron_color, edgecolors="none",
                zorder=3,
            )
        # Exon rectangles.
        for ex_start, ex_end in g.exons:
            ax.add_patch(plt.Rectangle(
                (ex_start, y_center - exon_height / 2),
                width=max(ex_end - ex_start, 1),
                height=exon_height,
                facecolor=gene_exon_color, edgecolor="none", zorder=4,
            ))
        # Label below the gene at its visible-window midpoint.
        label = f"{g.gene_symbol} →" if g.strand == "+" else f"← {g.gene_symbol}"
        x_lab = (max(g.start, win_lo) + min(g.end, win_hi)) / 2
        ax.text(
            x_lab, y_center - label_y_offset,
            label, ha="center", va="top",
            fontsize=7.5, color=gene_label_color,
            fontweight=gene_label_weight,
            zorder=5,
        )

    # Snug y-bounds so the axes box matches the drawing.
    ax.set_ylim(0, n_rows * row_height)
    ax.set_yticks([])
    ax.set_ylabel(track_label, fontsize=8)
    ax.grid(True, axis="x", linestyle=":", alpha=0.4)
    # Drop the top spine so the gene track flows visually into the plot above it.
    ax.spines["top"].set_visible(False)

    if lead_position is not None:
        ax.axvline(lead_position, color=LEAD_COLOR, linewidth=1.0,
                   linestyle="--", alpha=0.6, zorder=1)

    if not show_xticklabels:
        ax.tick_params(axis="x", labelbottom=False)


def _render_locuscompare_scatter(ax, pairs: list[HarmonisedRegionPair], lead_variant_id: str):
    """-log10(p_eQTL) vs -log10(p_GWAS) colored by r²."""
    bins_data: dict[str, list[HarmonisedRegionPair]] = {label: [] for *_, label in LD_R2_BINS}
    lead_pair: HarmonisedRegionPair | None = None
    for p in pairs:
        if p.variant_id == lead_variant_id:
            lead_pair = p
            continue
        _, label = _r2_color(p.r2_with_lead)
        bins_data[label].append(p)
    for _, _, color, label in reversed(LD_R2_BINS):
        members = bins_data.get(label, [])
        xs, ys = [], []
        for m in members:
            x = _safe_neg_log10(m.p_exposure)
            y = _safe_neg_log10(m.p_outcome)
            if x is None or y is None:
                continue
            xs.append(x)
            ys.append(y)
        if xs:
            ax.scatter(xs, ys, s=18.0, c=color, edgecolors="0.4",
                       linewidths=0.3, alpha=0.85, zorder=2)
    if lead_pair is not None:
        lx = _safe_neg_log10(lead_pair.p_exposure)
        ly = _safe_neg_log10(lead_pair.p_outcome)
        if lx is not None and ly is not None:
            ax.scatter([lx], [ly], s=140.0, marker="D",
                       facecolors=LEAD_COLOR, edgecolors="black",
                       linewidths=1.4, zorder=4)
            ax.annotate(lead_pair.variant_id, xy=(lx, ly),
                        xytext=(8, 8), textcoords="offset points",
                        fontsize=8, color="black")

    # Identity diagonal.
    all_x = [_safe_neg_log10(m.p_exposure) for m in pairs] + [
        _safe_neg_log10(m.p_outcome) for m in pairs
    ]
    finite_x = [v for v in all_x if v is not None]
    if finite_x:
        lo, hi = min(finite_x), max(finite_x)
        pad = 0.05 * max(hi - lo, 1.0)
        ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad],
                color="0.5", linestyle=":", linewidth=0.9, zorder=1,
                label="identity (y = x)")
    ax.set_xlabel("−log10 p (eQTL / exposure)", fontsize=9)
    ax.set_ylabel("−log10 p (GWAS / outcome)", fontsize=9)
    ax.set_title(
        "LocusCompare panel: clean diagonal supports H4; two clusters supports H3",
        fontsize=9,
    )
    ax.grid(True, linestyle=":", alpha=0.4)


# Said on the effect-size panel when it has nothing to plot and the caller gave no
# reason. Never bare axes: an empty panel reads as "no effect", which is a claim.
EFFECT_SIZE_EMPTY_FALLBACK_REASON = (
    "no variant joined across the two studies carries an effect size on both sides"
)


def _render_effect_size_scatter(
    ax,
    pairs: list[HarmonisedRegionPair],
    lead_variant_id: str,
    *,
    outcome_scale_label: str | None = None,
    unavailable_reason: str | None = None,
) -> int:
    """β_eQTL vs β_GWAS colored by r². Optional Wald-ratio slope through origin.

    Returns the number of points drawn. When that is zero the panel states why
    (`unavailable_reason`, else a generic sentence) instead of drawing empty axes.
    `outcome_scale_label` names the outcome axis scale when known (e.g. log odds
    ratio derived from the reported odds ratio).
    """
    bins_data: dict[str, list[HarmonisedRegionPair]] = {label: [] for *_, label in LD_R2_BINS}
    lead_pair: HarmonisedRegionPair | None = None
    for p in pairs:
        if p.variant_id == lead_variant_id:
            lead_pair = p
            continue
        _, label = _r2_color(p.r2_with_lead)
        bins_data[label].append(p)
    n_drawn = 0
    ax.axhline(0, color="0.85", linewidth=0.8, zorder=0)
    ax.axvline(0, color="0.85", linewidth=0.8, zorder=0)
    for _, _, color, label in reversed(LD_R2_BINS):
        members = bins_data.get(label, [])
        xs, ys = [], []
        for m in members:
            if m.beta_exposure is None or m.beta_outcome is None:
                continue
            xs.append(m.beta_exposure)
            ys.append(m.beta_outcome)
        if xs:
            n_drawn += len(xs)
            ax.scatter(xs, ys, s=18.0, c=color, edgecolors="0.4",
                       linewidths=0.3, alpha=0.85, zorder=2)
    if lead_pair is not None and lead_pair.beta_exposure is not None and lead_pair.beta_outcome is not None:
        n_drawn += 1
        ax.scatter([lead_pair.beta_exposure], [lead_pair.beta_outcome],
                   s=140.0, marker="D", facecolors=LEAD_COLOR,
                   edgecolors="black", linewidths=1.4, zorder=4)
        ax.annotate(lead_pair.variant_id,
                    xy=(lead_pair.beta_exposure, lead_pair.beta_outcome),
                    xytext=(8, 8), textcoords="offset points",
                    fontsize=8, color="black")
        if lead_pair.beta_exposure != 0:
            wr = lead_pair.beta_outcome / lead_pair.beta_exposure
            xs = [m.beta_exposure for m in pairs if m.beta_exposure is not None]
            if xs:
                lo, hi = min(xs), max(xs)
                pad = 0.05 * max(abs(lo), abs(hi), 1e-9)
                line_xs = [lo - pad, hi + pad]
                ax.plot(line_xs, [wr * x for x in line_xs],
                        color="black", linestyle="--", linewidth=0.9,
                        zorder=1, label=f"WR slope (lead) = {wr:+.3g}")
    ax.set_xlabel("β (eQTL / exposure)", fontsize=9)
    if outcome_scale_label:
        ax.set_ylabel(f"{outcome_scale_label}\n(GWAS / outcome)", fontsize=8)
    else:
        ax.set_ylabel("β (GWAS / outcome)", fontsize=9)
    ax.set_title("Effect-size scatter (LocusCompareR convention)", fontsize=9)
    if n_drawn == 0:
        reason = unavailable_reason or EFFECT_SIZE_EMPTY_FALLBACK_REASON
        ax.set_xticks([])
        ax.set_yticks([])
        ax.text(
            0.5, 0.5,
            textwrap.fill(f"No effect sizes to plot: {reason}.", width=46),
            transform=ax.transAxes, ha="center", va="center",
            fontsize=9, color="#2c3e50",
            bbox=dict(boxstyle="round,pad=0.5", facecolor="#fafafa",
                      edgecolor="#bfbfbf", linewidth=0.8),
        )
        return 0
    ax.grid(True, linestyle=":", alpha=0.4)
    if ax.get_legend_handles_labels()[0]:
        ax.legend(loc="best", fontsize=7, frameon=True)
    return n_drawn


def render_full_locuscompare(
    inp: RegionalLocusCompareInput,
    out_path: Path,
) -> Path:
    """Tier-2 renderer: 4-panel figure with caption disclosure.

    Layout follows the LocusZoom convention of GWAS-on-top, eQTL-below:
    - Two stacked manhattan tracks: GWAS outcome on top, eQTL exposure below.
      LD-coloured in LocusZoom palette; lead as purple diamond.
    - LocusCompare scatter (-log10 p_eQTL vs -log10 p_GWAS) with identity line.
    - Effect-size scatter (β_eQTL vs β_GWAS) with WR-slope-through-origin.
    - 6-field caption disclosure under the figure.

    Palindromic variants are excluded from the LD-coloured panels; the
    excluded count appears in the caveats line. Stratified-by-r²-bin
    downsampling kicks in past 5000 points. The manhattan tracks are forced
    to ±window_bp/2 around the lead so the lead is visually centered even
    when the data tails off (e.g. cis-eQTL window narrower than the GWAS
    window).
    """
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    pairs_excl_palindromic = [p for p in inp.pairs if not p.palindromic_excluded]
    n_palindromic_excluded = len(inp.pairs) - len(pairs_excl_palindromic)
    pairs = _stratified_downsample(pairs_excl_palindromic)

    # Compute window-centered xlim from the lead position so matplotlib does
    # not auto-snap to the data's actual extent.
    xlim_bp: tuple[int, int] | None = None
    lead_pos: int | None = next(
        (p.position for p in inp.pairs if p.variant_id == inp.lead_variant_id and p.position),
        None,
    )
    if lead_pos is not None and inp.window_bp:
        half = inp.window_bp // 2
        xlim_bp = (max(0, lead_pos - half), lead_pos + half)

    has_gene_track = bool(inp.gene_track)
    if has_gene_track:
        # Per-row gene-track height: ~0.45 inches per stacked row, with a
        # small floor so 1-2 rows still look right.
        n_gene_rows = max(_stack_genes_into_rows(inp.gene_track)) + 1
        gene_track_height = max(0.9, 0.45 * n_gene_rows)
        # Three stacked panels in the top section: GWAS, eQTL, genes. Heights
        # in inches so the figure-height calculation stays predictable.
        top_heights = [2.4, 2.4, gene_track_height]
    else:
        top_heights = [2.4, 2.4]

    bottom_height = 4.6  # square-ish scatter row
    fig_height = sum(top_heights) + bottom_height + 1.4  # title + caption padding
    fig = plt.figure(figsize=(13.0, fig_height))

    # Nested gridspec: outer 2-row split (top stack vs scatter row) with a
    # comfortable gap; inner top stack with a TIGHT gap so GWAS / eQTL /
    # genes read as one unified plot with three parts.
    outer = fig.add_gridspec(
        nrows=2,
        height_ratios=[sum(top_heights), bottom_height],
        hspace=0.22,
    )
    inner_top = outer[0].subgridspec(
        nrows=len(top_heights), ncols=1,
        height_ratios=top_heights,
        hspace=0.10,
    )
    inner_bottom = outer[1].subgridspec(
        nrows=1, ncols=2,
        width_ratios=[1.0, 1.0],
        wspace=0.28,
    )

    # LocusZoom convention: outcome (GWAS) on top, exposure (eQTL) below.
    ax_gwas = fig.add_subplot(inner_top[0])
    ax_eqtl = fig.add_subplot(inner_top[1], sharex=ax_gwas)
    if has_gene_track:
        ax_genes = fig.add_subplot(inner_top[2], sharex=ax_gwas)
    else:
        ax_genes = None
    ax_lc = fig.add_subplot(inner_bottom[0])
    ax_es = fig.add_subplot(inner_bottom[1])

    outcome_panel_title = inp.outcome_short_label or "GWAS outcome"
    exposure_panel_title = inp.exposure_short_label or "eQTL exposure"
    # Per-side full region beats joined-intersection for the LocusZoom
    # convention "bottom of plot fills the window". Falls back to joined
    # pairs when the orchestrator doesn't pass per-side extras.
    if inp.outcome_track_variants:
        _render_manhattan_panel_from_locus(
            ax_gwas, inp.outcome_track_variants,
            track_label="−log10 p",
            panel_title=outcome_panel_title,
            lead_variant_id=inp.lead_variant_id,
            r2_by_variant=inp.r2_by_variant,
            show_xticklabels=False,
            xlim_bp=xlim_bp,
        )
    else:
        _render_manhattan_panel(
            ax_gwas, pairs,
            use_p="outcome",
            track_label="−log10 p",
            lead_variant_id=inp.lead_variant_id,
            show_xticklabels=False,
            xlim_bp=xlim_bp,
        )
        ax_gwas.set_title(outcome_panel_title, fontsize=10, fontweight="bold", loc="left")
    if inp.exposure_track_variants:
        _render_manhattan_panel_from_locus(
            ax_eqtl, inp.exposure_track_variants,
            track_label="−log10 p",
            panel_title=exposure_panel_title,
            lead_variant_id=inp.lead_variant_id,
            r2_by_variant=inp.r2_by_variant,
            show_xticklabels=True,
            xlim_bp=xlim_bp,
        )
    else:
        _render_manhattan_panel(
            ax_eqtl, pairs,
            use_p="exposure",
            track_label="−log10 p",
            lead_variant_id=inp.lead_variant_id,
            show_xticklabels=True,
            xlim_bp=xlim_bp,
        )
        ax_eqtl.set_title(exposure_panel_title, fontsize=10, fontweight="bold", loc="left")
    chrom = inp.chromosome.lstrip("chr") if inp.chromosome else ""
    if has_gene_track and ax_genes is not None and xlim_bp is not None:
        # When a gene track is present, the eQTL panel hides its x-tick labels
        # so the gene track owns the genomic-position axis.
        ax_eqtl.tick_params(axis="x", labelbottom=False)
        render_gene_track(
            ax_genes, inp.gene_track,
            xlim_bp=xlim_bp,
            lead_position=lead_pos,
            track_label="Genes (GENCODE v39 protein_coding)",
            show_xticklabels=True,
            focal_gene_symbol=inp.focal_gene_symbol or None,
        )
        ax_genes.set_xlabel(
            f"genomic position, Mb (chromosome {chrom})" if chrom else "genomic position (Mb)",
            fontsize=9,
        )
    else:
        ax_eqtl.set_xlabel(
            f"genomic position, Mb (chromosome {chrom})" if chrom else "genomic position (Mb)",
            fontsize=9,
        )
    _render_locuscompare_scatter(ax_lc, pairs, inp.lead_variant_id)
    _render_effect_size_scatter(
        ax_es, pairs, inp.lead_variant_id,
        outcome_scale_label=inp.outcome_effect_scale_label,
        unavailable_reason=inp.effect_size_unavailable_reason,
    )

    title = inp.title or (
        f"Regional LocusCompare: {len(inp.pairs)} variants joined "
        f"(±{inp.window_bp // 2000} kb of {inp.lead_variant_id})"
    )
    fig.suptitle(title, fontsize=11)
    caption = _build_caption(inp, n_palindromic_excluded)
    fig.text(
        0.5, 0.005, caption, ha="center", va="bottom",
        fontsize=8, family="monospace",
    )
    fig.subplots_adjust(top=0.94, bottom=0.13)
    fig.savefig(out_path, dpi=160, bbox_inches="tight")
    plt.close(fig)
    return out_path


def _format_position_mb(x: float, _pos: int) -> str:
    """Genomic-position tick formatter — scaled to Mb with 3 decimals."""
    return f"{x / 1_000_000:.3f}"


def _safe_neg_log10(p: float | None) -> float | None:
    """Return -log10(p), guarding against p<=0 / p=None."""
    if p is None or p <= 0 or not math.isfinite(p):
        return None
    return -math.log10(p)


def _intersect_ids(result: WaldRatioResult) -> set[str]:
    return {p.variant_id for p in result.pairs if p.variant_id}


def _data_x_range(pairs: list[HarmonisedPair]) -> tuple[float, float]:
    xs = [p.beta_left for p in pairs]
    if not xs:
        return -1.0, 1.0
    lo, hi = min(xs), max(xs)
    pad = 0.1 * max(abs(lo), abs(hi), 1e-9)
    return lo - pad, hi + pad


def _auto_size_points(pairs: list[HarmonisedPair],
                      target_max_area: float = 260.0,
                      min_area: float = 24.0) -> list[float]:
    """Scale point areas so the largest PIP_product gets `target_max_area` px²
    and zero-PIP points still render at `min_area`. Linear in PIP_product so
    relative weights are preserved.
    """
    pips = [p.pip_product or 0.0 for p in pairs]
    max_pip = max(pips) if pips else 0.0
    if max_pip <= 0:
        return [min_area for _ in pairs]
    return [
        max(min_area, target_max_area * (pp / max_pip))
        for pp in pips
    ]


def render_diagnostic_plot(
    result: WaldRatioResult,
    out_path: Path,
    *,
    title: str | None = None,
    left_label: str = "β (exposure / left credible set)",
    right_label: str = "β (outcome / right credible set)",
    target_max_area: float = 260.0,
) -> Path:
    """Write a PNG to `out_path` and return it.

    Renders even when `result.call == UNDETERMINED` (so the user can audit
    why no lead was picked) — but marks the lead and the WR slope only when
    they exist.
    """
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(7.5, 6.0))
    pairs = result.pairs

    if not pairs:
        ax.text(
            0.5, 0.5,
            "No intersected variants between the two credible sets.",
            transform=ax.transAxes, ha="center", va="center",
        )
        fig.suptitle(title or "Wald-ratio diagnostic: no intersection")
        fig.savefig(out_path, dpi=160, bbox_inches="tight")
        plt.close(fig)
        return out_path

    palindromic = [p for p in pairs if p.palindromic]
    keepers = [p for p in pairs if not p.palindromic]

    sizes_keepers = _auto_size_points(keepers, target_max_area=target_max_area)
    sizes_palindromic = _auto_size_points(palindromic, target_max_area=target_max_area)

    ax.axhline(0, color="0.85", linewidth=0.8, zorder=0)
    ax.axvline(0, color="0.85", linewidth=0.8, zorder=0)

    if keepers:
        ax.scatter(
            [p.beta_left for p in keepers],
            [p.beta_right for p in keepers],
            s=sizes_keepers,
            c=[p.pip_product for p in keepers],
            cmap="viridis",
            edgecolors="0.3",
            linewidths=0.5,
            alpha=0.85,
            label=f"intersected ({len(keepers)})",
            zorder=2,
        )

    if palindromic:
        ax.scatter(
            [p.beta_left for p in palindromic],
            [p.beta_right for p in palindromic],
            s=sizes_palindromic,
            facecolors="none",
            edgecolors="crimson",
            linewidths=1.0,
            label=f"palindromic, excluded ({len(palindromic)})",
            zorder=3,
        )

    # Lead annotation + WR slope through origin (standard TwoSampleMR convention).
    if result.lead is not None and result.wr is not None:
        x_lo, x_hi = _data_x_range(keepers)
        xs = [x_lo, x_hi]
        ax.plot(
            xs, [result.wr * x for x in xs],
            color="black", linestyle="--", linewidth=1.0,
            label=f"WR slope = {result.wr:+.4g}",
            zorder=4,
        )

        # Vertical 95 % CI error bar at the lead — uses SE_outcome (the dominant
        # source of uncertainty for most coloc pairs). Horizontal bar uses SE_exp.
        if result.lead.se_right is not None:
            ax.errorbar(
                [result.lead.beta_left],
                [result.lead.beta_right],
                yerr=1.96 * result.lead.se_right,
                xerr=(1.96 * result.lead.se_left) if result.lead.se_left is not None else None,
                fmt="none",
                ecolor="0.25",
                elinewidth=0.8,
                capsize=3,
                zorder=5,
            )

        ax.scatter(
            [result.lead.beta_left],
            [result.lead.beta_right],
            s=sizes_keepers[keepers.index(result.lead)] + 80.0
            if result.lead in keepers else target_max_area + 80.0,
            facecolors="none",
            edgecolors="black",
            linewidths=1.6,
            label="lead (argmax PIPL × PIPR)",
            zorder=6,
        )
        ax.annotate(
            result.lead.variant_id,
            xy=(result.lead.beta_left, result.lead.beta_right),
            xytext=(8, 8), textcoords="offset points",
            fontsize=8, color="black",
        )

    ax.set_xlabel(left_label)
    ax.set_ylabel(right_label)

    if title is None:
        if result.wr is not None:
            title = (
                f"Lead-variant Wald-ratio = {result.wr:+.3f}"
                + (f" (95 % CI {result.ci_lo:+.3f} … {result.ci_hi:+.3f})"
                   if result.ci_lo is not None else "")
                + f"  •  call: {result.call.value}"
                + f"  •  n_intersected={result.n_intersected}, n_eligible={result.n_eligible}"
            )
        else:
            title = (
                f"Wald-ratio undetermined  •  call: {result.call.value}"
                f"  •  n_intersected={result.n_intersected}"
            )
    ax.set_title(title, fontsize=10)
    ax.legend(loc="best", fontsize=8, frameon=True)
    ax.grid(True, linestyle=":", alpha=0.4)

    fig.tight_layout()
    fig.savefig(out_path, dpi=160, bbox_inches="tight")
    plt.close(fig)
    return out_path


def _track_xlims(
    left: list[LocusVariant],
    right: list[LocusVariant],
    pad_frac: float = 0.05,
) -> tuple[int, int] | None:
    positions: list[int] = []
    for v in (*left, *right):
        if v.position is not None:
            positions.append(v.position)
    if not positions:
        return None
    lo, hi = min(positions), max(positions)
    span = max(hi - lo, 1)
    pad = max(int(round(span * pad_frac)), 1)
    return lo - pad, hi + pad


def _render_one_track(
    ax,
    locus: list[LocusVariant],
    intersected_ids: set[str],
    lead_id: str | None,
    *,
    track_label: str,
    color_main: str,
    show_xticklabels: bool = True,
):
    xs, ys, sizes, alphas = [], [], [], []
    inter_xs, inter_ys = [], []
    lead_xy: tuple[int, float] | None = None

    for v in locus:
        if v.position is None:
            continue
        nlp = _safe_neg_log10(v.p_value)
        if nlp is None:
            continue
        xs.append(v.position)
        ys.append(nlp)
        sizes.append(20.0 + 80.0 * (v.pip or 0.0))
        alphas.append(0.4 + 0.5 * (v.pip or 0.0))
        if v.variant_id and v.variant_id in intersected_ids:
            inter_xs.append(v.position)
            inter_ys.append(nlp)
        if v.variant_id == lead_id:
            lead_xy = (v.position, nlp)

    if xs:
        ax.scatter(xs, ys, s=sizes, c=color_main, alpha=0.55,
                   edgecolors="0.3", linewidths=0.4, zorder=2,
                   label=f"all CS variants ({len(xs)})")
    if inter_xs:
        ax.scatter(inter_xs, inter_ys, s=70.0,
                   facecolors="none", edgecolors="firebrick", linewidths=1.2,
                   zorder=3, label=f"intersected ({len(inter_xs)})")
    if lead_xy is not None:
        ax.scatter([lead_xy[0]], [lead_xy[1]], s=140.0,
                   facecolors="none", edgecolors="black", linewidths=1.8,
                   zorder=4, label="lead (argmax PIPL × PIPR)")
        ax.annotate(
            lead_id, xy=lead_xy, xytext=(8, 6), textcoords="offset points",
            fontsize=8, color="black",
        )
    ax.set_ylabel(track_label, fontsize=9)
    ax.grid(True, linestyle=":", alpha=0.4)
    ax.legend(loc="upper right", fontsize=7, frameon=True)
    if not show_xticklabels:
        ax.tick_params(axis="x", labelbottom=False)


def render_regional_locuscompare(
    left_locus: list[LocusVariant],
    right_locus: list[LocusVariant],
    result: WaldRatioResult,
    out_path: Path,
    *,
    title: str | None = None,
    left_track_label: str = "−log10 p (eQTL / left CS)",
    right_track_label: str = "−log10 p (GWAS / right CS)",
) -> Path:
    """Tier-1 LocusCompare: stacked manhattan tracks (eQTL on top, GWAS below)
    using ONLY credible-set members already returned by `credibleSet.locus`,
    plus a per-pair −log10(p) scatter at the intersected variants.

    Limitations vs. full LocusCompare:
    - Only credible-set members are plotted (5–100 variants typically), not all
      variants in a ±500 kb region. The full version is v2 — needs eQTL
      Catalogue / GTEx all-pairs + GWAS Catalog harmonised summary stats.
    - No LD r² coloring (no 1000G PLINK2 reference in v1).
    """
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    intersected_ids = _intersect_ids(result)
    lead_id = result.lead.variant_id if result.lead is not None else None

    fig, axes = plt.subplots(
        nrows=3, ncols=1,
        figsize=(8.5, 10.5),
        gridspec_kw={"height_ratios": [3, 3, 4], "hspace": 0.55},
    )
    ax_eqtl, ax_gwas, ax_lc = axes

    xlims = _track_xlims(left_locus, right_locus)
    chrom = next(
        (v.chromosome for v in (*left_locus, *right_locus) if v.chromosome),
        None,
    )

    _render_one_track(
        ax_eqtl, left_locus, intersected_ids, lead_id,
        track_label=left_track_label, color_main="steelblue",
        show_xticklabels=False,
    )
    _render_one_track(
        ax_gwas, right_locus, intersected_ids, lead_id,
        track_label=right_track_label, color_main="darkorange",
        show_xticklabels=True,
    )
    if xlims is not None:
        ax_eqtl.set_xlim(*xlims)
        ax_gwas.set_xlim(*xlims)
    for ax in (ax_eqtl, ax_gwas):
        ax.xaxis.set_major_formatter(FuncFormatter(_format_position_mb))
    ax_gwas.set_xlabel(
        f"genomic position, Mb (chromosome {chrom})"
        if chrom else "genomic position (Mb)"
    )

    # Bottom panel: per-pair −log10(p) scatter (THIS is the LocusCompare
    # diagonal — clean diagonal supports H4; two clusters supports H3).
    pair_xs, pair_ys, pair_sizes, pair_labels = [], [], [], []
    lead_xy_lc: tuple[float, float] | None = None
    for p in result.pairs:
        nlp_l = _safe_neg_log10(p.p_left)
        nlp_r = _safe_neg_log10(p.p_right)
        if nlp_l is None or nlp_r is None:
            continue
        pair_xs.append(nlp_l)
        pair_ys.append(nlp_r)
        pair_sizes.append(40.0 + 220.0 * (p.pip_product or 0.0) /
                          max((result.lead.pip_product if result.lead else 1.0), 1e-9))
        pair_labels.append(p.variant_id)
        if result.lead is not None and p.variant_id == result.lead.variant_id:
            lead_xy_lc = (nlp_l, nlp_r)
    if pair_xs:
        ax_lc.scatter(
            pair_xs, pair_ys, s=pair_sizes,
            c=[(p.pip_product or 0.0) for p in result.pairs
               if _safe_neg_log10(p.p_left) is not None and _safe_neg_log10(p.p_right) is not None],
            cmap="viridis", edgecolors="0.3", linewidths=0.5, alpha=0.85,
            label=f"intersected ({len(pair_xs)})",
            zorder=2,
        )
        if lead_xy_lc is not None:
            ax_lc.scatter([lead_xy_lc[0]], [lead_xy_lc[1]],
                          s=160.0, facecolors="none", edgecolors="black",
                          linewidths=1.8, zorder=3,
                          label="lead")
            ax_lc.annotate(
                result.lead.variant_id, xy=lead_xy_lc,
                xytext=(8, 8), textcoords="offset points",
                fontsize=8, color="black",
            )
        # Identity line — visual reference for "clean diagonal supports H4".
        lo = min(min(pair_xs), min(pair_ys))
        hi = max(max(pair_xs), max(pair_ys))
        pad = 0.05 * max(hi - lo, 1.0)
        ax_lc.plot([lo - pad, hi + pad], [lo - pad, hi + pad],
                   color="0.5", linestyle=":", linewidth=0.9,
                   zorder=1, label="identity (y = x)")
        ax_lc.set_xlim(lo - pad, hi + pad)
        ax_lc.set_ylim(lo - pad, hi + pad)
    else:
        ax_lc.text(0.5, 0.5,
                   "No intersected variants with p-values on both sides.",
                   transform=ax_lc.transAxes, ha="center", va="center")

    ax_lc.set_xlabel(left_track_label)
    ax_lc.set_ylabel(right_track_label)
    ax_lc.set_title(
        "LocusCompare panel: clean diagonal ⇒ H4; two clusters ⇒ H3",
        fontsize=9,
    )
    ax_lc.legend(loc="upper left", fontsize=7, frameon=True)
    ax_lc.grid(True, linestyle=":", alpha=0.4)

    if title is None:
        cs_l = len([v for v in left_locus if v.position is not None])
        cs_r = len([v for v in right_locus if v.position is not None])
        title = (
            f"Regional credible-set diagnostic  •  "
            f"left CS: {cs_l} variants  •  right CS: {cs_r} variants  •  "
            f"intersected: {result.n_intersected}"
        )
        if result.wr is not None:
            title += f"  •  WR={result.wr:+.3g}"
    fig.suptitle(title, fontsize=10)
    fig.subplots_adjust(top=0.94)
    fig.savefig(out_path, dpi=160, bbox_inches="tight")
    plt.close(fig)
    return out_path


__all__ = [
    "GeneTrackEntry",
    "HarmonisedRegionPair",
    "LD_R2_BINS",
    "LEAD_COLOR",
    "RegionalLocusCompareInput",
    "harmonise_regions_for_locuscompare",
    "render_diagnostic_plot",
    "render_full_locuscompare",
    "render_gene_track",
    "render_regional_locuscompare",
]
