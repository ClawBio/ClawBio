"""The three figures: outcome Manhattan per dataset, replication matrix, paired effects.

Static matplotlib PNGs. Colors are the reference categorical slots in fixed order (one
slot per outcome group, the same group the same color in every panel) and the blue <->
red diverging pair with a gray midpoint for signed z. Untested cells are hatched gray and
tested-null cells are plain light gray, so a missing test never reads as a null.
"""

from __future__ import annotations

import math
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
from matplotlib.patches import Patch, Rectangle  # noqa: E402

import snp_phewas_schema as S  # noqa: E402
from snp_phewas_catalog import get_dataset  # noqa: E402
from snp_phewas_replication import (  # noqa: E402
    NO_COMPARISON,
    V_LOW_CONFIDENCE,
    V_NOT_READ,
    V_NOT_TESTED,
    VERDICT_LABELS,
    _missing_kind,
    VERDICT_ORDER,
    ReplicationTable,
)

SERIES = ("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948")
INK, INK_2, MUTED = "#1a1a19", "#5f5e58", "#b9b8b0"
SURFACE = "#fcfcfb"
NULL_FILL, UNTESTED_FILL = "#f0efec", "#dcdbd5"
DIVERGING = LinearSegmentedColormap.from_list(
    "z", ["#184f95", "#6da7ec", NULL_FILL, "#ef8a89", "#a32a2a"])
GROUPS = ("disease", "measurement", "questionnaire", "medication")
GROUP_LABELS = {"disease": "disease", "measurement": "measurement",
                "questionnaire": "questionnaire answer", "medication": "medication use"}
NULL_POINT = "#c9c8c1"
# Signals take categorical slots; "neither significant" and "not tested" stay recessive.
_MISSING_CELL_TEXT = {V_NOT_TESTED: "not tested", V_NOT_READ: "not read",
                      V_LOW_CONFIDENCE: "low confidence"}

VERDICT_COLORS = {
    "replicates": SERIES[0], "concordant_one_significant": SERIES[2],
    "discordant": SERIES[7], "single_dataset_signal": SERIES[1],
    "direction_unknown": SERIES[6], "both_null": NULL_POINT, "not_tested_in_one": UNTESTED_FILL,
    "not_read_in_one": UNTESTED_FILL, "low_confidence_in_one": UNTESTED_FILL,
}
GENOME_WIDE = 5e-8


def outcome_group(row: S.OutcomeRow) -> str:
    """Coarse outcome group used to order and color the Manhattan x axis."""
    fam = get_dataset(row.dataset_id).family
    if fam == "pan_ukbb":
        tt = row.outcome_id.split("-", 1)[0]
        return {"phecode": "disease", "icd10": "disease", "biomarkers": "measurement",
                "continuous": "measurement", "categorical": "questionnaire",
                "prescriptions": "medication"}.get(tt, "disease")
    cat = (row.outcome_category_native or "").lower()
    if "drug purchase" in cat:
        return "medication"
    if "quantitative" in cat or row.outcome_kind == "continuous":
        return "measurement"
    return "disease"


def _style(ax) -> None:
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(MUTED)
    ax.tick_params(colors=INK_2, labelsize=8)
    ax.grid(axis="y", color="#ecebe6", linewidth=0.6)
    ax.set_axisbelow(True)


def manhattan(result: Any, out_png: str, *, n_labels: int = 10) -> str:
    """One stacked panel per dataset: outcomes on x (grouped), -log10 p on y. Labels go on
    the strongest outcomes above the dataset's phenome-wide line only."""
    ids = list(result.dataset_ids)
    fig, axes = plt.subplots(len(ids), 1, figsize=(12, 3.4 * len(ids)), squeeze=False,
                             facecolor=SURFACE)
    for ax, dsid in zip(axes[:, 0], ids):
        _style(ax)
        spec = get_dataset(dsid)
        rows = result.dataset_rows(dsid)
        tested = [r for r in rows if r.status in S.TESTED_STATUSES and r.neglog10_p is not None]
        order = sorted(rows, key=lambda r: (GROUPS.index(outcome_group(r)),
                                            r.outcome_category_native, r.outcome_id))
        xpos = {r.outcome_id: i for i, r in enumerate(order)}
        for g_i, g in enumerate(GROUPS):
            pts = [r for r in tested if outcome_group(r) == g and r.status == S.ST_OK]
            for marker, sel in (("^", lambda r: r.beta is not None and r.beta > 0),
                                ("v", lambda r: r.beta is not None and r.beta <= 0)):
                ps = [r for r in pts if sel(r)]
                ax.scatter([xpos[r.outcome_id] for r in ps], [r.neglog10_p for r in ps],
                           marker=marker, s=22, color=SERIES[g_i], linewidths=0,
                           label=None)
            open_pts = [r for r in pts if r.beta is None]
            ax.scatter([xpos[r.outcome_id] for r in open_pts], [r.neglog10_p for r in open_pts],
                       marker="o", s=22, facecolors="none", edgecolors=SERIES[g_i], linewidths=1)
        low = [r for r in tested if r.status == S.ST_LOW_CONFIDENCE]
        ax.scatter([xpos[r.outcome_id] for r in low], [r.neglog10_p for r in low],
                   marker="o", s=14, color=MUTED, linewidths=0)
        n_tested = len(tested)
        n_bonf = S.phenome_wide_n(rows)
        if n_bonf:
            bonf = -math.log10(0.05 / n_bonf)
            ax.axhline(bonf, color=INK_2, linewidth=1, linestyle="-")
            ax.text(len(order), bonf, f" 0.05/{n_bonf:,}", va="center", fontsize=7, color=INK_2)
        ax.axhline(-math.log10(GENOME_WIDE), color=MUTED, linewidth=1, linestyle="--")
        ax.text(len(order), -math.log10(GENOME_WIDE), " 5e-8 (reference)", va="center",
                fontsize=7, color=MUTED)
        line = -math.log10(0.05 / n_bonf) if n_bonf else float("inf")
        top = sorted((r for r in tested if r.status == S.ST_OK and r.neglog10_p >= line),
                     key=lambda r: -r.neglog10_p)[:n_labels]
        ymax = max([r.neglog10_p for r in tested] + [1.0])
        placed: list[tuple[float, float]] = []
        for r in top:  # skip a label that would sit on one already placed
            x, y = xpos[r.outcome_id], r.neglog10_p
            if any(abs(x - px) < 0.05 * len(order) and abs(y - py) < 0.05 * ymax
                   for px, py in placed):
                continue
            placed.append((x, y))
            ax.annotate(r.outcome_label[:38], (x, y), xytext=(3, 2),
                        textcoords="offset points", fontsize=6.5, color=INK)
        ax.set_xlim(-5, len(order) * 1.06)
        ax.set_ylabel("-log10 p", fontsize=9, color=INK_2)
        ax.set_xticks([])
        ax.set_title(f"{spec.label}: tested in {n_tested:,} of {len(rows):,} outcomes",
                     fontsize=10, color=INK, loc="left")
    handles = [Patch(color=SERIES[i], label=GROUP_LABELS[g]) for i, g in enumerate(GROUPS)]
    handles += [plt.Line2D([], [], marker="^", ls="", color=INK_2, label="effect allele raises"),
                plt.Line2D([], [], marker="v", ls="", color=INK_2, label="effect allele lowers"),
                plt.Line2D([], [], marker="o", ls="", mfc="none", color=INK_2, label="direction withheld"),
                plt.Line2D([], [], marker="o", ls="", color=MUTED, label="low confidence")]
    axes[0, 0].legend(handles=handles, fontsize=7, ncol=4, frameon=False, loc="upper left",
                      bbox_to_anchor=(0, 1.32))
    v = result.variant
    fig.suptitle(f"{v.rs_id or v.variant_id_grch38} ({v.variant_id_grch38}); effect allele "
                 f"{v.grch38.alt}", fontsize=11, color=INK, x=0.01, ha="left", y=1.0)
    fig.tight_layout()
    fig.savefig(out_png, dpi=150, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)
    return out_png


def replication_matrix(table: ReplicationTable, out_png: str, *, max_rows: int = 40) -> str:
    """Harmonized outcome x dataset, cell = signed z (clipped at |z| = 8), verdict strip."""
    ids = table.dataset_ids
    pair = table.pairs[0] if table.pairs else None
    rows = table.rows[:max_rows]
    h = max(2.5, 0.28 * len(rows) + 1.6)
    fig, ax = plt.subplots(figsize=(3.0 + 1.6 * len(ids) + 3.2, h), facecolor=SURFACE)
    ax.set_facecolor(SURFACE)
    for i, r in enumerate(rows):
        y = len(rows) - 1 - i
        for j, ds in enumerate(ids):
            c = r.cells.get(ds)
            if c is None or c.status not in (S.ST_OK, S.ST_API_CENSORED):
                ax.add_patch(Rectangle((j, y), 0.96, 0.9, facecolor=UNTESTED_FILL,
                                       hatch="////", edgecolor="#ffffff", linewidth=0))
                ax.text(j + 0.48, y + 0.45, _MISSING_CELL_TEXT[_missing_kind(c)], ha="center",
                        va="center", fontsize=6, color=INK_2)
                continue
            if c.status == S.ST_API_CENSORED or c.z is None:
                ax.add_patch(Rectangle((j, y), 0.96, 0.9, facecolor=NULL_FILL, linewidth=0))
                ax.text(j + 0.48, y + 0.45, "p>=0.05" if c.status == S.ST_API_CENSORED else
                        f"p={c.p_value:.2g}", ha="center", va="center", fontsize=6, color=INK_2)
                continue
            z = max(-8.0, min(8.0, c.z))
            ax.add_patch(Rectangle((j, y), 0.96, 0.9, facecolor=DIVERGING((z + 8) / 16),
                                   linewidth=0))
            ax.text(j + 0.48, y + 0.45, f"p={c.p_value:.2g}", ha="center", va="center",
                    fontsize=6, color=INK if abs(z) < 5 else "#ffffff")
        if pair:
            v = r.verdicts[pair]
            ax.add_patch(Rectangle((len(ids) + 0.1, y + 0.1), 0.25, 0.7,
                                   facecolor=VERDICT_COLORS[v], linewidth=0,
                                   hatch="////" if v in NO_COMPARISON else None))
            ax.text(len(ids) + 0.45, y + 0.45, VERDICT_LABELS[v], va="center", fontsize=6,
                    color=INK)
        ax.text(-0.05, y + 0.45, r.harmonized_label[:48], ha="right", va="center",
                fontsize=6.5, color=INK)
    ax.set_xlim(-0.05, len(ids) + 3.2)
    ax.set_ylim(-0.2, len(rows))
    ax.set_xticks([j + 0.48 for j in range(len(ids))])
    ax.set_xticklabels([get_dataset(d).label.replace(", ", "\n") for d in ids], fontsize=7,
                       color=INK_2)
    ax.xaxis.tick_top()
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    title = f"Replication across biobanks: {len(rows)} of {len(table.rows)} shared outcomes shown"
    if pair:
        title += (f"; significance at 0.05/{table.n_shared[pair]:,} shared outcomes"
                  if table.n_shared[pair] else "")
    ax.set_title(title, fontsize=9, color=INK, loc="left", pad=28)
    sm = plt.cm.ScalarMappable(cmap=DIVERGING, norm=plt.Normalize(-8, 8))
    cb = fig.colorbar(sm, ax=ax, fraction=0.025, pad=0.01)
    cb.set_label("signed z (beta / SE), clipped at 8", fontsize=7, color=INK_2)
    cb.ax.tick_params(labelsize=6)
    fig.tight_layout()
    fig.savefig(out_png, dpi=150, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)
    return out_png


def paired_effects(table: ReplicationTable, out_png: str) -> str:
    """Beta in dataset A vs dataset B with 95% CIs, one panel for binary outcomes (log-OR)
    and one for continuous outcomes in z units (effect scales differ across cohorts)."""
    if not table.pairs:
        raise ValueError("a paired effect plot needs two datasets")
    a, b = table.pairs[0]
    both = [r for r in table.rows
            if all(r.cells.get(d) is not None and r.cells[d].status == S.ST_OK
                   and r.cells[d].beta is not None for d in (a, b))]
    kinds = [k for k in ("binary", "continuous") if any(r.cells[a].outcome_kind == k for r in both)]
    kinds = kinds or ["binary"]
    fig, axes = plt.subplots(1, len(kinds), figsize=(5.8 * len(kinds), 5), facecolor=SURFACE,
                             squeeze=False)
    for ax, kind in zip(axes[0], kinds):
        _style(ax)
        pts = [r for r in both if r.cells[a].outcome_kind == kind]

        def xy(r):
            if kind == "binary":
                return r.cells[a].beta, r.cells[b].beta, 1.96 * r.cells[a].se, 1.96 * r.cells[b].se
            return r.cells[a].z, r.cells[b].z, 1.96, 1.96

        for v in VERDICT_ORDER:
            sel = [xy(r) for r in pts if r.verdicts[(a, b)] == v]
            if not sel:
                continue
            xs, ys, xe, ye = zip(*sel)
            label = f"{VERDICT_LABELS[v]} ({len(sel)})"
            if v == "both_null":  # recessive: no error bars, small, drawn first
                ax.scatter(xs, ys, s=9, color=NULL_POINT, linewidths=0, label=label, zorder=1)
            else:
                ax.errorbar(xs, ys, xerr=xe, yerr=ye, fmt="o", ms=5, color=VERDICT_COLORS[v],
                            ecolor=VERDICT_COLORS[v], elinewidth=0.9, label=label, zorder=3)
        lim = max([abs(c) for r in pts for c in xy(r)[:2]] + [0.1]) * 1.15
        ax.plot([-lim, lim], [-lim, lim], color=MUTED, linewidth=1, linestyle="--", zorder=0)
        ax.axhline(0, color=MUTED, linewidth=0.6)
        ax.axvline(0, color=MUTED, linewidth=0.6)
        unit = "log odds ratio" if kind == "binary" else "z (beta / SE)"
        ax.set_xlabel(f"{get_dataset(a).label}: {unit}", fontsize=8, color=INK_2)
        ax.set_ylabel(f"{get_dataset(b).label}: {unit}", fontsize=8, color=INK_2)
        title = f"{'Case/control' if kind == 'binary' else 'Continuous'} outcomes tested in both: {len(pts)}"
        if len(pts) >= 3:
            title += f"\neffect correlation r = {_corr(pts, a, b, kind):.2f} (n = {len(pts)})"
        ax.set_title(title, fontsize=9, color=INK, loc="left")
        if pts:
            ax.legend(fontsize=6.5, frameon=False, loc="upper left")
    fig.tight_layout()
    fig.savefig(out_png, dpi=150, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)
    return out_png


def _corr(pts, a, b, kind) -> float:
    xs = [(r.cells[a].beta if kind == "binary" else r.cells[a].z) for r in pts]
    ys = [(r.cells[b].beta if kind == "binary" else r.cells[b].z) for r in pts]
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    sx = math.sqrt(sum((x - mx) ** 2 for x in xs))
    sy = math.sqrt(sum((y - my) ** 2 for y in ys))
    return sxy / (sx * sy) if sx and sy else float("nan")
