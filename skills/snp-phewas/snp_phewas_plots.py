"""The three figures: outcome Manhattan per dataset, replication matrix, paired effects.

Static matplotlib PNGs. The Manhattan groups outcomes into body-system categories shared by
every biobank (alternating gray bands for the outcomes below the line, category names on the
x axis) and colors only the outcomes above the phenome-wide line, by direction, with the
blue <-> red pair the replication matrix uses for signed z. Untested cells are hatched gray
and tested-null cells are plain light gray, so a missing test never reads as a null.
"""

from __future__ import annotations

import math
import re
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.patheffects as pe  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.backends.backend_agg import FigureCanvasAgg  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402
from matplotlib.path import Path  # noqa: E402
from matplotlib.scale import FuncScale  # noqa: E402
from matplotlib.transforms import Bbox  # noqa: E402

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
GRID = "#ecebe6"
# Manhattan direction poles: the replication matrix's diverging pair (red = positive z).
RAISES, LOWERS = "#d03b3b", "#2a78d6"
NULL_TONES = ("#a6a59e", "#cfcec7")  # alternating category bands for outcomes below the line
COMPRESS_KNEE = 20.0  # -log10 p above which the Manhattan y axis is log-compressed

# Body-system categories shared by every biobank, in x order.
CATEGORIES = (
    "infectious", "neoplasms", "blood & immune", "endocrine & metabolic", "mental",
    "nervous system", "eye & ear", "circulatory", "respiratory", "digestive", "skin",
    "musculoskeletal", "genitourinary", "pregnancy & perinatal", "congenital", "symptoms",
    "injury & external", "measurements", "medications", "other")
_ICD_CHAPTER = {
    "I": "infectious", "II": "neoplasms", "III": "blood & immune",
    "IV": "endocrine & metabolic", "V": "mental", "VI": "nervous system", "VII": "eye & ear",
    "VIII": "eye & ear", "IX": "circulatory", "X": "respiratory", "XI": "digestive",
    "XII": "skin", "XIII": "musculoskeletal", "XIV": "genitourinary",
    "XV": "pregnancy & perinatal", "XVI": "pregnancy & perinatal", "XVII": "congenital",
    "XVIII": "symptoms", "XIX": "injury & external", "XX": "injury & external",
    "XXI": "other", "XXII": "other"}
_PHECODE = {
    "infectious diseases": "infectious", "neoplasms": "neoplasms",
    "hematopoietic": "blood & immune", "endocrine/metabolic": "endocrine & metabolic",
    "mental disorders": "mental", "neurological": "nervous system", "sense organs": "eye & ear",
    "circulatory system": "circulatory", "respiratory": "respiratory",
    "digestive": "digestive", "dermatologic": "skin", "musculoskeletal": "musculoskeletal",
    "genitourinary": "genitourinary", "pregnancy complications": "pregnancy & perinatal",
    "congenital anomalies": "congenital", "symptoms": "symptoms",
    "injuries & poisonings": "injury & external"}
# FinnGen's non-chapter endpoint groups, matched on the group name's leading words.
_FINNGEN_GROUP = (
    ("rheuma", "musculoskeletal"), ("psychiatric", "mental"), ("alcohol", "mental"),
    ("gastrointestinal", "digestive"), ("neurological", "nervous system"),
    ("diabetes", "endocrine & metabolic"), ("asthma", "respiratory"), ("copd", "respiratory"),
    ("interstitial lung", "respiratory"), ("cardiometabolic", "circulatory"),
    ("drug purchase", "medications"), ("quantitative", "measurements"),
    ("diseases marked as autimmune", "blood & immune"))
_ROMAN = re.compile(r"^(?:Chapter\s+)?([IVX]+)\s")


def outcome_category(row: S.OutcomeRow) -> str:
    """The body-system category a row is drawn under, from the source's own grouping."""
    fam = get_dataset(row.dataset_id).family
    native = (row.outcome_category_native or "").strip()
    if fam == "pan_ukbb":
        tt = row.outcome_id.split("-", 1)[0]
        if tt in ("biomarkers", "continuous"):
            return "measurements"
        if tt == "prescriptions":
            return "medications"
        if tt == "phecode":
            return _PHECODE.get(native.lower(), "other")
        if tt == "icd10":
            m = _ROMAN.match(native)
            return _ICD_CHAPTER.get(m.group(1), "other") if m else "other"
        if "Medications" in native:
            return "medications"
        return "other"
    if row.outcome_id.startswith("RX_"):
        return "medications"
    if row.outcome_kind == "continuous":
        return "measurements"
    m = _ROMAN.match(native)
    if m and m.group(1) in _ICD_CHAPTER:
        return _ICD_CHAPTER[m.group(1)]
    low = native.lower().removeprefix("comorbidities of ")
    for key, cat in _FINNGEN_GROUP:
        if low.startswith(key):
            return cat
    return "other"


def _short(label: str, n: int = 44) -> str:
    """A long outcome label cut from the MIDDLE, at word boundaries: the end is usually what
    tells siblings apart ("Medication for cholesterol, ...: Insulin" vs "...: None of the
    above"; "LDL direct, adjusted by medication (raw units)" vs "(rank-normalized)"). A
    "question: answer" label keeps the whole answer when it fits. Otherwise the front and
    back word counts are chosen to fill the width, the front preferred on a tie."""
    if len(label) <= n:
        return label
    head, sep, answer = label.rpartition(": ")
    if sep and len(answer) <= n - 12:
        return _fit_words(head.split(" "), n - len(answer) - 3, front_only=True) + "…: " + answer
    words = label.split(" ")
    best, best_key = None, None
    for k_back in range(1, len(words)):
        back = " ".join(words[-k_back:])
        room = n - len(back) - 3
        if room < 1:
            break
        front = _fit_words(words[:-k_back], room, front_only=True)
        if not front:
            continue
        # a front shorter than 40% of the width reads as a fragment; fill the rest
        key = (len(front) >= 0.4 * n, len(front) + len(back), len(front))
        if best_key is None or key > best_key:
            best, best_key = f"{front} … {back}", key
    best = best or label[: n - 1] + "…"
    open_, close = best.count("("), best.count(")")
    return best.replace(")", "", close - open_) if close > open_ else best


def _fit_words(words: list[str], room: int, *, front_only: bool) -> str:
    """The longest run of leading whole words within `room` characters, trailing
    punctuation trimmed."""
    out = ""
    for w in words:
        nxt = f"{out} {w}" if out else w
        if len(nxt) > room:
            break
        out = nxt
    return out.rstrip(",;:( ")


def _distinct_short(labels: list[str], n: int = 44) -> list[str]:
    """`_short` for one panel's labels, lengthened until no two read the same. Every full
    label in a dataset's catalog is distinct, so this always ends."""
    out = [_short(lab, n) for lab in labels]
    width = n
    while len(set(out)) < len(out) and width < max(map(len, labels), default=0):
        width += 12
        seen: dict[str, int] = {}
        for s in out:
            seen[s] = seen.get(s, 0) + 1
        out = [_short(lab, width) if seen[s] > 1 else s for lab, s in zip(labels, out)]
    return out


def _compressed_scale(ax, knee: float) -> None:
    """Linear up to `knee`, log-compressed above it: the region around the significance
    lines keeps its resolution when one outcome sits at p = 1e-170."""
    def fwd(y):
        y = np.asarray(y, dtype=float)
        return np.where(y <= knee, y, knee * (1 + np.log10(np.maximum(y, knee) / knee)))

    def inv(u):
        u = np.asarray(u, dtype=float)
        return np.where(u <= knee, u, knee * 10 ** (u / knee - 1))

    ax.set_yscale(FuncScale(ax, (fwd, inv)))


def _nice_ticks(top: float, knee: float | None) -> list[float]:
    if knee is None:
        step = next(s for s in (1, 2, 5, 10, 20, 25, 50, 100) if top / s <= 7)
        return [i * step for i in range(int(top / step) + 1)]
    ticks = [knee * f for f in (0, 0.25, 0.5, 0.75, 1)]
    for t in (30, 50, 100, 200, 300, 500, 1000):
        if knee < t <= top * 1.02:
            ticks.append(t)
    return ticks


def _place_labels(fig, ax, cands, avoid_px, renderer, floor_px) -> tuple[int, int]:
    """Greedy measured placement, strongest first: each label takes the first candidate
    spot whose measured text box stays inside the axes and above the null band, clears every
    label and marker already there, and whose leader line crosses no placed label. A label
    that fits nowhere is left off. Returns (placed, left off) so the panel can say so."""
    offsets = [(5, 0), (-5, 0), (5, 7), (5, -7), (-5, 7), (-5, -7), (0, 9),
               (18, 12), (18, -12), (-18, 12), (-18, -12), (30, 0), (-30, 0),
               (28, 22), (-28, 22), (28, -22), (-28, -22), (45, 10), (-45, 10),
               (45, -10), (-45, -10), (0, 26), (0, -22)]
    offsets += [(sx * dx, dy) for dx in (60, 90) for dy in (18, -18, 36, -36, 54, 0)
                for sx in (1, -1)]
    ax_box = ax.get_window_extent(renderer)
    placed, leaders, dropped = [], [], 0
    px_per_pt = fig.dpi / 72.0
    for x, y, text in cands:
        ok = False
        for dx, dy in offsets:
            ha = "left" if dx > 0 else "right" if dx < 0 else "center"
            va = "center" if dy == 0 else "bottom" if dy > 0 else "top"
            far = (dx, dy) not in ((5, 0), (-5, 0))
            t = ax.annotate(
                text, (x, y), xytext=(dx, dy), textcoords="offset points", ha=ha, va=va,
                fontsize=7, color=INK, zorder=6,
                path_effects=[pe.withStroke(linewidth=2.5, foreground=SURFACE)])
            bb = t.get_window_extent(renderer).padded(1.5)
            # inside the axes and clear of the null band under the line
            inside = (bb.x0 >= ax_box.x0 and bb.x1 <= ax_box.x1 and bb.y0 >= floor_px
                      and bb.y1 <= ax_box.y1)
            px0, py0 = ax.transData.transform((x, y))
            seg = Path([(px0, py0), (px0 + dx * px_per_pt * 0.92, py0 + dy * px_per_pt * 0.92)])
            clear = (inside and not any(bb.overlaps(p) for p in placed)
                     and not any(bb.overlaps(m) for m in avoid_px)
                     and not (far and any(seg.intersects_bbox(p) for p in placed))
                     and not any(lead.intersects_bbox(bb) for lead in leaders))
            if clear:
                placed.append(bb)
                if far:  # leader from the marker to the label's anchor
                    leaders.append(seg)
                    ax.annotate("", (x, y), xytext=(dx * 0.92, dy * 0.92),
                                textcoords="offset points", zorder=5,
                                arrowprops=dict(arrowstyle="-", color=INK_2, lw=0.55,
                                                alpha=0.75, shrinkA=0, shrinkB=3.5))
                ok = True
                break
            t.remove()
        dropped += not ok
    return len(placed), dropped


def _style(ax) -> None:
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(MUTED)
    ax.tick_params(colors=INK_2, labelsize=8)
    ax.grid(axis="y", color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)


def manhattan(result: Any, out_png: str, *, stats: dict | None = None) -> str:
    """One stacked panel per dataset: outcomes read on x, grouped by body-system category;
    -log10 p on y, log-compressed above COMPRESS_KNEE when the top hit is far above it.
    Every outcome at or above the dataset's phenome-wide line is labeled by measured
    placement; one that fits nowhere is left off and the panel subtitle says how many.
    `stats`, when given, is filled per dataset with n_above_line / n_labeled / n_unlabeled."""
    ids = list(result.dataset_ids)
    panel_h = 3.7
    fig = plt.figure(figsize=(12, panel_h * len(ids) + 0.9), facecolor=SURFACE)
    top_pad = 0.9 / (panel_h * len(ids) + 0.9)
    gs = fig.add_gridspec(len(ids), 1, left=0.06, right=0.985, top=1 - top_pad - 0.02,
                          bottom=0.09 / len(ids), hspace=0.62)
    # An explicit raster canvas to measure label text on: a figure's own canvas need not
    # have a renderer (matplotlib 3.11 detaches it on close).
    canvas = FigureCanvasAgg(fig)
    renderer = canvas.get_renderer()
    for i, dsid in enumerate(ids):
        ax = fig.add_subplot(gs[i, 0])
        _style(ax)
        ax.tick_params(length=3)
        spec = get_dataset(dsid)
        rows = result.dataset_rows(dsid)
        drawn = [r for r in rows if r.status in S.TESTED_STATUSES and r.neglog10_p is not None]
        cat = {id(r): outcome_category(r) for r in drawn}
        drawn.sort(key=lambda r: (CATEGORIES.index(cat[id(r)]), r.outcome_label.lower()))
        # each category gets at least a minimum band, plus a gap, so small ones stay legible
        min_w, gap = max(4.0, 0.03 * len(drawn)), max(2.0, 0.006 * len(drawn))
        x, pos = {}, 0.0
        for c in CATEGORIES:
            members = [r for r in drawn if cat[id(r)] == c]
            if not members:
                continue
            w = max(len(members), min_w)
            start = pos + (w - len(members)) / 2
            for j, r in enumerate(members):
                x[id(r)] = start + j
            pos += w + gap
        span = max(pos - gap, 1.0)
        n_bonf = S.phenome_wide_n(rows)
        line = -math.log10(0.05 / n_bonf) if n_bonf else float("inf")
        sig = [r for r in drawn if r.status == S.ST_OK and r.neglog10_p >= line]
        sig_ids = {id(r) for r in sig}
        ymax = max([r.neglog10_p for r in drawn] + [line if n_bonf else 0, 8.0])
        knee = COMPRESS_KNEE if ymax > 2 * COMPRESS_KNEE else None
        if knee:
            _compressed_scale(ax, knee)
        # category bands: alternating tones for the nulls, labels at band centers
        centers, labels = [], []
        used = [c for c in CATEGORIES if any(cat[id(r)] == c for r in drawn)]
        for k, c in enumerate(used):
            members = [r for r in drawn if cat[id(r)] == c]
            xs = [x[id(r)] for r in members]
            centers.append((min(xs) + max(xs)) / 2)
            labels.append(c)
            null = [r for r in members if id(r) not in sig_ids]
            lowc = [r for r in null if r.status == S.ST_LOW_CONFIDENCE]
            reg = [r for r in null if r.status != S.ST_LOW_CONFIDENCE]
            ax.scatter([x[id(r)] for r in reg], [r.neglog10_p for r in reg], s=9,
                       color=NULL_TONES[k % 2], linewidths=0, zorder=2, rasterized=True)
            ax.scatter([x[id(r)] for r in lowc], [r.neglog10_p for r in lowc], s=9,
                       facecolors="none", edgecolors=NULL_TONES[k % 2], linewidths=0.6,
                       zorder=2, rasterized=True)
        for r in sig:
            if r.beta is None or r.beta == 0 or r.direction_withheld_reason:
                ax.scatter(x[id(r)], r.neglog10_p, s=34, facecolors=SURFACE, edgecolors=INK_2,
                           linewidths=1.1, zorder=4)
            else:
                up = r.beta > 0
                ax.scatter(x[id(r)], r.neglog10_p, marker="^" if up else "v", s=40,
                           color=RAISES if up else LOWERS, edgecolors=SURFACE,
                           linewidths=0.6, zorder=4)
        if n_bonf:
            ax.axhline(line, color=INK_2, linewidth=0.9, zorder=3)
        ax.axhline(-math.log10(GENOME_WIDE), color=MUTED, linewidth=0.9, linestyle=(0, (4, 3)),
                   zorder=3)
        ax.set_xlim(-gap, span + gap)
        top = (knee * (1 + math.log10(ymax * 1.25 / knee)) if knee else ymax * 1.12)
        ax.set_ylim(0, knee * 10 ** (top / knee - 1) if knee else top)
        ax.set_yticks(_nice_ticks(ymax * 1.2, knee))
        ax.set_yticklabels([f"{t:g}" for t in _nice_ticks(ymax * 1.2, knee)])
        if knee:  # the compressed band is shaded; the subtitle says what the shading means
            ax.axhspan(knee, ax.get_ylim()[1], color="#f3f2ee", zorder=0, linewidth=0)
        ax.set_xticks(centers)
        ax.set_xticklabels(labels, rotation=35, ha="right", rotation_mode="anchor",
                           fontsize=7.5, color=INK_2)
        ax.tick_params(axis="x", length=0, pad=4)
        ax.set_ylabel("−log₁₀ p", fontsize=9, color=INK_2)
        # labels: every outcome above the line, measured placement
        canvas.draw()
        disp = ax.transData
        avoid = []
        for r in sig:
            px, py = disp.transform((x[id(r)], r.neglog10_p))
            avoid.append(Bbox.from_extents(px - 4, py - 4, px + 4, py + 4))
        ordered = sorted(sig, key=lambda r: -r.neglog10_p)
        texts = _distinct_short([S.display_label(r.outcome_label) for r in ordered])
        cands = [(x[id(r)], r.neglog10_p, t) for r, t in zip(ordered, texts)]
        floor_px = disp.transform((0, line if n_bonf else 0))[1] - 6
        n_lab, n_drop = _place_labels(fig, ax, cands, avoid, renderer, floor_px)
        if stats is not None:
            stats[dsid] = {"n_above_line": len(sig), "n_labeled": n_lab,
                           "n_unlabeled": n_drop}
        tested = sum(1 for r in rows if r.status in S.TESTED_STATUSES)
        ax.set_title(spec.label, fontsize=10.5, color=INK, loc="left", fontweight="bold",
                     pad=16)
        n_not_looked = sum(1 for r in rows if r.status == S.ST_NOT_QUERIED)
        if n_not_looked:  # a narrowed read says so first, and the line is not phenome-wide
            why = (" (those matched to the other biobank)"
                   if (getattr(result, "provenance", None) or {}).get(dsid, {}).get("scope")
                   == "harmonized" else "")
            sub = (f"{len(rows) - n_not_looked:,} of {len(rows):,} outcomes looked up"
                   f"{why}  ·  {tested:,} have a value  ·  "
                   f"{len(sig)} above the line p < 0.05 / {n_bonf:,}")
        else:
            sub = (f"{tested:,} of {len(rows):,} outcomes have a value  ·  {len(sig)} above "
                   f"the phenome-wide line p < 0.05 / {n_bonf:,}")
        if n_drop:
            sub += f"  ·  {n_drop} not labeled for space (see table)"
        if knee:
            sub += f"  ·  shaded: y axis compressed above {knee:g}"
        ax.text(0, 1.02, sub, transform=ax.transAxes, fontsize=8, color=INK_2, va="bottom")
    v = result.variant
    fig.text(0.06, 1 - 0.25 * top_pad, f"{v.rs_id or v.variant_id_grch38}  ·  "
             f"{v.variant_id_grch38}  ·  effect allele {v.grch38.alt}", fontsize=12, color=INK,
             fontweight="bold", va="center")
    handles = [
        Line2D([], [], marker="^", ls="", color=RAISES, ms=7, label="effect allele raises"),
        Line2D([], [], marker="v", ls="", color=LOWERS, ms=7, label="effect allele lowers"),
        Line2D([], [], marker="o", ls="", mfc=SURFACE, mec=INK_2, ms=6, label="no direction (withheld or beta 0)"),
        Line2D([], [], marker="o", ls="", color=NULL_TONES[0], ms=4, label="below the line"),
        Line2D([], [], marker="o", ls="", mfc="none", mec=NULL_TONES[0], ms=4,
               label="low confidence"),
        Line2D([], [], color=INK_2, lw=0.9, label="significance line (Bonferroni)"),
        Line2D([], [], color=MUTED, lw=0.9, ls=(0, (4, 3)), label="5e-8 (reference)"),
    ]
    fig.legend(handles=handles, loc="center left", bbox_to_anchor=(0.055, 1 - 0.62 * top_pad),
               ncol=7, fontsize=7.5, frameon=False, handletextpad=0.4, columnspacing=1.4,
               labelcolor=INK_2)
    fig.savefig(out_png, dpi=150, facecolor=SURFACE, bbox_inches="tight", pad_inches=0.15)
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
            if c.status == S.ST_API_CENSORED:
                ax.add_patch(Rectangle((j, y), 0.96, 0.9, facecolor=NULL_FILL, linewidth=0))
                ax.text(j + 0.48, y + 0.45, "p>=0.05", ha="center", va="center", fontsize=6,
                        color=INK_2)
                continue
            if c.z is None:  # a value whose direction is withheld is not a tested null
                ax.add_patch(Rectangle((j, y), 0.96, 0.9, facecolor=SURFACE, edgecolor=INK_2,
                                       linewidth=0.8))
                ax.text(j + 0.48, y + 0.45, f"p={c.p_value:.2g}, no direction", ha="center",
                        va="center", fontsize=6, color=INK)
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
