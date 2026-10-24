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
from functools import lru_cache
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
from snp_phewas_catalog import UnknownDatasetError, get_dataset  # noqa: E402
from snp_phewas_replication import (  # noqa: E402
    NO_COMPARISON,
    T_ALL_THREE,
    T_DIRECTION_UNKNOWN,
    T_DISCORDANT,
    T_NULL,
    T_ONE,
    T_TWO,
    THREE_WAY_LABELS,
    THREE_WAY_ORDER,
    V_CONCORDANT_ONE,
    V_DIRECTION_UNKNOWN,
    V_DISCORDANT,
    V_LOW_CONFIDENCE,
    V_NO_COUNTERPART,
    V_NOT_READ,
    V_NOT_TESTED,
    V_REPLICATES,
    V_SINGLE,
    VERDICT_LABELS,
    VERDICT_ORDER,
    ReplicationTable,
    _missing_kind,
    cell_signal,
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


@lru_cache(maxsize=1)
def _pan_phecode_categories() -> tuple[tuple[float, str, str], ...]:
    """(numeric code, normalized code, category) of every Pan-UKBB phecode outcome, sorted:
    the body-system grouping an MVP phecode borrows (MVP's GWAS Catalog records carry no
    category of their own; the two codings are the same phecode 1.2 family)."""
    import csv
    import gzip

    from snp_phewas_crosswalk import normalize_phecode

    spec = get_dataset("pan_ukbb_eur")
    out: dict[str, str] = {}
    with gzip.open(spec.outcome_catalog_path, "rt", encoding="utf-8", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["trait_type"] == "phecode":
                out.setdefault(normalize_phecode(r["phenocode"]),
                               _PHECODE.get((r["category"] or "").strip().lower(), "other"))
    return tuple(sorted((float(c), c, cat) for c, cat in out.items()))


def phecode_category(code: str) -> str:
    """The category of a phecode by the SAME code's Pan-UKBB grouping. A code Pan-UKBB lacks
    takes the category its numeric neighbors share (phecode chapters are contiguous code
    ranges), else the category of a Pan-UKBB code with the same integer part, else other."""
    from bisect import bisect_left

    from snp_phewas_crosswalk import normalize_phecode

    table = _pan_phecode_categories()
    c = normalize_phecode(code)
    x = float(c)
    i = bisect_left(table, (x,))
    if i < len(table) and table[i][1] == c:
        return table[i][2]
    lo = table[i - 1][2] if i > 0 else None
    hi = table[i][2] if i < len(table) else None
    if lo is not None and lo == hi:
        return lo
    whole = c.split(".", 1)[0]
    same = [cat for _x, k, cat in table if k.split(".", 1)[0] == whole]
    return same[0] if same else "other"


def outcome_category(row: S.OutcomeRow) -> str:
    """The body-system category a row is drawn under, from the source's own grouping."""
    fam = get_dataset(row.dataset_id).family
    native = (row.outcome_category_native or "").strip()
    if fam == "mvp":
        if native.startswith("PheCode "):
            return phecode_category(native.split(" ", 1)[1])
        return "measurements" if row.outcome_kind == "continuous" else "other"
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
                # a patch, not only a halo: a halo covers the glyphs and leaves the gaps
                # between words open, so a line under a label read as "I21.Acute.myocardial"
                bbox=dict(boxstyle="square,pad=0.12", fc=SURFACE, ec="none", alpha=0.85),
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
        n_unreadable = sum(1 for r in rows if r.status == S.ST_NOT_READABLE)
        n_read = len(rows) - n_unreadable  # the outcomes this release can read at all
        rd = "readable " if n_unreadable else ""
        if n_not_looked:  # a narrowed read says so first, and the line is not phenome-wide
            have = "" if tested == n_read - n_not_looked else f"  ·  {tested:,} have a value"
            # why it was narrowed, only when that is the reason (a requested subset is not)
            why = (" (matched to another biobank)"
                   if (getattr(result, "provenance", None) or {}).get(dsid, {}).get("scope")
                   == "harmonized" else "")
            sub = (f"{n_read - n_not_looked:,} of {n_read:,} {rd}outcomes looked up{why}"
                   f"{have}  ·  {len(sig)} above the line p < 0.05 / {n_bonf:,}")
        else:
            sub = (f"{tested:,} of {n_read:,} {rd}outcomes have a value  ·  {len(sig)} "
                   f"above the phenome-wide line p < 0.05 / {n_bonf:,}")
        if n_unreadable:
            sub += f"  ·  {n_unreadable:,} with no readable file"
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


THREE_WAY_COLORS = {
    T_ALL_THREE: SERIES[0], T_TWO: "#7fb0ea", T_DISCORDANT: SERIES[7],
    T_DIRECTION_UNKNOWN: SERIES[6], T_ONE: SERIES[1], T_NULL: NULL_POINT,
    V_LOW_CONFIDENCE: UNTESTED_FILL, V_NOT_READ: UNTESTED_FILL, V_NOT_TESTED: UNTESTED_FILL,
}


def _fmt_p(nlp: float | None) -> str:
    """p from -log10 p, readable at any size (10 ** -586 is 0.0 as a float)."""
    if nlp is None:
        return "n/a"
    if nlp < 3:
        return f"{10 ** -nlp:.2g}"
    exp = math.ceil(nlp)
    mant = 10 ** (exp - nlp)
    if mant >= 9.95:
        mant, exp = 1.0, exp - 1
    return f"{mant:.2g}e-{exp}"


def _short_name(dsid: str) -> str:
    try:
        return get_dataset(dsid).short_label
    except UnknownDatasetError:  # a table built outside the catalog (tests, a caller's own)
        return dsid


def _family(dsid: str) -> str:
    try:
        return get_dataset(dsid).family
    except UnknownDatasetError:
        return dsid


def _long_name(dsid: str) -> str:
    try:
        return get_dataset(dsid).label
    except UnknownDatasetError:
        return dsid


def select_matrix_rows(table: ReplicationTable, max_rows: int = 40, hard_cap: int = 80
                       ) -> tuple[list, int, int]:
    """(rows to draw, n significant, n null-everywhere not drawn).

    Every outcome significant in at least one dataset (at the threshold of a pair it is in)
    comes first, strongest first, up to `hard_cap`; the cap is then filled with the
    strongest of the rest. A figure of mostly tested-in-both-significant-in-neither rows is
    noise: on rs12740374 (two biobanks, shared-outcome read) 31 of the 40 rows drawn by the
    old top-40-by-p rule were null in both."""
    sig = [r for r in table.rows if r.significant_anywhere(table.alpha)]
    rest = [r for r in table.rows if not r.significant_anywhere(table.alpha)]
    drawn = sig[:hard_cap]
    drawn += rest[:max(0, max_rows - len(drawn))]
    drawn_ids = {id(r) for r in drawn}

    def null_everywhere(r) -> bool:
        tested = [c for c in r.cells.values() if c is not None]
        return all(cell_signal(c, 1.0).tested for c in tested)

    n_null_hidden = sum(1 for r in rest if id(r) not in drawn_ids and null_everywhere(r))
    return drawn, len(sig), n_null_hidden


def pair_strips(table: ReplicationTable) -> list[tuple[str, Any]]:
    """(label, key) per verdict strip of the matrix: one per pair, then "all three".

    One line per pair, named by biobank ("Pan-UKBB × MVP"): a label split after the "×" left
    "Pan-UKBB EUR ×" alone on a line, read as the title of the column under it. The ancestry
    is in the dataset headings; it is kept only where two datasets share a biobank (the
    catalog's `family`, the same test `replication_lines` uses)."""
    ids = table.dataset_ids
    fam = {d: _short_name(d).split(" ")[0] for d in ids}
    if len({_family(d) for d in ids}) < len(ids):  # two datasets of one biobank
        fam = {d: _short_name(d) for d in ids}
    strips: list[tuple[str, Any]] = [(f"{fam[a]} × {fam[b]}", p) for p in table.pairs
                                     for a, b in [p]]
    if table.n_shared_3 is not None:
        strips.append(("all three", "three_way"))
    return strips


MATRIX_DPI = 150


def replication_matrix(table: ReplicationTable, out_png: str, *, max_rows: int = 40) -> str:
    """Harmonized outcome x dataset, cell = signed z (clipped at |z| = 8) where the outcome
    is significant in that dataset, plain light gray where it was tested and is not, hatched
    where it was not tested or not read, outlined where its direction is withheld, blank
    where the dataset has no such outcome. One narrow verdict strip per dataset pair, and a
    three-way strip with three datasets; a legend spells the verdicts out once."""
    ids = table.dataset_ids
    pairs = table.pairs
    rows, n_sig, n_null_hidden = select_matrix_rows(table, max_rows)
    strips = pair_strips(table)
    n_cols = len(ids)
    strip_x0 = n_cols + 0.25
    strip_w, strip_gap = 0.34, 0.12
    strips_end = strip_x0 + len(strips) * (strip_w + strip_gap)
    h = max(3.2, 0.26 * len(rows) + 2.6)
    width = 3.3 + 1.25 * n_cols + 0.5 * len(strips) + 0.9
    # created at the save resolution: `_pair_heading` measures text, and text does not scale
    # exactly with dpi, so a measurement at 100 dpi let labels cross the bracket at 150
    fig, ax = plt.subplots(figsize=(width, h), facecolor=SURFACE, dpi=MATRIX_DPI)
    ax.set_facecolor(SURFACE)
    used_pair, used_three = set(), set()

    def sig_in(ds: str, r) -> bool:
        return any(ds in p and r.verdicts.get(p) != V_NO_COUNTERPART
                   and cell_signal(r.cells.get(ds), table.alpha[p]).significant for p in pairs)

    row_labels = _distinct_short([S.display_label(r.harmonized_label) for r in rows], 46)
    for i, r in enumerate(rows):
        y = len(rows) - 1 - i
        for j, ds in enumerate(ids):
            c = r.cells.get(ds)
            if c is None:  # the dataset has no such outcome: blank, never hatched
                continue
            if c.status not in (S.ST_OK, S.ST_API_CENSORED):
                ax.add_patch(Rectangle((j, y), 0.96, 0.9, facecolor=UNTESTED_FILL,
                                       hatch="////", edgecolor="#ffffff", linewidth=0))
                ax.text(j + 0.48, y + 0.45, _MISSING_CELL_TEXT[_missing_kind(c)], ha="center",
                        va="center", fontsize=6, color=INK_2)
                continue
            if c.status == S.ST_API_CENSORED:
                ax.add_patch(Rectangle((j, y), 0.96, 0.9, facecolor=NULL_FILL, linewidth=0))
                ax.text(j + 0.48, y + 0.45, "p≥0.05", ha="center", va="center", fontsize=6,
                        color=INK_2)
                continue
            if c.z is None or c.direction_withheld_reason:  # withheld is not a tested null
                ax.add_patch(Rectangle((j + 0.02, y + 0.02), 0.92, 0.86, facecolor=SURFACE,
                                       edgecolor=INK_2, linewidth=0.8))
                ax.text(j + 0.48, y + 0.45, f"p={_fmt_p(c.neglog10_p)}\nno direction", ha="center",
                        va="center", fontsize=5.5, color=INK, linespacing=0.95)
                continue
            if not sig_in(ds, r):
                ax.add_patch(Rectangle((j, y), 0.96, 0.9, facecolor=NULL_FILL, linewidth=0))
                ax.text(j + 0.48, y + 0.45, f"p={_fmt_p(c.neglog10_p)}", ha="center",
                        va="center", fontsize=6, color=INK_2)
                continue
            z = max(-8.0, min(8.0, c.z))
            ax.add_patch(Rectangle((j, y), 0.96, 0.9, facecolor=DIVERGING((z + 8) / 16),
                                   linewidth=0))
            ax.text(j + 0.48, y + 0.45, f"p={_fmt_p(c.neglog10_p)}", ha="center",
                    va="center", fontsize=6, color=INK if abs(z) < 5 else "#ffffff")
        for k, (_name, key) in enumerate(strips):
            x = strip_x0 + k * (strip_w + strip_gap)
            if key == "three_way":
                v = r.three_way
                if v is None:
                    continue
                used_three.add(v)
                color, hatch = THREE_WAY_COLORS[v], "////" if v in NO_COMPARISON else None
            else:
                v = r.verdicts[key]
                if v == V_NO_COUNTERPART:
                    continue
                used_pair.add(v)
                color, hatch = VERDICT_COLORS[v], "////" if v in NO_COMPARISON else None
            ax.add_patch(Rectangle((x, y + 0.08), strip_w, 0.74, facecolor=color, hatch=hatch,
                                   edgecolor="#ffffff" if hatch else color, linewidth=0))
        ax.text(-0.06, y + 0.45, row_labels[i], ha="right", va="center",
                fontsize=6.5, color=INK)
    ax.set_xlim(-0.06, strips_end + 0.9)
    ax.set_ylim(-0.2, len(rows) + 5.0)  # room for the rotated strip headers and their heading
    labels = [ax.text(strip_x0 + k * (strip_w + strip_gap) + strip_w / 2, len(rows) + 0.15,
                      name, ha="center", va="bottom", fontsize=6, color=INK_2, rotation=90)
              for k, (name, _key) in enumerate(strips)]
    for j, d in enumerate(ids):
        ax.text(j + 0.48, len(rows) + 0.15, _long_name(d).replace(", ", "\n"),
                ha="center", va="bottom", fontsize=6.5, color=INK_2)
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    n_hidden = len(table.rows) - len(rows)
    title = (f"Replication across biobanks: {len(rows)} of {len(table.rows):,} shared "
             f"outcomes shown, the {min(n_sig, len(rows))} significant in at least one "
             f"biobank first")
    if n_hidden:
        title += (f"\n{n_hidden:,} not drawn ({n_null_hidden:,} of them tested and significant "
                  f"nowhere); all are in replication.tsv")
    thr = "; ".join(f"{_short_name(a)} × {_short_name(b)} 0.05/{table.n_shared[(a, b)]:,}"
                    for a, b in pairs if table.n_shared[(a, b)])
    if table.n_shared_3:
        # fewer outcomes are shared by all three than by any pair, so this threshold is the
        # loosest: an "all three" verdict can stand where a pair strip does not
        thr += (f"; all three 0.05/{table.n_shared_3:,} (fewer shared outcomes, so a looser "
                f"line than the pairs')")
    ax.set_title(title + ("\nsignificance per pair: " + thr if thr else ""), fontsize=8.5,
                 color=INK, loc="left", pad=6)
    # the colorbar sits beside the strips, inside the figure, not at the far right edge
    cax = ax.inset_axes([strips_end + 0.25, max(0.0, len(rows) - 6.5), 0.18,
                         min(6.0, max(1.5, len(rows) - 0.5))], transform=ax.transData)
    sm = plt.cm.ScalarMappable(cmap=DIVERGING, norm=plt.Normalize(-8, 8))
    cb = fig.colorbar(sm, cax=cax)
    cb.set_label("signed z, clipped at 8\n(significant cells only)", fontsize=6, color=INK_2)
    cb.ax.tick_params(labelsize=5.5)
    handles = [Rectangle((0, 0), 1, 1, facecolor=VERDICT_COLORS[v],
                         hatch="////" if v in NO_COMPARISON else None, edgecolor="#ffffff",
                         label=VERDICT_LABELS[v]) for v in VERDICT_ORDER if v in used_pair]
    handles += [Rectangle((0, 0), 1, 1, facecolor=THREE_WAY_COLORS[v],
                          hatch="////" if v in NO_COMPARISON else None, edgecolor="#ffffff",
                          label=f"all three: {THREE_WAY_LABELS[v]}")
                for v in THREE_WAY_ORDER if v in used_three]
    if handles:
        fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.02, 0.0), ncol=2,
                   fontsize=6, frameon=False, handlelength=1.2, columnspacing=1.2)
    fig.tight_layout()
    if labels:  # the strips are one group: a heading over a bracket says so
        _pair_heading(fig, ax, labels, len(rows) + 0.15, strip_x0, strips_end - strip_gap)
    fig.savefig(out_png, dpi=MATRIX_DPI, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)
    return out_png


def _pair_heading(fig, ax, labels: list, base: float, x0: float, x1: float) -> None:
    """A bracket and "replication by pair" over the verdict-strip labels, placed AFTER the
    final layout: the labels' height is fixed in points, so their height in data units depends
    on the axes' final size, and a position measured before `tight_layout` put the heading on
    the title. The y range is extended to hold label + gap + heading, then nothing
    rescales before the save."""
    renderer = FigureCanvasAgg(fig).get_renderer()
    height_px = ax.get_window_extent(renderer).height
    pt = fig.dpi / 72.0
    label_px = max(t.get_window_extent(renderer).height for t in labels)
    gap_px, head_px = 4 * pt, 6.5 * pt * 1.4
    y_lo, y_hi = ax.get_ylim()
    k = (label_px + gap_px + head_px + 2 * pt) / height_px
    if k < 0.9:  # the data top that leaves exactly that many pixels above `base`, and room
        # for the two-line dataset headings that sit at the same base
        y_hi = max(base + 1.0, (base - y_lo * k) / (1 - k))
        ax.set_ylim(y_lo, y_hi)
    per_px = (y_hi - y_lo) / height_px
    y_line = base + (label_px + gap_px) * per_px
    ax.plot([x0, x1], [y_line, y_line], color=INK_2, lw=0.6, solid_capstyle="butt")
    ax.text((x0 + x1) / 2, y_line + 2 * pt * per_px, "replication by pair", ha="center",
            va="bottom", fontsize=6.5, color=INK_2)


_SIGNIFICANT_VERDICTS = (V_REPLICATES, V_CONCORDANT_ONE, V_DISCORDANT, V_SINGLE,
                         V_DIRECTION_UNKNOWN)


def paired_effects(table: ReplicationTable, out_png: str) -> str:
    """Beta in one dataset against another with 95% CIs: one row of panels per dataset pair,
    one panel per outcome kind present (binary in log-OR, continuous in z units, since effect
    scales differ across cohorts). Each panel gives the effect correlation over every
    outcome tested in both AND over the outcomes significant in either, each with its n:
    over hundreds of shared outcomes the first is dominated by nulls. Two datasets give one
    row, as before."""
    if not table.pairs:
        raise ValueError("a paired effect plot needs two datasets")
    pairs = [p for p in table.pairs if table.n_shared[p]] or table.pairs[:1]
    pts_by: dict[tuple[str, str], list] = {}
    for a, b in pairs:
        pts_by[(a, b)] = [r for r in table.rows
                          if r.verdicts.get((a, b)) != V_NO_COUNTERPART
                          and all(r.cells.get(d) is not None and r.cells[d].status == S.ST_OK
                                  and r.cells[d].beta is not None for d in (a, b))]
    # One panel per (pair, kind) that has outcomes; an empty combination is not drawn (with
    # three biobanks only one pair has continuous outcomes, and empty panels filled half the
    # figure). Nothing at all still gives one panel that says so.
    panels = [((a, b), k) for a, b in pairs for k in ("binary", "continuous")
              if any(r.cells[a].outcome_kind == k for r in pts_by[(a, b)])]
    panels = panels or [(pairs[0], "binary")]
    ncol = 1 if len(panels) == 1 else 2
    nrow = math.ceil(len(panels) / ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(5.8 * ncol, 5 * nrow), facecolor=SURFACE,
                             squeeze=False)
    for ax in axes.flat[len(panels):]:
        ax.set_visible(False)
    for ax, ((a, b), kind) in zip(axes.flat, panels):
        _style(ax)
        pts = [r for r in pts_by[(a, b)] if r.cells[a].outcome_kind == kind]
        if not pts:  # say it rather than draw an empty grid
            ax.set_title(_paired_title(pts, a, b, kind), fontsize=8.5, color=INK, loc="left")
            ax.text(0.5, 0.5, "no outcome of this kind was tested in both", ha="center",
                    va="center", fontsize=8, color=INK_2, transform=ax.transAxes)
            ax.set_xticks([])
            ax.set_yticks([])
            ax.grid(False)
            continue

        def xy(r):
            if kind == "binary":
                return (r.cells[a].beta, r.cells[b].beta, 1.96 * r.cells[a].se,
                        1.96 * r.cells[b].se)
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
        ax.set_xlabel(f"{_long_name(a)}: {unit}", fontsize=8, color=INK_2)
        ax.set_ylabel(f"{_long_name(b)}: {unit}", fontsize=8, color=INK_2)
        ax.set_title(_paired_title(pts, a, b, kind), fontsize=8.5, color=INK, loc="left")
        if pts:
            ax.legend(fontsize=6.5, frameon=False, loc="upper left")
    fig.tight_layout()
    fig.savefig(out_png, dpi=150, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)
    return out_png


def paired_correlations(pts, a: str, b: str, kind: str) -> dict[str, Any]:
    """{n_all, r_all, n_sig, r_sig}: the effect correlation over every outcome tested in
    both, and over the outcomes significant in either (r is None below three points)."""
    sig = [r for r in pts if r.verdicts[(a, b)] in _SIGNIFICANT_VERDICTS]
    return {"n_all": len(pts), "r_all": _corr(pts, a, b, kind) if len(pts) >= 3 else None,
            "n_sig": len(sig), "r_sig": _corr(sig, a, b, kind) if len(sig) >= 3 else None}


def _paired_title(pts, a, b, kind) -> str:
    c = paired_correlations(pts, a, b, kind)
    head = (f"{_short_name(a)} × {_short_name(b)}, "
            f"{'case/control' if kind == 'binary' else 'continuous'} outcomes tested in both: "
            f"{c['n_all']}")
    parts = []
    if c["r_all"] is not None:
        parts.append(f"r = {c['r_all']:.2f} over all {c['n_all']}")
    if c["r_sig"] is not None:
        parts.append(f"r = {c['r_sig']:.2f} over the {c['n_sig']} significant in either")
    elif c["n_all"]:
        parts.append(f"{c['n_sig']} significant in either (too few for r)")
    return head + ("\n" + "; ".join(parts) if parts else "")


def _corr(pts, a, b, kind) -> float:
    xs = [(r.cells[a].beta if kind == "binary" else r.cells[a].z) for r in pts]
    ys = [(r.cells[b].beta if kind == "binary" else r.cells[b].z) for r in pts]
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    sx = math.sqrt(sum((x - mx) ** 2 for x in xs))
    sy = math.sqrt(sum((y - my) ** 2 for y in ys))
    return sxy / (sx * sy) if sx and sy else float("nan")
