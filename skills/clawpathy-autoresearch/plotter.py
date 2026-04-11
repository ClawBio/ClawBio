"""Karpathy-style progress plot for clawpathy-autoresearch.

Generates a scatter + step-line chart showing experiment iterations,
with kept improvements highlighted and annotated.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


@dataclass
class ExperimentRecord:
    """A single experiment iteration."""

    experiment: int
    score: float
    kept: bool
    label: str


def plot_progress(
    history: list[ExperimentRecord],
    output_path: Path,
    title: str | None = None,
) -> Path:
    """Generate a Karpathy-style progress plot.

    Matches the aesthetic of github.com/karpathy/autoresearch:
    - Grey dots for discarded experiments (scattered widely)
    - Green dots + step line for kept improvements
    - Rotated italic labels on kept points
    - Clean, minimal aesthetic
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    kept = [r for r in history if r.kept]
    discarded = [r for r in history if not r.kept]

    n_total = len(history)
    n_kept = len(kept)

    if title is None:
        title = f"Autoresearch Progress: {n_total} Experiments, {n_kept} Kept Improvements"

    fig, ax = plt.subplots(figsize=(18, 9))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#cccccc")
    ax.spines["bottom"].set_color("#cccccc")
    ax.tick_params(colors="#888888", labelsize=11)

    # Discarded: visible grey dots showing the failed attempts
    if discarded:
        ax.scatter(
            [r.experiment for r in discarded],
            [r.score for r in discarded],
            c="#cccccc",
            s=35,
            alpha=0.5,
            zorder=2,
            label="Discarded",
        )

    if kept:
        kept_x = [r.experiment for r in kept]
        kept_y = [r.score for r in kept]
        max_x = max(r.experiment for r in history)

        # Running best step line
        running_best_x = []
        running_best_y = []
        current_best = kept_y[0]
        for i, (x, y) in enumerate(zip(kept_x, kept_y)):
            if i > 0:
                running_best_x.append(x)
                running_best_y.append(current_best)
            running_best_x.append(x)
            running_best_y.append(y)
            current_best = y

        if kept_x[-1] < max_x:
            running_best_x.append(max_x)
            running_best_y.append(current_best)

        ax.plot(
            running_best_x,
            running_best_y,
            c="#2ecc71",
            linewidth=2.5,
            zorder=3,
            solid_capstyle="round",
            label="Running best",
        )
        ax.scatter(
            kept_x,
            kept_y,
            c="#2ecc71",
            s=80,
            zorder=4,
            edgecolors="white",
            linewidths=1.5,
            label="Kept",
        )

        # Labels on kept points with connector lines.
        # Skip near-duplicate scores (within 1% of previous kept) to reduce
        # clutter, and spread labels vertically to avoid overlap.
        labelled_scores: list[float] = []
        label_items: list[tuple[int, ExperimentRecord]] = []
        for idx, r in enumerate(kept):
            # Skip if score is within 1% of an already-labelled score
            if any(abs(r.score - s) / max(s, 0.001) < 0.01 for s in labelled_scores):
                continue
            labelled_scores.append(r.score)
            label_items.append((idx, r))

        # Place labels with greedy vertical de-overlap.
        # Compute desired y in data coords, then push apart any that
        # are too close (within min_gap in data-space).
        y_range = ax.get_ylim()
        data_height = y_range[1] - y_range[0]
        min_gap = data_height * 0.045  # minimum vertical gap between labels

        # Build list of (desired_y, experiment_x, label_text)
        raw_placements = []
        for seq, (idx, r) in enumerate(label_items):
            label_text = r.label if len(r.label) <= 42 else r.label[:39] + "..."
            # Start label slightly above the point
            desired_y = r.score + min_gap * 0.5
            raw_placements.append((desired_y, r.experiment, r.score, label_text))

        # Sort by desired_y so we can push overlaps upward
        raw_placements.sort(key=lambda t: t[0])

        # Greedy push: ensure each label is at least min_gap above the previous
        placed_y = []
        for desired_y, x_data, y_data, text in raw_placements:
            y = desired_y
            if placed_y and y < placed_y[-1] + min_gap:
                y = placed_y[-1] + min_gap
            placed_y.append(y)

        # Draw annotations
        for i, (_, x_data, y_data, text) in enumerate(raw_placements):
            label_y = placed_y[i]
            ax.annotate(
                text,
                (x_data, y_data),
                xytext=(x_data + 0.6, label_y),
                fontsize=10,
                fontweight="medium",
                color="#2b8a4e",
                alpha=0.9,
                ha="left",
                va="center",
                arrowprops=dict(
                    arrowstyle="-",
                    color="#bbddbb",
                    lw=0.7,
                    shrinkA=0,
                    shrinkB=3,
                ),
            )

        # Extend y-axis top to make room for pushed-up labels
        if placed_y:
            top_label = placed_y[-1] + min_gap
            current_top = ax.get_ylim()[1]
            if top_label > current_top:
                ax.set_ylim(ax.get_ylim()[0], top_label)

    ax.set_xlabel("Experiment #", fontsize=13, color="#555555", labelpad=10)
    ax.set_ylabel("Mean Reproduction Error (lower is better)", fontsize=13, color="#555555", labelpad=10)
    ax.set_title(title, fontsize=15, fontweight="bold", color="#333333", pad=20)

    legend = ax.legend(
        loc="upper right",
        framealpha=0.9,
        edgecolor="#eeeeee",
        fontsize=10,
        labelcolor="#666666",
    )
    legend.get_frame().set_linewidth(0.5)

    ax.grid(True, alpha=0.1, color="#cccccc", linewidth=0.5)
    ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))

    plt.tight_layout()
    fig.savefig(output_path, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    return output_path
