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
    - Very faint grey dots for discarded experiments
    - Green dots + step line for kept improvements
    - Rotated italic labels on kept points
    - Clean, minimal aesthetic with lots of white space
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    kept = [r for r in history if r.kept]
    discarded = [r for r in history if not r.kept]

    n_total = len(history)
    n_kept = len(kept)

    if title is None:
        title = f"Autoresearch Progress: {n_total} Experiments, {n_kept} Kept Improvements"

    fig, ax = plt.subplots(figsize=(16, 8))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    # Remove top and right spines for clean look
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#cccccc")
    ax.spines["bottom"].set_color("#cccccc")
    ax.tick_params(colors="#888888", labelsize=11)

    # Discarded: very faint dots (nearly invisible like Karpathy's)
    if discarded:
        ax.scatter(
            [r.experiment for r in discarded],
            [r.score for r in discarded],
            c="#d4d4d4",
            s=25,
            alpha=0.35,
            zorder=2,
            label="Discarded",
        )

    # Kept: green dots + step line
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

        # Extend step line to the right edge
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

        # Rotated italic labels like Karpathy's
        # Alternate offset direction to reduce overlap
        for idx, r in enumerate(kept):
            y_offset = -10 if idx % 2 == 0 else 12
            ax.annotate(
                r.label,
                (r.experiment, r.score),
                textcoords="offset points",
                xytext=(10, y_offset),
                fontsize=7.5,
                fontstyle="italic",
                color="#5a9e6f",
                alpha=0.85,
                rotation=30,
                rotation_mode="anchor",
                ha="left",
                va="top" if y_offset < 0 else "bottom",
            )

    ax.set_xlabel("Experiment #", fontsize=13, color="#555555", labelpad=10)
    ax.set_ylabel("Mean Reproduction Error (lower is better)", fontsize=13, color="#555555", labelpad=10)
    ax.set_title(title, fontsize=15, fontweight="bold", color="#333333", pad=20)

    # Legend: top-right, minimal frame
    legend = ax.legend(
        loc="upper right",
        framealpha=0.9,
        edgecolor="#eeeeee",
        fontsize=10,
        labelcolor="#666666",
    )
    legend.get_frame().set_linewidth(0.5)

    # Very subtle grid
    ax.grid(True, alpha=0.1, color="#cccccc", linewidth=0.5)
    ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))

    plt.tight_layout()
    fig.savefig(output_path, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    return output_path
