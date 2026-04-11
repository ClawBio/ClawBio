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

    - Grey dots: discarded experiments
    - Green dots + step line: kept improvements (running best)
    - Labels on kept points showing what changed
    - Title: "N Experiments, M Kept Improvements"
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    kept = [r for r in history if r.kept]
    discarded = [r for r in history if not r.kept]

    n_total = len(history)
    n_kept = len(kept)

    if title is None:
        title = f"Autoresearch Progress: {n_total} Experiments, {n_kept} Kept Improvements"

    fig, ax = plt.subplots(figsize=(14, 7))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    # Discarded: grey dots
    if discarded:
        ax.scatter(
            [r.experiment for r in discarded],
            [r.score for r in discarded],
            c="#cccccc",
            s=30,
            alpha=0.5,
            zorder=2,
            label="Discarded",
        )

    # Kept: green dots + step line
    if kept:
        kept_x = [r.experiment for r in kept]
        kept_y = [r.score for r in kept]

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

        ax.plot(
            running_best_x,
            running_best_y,
            c="#2ecc71",
            linewidth=2,
            zorder=3,
            label="Running best",
        )
        ax.scatter(
            kept_x,
            kept_y,
            c="#2ecc71",
            s=60,
            zorder=4,
            edgecolors="white",
            linewidths=0.5,
            label="Kept",
        )

        # Annotate kept points
        for r in kept:
            ax.annotate(
                r.label,
                (r.experiment, r.score),
                textcoords="offset points",
                xytext=(8, -12),
                fontsize=7,
                color="#666666",
                alpha=0.8,
            )

    ax.set_xlabel("Experiment #", fontsize=12)
    ax.set_ylabel("Reproduction Score (higher is better)", fontsize=12)
    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.legend(loc="upper left", framealpha=0.9)
    ax.grid(True, alpha=0.2)
    ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))

    plt.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)

    return output_path
