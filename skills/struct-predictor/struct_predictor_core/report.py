"""
report.py — Markdown report, matplotlib figures, result.json, reproducibility bundle.
"""
from __future__ import annotations

import importlib.metadata
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


DISCLAIMER = (
    "ClawBio is a research and educational tool. It is not a medical device "
    "and does not provide clinical diagnoses. Consult a healthcare "
    "professional before making any medical decisions."
)

# pLDDT confidence bands (AlphaFold2 / Boltz convention)
BANDS = [
    ("Very high", 90.0, 100.0, "#0053D6"),
    ("High",      70.0,  90.0, "#65CBF3"),
    ("Low",       50.0,  70.0, "#FFDB13"),
    ("Very low",   0.0,  50.0, "#FF7D45"),
]


def generate_report(
    output_dir: Path,
    sequences_info: list[dict],
    plddt: np.ndarray,
    pae: np.ndarray,
    chain_boundaries: list[dict],
    cif_path: Path,
    cmd: str,
    input_label: str,
    demo: bool,
) -> None:
    """Write all output artefacts to output_dir."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "figures").mkdir(exist_ok=True)
    (output_dir / "reproducibility").mkdir(exist_ok=True)

    band_bd = _confidence_band_breakdown(plddt)
    mean_plddt = float(np.mean(plddt))
    min_plddt  = float(np.min(plddt))

    md = _build_markdown(
        sequences_info=sequences_info,
        plddt=plddt,
        chain_boundaries=chain_boundaries,
        band_bd=band_bd,
        mean_plddt=mean_plddt,
        min_plddt=min_plddt,
        input_label=input_label,
        demo=demo,
    )
    (output_dir / "report.md").write_text(md)

    result = {
        "skill": "struct-predictor",
        "input": input_label,
        "demo": demo,
        "n_residues": len(plddt),
        "n_chains": len(chain_boundaries),
        "mean_plddt": round(mean_plddt, 2),
        "min_plddt": round(min_plddt, 2),
        "band_breakdown": band_bd,
        "sequences": [
            {"name": s["name"], "length": len(s["sequence"])}
            for s in sequences_info
        ],
    }
    (output_dir / "result.json").write_text(json.dumps(result, indent=2))

    _plot_plddt(output_dir / "figures" / "plddt.png", plddt, chain_boundaries)
    _plot_pae(output_dir / "figures" / "pae.png", pae, chain_boundaries)

    from struct_predictor_core.viewer import generate_viewer_html
    generate_viewer_html(
        output_path=output_dir / "viewer.html",
        cif_path=cif_path,
        plddt=plddt,
        chain_boundaries=chain_boundaries,
    )

    _write_reproducibility(output_dir / "reproducibility", cmd)


def _confidence_band_breakdown(plddt: np.ndarray) -> dict[str, float]:
    """Return percentage of residues in each pLDDT confidence band.

    Uses a first-match strategy (highest band first) so every value,
    including exactly 100.0, is assigned to exactly one band and
    percentages always sum to 100.
    """
    n = len(plddt)
    bd = {}
    unassigned = np.ones(n, dtype=bool)
    for label, lo, _hi, _ in BANDS:
        mask = unassigned & (plddt >= lo)
        bd[label] = round(100.0 * int(np.sum(mask)) / n, 1) if n > 0 else 0.0
        unassigned &= ~mask
    return bd


def _build_markdown(
    sequences_info: list[dict],
    plddt: np.ndarray,
    chain_boundaries: list[dict],
    band_bd: dict,
    mean_plddt: float,
    min_plddt: float,
    input_label: str,
    demo: bool,
) -> str:
    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    n = len(plddt)

    demo_banner = "\n> **Demo output** — Trp-cage miniprotein (1L2Y). No input file was provided.\n" if demo else ""

    lines = [
        "# Struct Predictor Report",
        "",
        f"**Date**: {now}  ",
        f"**Skill**: struct-predictor  ",
        f"**Engine**: Boltz-2  ",
        f"**Input**: {input_label}  ",
        demo_banner,
        "---",
        "",
        "## Sequence Summary",
        "",
        "| Chain | Name | Length |",
        "|-------|------|--------|",
    ]
    for s in sequences_info:
        lines.append(f"| {s.get('chain_id', '?')} | {s['name']} | {len(s['sequence'])} |")

    lines += [
        "",
        "---",
        "",
        "## Confidence Summary",
        "",
        "| Metric | Value |",
        "|--------|-------|",
        f"| Residues | {n} |",
        f"| Mean pLDDT | {mean_plddt:.1f} |",
        f"| Min pLDDT | {min_plddt:.1f} |",
        "",
        "### pLDDT Band Breakdown",
        "",
        "| Band | Range | % Residues |",
        "|------|-------|-----------|",
    ]
    for label, lo, hi, _ in BANDS:
        lines.append(f"| {label} | {lo:.0f}–{hi:.0f} | {band_bd[label]:.1f}% |")

    lines += [
        "",
        "### Chain Boundaries",
        "",
        "| Chain | Start residue | End residue | Length |",
        "|-------|--------------|------------|--------|",
    ]
    for cb in chain_boundaries:
        length = cb["end"] - cb["start"] + 1
        lines.append(f"| {cb['chain_id']} | {cb['start']+1} | {cb['end']+1} | {length} |")

    lines += [
        "",
        "---",
        "",
        "## Figures",
        "",
        "**[Open 3D Viewer](viewer.html)** — interactive structure coloured by pLDDT (open in browser)",
        "",
        "![pLDDT per residue](figures/plddt.png)",
        "",
        "![PAE heatmap](figures/pae.png)",
        "",
        "---",
        "",
        "## Interpretation",
        "",
        "- **Very high pLDDT (≥ 90)**: backbone prediction is expected to be accurate within ~0.5 Å.",
        "- **High pLDDT (70–90)**: generally reliable, small side-chain errors possible.",
        "- **Low pLDDT (50–70)**: disordered or poorly predicted region; interpret with caution.",
        "- **Very low pLDDT (< 50)**: likely intrinsically disordered; structure should not be used.",
        "- **PAE**: low values (dark) indicate confident relative positioning; high values indicate uncertain domain/chain packing.",
        "",
        "---",
        "",
        "## Citations",
        "",
        "- Passaro S et al. (2025) *Boltz-2: Towards Accurate and Efficient Binding Affinity Prediction*. bioRxiv. doi:10.1101/2025.06.14.659707",
        "- Wohlwend J et al. (2024) *Boltz-1: Democratizing Biomolecular Interaction Modeling*. bioRxiv. doi:10.1101/2024.11.19.624167",
        "- Jumper J et al. (2021) *AlphaFold2 pLDDT definition*. Nature. doi:10.1038/s41586-021-03819-2",
        "",
        "---",
        "",
        f"> {DISCLAIMER}",
    ]
    return "\n".join(lines)


def _plot_plddt(out_path: Path, plddt: np.ndarray, chain_boundaries: list[dict]) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.patches as mpatches
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(max(6, len(plddt) / 10), 3))
    x = np.arange(1, len(plddt) + 1)

    for label, lo, hi, colour in BANDS:
        mask = (plddt >= lo) & (plddt < hi)
        if np.any(mask):
            ax.scatter(x[mask], plddt[mask], c=colour, s=8, label=label, zorder=3)
        ax.axhspan(lo, hi, alpha=0.08, color=colour)

    ax.plot(x, plddt, color="black", linewidth=0.8, zorder=2)

    for cb in chain_boundaries[1:]:
        ax.axvline(cb["start"] + 1, color="gray", linestyle="--", linewidth=0.8, alpha=0.6)

    ax.set_xlim(0.5, len(plddt) + 0.5)
    ax.set_ylim(0, 100)
    ax.set_xlabel("Residue")
    ax.set_ylabel("pLDDT score")
    ax.set_title("Per-residue pLDDT confidence")
    handles = [mpatches.Patch(color=c, label=lbl) for lbl, lo, hi, c in BANDS]
    ax.legend(handles=handles, fontsize=7, loc="lower right")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def _plot_pae(out_path: Path, pae: np.ndarray, chain_boundaries: list[dict]) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(5, 4))
    vmax = max(float(pae.max()), 1.0)

    im = ax.imshow(pae, cmap="Greens_r", vmin=0, vmax=vmax)
    fig.colorbar(im, ax=ax, label="PAE (Å)")
    ax.set_xticks([])
    ax.set_yticks([])

    for cb in chain_boundaries[1:]:
        ax.axhline(cb["start"] - 0.5, color="black", linewidth=1.0)
        ax.axvline(cb["start"] - 0.5, color="black", linewidth=1.0)

    ax.set_title("Predicted Aligned Error (PAE)")
    ax.set_xlabel("Scored residue")
    ax.set_ylabel("Aligned residue")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def _write_reproducibility(repro_dir: Path, cmd: str) -> None:
    repro_dir.mkdir(parents=True, exist_ok=True)
    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")

    (repro_dir / "commands.sh").write_text(
        f"#!/bin/bash\n# Generated by struct-predictor on {now}\n\n{cmd}\n"
    )

    try:
        version_str = importlib.metadata.version("boltz")
    except importlib.metadata.PackageNotFoundError:
        version_str = "boltz not installed"

    (repro_dir / "environment.txt").write_text(
        f"# Generated by struct-predictor on {now}\nboltz: {version_str}\n"
    )
