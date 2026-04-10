#!/usr/bin/env python3
"""Treatment effect analysis on organ-specific biological age using Filbin COVID data.

Runs the proteomics clock on Day 0 and Day 7, then compares delta-ages
across COVID severity groups. This demonstrates longitudinal use of the
proteomics clock skill.

Prerequisites:
    1. Download Filbin data from https://doi.org/10.17632/nf853r8xsj.2
    2. Preprocess: python examples/fetch_filbin.py --data-dir <dir> --output /tmp/filbin_prepared
    3. Run this script: python examples/treatment_effect_covid.py --data-dir /tmp/filbin_prepared --output /tmp/filbin_treatment
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

# Add skill directory to path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from proteomics_clock import run_analysis


ORGANS_OF_INTEREST = ["Heart", "Brain", "Kidney", "Lung", "Immune", "Organismal"]


def run_clock_on_timepoint(input_path: Path, output_dir: Path, organs: list[str]) -> pd.DataFrame:
    """Run proteomics clock on a single timepoint and return predictions."""
    result = run_analysis(
        output_dir=output_dir,
        organs=organs,
        generation="both",
        fold=1,
        convert_mortality_to_years=True,
        age_column=None,
        verbose=False,
        input_path=input_path,
    )
    gen1 = pd.read_csv(output_dir / "tables" / "predictions_gen1.csv")
    return gen1


def compute_delta_ages(
    preds_t0: pd.DataFrame,
    preds_t1: pd.DataFrame,
    meta: pd.DataFrame,
    severity_col: str,
) -> pd.DataFrame:
    """Compute delta (T1 - T0) biological age per subject, merge with severity."""
    # Extract subject_id from sample_id (remove _D suffix)
    preds_t0 = preds_t0.copy()
    preds_t1 = preds_t1.copy()
    preds_t0["subject_id"] = preds_t0["sample_id"].astype(str).str.replace(r"_D\d+", "", regex=True)
    preds_t1["subject_id"] = preds_t1["sample_id"].astype(str).str.replace(r"_D\d+", "", regex=True)

    organ_cols = [c for c in preds_t0.columns if c not in ("sample_id", "subject_id")]
    common_subjects = set(preds_t0["subject_id"]) & set(preds_t1["subject_id"])

    t0 = preds_t0[preds_t0["subject_id"].isin(common_subjects)].set_index("subject_id")[organ_cols]
    t1 = preds_t1[preds_t1["subject_id"].isin(common_subjects)].set_index("subject_id")[organ_cols]
    t0 = t0.loc[t0.index.sort_values()]
    t1 = t1.loc[t1.index.sort_values()]

    common_idx = t0.index.intersection(t1.index)
    delta = t1.loc[common_idx] - t0.loc[common_idx]
    delta = delta.reset_index()

    # Merge with severity
    meta = meta.copy()
    meta["subject_id"] = meta["subject_id"].astype(str)
    delta = delta.merge(meta[["subject_id", severity_col]].drop_duplicates(), on="subject_id", how="left")

    return delta


def plot_delta_boxplots(delta: pd.DataFrame, severity_col: str, organs: list[str], outpath: Path) -> None:
    """Boxplot of delta-age by severity group for each organ."""
    groups = sorted(delta[severity_col].dropna().unique())
    n_organs = len(organs)
    fig, axes = plt.subplots(1, n_organs, figsize=(4 * n_organs, 5), sharey=False)
    if n_organs == 1:
        axes = [axes]

    for ax, organ in zip(axes, organs):
        if organ not in delta.columns:
            ax.set_title(organ)
            ax.text(0.5, 0.5, "N/A", ha="center", va="center", transform=ax.transAxes)
            continue

        data_by_group = []
        labels = []
        for g in groups:
            vals = delta.loc[delta[severity_col] == g, organ].dropna()
            data_by_group.append(vals.values)
            labels.append(f"{g}\n(n={len(vals)})")

        ax.boxplot(data_by_group, tick_labels=labels)
        ax.axhline(0, color="gray", linestyle="--", alpha=0.5)
        ax.set_title(organ)
        ax.set_ylabel("Delta biological age (years)")

        # t-test if 2 groups
        if len(data_by_group) == 2 and len(data_by_group[0]) > 1 and len(data_by_group[1]) > 1:
            t_stat, p_val = stats.ttest_ind(data_by_group[0], data_by_group[1], equal_var=False)
            ax.set_xlabel(f"p={p_val:.3g}")

    fig.suptitle("Change in Organ Biological Age (Day 7 - Day 0) by COVID Severity", fontsize=12)
    plt.tight_layout()
    plt.savefig(outpath, dpi=180)
    plt.close()


def plot_delta_heatmap(delta: pd.DataFrame, severity_col: str, organs: list[str], outpath: Path) -> None:
    """Heatmap of mean delta-age per severity group x organ."""
    import seaborn as sns

    available_organs = [o for o in organs if o in delta.columns]
    groups = sorted(delta[severity_col].dropna().unique())
    means = delta.groupby(severity_col)[available_organs].mean().loc[groups]

    fig, ax = plt.subplots(figsize=(max(6, len(available_organs)), max(3, 0.5 * len(groups))))
    sns.heatmap(means, cmap="RdBu_r", center=0, annot=True, fmt=".2f", ax=ax)
    ax.set_title("Mean Delta Biological Age (Day 7 - Day 0)")
    ax.set_ylabel("COVID Severity")
    ax.set_xlabel("Organ Clock")
    plt.tight_layout()
    plt.savefig(outpath, dpi=180)
    plt.close()


def write_treatment_report(
    output_dir: Path,
    delta: pd.DataFrame,
    severity_col: str,
    organs: list[str],
    n_t0: int,
    n_t7: int,
    n_paired: int,
) -> None:
    """Write markdown report summarising treatment effect analysis."""
    groups = sorted(delta[severity_col].dropna().unique())
    available_organs = [o for o in organs if o in delta.columns]

    rows = []
    for organ in available_organs:
        for g in groups:
            vals = delta.loc[delta[severity_col] == g, organ].dropna()
            if len(vals) > 0:
                rows.append(f"| {organ} | {g} | {len(vals)} | {vals.mean():.2f} | {vals.std(ddof=1):.2f} |")

    summary_block = "\n".join(rows) if rows else "| - | - | 0 | - | - |"

    report = f"""# COVID-19 Organ Aging Treatment Effect Analysis

**Dataset**: Filbin et al. (2021) — Longitudinal COVID-19 Olink proteomics
**Timepoints**: Day 0 → Day 7
**Severity column**: `{severity_col}`
**Samples**: {n_t0} (Day 0), {n_t7} (Day 7), {n_paired} paired subjects

## Delta Biological Age by Severity Group

| Organ | Severity | N | Mean Delta (years) | Std |
|---|---|---:|---:|---:|
{summary_block}

## Figures

- Delta boxplots: `figures/delta_boxplots.png`
- Delta heatmap: `figures/delta_heatmap.png`

## Interpretation

Positive delta = biological age increased over 7 days (accelerated aging).
Negative delta = biological age decreased (recovery/improvement).

Compare severity groups to assess whether more severe COVID leads to
greater organ-specific biological age acceleration.

## Disclaimer

ClawBio is a research and educational tool. It is not a medical device and
does not provide clinical diagnoses. Consult a healthcare professional
before making any medical decisions.
"""
    (output_dir / "treatment_report.md").write_text(report)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Treatment effect analysis on organ biological age (Filbin COVID data)"
    )
    parser.add_argument(
        "--data-dir", required=True,
        help="Directory with preprocessed Filbin CSVs (from fetch_filbin.py)"
    )
    parser.add_argument("--output", required=True, help="Output directory")
    parser.add_argument(
        "--severity-col", default="Acuity",
        help="Column name for severity grouping (default: Acuity)"
    )
    parser.add_argument(
        "--organs", default=",".join(ORGANS_OF_INTEREST),
        help="Comma-separated organs to analyze"
    )
    args = parser.parse_args()

    data_dir = Path(args.data_dir)
    output_dir = Path(args.output)
    organs = [o.strip() for o in args.organs.split(",")]

    # Find timepoint files
    t0_path = data_dir / "filbin_olink_D0.csv.gz"
    t7_path = data_dir / "filbin_olink_D7.csv.gz"
    meta_path = data_dir / "filbin_clinical_metadata.csv"

    for p in [t0_path, t7_path, meta_path]:
        if not p.exists():
            print(f"ERROR: Missing {p}. Run fetch_filbin.py first.")
            sys.exit(1)

    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "figures").mkdir(exist_ok=True)

    # Run clock on each timepoint
    print("Running proteomics clock on Day 0 ...")
    clock_t0_dir = output_dir / "clock_D0"
    preds_t0 = run_clock_on_timepoint(t0_path, clock_t0_dir, organs)

    print("Running proteomics clock on Day 7 ...")
    clock_t7_dir = output_dir / "clock_D7"
    preds_t7 = run_clock_on_timepoint(t7_path, clock_t7_dir, organs)

    # Load metadata and compute deltas
    meta = pd.read_csv(meta_path)
    severity_col = args.severity_col
    if severity_col not in meta.columns:
        print(f"WARNING: '{severity_col}' not in metadata. Available: {list(meta.columns)}")
        severity_col = meta.columns[1]
        print(f"Falling back to: {severity_col}")

    delta = compute_delta_ages(preds_t0, preds_t7, meta, severity_col)
    n_paired = delta["subject_id"].nunique()

    print(f"Paired subjects: {n_paired}")
    print(f"Severity groups: {sorted(delta[severity_col].dropna().unique())}")

    # Plots
    available_organs = [o for o in organs if o in delta.columns]
    plot_delta_boxplots(delta, severity_col, available_organs, output_dir / "figures" / "delta_boxplots.png")
    plot_delta_heatmap(delta, severity_col, available_organs, output_dir / "figures" / "delta_heatmap.png")

    # Save tables
    delta.to_csv(output_dir / "delta_ages.csv", index=False)

    # Report
    write_treatment_report(
        output_dir=output_dir,
        delta=delta,
        severity_col=severity_col,
        organs=organs,
        n_t0=len(preds_t0),
        n_t7=len(preds_t7),
        n_paired=n_paired,
    )

    print(f"\nDone. Report: {output_dir / 'treatment_report.md'}")
    print(json.dumps({
        "n_t0": len(preds_t0),
        "n_t7": len(preds_t7),
        "n_paired": n_paired,
        "organs": available_organs,
        "output_dir": str(output_dir),
    }, indent=2))


if __name__ == "__main__":
    main()
