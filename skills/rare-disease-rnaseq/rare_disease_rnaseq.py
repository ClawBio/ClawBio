#!/usr/bin/env python3
"""Blood RNA-seq expression-outlier detection for rare-disease diagnostics.

Reproduces the diagnostic principle of the Genomics England NGRL paper
(medRxiv 2026.03.19.26348811): per-gene outlier scoring of case samples
against a control reference panel, filtered by a haploinsufficient
disease-gene panel.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from clawbio.common.report import write_result_json


SKILL = "rare-disease-rnaseq"
VERSION = "0.1.0"
DEFAULT_Z_THRESHOLD = 3.0
DEFAULT_PANEL = Path(__file__).parent / "data" / "disease_panel.csv"

DISCLAIMER = (
    "ClawBio is a research and educational tool. "
    "It is not a medical device and does not provide clinical diagnoses. "
    "Consult a healthcare professional before making any medical decisions."
)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--counts", type=Path, help="Counts matrix (csv/tsv): genes x samples")
    p.add_argument("--cases", type=Path, help="File with one case sample ID per line")
    p.add_argument("--controls", type=Path, help="File with one control sample ID per line")
    p.add_argument("--panel", type=Path, default=DEFAULT_PANEL,
                   help="Disease-gene panel CSV (gene,mechanism,phenotype)")
    p.add_argument("--z-threshold", type=float, default=DEFAULT_Z_THRESHOLD,
                   help=f"|z| threshold for outlier call (default {DEFAULT_Z_THRESHOLD})")
    p.add_argument("--output", type=Path, default=Path("rdoutlier_report"))
    p.add_argument("--demo", action="store_true",
                   help="Generate synthetic data and run end-to-end")
    p.add_argument("--seed", type=int, default=42)
    return p.parse_args()


def _sep_for(path: Path) -> str:
    return "\t" if path.suffix.lower() == ".tsv" else ","


def generate_demo(seed: int, output_dir: Path) -> tuple[Path, Path, Path]:
    """Synthetic blood RNA-seq: 100 controls + 2 cases with injected outliers.

    Case-001: FBN1 expression collapsed to ~10% of mean (Marfan-like LoF)
    Case-002: NF1 expression elevated 4x mean (atypical up-outlier)
    """
    rng = np.random.default_rng(seed)
    panel_df = pd.read_csv(DEFAULT_PANEL)
    panel_genes = panel_df["gene"].tolist()
    n_panel = len(panel_genes)
    n_background = 150
    background_genes = [f"BG{i:04d}" for i in range(n_background)]
    genes = panel_genes + background_genes
    n_genes = len(genes)
    n_controls = 100
    n_cases = 2
    control_ids = [f"CTRL-{i:03d}" for i in range(n_controls)]
    case_ids = ["CASE-001", "CASE-002"]
    samples = control_ids + case_ids

    # Per-gene mean expression (counts) drawn from a realistic blood-RNAseq
    # distribution: log-normal-ish, range ~50 to ~5000.
    gene_means = np.exp(rng.normal(loc=6.5, scale=1.0, size=n_genes))
    # Per-gene dispersion; negative-binomial-like noise
    counts = np.zeros((n_genes, len(samples)), dtype=int)
    for j, _ in enumerate(samples):
        # Per-sample size factor (library size variation)
        sf = rng.lognormal(mean=0.0, sigma=0.15)
        for i in range(n_genes):
            mu = gene_means[i] * sf
            # Use Poisson + extra noise as a stand-in for NB
            counts[i, j] = int(rng.poisson(mu) + rng.normal(0, mu * 0.10))
    counts = np.clip(counts, 0, None)

    # Inject outliers
    fbn1_idx = genes.index("FBN1")
    nf1_idx = genes.index("NF1")
    case1_col = samples.index("CASE-001")
    case2_col = samples.index("CASE-002")
    counts[fbn1_idx, case1_col] = int(gene_means[fbn1_idx] * 0.10)
    counts[nf1_idx, case2_col] = int(gene_means[nf1_idx] * 4.0)

    counts_df = pd.DataFrame(counts, index=genes, columns=samples)
    counts_df.index.name = "gene"

    output_dir.mkdir(parents=True, exist_ok=True)
    counts_path = output_dir / "demo_counts.csv"
    cases_path = output_dir / "demo_cases.txt"
    controls_path = output_dir / "demo_controls.txt"
    counts_df.to_csv(counts_path)
    cases_path.write_text("\n".join(case_ids) + "\n")
    controls_path.write_text("\n".join(control_ids) + "\n")
    return counts_path, cases_path, controls_path


def cpm_log2(counts: pd.DataFrame) -> pd.DataFrame:
    lib = counts.sum(axis=0).replace(0, 1)
    cpm = counts.divide(lib, axis=1) * 1e6
    return np.log2(cpm + 1.0)


def per_gene_robust_stats(expr: pd.DataFrame, control_ids: list[str]) -> pd.DataFrame:
    ctrl = expr[control_ids]
    median = ctrl.median(axis=1)
    # Median absolute deviation
    mad = (ctrl.subtract(median, axis=0)).abs().median(axis=1)
    return pd.DataFrame({"median": median, "mad": mad})


def call_outliers(
    expr: pd.DataFrame,
    case_ids: list[str],
    stats: pd.DataFrame,
    panel: pd.DataFrame,
    z_threshold: float,
) -> pd.DataFrame:
    panel_lookup = panel.set_index("gene")[["mechanism", "phenotype"]]
    rows = []
    # Modified z-score: 0.6745 * (x - median) / MAD; skip zero-MAD genes
    for gene in expr.index:
        m = stats.at[gene, "median"]
        d = stats.at[gene, "mad"]
        if d <= 0 or np.isnan(d):
            continue
        for case in case_ids:
            x = expr.at[gene, case]
            z = 0.6745 * (x - m) / d
            if abs(z) >= z_threshold:
                in_panel = gene in panel_lookup.index
                rows.append({
                    "case": case,
                    "gene": gene,
                    "z_score": round(float(z), 2),
                    "direction": "down" if z < 0 else "up",
                    "log2_cpm_case": round(float(x), 2),
                    "log2_cpm_control_median": round(float(m), 2),
                    "in_panel": in_panel,
                    "mechanism": panel_lookup.at[gene, "mechanism"] if in_panel else "",
                    "phenotype": panel_lookup.at[gene, "phenotype"] if in_panel else "",
                })
    columns = ["case", "gene", "z_score", "direction", "log2_cpm_case",
               "log2_cpm_control_median", "in_panel", "mechanism", "phenotype"]
    df = pd.DataFrame(rows, columns=columns)
    if df.empty:
        return df
    df = df.sort_values(["case", "in_panel", "z_score"],
                        key=lambda c: c.abs() if c.name == "z_score" else c,
                        ascending=[True, False, False])
    return df.reset_index(drop=True)


def plot_heatmap(
    expr: pd.DataFrame, case_ids: list[str], stats: pd.DataFrame,
    outliers: pd.DataFrame, output_path: Path, top_n: int = 20,
) -> None:
    if outliers.empty:
        return
    top_genes = outliers.head(top_n)["gene"].unique().tolist()
    z_matrix = []
    for gene in top_genes:
        m = stats.at[gene, "median"]
        d = stats.at[gene, "mad"] or 1.0
        row = [0.6745 * (expr.at[gene, c] - m) / d for c in case_ids]
        z_matrix.append(row)
    z_arr = np.array(z_matrix)
    fig, ax = plt.subplots(figsize=(max(4, 0.6 * len(case_ids)), max(4, 0.4 * len(top_genes))))
    vmax = max(3.0, np.nanmax(np.abs(z_arr)))
    im = ax.imshow(z_arr, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    ax.set_xticks(range(len(case_ids)))
    ax.set_xticklabels(case_ids, rotation=45, ha="right")
    ax.set_yticks(range(len(top_genes)))
    ax.set_yticklabels(top_genes)
    ax.set_title("Per-case modified z-score (top outlier genes)")
    fig.colorbar(im, ax=ax, label="z-score")
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


def write_report(
    output_dir: Path, outliers: pd.DataFrame, case_ids: list[str],
    n_controls: int, n_genes_screened: int, z_threshold: float,
) -> None:
    lines = [
        "# Rare-Disease Blood RNA-seq Outlier Report",
        "",
        f"Generated: {datetime.now(timezone.utc).isoformat()}",
        f"Skill: `{SKILL}` v{VERSION}",
        "",
        "## Cohort summary",
        "",
        f"- Cases analysed: **{len(case_ids)}**",
        f"- Control reference panel: **{n_controls}** samples",
        f"- Genes screened: **{n_genes_screened}**",
        f"- Outlier threshold: **|z| >= {z_threshold}**",
        f"- Outlier events called: **{len(outliers)}**",
        "",
        "## Per-case candidate diagnoses",
        "",
    ]
    if outliers.empty:
        lines.append("_No outliers above threshold._")
    else:
        for case in case_ids:
            sub = outliers[outliers["case"] == case]
            lines.append(f"### {case}")
            lines.append("")
            if sub.empty:
                lines.append("_No outliers above threshold._")
                lines.append("")
                continue
            in_panel = sub[sub["in_panel"]]
            off_panel = sub[~sub["in_panel"]]
            lines.append(f"- Total outlier genes: {len(sub)}")
            lines.append(f"- In disease panel: {len(in_panel)}")
            lines.append(f"- Off-panel (informational): {len(off_panel)}")
            lines.append("")
            if not in_panel.empty:
                lines.append("**Disease-panel outliers (clinical priority):**")
                lines.append("")
                lines.append("| Gene | Direction | z-score | Mechanism | Associated phenotype |")
                lines.append("|------|-----------|---------|-----------|----------------------|")
                for _, row in in_panel.iterrows():
                    lines.append(
                        f"| {row['gene']} | {row['direction']} | {row['z_score']} | "
                        f"{row['mechanism']} | {row['phenotype']} |"
                    )
                lines.append("")
                lines.append("**Suggested next step:** correlate flagged outlier(s) with phenotype, "
                             "review WGS variant calls in the same gene, present at MDT.")
                lines.append("")
    lines += [
        "## Method",
        "",
        "Per-gene robust outlier scoring on log2(CPM+1). For each gene, "
        "modified z-score = 0.6745 * (x - median) / MAD computed on the control "
        "reference panel. Outliers are events with |z| above threshold; clinical "
        "priority is given to genes in the dosage-sensitive disease panel.",
        "",
        "**Production swap:** clinical-grade calls use the DROP pipeline "
        "(gagneurlab/drop) with OUTRIDER (denoising autoencoder, NB-distributed) "
        "and FRASER2 (splicing outliers). The skill's I/O contract is unchanged.",
        "",
        "---",
        "",
        f"*{DISCLAIMER}*",
        "",
    ]
    (output_dir / "report.md").write_text("\n".join(lines))


def write_reproducibility(output_dir: Path, args: argparse.Namespace,
                          input_checksum: str) -> None:
    repro = output_dir / "reproducibility"
    repro.mkdir(parents=True, exist_ok=True)
    cmd = " ".join([sys.executable, "rare_disease_rnaseq.py"] + sys.argv[1:])
    (repro / "commands.sh").write_text(f"#!/usr/bin/env bash\nset -euo pipefail\n{cmd}\n")
    (repro / "environment.yml").write_text(
        "name: clawbio-rdoutlier\n"
        "dependencies:\n"
        "  - python>=3.10\n"
        "  - pandas\n"
        "  - numpy\n"
        "  - matplotlib\n"
    )
    (repro / "checksums.sha256").write_text(f"{input_checksum}  counts_input\n")


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    args = parse_args()
    output_dir = args.output.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "figures").mkdir(exist_ok=True)
    (output_dir / "tables").mkdir(exist_ok=True)

    if args.demo:
        counts_path, cases_path, controls_path = generate_demo(args.seed, output_dir / "demo_data")
    else:
        if not (args.counts and args.cases and args.controls):
            print("ERROR: --counts, --cases and --controls required (or use --demo)",
                  file=sys.stderr)
            return 2
        counts_path, cases_path, controls_path = args.counts, args.cases, args.controls

    counts = pd.read_csv(counts_path, sep=_sep_for(counts_path), index_col=0)
    case_ids = [s.strip() for s in cases_path.read_text().splitlines() if s.strip()]
    control_ids = [s.strip() for s in controls_path.read_text().splitlines() if s.strip()]

    missing = [s for s in case_ids + control_ids if s not in counts.columns]
    if missing:
        print(f"ERROR: sample IDs not in counts matrix: {missing[:5]}", file=sys.stderr)
        return 3

    panel = pd.read_csv(args.panel)
    expr = cpm_log2(counts)
    stats = per_gene_robust_stats(expr, control_ids)
    outliers = call_outliers(expr, case_ids, stats, panel, args.z_threshold)

    outliers.to_csv(output_dir / "tables" / "outlier_calls.csv", index=False)
    stats.to_csv(output_dir / "tables" / "per_gene_stats.csv")

    plot_heatmap(expr, case_ids, stats, outliers,
                 output_dir / "figures" / "case_outlier_heatmap.png")

    write_report(output_dir, outliers, case_ids,
                 n_controls=len(control_ids),
                 n_genes_screened=int((stats["mad"] > 0).sum()),
                 z_threshold=args.z_threshold)

    checksum = sha256_of(counts_path)
    write_reproducibility(output_dir, args, checksum)

    panel_outliers = outliers[outliers["in_panel"]] if not outliers.empty else outliers
    summary = {
        "cases": len(case_ids),
        "controls": len(control_ids),
        "genes_screened": int((stats["mad"] > 0).sum()),
        "z_threshold": args.z_threshold,
        "outlier_events": int(len(outliers)),
        "panel_outlier_events": int(len(panel_outliers)),
        "cases_with_panel_hit": int(panel_outliers["case"].nunique()) if not panel_outliers.empty else 0,
    }
    write_result_json(output_dir, SKILL, VERSION, summary,
                      data={"outliers": outliers.to_dict(orient="records")},
                      input_checksum=checksum)

    print(json.dumps(summary, indent=2))
    print(f"\nReport: {output_dir / 'report.md'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
