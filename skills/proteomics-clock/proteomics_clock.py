#!/usr/bin/env python3
"""Compute organ-specific biological age from Olink proteomic data.

Implements the organ aging clocks from Goeminne et al. (2025) Cell Metabolism.
Coefficients are downloaded at runtime from the organAging GitHub repository.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests


DISCLAIMER = (
    "ClawBio is a research and educational tool. "
    "It is not a medical device and does not provide clinical diagnoses. "
    "Consult a healthcare professional before making any medical decisions."
)

CITATION = (
    "Goeminne LJE, Vladimirova A, Eames A, Tyshkovskiy A, Argentieri MA, "
    "Ying K, Moqri M, Gladyshev VN (2025). Plasma protein-based organ-specific "
    "aging and mortality models unveil diseases as accelerated aging of organismal "
    "systems. Cell Metabolism 37(1):205-222.e6. "
    "DOI: 10.1016/j.cmet.2024.10.005"
)

# Pinned to organAging commit 5147b03 for reproducibility.
# To update: change this SHA and clear the local cache (~/.cache/clawbio/proteomics-clock/).
ORGANAGING_COMMIT = "5147b0301ec7f4abdb10ef650d04f47454ddc8fd"
GITHUB_RAW_BASE = (
    f"https://raw.githubusercontent.com/ludgergoeminne/organAging/{ORGANAGING_COMMIT}/data/output_Python"
)

ORGAN_PROTEINS_URL = f"{GITHUB_RAW_BASE}/GTEx_4x_FC_genes.json"

COEF_URL_TEMPLATES = {
    "gen1": f"{GITHUB_RAW_BASE}/instance_0/chronological_models/{{organ}}_coefs_GTEx_4x_FC.csv",
    "gen2": f"{GITHUB_RAW_BASE}/instance_0/mortality_based_models/{{organ}}_mortality_coefs_GTEx_4x_FC.csv",
}

# Gompertz conversion constants (from organAging README)
GOMPERTZ_INTERCEPT = -9.94613787413831
GOMPERTZ_SLOPE = 0.0897860500778604
GOMPERTZ_AVG_HAZARD = -4.801912

ALL_ORGANS = [
    "Adipose", "Adrenal", "Artery", "Brain", "Conventional", "Esophagus",
    "Female", "Heart", "Immune", "Intestine", "Kidney", "Liver", "Lung",
    "Male", "Multi-organ", "Muscle", "Organismal", "Pancreas", "Pituitary",
    "Salivary", "Skin", "Stomach", "Thyroid",
]

# Bladder has 0 proteins — always excluded
DEFAULT_ORGANS = list(ALL_ORGANS)


def _cache_dir() -> Path:
    base = os.environ.get("CLAWBIO_CACHE")
    if base:
        d = Path(base) / "proteomics-clock"
    else:
        d = Path.home() / ".cache" / "clawbio" / "proteomics-clock"
    d.mkdir(parents=True, exist_ok=True)
    return d


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def _get_skill_data_dir() -> Path:
    return Path(__file__).resolve().parent / "data"


# ---------------------------------------------------------------------------
# Parsing helpers
# ---------------------------------------------------------------------------


def parse_organ_list(organs: str | None) -> list[str]:
    if organs is None:
        return list(DEFAULT_ORGANS)
    values = [item.strip() for item in organs.split(",") if item.strip()]
    if not values:
        raise ValueError("Organ list is empty")
    return values


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------


def load_input(path: Path) -> pd.DataFrame:
    suffixes = [s.lower() for s in path.suffixes]
    if suffixes[-2:] == [".csv", ".gz"]:
        return pd.read_csv(path, compression="gzip")
    if suffixes[-2:] == [".tsv", ".gz"]:
        return pd.read_csv(path, sep="\t", compression="gzip")
    if suffixes and suffixes[-1] == ".csv":
        return pd.read_csv(path)
    if suffixes and suffixes[-1] == ".tsv":
        return pd.read_csv(path, sep="\t")
    raise ValueError(
        f"Unsupported format: {path.name}. Use .csv, .tsv, .csv.gz, or .tsv.gz"
    )


# ---------------------------------------------------------------------------
# Coefficient downloading and caching
# ---------------------------------------------------------------------------


def _download_text(url: str, cache_name: str) -> str:
    cached = _cache_dir() / cache_name
    if cached.exists():
        return cached.read_text()
    resp = requests.get(url, timeout=30)
    resp.raise_for_status()
    cached.write_text(resp.text)
    return resp.text


def download_organ_proteins() -> dict[str, list[str]]:
    text = _download_text(ORGAN_PROTEINS_URL, "GTEx_4x_FC_genes.json")
    return json.loads(text)


def download_coefficients(
    organs: list[str], generation: str, fold: int = 1
) -> dict[str, dict[str, float]]:
    """Download and parse coefficients for the given organs and generation.

    Returns {organ: {protein_or_Intercept: coefficient}}.
    """
    gens = ["gen1", "gen2"] if generation == "both" else [generation]
    result: dict[str, dict[str, float]] = {}

    for gen in gens:
        url_template = COEF_URL_TEMPLATES[gen]
        for organ in organs:
            key = f"{organ}_{gen}"
            cache_name = f"{organ}_{gen}_fold{fold}.csv"
            url = url_template.format(organ=organ)

            try:
                text = _download_text(url, cache_name)
            except requests.HTTPError:
                continue

            reader = csv.reader(io.StringIO(text))
            header = next(reader)
            rows = list(reader)

            if fold < 1 or fold > len(rows):
                fold_idx = 0
            else:
                fold_idx = fold - 1

            row = rows[fold_idx]
            coefs = {}
            for col_name, val in zip(header, row):
                try:
                    coefs[col_name] = float(val)
                except ValueError:
                    continue

            result[key] = coefs

    return result


# ---------------------------------------------------------------------------
# Prediction
# ---------------------------------------------------------------------------


def predict_organ_ages(
    df: pd.DataFrame,
    coefficients: dict[str, dict[str, float]],
    organs: list[str],
    generation: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Predict organ ages for given generation.

    Returns (predictions_df, missing_proteins_df).
    """
    predictions = pd.DataFrame({"sample_id": df["sample_id"] if "sample_id" in df.columns else range(len(df))})
    missing_rows: list[dict[str, str]] = []

    for organ in organs:
        key = f"{organ}_{generation}"
        if key in coefficients:
            coefs = coefficients[key]
        elif organ in coefficients:
            coefs = coefficients[organ]
        else:
            continue
        intercept = coefs.get("Intercept", 0.0) if generation == "gen1" else 0.0
        protein_coefs = {k: v for k, v in coefs.items() if k != "Intercept"}

        available = [p for p in protein_coefs if p in df.columns]
        missing = [p for p in protein_coefs if p not in df.columns]

        for protein in missing:
            missing_rows.append({"organ": organ, "protein": protein, "generation": generation})

        if not available:
            predictions[organ] = np.nan
            continue

        vals = df[available].values.astype(float)
        weights = np.array([protein_coefs[p] for p in available])
        predictions[organ] = intercept + vals @ weights

    missing_df = pd.DataFrame(missing_rows, columns=["organ", "protein", "generation"])
    return predictions, missing_df


def mortality_to_years(hazard_values: np.ndarray) -> np.ndarray:
    """Convert relative log(mortality hazard) to age in years via Gompertz."""
    return (
        (-GOMPERTZ_AVG_HAZARD + hazard_values) / GOMPERTZ_SLOPE
        - GOMPERTZ_INTERCEPT
    )


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------


def plot_organ_distributions(predictions: pd.DataFrame, outpath: Path) -> None:
    numeric = predictions.select_dtypes(include=["number"])
    numeric = numeric.dropna(axis=1, how="all")
    if numeric.empty:
        return
    fig, ax = plt.subplots(figsize=(max(8, 0.8 * len(numeric.columns)), 6))
    ax.boxplot(
        [numeric[col].dropna().values for col in numeric.columns],
        tick_labels=numeric.columns,
    )
    ax.set_xticklabels(numeric.columns, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Predicted biological age (years)")
    ax.set_title("Organ-Specific Proteomic Age Predictions")
    plt.tight_layout()
    plt.savefig(outpath, dpi=180)
    plt.close()


def plot_organ_correlation(predictions: pd.DataFrame, outpath: Path) -> None:
    import seaborn as sns
    numeric = predictions.select_dtypes(include=["number"]).dropna(axis=1, how="all")
    if numeric.shape[1] < 2:
        return
    corr = numeric.corr("pearson")
    fig, ax = plt.subplots(figsize=(max(6, 0.5 * len(corr.columns)), max(5, 0.4 * len(corr.columns))))
    sns.heatmap(corr, vmin=-1, vmax=1, cmap="RdBu_r", annot=corr.shape[0] <= 12, fmt=".2f", ax=ax)
    ax.set_title("Organ Clock Correlation (Pearson)")
    plt.tight_layout()
    plt.savefig(outpath, dpi=180)
    plt.close()


def plot_organ_heatmap(predictions: pd.DataFrame, sample_ids: pd.Series, outpath: Path) -> None:
    import seaborn as sns
    numeric = predictions.select_dtypes(include=["number"]).dropna(axis=1, how="all")
    if numeric.empty:
        return
    fig, ax = plt.subplots(figsize=(max(8, 0.6 * len(numeric.columns)), max(4, 0.3 * len(numeric))))
    sns.heatmap(numeric.values, xticklabels=numeric.columns, yticklabels=sample_ids.values,
                cmap="YlOrRd", ax=ax)
    ax.set_title("Biological Age by Organ and Sample")
    ax.set_xlabel("Organ Clock")
    ax.set_ylabel("Sample")
    plt.xticks(rotation=45, ha="right", fontsize=8)
    plt.yticks(fontsize=7)
    plt.tight_layout()
    plt.savefig(outpath, dpi=180)
    plt.close()


# ---------------------------------------------------------------------------
# Report generation
# ---------------------------------------------------------------------------


def write_reproducibility(output_dir: Path, input_desc: str, command_args: list[str]) -> None:
    repro_dir = output_dir / "reproducibility"
    repro_dir.mkdir(parents=True, exist_ok=True)

    command_text = "python skills/proteomics-clock/proteomics_clock.py " + " ".join(command_args) + "\n"
    (repro_dir / "commands.sh").write_text(command_text)

    env_text = """name: clawbio-proteomics-clock
channels:
  - conda-forge
dependencies:
  - python>=3.11
  - pandas
  - numpy
  - matplotlib
  - seaborn
  - requests
"""
    (repro_dir / "environment.yml").write_text(env_text)

    checksums = []
    input_path = Path(input_desc)
    if input_path.exists():
        checksums.append(f"{_sha256(input_path)}  {input_path.name}")
    for subdir in ["tables", "figures"]:
        for path in sorted((output_dir / subdir).glob("*")):
            if path.is_file():
                checksums.append(f"{_sha256(path)}  {subdir}/{path.name}")
    (repro_dir / "checksums.sha256").write_text("\n".join(checksums) + "\n")


def write_report(
    output_dir: Path,
    input_desc: str,
    organs: list[str],
    generation: str,
    fold: int,
    gen1_preds: pd.DataFrame | None,
    gen2_preds: pd.DataFrame | None,
    missing_df: pd.DataFrame,
    metadata_json: dict,
) -> None:
    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    n_samples = gen1_preds.shape[0] if gen1_preds is not None else (gen2_preds.shape[0] if gen2_preds is not None else 0)

    summary_rows = []
    for preds, gen_label in [(gen1_preds, "gen1"), (gen2_preds, "gen2")]:
        if preds is None:
            continue
        for col in preds.columns:
            if col == "sample_id":
                continue
            series = pd.to_numeric(preds[col], errors="coerce")
            valid = series.dropna()
            if valid.empty:
                continue
            summary_rows.append(
                f"| {col} | {gen_label} | {valid.count()} | {valid.mean():.2f} | {valid.std(ddof=1):.2f} |"
            )

    summary_block = "\n".join(summary_rows) if summary_rows else "| (none) | - | 0 | n/a | n/a |"

    missing_counts = {}
    if not missing_df.empty:
        for (organ, gen), group in missing_df.groupby(["organ", "generation"]):
            missing_counts[f"{organ} ({gen})"] = len(group)
    missing_block = (
        "No missing proteins were detected."
        if not missing_counts
        else "\n".join([f"- {k}: {v} missing" for k, v in sorted(missing_counts.items())])
    )

    report = f"""# ClawBio Proteomics Clock Report

**Date**: {now}
**Input**: `{input_desc}`
**Samples**: {n_samples}
**Organs requested**: {", ".join(organs)}
**Generation**: {generation}
**Fold**: {fold}

## Method

Organ-specific proteomic aging clocks from Goeminne et al. (2025).
Predictions are computed as weighted sums of Olink NPX protein values
using elastic net coefficients trained on UK Biobank data.

- **Gen1 (chronological)**: predicted age = intercept + sum(NPX_i * coef_i)
- **Gen2 (mortality-based)**: predicted log-hazard = sum(NPX_i * coef_i), converted to years via Gompertz transform

## Outputs

- Gen1 predictions: `tables/predictions_gen1.csv`
- Gen2 predictions: `tables/predictions_gen2.csv`
- Summary: `tables/prediction_summary.csv`
- Missing proteins: `tables/missing_proteins.csv`
- Metadata: `tables/clock_metadata.json`

## Prediction Summary

| Organ | Generation | N | Mean | Std |
|---|---|---:|---:|---:|
{summary_block}

## Missing Proteins

{missing_block}

## Reproducibility

- Commands: `reproducibility/commands.sh`
- Environment: `reproducibility/environment.yml`
- Checksums: `reproducibility/checksums.sha256`

## Citation

{CITATION}

## Disclaimer

{DISCLAIMER}
"""
    (output_dir / "report.md").write_text(report)


# ---------------------------------------------------------------------------
# Main analysis pipeline
# ---------------------------------------------------------------------------


def run_analysis(
    output_dir: Path,
    organs: list[str],
    generation: str,
    fold: int,
    convert_mortality_to_years: bool,
    age_column: str | None,
    verbose: bool,
    input_path: Path | None = None,
) -> dict:
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(
            f"Output directory '{output_dir}' is not empty. "
            "Choose a new --output path to avoid overwriting."
        )

    if input_path is None or not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "tables").mkdir(exist_ok=True)
    (output_dir / "figures").mkdir(exist_ok=True)

    if verbose:
        print(f"Loading input: {input_path}")
    df = load_input(input_path)
    if df.empty:
        raise ValueError("Input data is empty")

    if verbose:
        print(f"Downloading coefficients for {len(organs)} organs, generation={generation}, fold={fold}")
    coefficients = download_coefficients(organs, generation, fold)
    if not coefficients:
        raise RuntimeError("No coefficients could be downloaded. Check network and organ names.")

    gen1_preds = None
    gen2_preds = None
    all_missing = pd.DataFrame(columns=["organ", "protein", "generation"])

    if generation in ("gen1", "both"):
        if verbose:
            print("Predicting gen1 (chronological) ages...")
        gen1_preds, gen1_missing = predict_organ_ages(df, coefficients, organs, "gen1")
        all_missing = pd.concat([all_missing, gen1_missing], ignore_index=True)

    if generation in ("gen2", "both"):
        if verbose:
            print("Predicting gen2 (mortality-based) ages...")
        gen2_preds, gen2_missing = predict_organ_ages(df, coefficients, organs, "gen2")
        all_missing = pd.concat([all_missing, gen2_missing], ignore_index=True)

        if convert_mortality_to_years and gen2_preds is not None:
            for col in gen2_preds.columns:
                if col == "sample_id":
                    continue
                vals = gen2_preds[col].values.astype(float)
                gen2_preds[col] = mortality_to_years(vals)

    # Save tables
    if gen1_preds is not None:
        gen1_preds.to_csv(output_dir / "tables" / "predictions_gen1.csv", index=False)
    else:
        pd.DataFrame().to_csv(output_dir / "tables" / "predictions_gen1.csv", index=False)

    if gen2_preds is not None:
        gen2_preds.to_csv(output_dir / "tables" / "predictions_gen2.csv", index=False)
    else:
        pd.DataFrame().to_csv(output_dir / "tables" / "predictions_gen2.csv", index=False)

    # Build summary
    summary_parts = []
    for preds, gen_label in [(gen1_preds, "gen1"), (gen2_preds, "gen2")]:
        if preds is None:
            continue
        numeric = preds.select_dtypes(include=["number"])
        if numeric.empty:
            continue
        desc = numeric.describe().T.reset_index().rename(columns={"index": "organ"})
        desc["generation"] = gen_label
        summary_parts.append(desc)
    if summary_parts:
        summary = pd.concat(summary_parts, ignore_index=True)
    else:
        summary = pd.DataFrame()
    summary.to_csv(output_dir / "tables" / "prediction_summary.csv", index=False)
    all_missing.to_csv(output_dir / "tables" / "missing_proteins.csv", index=False)

    # Clock metadata
    metadata = {
        "organs_requested": organs,
        "generation": generation,
        "fold": fold,
        "n_coefficients_loaded": len(coefficients),
        "convert_mortality_to_years": convert_mortality_to_years,
    }
    (output_dir / "tables" / "clock_metadata.json").write_text(json.dumps(metadata, indent=2))

    # Figures — use gen1 preferentially, fall back to gen2
    plot_preds = gen1_preds if gen1_preds is not None else gen2_preds
    if plot_preds is not None:
        sample_ids = plot_preds["sample_id"] if "sample_id" in plot_preds.columns else pd.Series(range(len(plot_preds)))
        numeric_preds = plot_preds.drop(columns=["sample_id"], errors="ignore")
        plot_organ_distributions(numeric_preds, output_dir / "figures" / "organ_distributions.png")
        plot_organ_correlation(numeric_preds, output_dir / "figures" / "organ_correlation.png")
        plot_organ_heatmap(numeric_preds, sample_ids, output_dir / "figures" / "organ_heatmap.png")

    # Report
    write_report(
        output_dir=output_dir,
        input_desc=str(input_path),
        organs=organs,
        generation=generation,
        fold=fold,
        gen1_preds=gen1_preds,
        gen2_preds=gen2_preds,
        missing_df=all_missing,
        metadata_json=metadata,
    )

    clocks_found = []
    for preds in [gen1_preds, gen2_preds]:
        if preds is not None:
            clocks_found.extend([c for c in preds.columns if c != "sample_id"])

    return {
        "output_dir": str(output_dir),
        "input": str(input_path),
        "n_samples": int(df.shape[0]),
        "clocks_found": clocks_found,
    }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Compute organ-specific biological age from Olink proteomic data "
        "(Goeminne et al. 2025)"
    )
    parser.add_argument("--input", help="Olink NPX file (.csv/.tsv/.csv.gz/.tsv.gz)")
    parser.add_argument("--output", required=True, help="Output directory")
    parser.add_argument(
        "--organs",
        help="Comma-separated organ list; defaults to all 23 organs (excluding Bladder)",
    )
    parser.add_argument(
        "--generation",
        default="both",
        choices=["gen1", "gen2", "both"],
        help="Model generation: gen1 (chronological), gen2 (mortality), or both (default: both)",
    )
    parser.add_argument(
        "--fold",
        type=int,
        default=1,
        help="Cross-validation fold to use (1-5, default: 1)",
    )
    parser.add_argument(
        "--no-convert-mortality",
        action="store_true",
        help="Keep gen2 output as log-hazard ratios instead of converting to years",
    )
    parser.add_argument(
        "--age-column",
        default="age",
        help="Column name for chronological age (for residual calculation)",
    )
    parser.add_argument("--demo", action="store_true", help="Use bundled synthetic demo data")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose output")
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    if args.demo:
        input_path = _get_skill_data_dir() / "demo_olink_npx.csv.gz"
    else:
        if not args.input:
            parser.error("Provide --input or use --demo")
        input_path = Path(args.input)

    organs = parse_organ_list(args.organs)

    result = run_analysis(
        output_dir=Path(args.output),
        organs=organs,
        generation=args.generation,
        fold=args.fold,
        convert_mortality_to_years=not args.no_convert_mortality,
        age_column=args.age_column,
        verbose=args.verbose,
        input_path=input_path,
    )

    # Reproducibility
    command_args = []
    if args.demo:
        command_args.append("--demo")
    else:
        command_args.extend(["--input", str(input_path)])
    command_args.extend(["--output", str(args.output)])
    command_args.extend(["--organs", ",".join(organs)])
    command_args.extend(["--generation", args.generation])
    command_args.extend(["--fold", str(args.fold)])
    if args.no_convert_mortality:
        command_args.append("--no-convert-mortality")
    if args.verbose:
        command_args.append("--verbose")

    write_reproducibility(
        output_dir=Path(args.output),
        input_desc=result["input"],
        command_args=command_args,
    )

    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
