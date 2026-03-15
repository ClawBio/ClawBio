#!/usr/bin/env python3

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyaging as pya


DISCLAIMER = (
    "ClawBio is a research and educational tool. "
    "It is not a medical device and does not provide clinical diagnoses. "
    "Consult a healthcare professional before making any medical decisions."
)

DEFAULT_CLOCKS = [
    "Horvath2013",
    "AltumAge",
    "PCGrimAge",
    "GrimAge2",
    "DunedinPACE",
]


def parse_clock_list(clocks: str | None) -> list[str]:
    if not clocks:
        return list(DEFAULT_CLOCKS)
    values = [item.strip() for item in clocks.split(",") if item.strip()]
    if not values:
        raise ValueError("Clock list is empty")
    return values


def parse_metadata_cols(metadata_cols: str | None) -> list[str]:
    if not metadata_cols:
        return ["gender", "tissue_type", "dataset"]
    return [item.strip() for item in metadata_cols.split(",") if item.strip()]


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def _load_dataframe(path: Path) -> pd.DataFrame:
    suffix = path.suffix.lower()
    if suffix in {".pkl", ".pickle"}:
        return pd.read_pickle(path)
    if suffix == ".csv":
        return pd.read_csv(path)
    if suffix == ".tsv":
        return pd.read_csv(path, sep="\t")
    raise ValueError("Unsupported input format. Use .pkl, .pickle, .csv, or .tsv")


def _get_repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _resolve_geo_pickle(geo_id: str, verbose: bool) -> Path:
    pya.data.download_example_data(geo_id, verbose=verbose)

    candidates = [
        Path.cwd() / "pyaging_data" / f"{geo_id}.pkl",
        _get_repo_root() / "pyaging_data" / f"{geo_id}.pkl",
        Path(__file__).resolve().parent / "pyaging_data" / f"{geo_id}.pkl",
    ]

    for path in candidates:
        if path.exists():
            return path

    raise FileNotFoundError(
        f"Could not locate downloaded pickle for {geo_id}. Expected one of: "
        + ", ".join(str(path) for path in candidates)
    )


def _ensure_female_column(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    cols_lower = {col.lower(): col for col in out.columns}
    if "female" in cols_lower:
        return out

    if "gender" in cols_lower:
        gender_col = cols_lower["gender"]
        out["female"] = out[gender_col].astype(str).str.upper().eq("F").astype(int)
    elif "sex" in cols_lower:
        sex_col = cols_lower["sex"]
        values = out[sex_col].astype(str).str.upper().str.strip()
        female_mask = values.isin({"F", "FEMALE", "XX", "1"})
        out["female"] = female_mask.astype(int)

    return out


def _aggregate_epicv2_probes(df: pd.DataFrame, verbose: bool) -> pd.DataFrame:
    if hasattr(pya, "pp") and hasattr(pya.pp, "epicv2_probe_aggregation"):
        return pya.pp.epicv2_probe_aggregation(df, verbose=verbose)
    return df


def _df_to_adata(df: pd.DataFrame, metadata_cols: list[str], imputer_strategy: str, verbose: bool):
    available_metadata_cols = [col for col in metadata_cols if col in df.columns]

    if hasattr(pya, "pp") and hasattr(pya.pp, "df_to_adata"):
        return pya.pp.df_to_adata(
            df,
            metadata_cols=available_metadata_cols,
            imputer_strategy=imputer_strategy,
            verbose=verbose,
        )

    if hasattr(pya, "preprocess") and hasattr(pya.preprocess, "df_to_adata"):
        return pya.preprocess.df_to_adata(
            df,
            metadata_cols=available_metadata_cols,
            imputer_strategy=imputer_strategy,
            verbose=verbose,
        )

    raise RuntimeError("PyAging does not expose df_to_adata in pp or preprocess namespace")


def _predict_age(adata, clocks: list[str], verbose: bool) -> None:
    pya.pred.predict_age(adata, clocks, verbose=verbose)


def _find_clock_columns(obs: pd.DataFrame, clocks: list[str]) -> list[str]:
    by_lower = {col.lower(): col for col in obs.columns}
    found = []
    for clock in clocks:
        col = by_lower.get(clock.lower())
        if col:
            found.append(col)
    return found


def _collect_missing_features(adata, clocks: list[str]) -> pd.DataFrame:
    rows: list[dict[str, str]] = []
    uns = getattr(adata, "uns", {})
    for clock in clocks:
        key = f"{clock.lower()}_missing_features"
        value = uns.get(key, [])
        if isinstance(value, (list, tuple, np.ndarray)):
            for feature in value:
                rows.append({"clock": clock, "feature": str(feature)})
        elif value:
            rows.append({"clock": clock, "feature": str(value)})

    return pd.DataFrame(rows, columns=["clock", "feature"])


def _collect_clock_metadata(adata, clocks: list[str]) -> dict[str, dict]:
    out: dict[str, dict] = {}
    uns = getattr(adata, "uns", {})
    for clock in clocks:
        key = f"{clock.lower()}_metadata"
        value = uns.get(key)
        if isinstance(value, dict):
            out[clock] = value
    return out


def plot_clock_distributions(predictions: pd.DataFrame, outpath: Path) -> None:
    numeric = predictions.select_dtypes(include=["number"])
    if numeric.empty:
        return
    plt.figure(figsize=(max(7, 1.5 * len(numeric.columns)), 5))
    plt.boxplot([numeric[col].dropna().values for col in numeric.columns], tick_labels=numeric.columns)
    plt.xticks(rotation=35, ha="right")
    plt.ylabel("Predicted value")
    plt.title("Methylation Clock Predictions")
    plt.tight_layout()
    plt.savefig(outpath, dpi=180)
    plt.close()


def plot_clock_correlation(predictions: pd.DataFrame, outpath: Path) -> None:
    numeric = predictions.select_dtypes(include=["number"])
    if numeric.shape[1] < 2:
        return

    corr = numeric.corr("pearson")
    fig, ax = plt.subplots(figsize=(6, 5))
    im = ax.imshow(corr.values, vmin=-1, vmax=1)
    ax.set_xticks(range(len(corr.columns)))
    ax.set_xticklabels(corr.columns, rotation=35, ha="right")
    ax.set_yticks(range(len(corr.index)))
    ax.set_yticklabels(corr.index)
    ax.set_title("Clock Correlation (Pearson)")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    plt.tight_layout()
    plt.savefig(outpath, dpi=180)
    plt.close()


def write_reproducibility(output_dir: Path, input_desc: str, clocks: list[str], command_args: list[str]) -> None:
    repro_dir = output_dir / "reproducibility"
    repro_dir.mkdir(parents=True, exist_ok=True)

    command_text = "python skills/methylation-clock/methylation_clock.py " + " ".join(command_args) + "\n"
    (repro_dir / "commands.sh").write_text(command_text)

    env_text = """name: clawbio-methylation-clock
channels:
  - conda-forge
dependencies:
  - python>=3.11
  - pandas
  - numpy
  - matplotlib
  - pyaging
"""
    (repro_dir / "environment.yml").write_text(env_text)

    checksums = []
    input_path = Path(input_desc)
    if input_path.exists():
        checksums.append(f"{_sha256(input_path)}  {input_path.name}")
    for path in sorted((output_dir / "tables").glob("*")):
        if path.is_file():
            checksums.append(f"{_sha256(path)}  tables/{path.name}")
    for path in sorted((output_dir / "figures").glob("*")):
        if path.is_file():
            checksums.append(f"{_sha256(path)}  figures/{path.name}")

    (repro_dir / "checksums.sha256").write_text("\n".join(checksums) + "\n")


def write_report(
    output_dir: Path,
    input_desc: str,
    clocks: list[str],
    predictions: pd.DataFrame,
    missing_df: pd.DataFrame,
    metadata_json: dict[str, dict],
) -> None:
    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")

    clock_cols = [col for col in predictions.columns if col != "sample_id"]
    summary_rows = []
    for clock in clock_cols:
        series = pd.to_numeric(predictions[clock], errors="coerce")
        summary_rows.append(
            f"| {clock} | {series.notna().sum()} | {series.mean():.3f} | {series.std(ddof=1):.3f} |"
        )

    missing_counts = (
        missing_df.groupby("clock").size().to_dict() if not missing_df.empty else {}
    )

    report = f"""# ClawBio Methylation Clock Report

**Date**: {now}
**Input**: `{input_desc}`
**Samples**: {predictions.shape[0]}
**Clocks requested**: {", ".join(clocks)}

## Outputs

- Predictions table: `tables/predictions.csv`
- Summary table: `tables/prediction_summary.csv`
- Missing-features table: `tables/missing_features.csv`
- Clock metadata: `tables/clock_metadata.json`
- Distribution figure: `figures/clock_distributions.png`
- Correlation figure: `figures/clock_correlation.png`

## Prediction Summary

| Clock | Non-missing | Mean | Std |
|---|---:|---:|---:|
{"\n".join(summary_rows) if summary_rows else "| (none) | 0 | n/a | n/a |"}

## Missing Features

{"No missing clock features were reported by PyAging." if not missing_counts else "\n".join([f"- {clock}: {count}" for clock, count in sorted(missing_counts.items())])}

## Reproducibility

- Commands: `reproducibility/commands.sh`
- Environment: `reproducibility/environment.yml`
- Checksums: `reproducibility/checksums.sha256`

## Disclaimer

{DISCLAIMER}
"""
    (output_dir / "report.md").write_text(report)


def run_analysis(
    output_dir: Path,
    clocks: list[str],
    metadata_cols: list[str],
    imputer_strategy: str,
    verbose: bool,
    input_path: Path | None = None,
    geo_id: str | None = None,
    skip_epicv2_aggregation: bool = False,
) -> dict:
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(
            f"Output directory '{output_dir}' is not empty. "
            "Choose a new --output path to avoid overwriting existing reports."
        )

    if (input_path is None) == (geo_id is None):
        raise ValueError("Provide exactly one of input_path or geo_id")

    if geo_id:
        resolved_input = _resolve_geo_pickle(geo_id, verbose=verbose)
    else:
        resolved_input = input_path
        if resolved_input is None or not resolved_input.exists():
            raise FileNotFoundError(f"Input file not found: {resolved_input}")

    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "tables").mkdir(exist_ok=True)
    (output_dir / "figures").mkdir(exist_ok=True)

    df = _load_dataframe(resolved_input)
    if not isinstance(df, pd.DataFrame) or df.empty:
        raise ValueError("Loaded methylation data is empty or invalid")

    df = _ensure_female_column(df)
    if not skip_epicv2_aggregation:
        df = _aggregate_epicv2_probes(df, verbose=verbose)

    adata = _df_to_adata(
        df=df,
        metadata_cols=metadata_cols,
        imputer_strategy=imputer_strategy,
        verbose=verbose,
    )
    _predict_age(adata=adata, clocks=clocks, verbose=verbose)

    obs = adata.obs.copy()
    obs = obs.reset_index().rename(columns={"index": "sample_id"})
    clock_cols = _find_clock_columns(obs, clocks)
    if not clock_cols:
        raise RuntimeError("No requested clock outputs were found in adata.obs")

    prediction_cols = ["sample_id"] + clock_cols
    predictions = obs[prediction_cols].copy()
    summary = predictions[clock_cols].apply(pd.to_numeric, errors="coerce").describe().T
    summary = summary.reset_index().rename(columns={"index": "clock"})

    missing_df = _collect_missing_features(adata, clocks)
    metadata_json = _collect_clock_metadata(adata, clocks)

    predictions.to_csv(output_dir / "tables" / "predictions.csv", index=False)
    summary.to_csv(output_dir / "tables" / "prediction_summary.csv", index=False)
    missing_df.to_csv(output_dir / "tables" / "missing_features.csv", index=False)
    (output_dir / "tables" / "clock_metadata.json").write_text(json.dumps(metadata_json, indent=2))

    plot_clock_distributions(predictions[clock_cols], output_dir / "figures" / "clock_distributions.png")
    plot_clock_correlation(predictions[clock_cols], output_dir / "figures" / "clock_correlation.png")

    write_report(
        output_dir=output_dir,
        input_desc=str(resolved_input),
        clocks=clocks,
        predictions=predictions,
        missing_df=missing_df,
        metadata_json=metadata_json,
    )

    return {
        "output_dir": str(output_dir),
        "input": str(resolved_input),
        "n_samples": int(predictions.shape[0]),
        "clocks_found": clock_cols,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Compute epigenetic age with PyAging methylation clocks")
    parser.add_argument("--input", help="Local methylation file (.pkl/.pickle/.csv/.tsv)")
    parser.add_argument("--geo-id", help="GEO accession to download via PyAging (for example: GSE139307)")
    parser.add_argument("--output", required=True, help="Output directory")
    parser.add_argument("--clocks", help="Comma-separated clocks; defaults to core 5 clocks")
    parser.add_argument(
        "--metadata-cols",
        default="gender,tissue_type,dataset",
        help="Comma-separated metadata columns to preserve in AnnData obs",
    )
    parser.add_argument(
        "--imputer-strategy",
        default="knn",
        choices=["knn", "mean", "median", "most_frequent"],
        help="Imputer strategy for pya.pp.df_to_adata",
    )
    parser.add_argument("--skip-epicv2-aggregation", action="store_true", help="Skip EPICv2 duplicate-probe aggregation")
    parser.add_argument("--demo", action="store_true", help="Use bundled small demo pickle")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose PyAging logging")
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    if args.demo:
        demo_path = _get_repo_root() / "pyaging_data" / "GSE139307_small.pkl"
        input_path = demo_path
        geo_id = None
    else:
        input_path = Path(args.input) if args.input else None
        geo_id = args.geo_id

    if (input_path is None and geo_id is None) or (input_path is not None and geo_id is not None):
        parser.error("Provide exactly one input source: --input or --geo-id (or use --demo)")

    clocks = parse_clock_list(args.clocks)
    metadata_cols = parse_metadata_cols(args.metadata_cols)

    result = run_analysis(
        output_dir=Path(args.output),
        clocks=clocks,
        metadata_cols=metadata_cols,
        imputer_strategy=args.imputer_strategy,
        verbose=args.verbose,
        input_path=input_path,
        geo_id=geo_id,
        skip_epicv2_aggregation=args.skip_epicv2_aggregation,
    )

    command_args = []
    if args.demo:
        command_args.extend(["--demo"])
    elif input_path is not None:
        command_args.extend(["--input", str(input_path)])
    else:
        command_args.extend(["--geo-id", str(geo_id)])
    command_args.extend(["--output", str(args.output), "--clocks", ",".join(clocks)])
    command_args.extend(["--metadata-cols", ",".join(metadata_cols)])
    command_args.extend(["--imputer-strategy", args.imputer_strategy])
    if args.skip_epicv2_aggregation:
        command_args.append("--skip-epicv2-aggregation")
    if args.verbose:
        command_args.append("--verbose")

    write_reproducibility(
        output_dir=Path(args.output),
        input_desc=result["input"],
        clocks=clocks,
        command_args=command_args,
    )

    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
