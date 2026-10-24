#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import re
import sys
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version as package_version
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
try:
    import pyaging as pya
except ModuleNotFoundError as exc:
    if exc.name != "pyaging":
        raise
    pya = None

_PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from clawbio.common.checksums import sha256_file as _sha256  # noqa: E402
from clawbio.common.reproducibility import (  # noqa: E402
    ReproCommand,
    ReproPath,
    write_checksums,
    write_environment_yml,
    write_portable_commands_sh,
)


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
PROBE_NAME = re.compile(r"^(?:cg|ch|rs|nv)[0-9A-Za-z.\-]*(?:_|$)")


def parse_clock_list(clocks: str | None) -> list[str]:
    if not clocks:
        return list(DEFAULT_CLOCKS)
    values = [item.strip() for item in clocks.split(",") if item.strip()]
    if not values:
        raise ValueError("Clock list is empty")
    return values


def parse_metadata_cols(metadata_cols: str | None) -> list[str]:
    if not metadata_cols:
        return ["gender", "sex", "tissue_type", "dataset"]
    return [item.strip() for item in metadata_cols.split(",") if item.strip()]


def _load_dataframe(path: Path) -> pd.DataFrame:
    suffixes = [s.lower() for s in path.suffixes]
    if suffixes and suffixes[-1] in {".pkl", ".pickle"}:
        import warnings
        warnings.warn(
            f"Loading pickle file {path.name}. Untrusted pickles can execute arbitrary code.",
            stacklevel=2,
        )
        df = pd.read_pickle(path)
    elif suffixes[-2:] == [".csv", ".gz"]:
        df = pd.read_csv(path, compression="gzip")
    elif suffixes[-2:] == [".tsv", ".gz"]:
        df = pd.read_csv(path, sep="\t", compression="gzip")
    elif suffixes and suffixes[-1] == ".csv":
        df = pd.read_csv(path)
    elif suffixes and suffixes[-1] == ".tsv":
        df = pd.read_csv(path, sep="\t")
    else:
        raise ValueError("Unsupported input format. Use .pkl, .pickle, .csv, .tsv, .csv.gz, or .tsv.gz")

    if "sample_id" in df.columns:
        sample_ids = df["sample_id"]
        if (sample_ids.isna().any()
                or sample_ids.astype(str).str.strip().eq("").any()
                or not sample_ids.is_unique):
            raise ValueError("sample_id values must be present and unique")
        df = df.set_index("sample_id")
    if (df.index.hasnans
            or (df.index.astype(str).str.strip() == "").any()
            or not df.index.is_unique):
        raise ValueError("Sample identifiers must be present and unique")
    return df


def _validate_methylation_frame(df: pd.DataFrame, metadata_cols: list[str]) -> dict:
    """Require numeric beta values for at least one recognized methylation probe."""
    probe_cols = [col for col in df.columns if PROBE_NAME.match(str(col))]
    if not probe_cols:
        raise ValueError("No methylation probe columns found (expected cg/ch/rs/nv IDs)")
    probe_set = set(probe_cols)
    metadata_set = set(metadata_cols)
    metadata_probes = probe_set & metadata_set
    if metadata_probes:
        raise ValueError(
            f"--metadata-cols cannot include methylation probe columns: {sorted(metadata_probes)[:3]}"
        )
    nonnumeric = [
        col for col, dtype in df.dtypes.items()
        if col in probe_set and not pd.api.types.is_numeric_dtype(dtype)
    ]
    if nonnumeric:
        raise ValueError(f"Methylation probe columns must be numeric: {nonnumeric[:3]}")
    unknown_metadata = [
        col for col, dtype in df.dtypes.items()
        if col not in probe_set and col not in metadata_set
        and not pd.api.types.is_numeric_dtype(dtype)
    ]
    if unknown_metadata:
        raise ValueError(
            f"Nonnumeric columns must be listed in --metadata-cols: {unknown_metadata[:3]}"
        )
    values = df.loc[:, probe_cols].to_numpy(dtype=float, copy=False)
    if np.isinf(values).any() or ((values < 0) | (values > 1)).any():
        raise ValueError(
            "Methylation beta values outside [0, 1] or infinite; "
            "if your file contains percentages on [0, 100], divide by 100 before input"
        )
    observed = values[~np.isnan(values)]
    if observed.size == 0:
        raise ValueError("All methylation probe beta values are missing")
    return {
        "n_samples": int(len(df)),
        "n_probe_columns": int(len(probe_cols)),
        "n_observed_values": int(observed.size),
        "n_missing_values": int(values.size - observed.size),
        "beta_min": float(observed.min()),
        "beta_max": float(observed.max()),
    }


def _validate_clock_covariates(
    df: pd.DataFrame, clocks: list[str], metadata_cols: list[str] | None = None
) -> None:
    """Require PyAging's additional model inputs for GrimAge-family defaults."""
    if not {clock.lower() for clock in clocks} & {"pcgrimage", "grimage2"}:
        return
    hidden_covariates = {"age", "female"} & set(metadata_cols or [])
    if hidden_covariates:
        raise ValueError(
            f"--metadata-cols cannot include required model input columns: {sorted(hidden_covariates)}"
        )
    if "age" not in df.columns:
        raise ValueError("PCGrimAge and GrimAge2 require an age column in years")
    if not pd.api.types.is_numeric_dtype(df["age"]):
        raise ValueError("age must contain finite numeric years for PCGrimAge and GrimAge2")
    age = df["age"].to_numpy(dtype=float)
    if not np.isfinite(age).all() or (age < 0).any():
        raise ValueError("age must contain finite nonnegative years for PCGrimAge and GrimAge2")
    if "female" not in df.columns:
        raise ValueError("PCGrimAge and GrimAge2 require female (1 or 0), or known gender/sex")
    if not pd.api.types.is_numeric_dtype(df["female"]):
        raise ValueError("female must contain 1 or 0 for PCGrimAge and GrimAge2")
    female = df["female"].to_numpy(dtype=float)
    if not np.isin(female, [0, 1]).all():
        raise ValueError("female must contain 1 or 0 for every sample in PCGrimAge and GrimAge2")


def _get_repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _get_skill_data_dir() -> Path:
    return Path(__file__).resolve().parent / "data"


def _resolve_geo_pickle(geo_id: str, verbose: bool, output_dir: Path) -> Path:
    input_dir = output_dir / "input"
    input_dir.mkdir(parents=True, exist_ok=True)
    downloaded = Path(pya.data.download_example_data(geo_id, dir=str(input_dir), verbose=verbose))
    if not downloaded.is_file():
        raise FileNotFoundError(f"PyAging did not create the GEO input file: {downloaded}")
    return downloaded


def _ensure_female_column(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    cols_lower = {col.lower(): col for col in out.columns}
    if "female" in cols_lower:
        return out

    if "gender" in cols_lower:
        gender_col = cols_lower["gender"]
        values = out[gender_col].astype(str).str.upper().str.strip()
    elif "sex" in cols_lower:
        sex_col = cols_lower["sex"]
        values = out[sex_col].astype(str).str.upper().str.strip()
    else:
        return out

    female_mask = values.isin({"F", "FEMALE", "XX", "1"})
    male_mask = values.isin({"M", "MALE", "XY", "0"})
    out["female"] = np.where(female_mask, 1.0, np.where(male_mask, 0.0, np.nan))

    return out


def _aggregate_epicv2_probes(df: pd.DataFrame, verbose: bool) -> pd.DataFrame:
    if not (hasattr(pya, "pp") and hasattr(pya.pp, "epicv2_probe_aggregation")):
        return df
    # Newer PyAging (>=0.5) strips every column name at the first "_", which
    # mangles metadata columns such as "tissue_type". Aggregate probe columns only.
    is_probe = df.columns.str.match(r"^(cg|ch|rs|nv)[0-9A-Za-z.\-]*(_|$)")
    if is_probe.all():
        return pya.pp.epicv2_probe_aggregation(df, verbose=verbose)
    probes = pya.pp.epicv2_probe_aggregation(df.loc[:, is_probe], verbose=verbose)
    return pd.concat([df.loc[:, ~is_probe], probes], axis=1)


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


def plot_clock_distributions(
    predictions: pd.DataFrame,
    outpath: Path,
    clock_metadata: dict[str, dict] | None = None,
) -> None:
    numeric = predictions.select_dtypes(include=["number"])
    if numeric.empty:
        return
    metadata = {name.lower(): value for name, value in (clock_metadata or {}).items()}
    ncols = min(2, len(numeric.columns))
    nrows = (len(numeric.columns) + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 3.5 * nrows), squeeze=False)
    for ax, clock in zip(axes.flat, numeric.columns):
        values = numeric[clock].dropna().to_numpy()
        if len(values) >= 5:
            ax.boxplot(values, positions=[1], widths=0.35)
        ax.scatter(np.ones(len(values)), values, alpha=0.65, s=28, zorder=3)
        unit = metadata.get(clock.lower(), {}).get("unit", "predicted value")
        if isinstance(unit, list):
            unit = ", ".join(map(str, unit)) or "predicted value"
        ax.set_title(clock)
        ax.set_ylabel(str(unit))
        ax.set_xlim(0.5, 1.5)
        ax.set_xticks([])
        ax.grid(axis="y", alpha=0.2)
    for ax in list(axes.flat)[len(numeric.columns):]:
        ax.set_visible(False)
    fig.suptitle("Methylation clock predictions by clock")
    fig.tight_layout()
    fig.savefig(outpath, dpi=180)
    plt.close(fig)


def plot_clock_correlation(predictions: pd.DataFrame, outpath: Path) -> None:
    numeric = predictions.select_dtypes(include=["number"])
    if numeric.shape[1] < 2:
        return
    complete = numeric.dropna()
    if len(complete) < 10:
        return
    corr = complete.corr("pearson")
    fig, ax = plt.subplots(figsize=(6, 5))
    im = ax.imshow(corr.values, vmin=-1, vmax=1)
    ax.set_xticks(range(len(corr.columns)))
    ax.set_xticklabels(corr.columns, rotation=35, ha="right")
    ax.set_yticks(range(len(corr.index)))
    ax.set_yticklabels(corr.index)
    ax.set_title(f"Exploratory clock correlation (Pearson, n={len(complete)})")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    plt.tight_layout()
    plt.savefig(outpath, dpi=180)
    plt.close()


def repro_command_for_bundle(
    output_dir: Path,
    *,
    input_path: "Path | None",
    geo_id: "str | None",
    demo: bool,
    clocks: list[str],
    metadata_cols: list[str],
    imputer_strategy: str,
    skip_epicv2_aggregation: bool,
    verbose: bool = False,
) -> ReproCommand:
    args: list = []
    if demo:
        args.append("--demo")
    elif input_path is not None:
        args += ["--input", ReproPath(input_path, anchor="auto")]
    else:
        args += ["--geo-id", geo_id]
    args += [
        "--output", ReproPath(output_dir / "replay", anchor="output_dir"),
        "--clocks", ",".join(clocks),
        "--metadata-cols", ",".join(metadata_cols),
        "--imputer-strategy", imputer_strategy,
    ]
    if skip_epicv2_aggregation:
        args.append("--skip-epicv2-aggregation")
    if verbose:
        args.append("--verbose")
    return ReproCommand(
        script_path=Path("skills/methylation-clock/methylation_clock.py"),
        args=args,
        comment="Replay this ClawBio methylation-clock run",
    )


def write_reproducibility(
    output_dir: Path,
    *,
    input_path: "Path | None",
    geo_id: "str | None",
    demo: bool,
    clocks: list[str],
    metadata_cols: list[str],
    imputer_strategy: str,
    skip_epicv2_aggregation: bool,
    verbose: bool = False,
) -> None:
    repro_dir = output_dir / "reproducibility"
    repro_dir.mkdir(parents=True, exist_ok=True)

    pip_deps = []
    for name in ("pandas", "numpy", "matplotlib", "pyaging"):
        try:
            pip_deps.append(f"{name}=={package_version(name)}")
        except PackageNotFoundError:
            pip_deps.append(name)
    write_environment_yml(
        output_dir,
        env_name="clawbio-methylation-clock",
        pip_deps=pip_deps,
        python_version=f"{sys.version_info.major}.{sys.version_info.minor}",
    )

    write_portable_commands_sh(
        output_dir,
        repro_command_for_bundle(
            output_dir,
            input_path=input_path,
            geo_id=geo_id,
            demo=demo,
            clocks=clocks,
            metadata_cols=metadata_cols,
            imputer_strategy=imputer_strategy,
            skip_epicv2_aggregation=skip_epicv2_aggregation,
            verbose=verbose,
        ),
        repo_root=_PROJECT_ROOT,
    )

    artifacts = [output_dir / "report.md", output_dir / "result.json"]
    artifacts.extend(sorted((output_dir / "tables").glob("*")))
    artifacts.extend(sorted((output_dir / "figures").glob("*")))
    write_checksums(artifacts, output_dir, anchor=output_dir)


def write_report(
    output_dir: Path,
    input_desc: str,
    clocks: list[str],
    predictions: pd.DataFrame,
    missing_df: pd.DataFrame,
    metadata_json: dict[str, dict],
    input_qc: dict | None = None,
    synthetic_demo: bool = False,
) -> None:
    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")

    clock_cols = [col for col in predictions.columns if col != "sample_id"]
    summary_rows = []
    metadata_by_lower = {name.lower(): value for name, value in metadata_json.items()}
    for clock in clock_cols:
        series = pd.to_numeric(predictions[clock], errors="coerce")
        unit = metadata_by_lower.get(clock.lower(), {}).get("unit", "unspecified")
        if isinstance(unit, list):
            unit = ", ".join(str(item) for item in unit) or "unspecified"
        summary_rows.append(
            f"| {clock} | {unit} | {series.notna().sum()} | {series.mean():.3f} | {series.std(ddof=1):.3f} |"
        )

    missing_counts = (
        missing_df.groupby("clock").size().to_dict() if not missing_df.empty else {}
    )

    summary_block = "\n".join(summary_rows) if summary_rows else "| (none) | 0 | n/a | n/a |"
    missing_block = (
        "No missing clock features were reported by PyAging."
        if not missing_counts
        else "\n".join([f"- {clock}: {count}" for clock, count in sorted(missing_counts.items())])
    )
    figure_lines = [
        f"- {label}: `figures/{filename}`"
        for label, filename in (
            ("Distribution figure", "clock_distributions.png"),
            ("Correlation figure", "clock_correlation.png"),
        )
        if (output_dir / "figures" / filename).exists()
    ]
    figures_block = "\n".join(figure_lines) if figure_lines else "- No figures generated"
    qc_block = ""
    qc_output_line = ""
    if input_qc is not None:
        qc_output_line = "- Input QC: `tables/input_qc.json`"
        qc_block = f"""## Input QC

- Methylation probe columns before aggregation: {input_qc['n_probe_columns']}
- Methylation probe columns after aggregation: {input_qc['n_probe_columns_after_aggregation']}
- Observed beta values: {input_qc['n_observed_values']}
- Missing beta values: {input_qc['n_missing_values']}
- Observed beta range: {input_qc['beta_min']:.3f} to {input_qc['beta_max']:.3f}

"""

    demo_block = (
        "## Synthetic demo\n\nThis input contains simulated beta values with no biological "
        "age or sex signal. Predictions demonstrate software operation only.\n\n"
        if synthetic_demo else ""
    )

    report = f"""# ClawBio Methylation Clock Report

**Date**: {now}
**Input**: `{input_desc}`
**Samples**: {predictions.shape[0]}
**Clocks requested**: {", ".join(clocks)}

{demo_block}## Outputs

- Predictions table: `tables/predictions.csv`
- Summary table: `tables/prediction_summary.csv`
- Missing-features table: `tables/missing_features.csv`
- Clock metadata: `tables/clock_metadata.json`
{qc_output_line}
{figures_block}

{qc_block}## Prediction Summary

| Clock | Unit | Non-missing | Mean | Std |
|---|---|---:|---:|---:|
{summary_block}

Clock outputs can have different units and targets; compare them using each clock's metadata.

## Missing Features

{missing_block}

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
    demo: bool = False,
) -> dict:
    if pya is None:
        raise RuntimeError("PyAging is required for methylation-clock; install it with 'pip install pyaging'")
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(
            f"Output directory '{output_dir}' is not empty. "
            "Choose a new --output path to avoid overwriting existing reports."
        )

    if (input_path is None) == (geo_id is None):
        raise ValueError("Provide exactly one of input_path or geo_id")

    if geo_id:
        resolved_input = _resolve_geo_pickle(geo_id, verbose=verbose, output_dir=output_dir)
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
    input_qc = _validate_methylation_frame(df, metadata_cols)
    _validate_clock_covariates(df, clocks, metadata_cols)
    if not skip_epicv2_aggregation:
        df = _aggregate_epicv2_probes(df, verbose=verbose)
    input_qc["n_probe_columns_after_aggregation"] = sum(
        bool(PROBE_NAME.match(str(col))) for col in df.columns
    )

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
    (output_dir / "tables" / "input_qc.json").write_text(json.dumps(input_qc, indent=2))

    plot_clock_distributions(
        predictions[clock_cols], output_dir / "figures" / "clock_distributions.png", metadata_json
    )
    plot_clock_correlation(predictions[clock_cols], output_dir / "figures" / "clock_correlation.png")

    write_report(
        output_dir=output_dir,
        input_desc=str(resolved_input),
        clocks=clocks,
        predictions=predictions,
        missing_df=missing_df,
        metadata_json=metadata_json,
        input_qc=input_qc,
        synthetic_demo=demo,
    )

    result = {
        "output_dir": str(output_dir),
        "input": str(resolved_input),
        "n_samples": int(predictions.shape[0]),
        "clocks_found": clock_cols,
        "input_qc": input_qc,
        "synthetic_demo": demo,
    }
    (output_dir / "result.json").write_text(json.dumps(result, indent=2))
    return result


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Compute epigenetic age with PyAging methylation clocks")
    parser.add_argument("--input", help="Local methylation file (.pkl/.pickle/.csv/.tsv/.csv.gz/.tsv.gz)")
    parser.add_argument("--geo-id", help="GEO accession to download via PyAging (for example: GSE139307)")
    parser.add_argument("--output", required=True, help="Output directory")
    parser.add_argument("--clocks", help="Comma-separated clocks; defaults to core 5 clocks")
    parser.add_argument(
        "--metadata-cols",
        default="gender,sex,tissue_type,dataset",
        help="Comma-separated metadata columns to preserve in AnnData obs",
    )
    parser.add_argument(
        "--imputer-strategy",
        default="knn",
        choices=["knn", "mean", "median", "most_frequent"],
        help="Imputer strategy for pya.pp.df_to_adata",
    )
    parser.add_argument("--skip-epicv2-aggregation", action="store_true", help="Skip EPICv2 duplicate-probe aggregation")
    parser.add_argument("--demo", action="store_true", help="Use bundled two-sample synthetic beta-value demo")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose PyAging logging")
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    if args.demo:
        demo_path = _get_skill_data_dir() / "synthetic_beta_demo.csv.gz"
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
        demo=args.demo,
    )

    write_reproducibility(
        output_dir=Path(args.output),
        input_path=input_path,
        geo_id=geo_id,
        demo=args.demo,
        clocks=clocks,
        metadata_cols=metadata_cols,
        imputer_strategy=args.imputer_strategy,
        skip_epicv2_aggregation=args.skip_epicv2_aggregation,
        verbose=args.verbose,
    )

    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
