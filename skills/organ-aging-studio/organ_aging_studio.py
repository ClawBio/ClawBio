#!/usr/bin/env python3
"""Organ Aging Studio — interactive Goeminne proteomic clock with protein filters.

Wraps the published organAging coefficients (protein NPX × weight + intercept).
Designed for demos, web UIs, and hackathon show-and-tell.
"""

from __future__ import annotations

import argparse
import json
import shlex
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

_SKILL_DIR = Path(__file__).resolve().parent
_CLAWBIO_ROOT = _SKILL_DIR.parent.parent
_PROTEOMICS_CLOCK_DIR = _SKILL_DIR.parent / "proteomics-clock"

for _path in (_CLAWBIO_ROOT, _PROTEOMICS_CLOCK_DIR):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

from clawbio.common.reproducibility import (  # noqa: E402
    write_checksums,
    write_commands_sh,
    write_environment_yml,
)

from proteomics_clock import (  # noqa: E402
    DISCLAIMER,
    CITATION,
    assess_input_npx_scale,
    download_coefficients,
    download_organ_proteins,
    load_input,
    mortality_to_years,
    parse_organ_list,
)

DEMO_ORGANS = ["Heart", "Brain", "Liver", "Immune", "Organismal"]
DEMO_DATA = _PROTEOMICS_CLOCK_DIR / "data" / "demo_olink_npx.csv.gz"


def filter_protein_coefs(
    protein_coefs: dict[str, float],
    available_columns: set[str],
    *,
    top_n: int | None = None,
    min_abs_coef: float = 0.0,
) -> dict[str, float]:
    """Keep coefficients for proteins present in the sample, ranked by |coef|."""
    filtered = {
        protein: coef
        for protein, coef in protein_coefs.items()
        if protein in available_columns and abs(coef) >= min_abs_coef
    }
    if top_n is not None and top_n > 0:
        ranked = sorted(filtered.items(), key=lambda item: abs(item[1]), reverse=True)
        filtered = dict(ranked[:top_n])
    return filtered


def compute_organ_age(
    row: pd.Series,
    coefs: dict[str, float],
    *,
    generation: str,
    top_n: int | None = None,
    min_abs_coef: float = 0.0,
) -> dict[str, Any]:
    """Return predicted age plus per-protein contribution breakdown."""
    intercept = coefs.get("Intercept", 0.0) if generation == "gen1" else 0.0
    protein_coefs = {k: v for k, v in coefs.items() if k != "Intercept"}
    available = set(row.index.astype(str))
    active = filter_protein_coefs(
        protein_coefs, available, top_n=top_n, min_abs_coef=min_abs_coef
    )

    contributions: list[dict[str, Any]] = []
    total = float(intercept)
    for protein, coef in sorted(active.items(), key=lambda item: abs(item[1]), reverse=True):
        npx = float(row[protein])
        contrib = npx * coef
        total += contrib
        contributions.append(
            {
                "protein": protein,
                "npx": round(npx, 4),
                "coefficient": round(coef, 6),
                "contribution": round(contrib, 4),
            }
        )

    predicted = total
    if generation == "gen2":
        predicted = float(mortality_to_years(np.array([total]))[0])

    chrono = row.get("age")
    chrono_age = float(chrono) if chrono is not None and not pd.isna(chrono) else None
    delta = round(predicted - chrono_age, 2) if chrono_age is not None else None

    return {
        "predicted_age_years": round(predicted, 2),
        "chronological_age_years": chrono_age,
        "age_delta_years": delta,
        "raw_delta_years": delta,
        "intercept": round(intercept, 4),
        "n_proteins_used": len(contributions),
        "n_proteins_available": len(active),
        "n_proteins_in_model": len(protein_coefs),
        "is_partial_prediction": len(active) < len(protein_coefs),
        "contributions": contributions,
    }


def run_studio(
    *,
    input_path: Path,
    output_dir: Path,
    organs: list[str],
    generation: str = "gen1",
    fold: int = 1,
    sample_id: str | None = None,
    top_n: int | None = None,
    min_abs_coef: float = 0.0,
    demo: bool = False,
) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    df = load_input(input_path)
    if "sample_id" not in df.columns:
        df = df.copy()
        df.insert(0, "sample_id", [f"sample_{i}" for i in range(len(df))])

    if sample_id is not None:
        subset = df[df["sample_id"].astype(str) == str(sample_id)]
        if subset.empty:
            raise ValueError(f"sample_id '{sample_id}' not found in input")
        df = subset

    input_warnings = assess_input_npx_scale(df)

    coefficients = download_coefficients(organs, generation, fold)
    organ_protein_map = download_organ_proteins()

    results: dict[str, Any] = {
        "input": str(input_path),
        "generation": generation,
        "fold": fold,
        "organs": organs,
        "filters": {"top_n": top_n, "min_abs_coef": min_abs_coef},
        "input_sanity_warnings": input_warnings,
        "samples": [],
    }

    tables_dir = output_dir / "tables"
    tables_dir.mkdir(exist_ok=True)
    all_contributions: list[pd.DataFrame] = []

    for _, row in df.iterrows():
        sid = str(row["sample_id"])
        sample_result: dict[str, Any] = {"sample_id": sid, "organs": {}}
        for organ in organs:
            key = f"{organ}_{generation}"
            coefs = coefficients.get(key)
            if coefs is None:
                continue
            organ_result = compute_organ_age(
                row,
                coefs,
                generation=generation,
                top_n=top_n,
                min_abs_coef=min_abs_coef,
            )
            organ_result["organ"] = organ
            organ_result["model_proteins_total"] = len(
                [p for p in organ_protein_map.get(organ, []) if p in coefs]
            )
            sample_result["organs"][organ] = organ_result

            contrib_df = pd.DataFrame(organ_result["contributions"])
            if not contrib_df.empty:
                contrib_df.insert(0, "organ", organ)
                contrib_df.insert(0, "sample_id", sid)
                all_contributions.append(contrib_df)

        results["samples"].append(sample_result)

    if all_contributions:
        pd.concat(all_contributions, ignore_index=True).to_csv(
            tables_dir / "protein_contributions.csv", index=False
        )

    (output_dir / "result.json").write_text(json.dumps(results, indent=2))

    report = build_report(results)
    (output_dir / "report.md").write_text(report)

    cmd = ["python", "skills/organ-aging-studio/organ_aging_studio.py"]
    cmd += ["--demo"] if demo else ["--input", str(input_path)]
    cmd += [
        "--output", str(output_dir),
        "--organs", ",".join(organs),
        "--generation", generation,
        "--fold", str(fold),
        "--min-abs-coef", str(min_abs_coef),
    ]
    if sample_id is not None:
        cmd += ["--sample-id", str(sample_id)]
    if top_n is not None:
        cmd += ["--top-n", str(top_n)]
    write_commands_sh(output_dir, shlex.join(cmd))
    env_path = write_environment_yml(
        output_dir,
        "clawbio-organ-aging-studio",
        pip_deps=["pandas>=2.0", "numpy>=1.24", "requests>=2.28", "matplotlib"],
    )
    write_checksums(
        [
            output_dir / "report.md",
            output_dir / "result.json",
            tables_dir / "protein_contributions.csv",
            env_path,
        ],
        output_dir,
        anchor=output_dir,
    )

    return results


def build_report(results: dict[str, Any]) -> str:
    filters = results.get("filters", {})
    filters_active = bool(filters.get("top_n") or (filters.get("min_abs_coef") or 0.0) > 0)
    input_warnings = results.get("input_sanity_warnings", [])
    lines = [
        "# Organ Aging Studio Report",
        "",
        f"**Generated**: {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}",
        f"**Model**: Goeminne et al. (2025) organ-specific proteomic clocks",
        f"**Generation**: {results['generation']}",
        "",
        CITATION,
        "",
        "## Summary",
        "",
    ]

    if filters_active:
        lines.extend(
            [
                "> Interpretation note: `--top-n` and/or `--min-abs-coef` were used, so these are filtered partial-model predictions rather than the validated full published clock.",
                "",
            ]
        )

    for sample in results["samples"]:
        sid = sample["sample_id"]
        lines.append(f"### Sample `{sid}`")
        lines.append("")
        lines.append("| Organ | Predicted age | Chronological | Raw delta | Proteins used |")
        lines.append("|-------|---------------|---------------|-----------|---------------|")
        for organ, data in sample["organs"].items():
            chrono = data["chronological_age_years"]
            chrono_s = f"{chrono:.1f}" if chrono is not None else "—"
            delta = data["age_delta_years"]
            delta_s = f"{delta:+.1f}" if delta is not None else "—"
            lines.append(
                f"| {organ} | {data['predicted_age_years']:.1f} | {chrono_s} | {delta_s} | "
                f"{data['n_proteins_used']}/{data['n_proteins_in_model']} |"
            )
        lines.append("")

    lines.extend(
        [
            "## Input Sanity",
            "",
        ]
    )

    if input_warnings:
        lines.extend([f"- {warning}" for warning in input_warnings])
    else:
        lines.append("- No obvious unit issues detected.")

    lines.extend(
        [
            "## How prediction works",
            "",
            "```",
            "predicted_age = intercept + Σ (protein_NPX × coefficient)",
            "```",
            "",
            "Coefficients are downloaded from the pinned [organAging](https://github.com/ludgergoeminne/organAging) repository.",
            "Delta values are the raw predicted-minus-chronological gap, not age-residualised acceleration.",
            "",
            DISCLAIMER,
        ]
    )
    return "\n".join(lines)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Organ Aging Studio — proteomic clock with protein-level breakdown"
    )
    parser.add_argument("--input", help="Olink NPX table (.csv/.tsv/.gz)")
    parser.add_argument("--output", required=True, help="Output directory")
    parser.add_argument(
        "--organs",
        default=",".join(DEMO_ORGANS),
        help=f"Comma-separated organs (default: {','.join(DEMO_ORGANS)})",
    )
    parser.add_argument(
        "--generation",
        default="gen1",
        choices=["gen1", "gen2"],
        help="gen1 = chronological age; gen2 = mortality hazard → years",
    )
    parser.add_argument("--fold", type=int, default=1, help="CV fold 1-5")
    parser.add_argument("--sample-id", help="Analyse one sample only")
    parser.add_argument(
        "--top-n",
        type=int,
        help="Use only the top N proteins by |coefficient| (among those in the sample)",
    )
    parser.add_argument(
        "--min-abs-coef",
        type=float,
        default=0.0,
        help="Drop proteins with |coefficient| below this threshold",
    )
    parser.add_argument("--demo", action="store_true", help="Use bundled Olink demo data")
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    if not args.demo and not args.input:
        parser.error("Provide --input or use --demo")
    input_path = DEMO_DATA if args.demo else Path(args.input)

    organs = parse_organ_list(args.organs)
    result = run_studio(
        input_path=input_path,
        output_dir=Path(args.output),
        organs=organs,
        generation=args.generation,
        fold=args.fold,
        sample_id=args.sample_id,
        top_n=args.top_n,
        min_abs_coef=args.min_abs_coef,
        demo=args.demo,
    )
    n_samples = len(result["samples"])
    print(f"✓ Organ Aging Studio complete — {n_samples} sample(s) → {args.output}")


if __name__ == "__main__":
    main()
