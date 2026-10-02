#!/usr/bin/env python3
"""
Ontology Annotator — map messy metadata column values (tissue, cell type,
disease, trait) to standard ontology IDs (UBERON, CL, MONDO, EFO) via the
EBI OLS4 search API.

Usage:
    python ontology_annotator.py --input samples.csv --output report/
    python ontology_annotator.py --input samples.csv --output report/ \
        --columns tissue:uberon,cell_type:cl,disease:mondo,trait:efo --min-similarity 0.75
    python ontology_annotator.py --input adata.h5ad --output report/
    python ontology_annotator.py --demo --output /tmp/ontology_demo
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

_PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

SKILL_DIR = Path(__file__).resolve().parent
if str(SKILL_DIR) not in sys.path:
    sys.path.insert(0, str(SKILL_DIR))

from clawbio.common.report import write_result_json  # noqa: E402
from clawbio.common.reproducibility import (  # noqa: E402
    ReproCommand,
    ReproPath,
    write_checksums,
    write_environment_yml,
    write_portable_commands_sh,
)

from ontology_annotator_core.io_readers import (  # noqa: E402
    infer_columns,
    load_table,
    parse_columns_arg,
)
from ontology_annotator_core.ols4_client import OLS4Client, load_fixtures  # noqa: E402
from ontology_annotator_core.report import generate_markdown  # noqa: E402
from ontology_annotator_core.scoring import best_candidate, rank_candidates  # noqa: E402

DISCLAIMER = (
    "ClawBio is a research and educational tool. It is not a medical device "
    "and does not provide clinical diagnoses. Consult a healthcare "
    "professional before making any medical decisions."
)

DEMO_INPUT_PATH = SKILL_DIR / "examples" / "demo_input.csv"
DEMO_FIXTURES_PATH = SKILL_DIR / "data" / "ols4_fixtures.json"
DEFAULT_CACHE_DIR = Path.home() / ".clawbio" / "ontology_annotator_cache"
DEFAULT_MIN_SIMILARITY = 0.8
CATALOG_PATH = _PROJECT_ROOT / "skills" / "catalog.json"


# ---------------------------------------------------------------------------
# Core pipeline
# ---------------------------------------------------------------------------


def annotate_table(
    df,
    columns: dict[str, str],
    client: OLS4Client,
    min_similarity: float,
) -> tuple[Any, dict[str, dict[str, int]], list[dict]]:
    """Annotate `df` in place-ish (returns a copy) for every (col, ontology)
    in `columns`. Caches OLS4 lookups per unique normalised value within a
    single run via the client's own disk/fixture cache, plus an in-process
    memo so repeated values in the same table never re-query at all.

    Returns (annotated_df, per_column_summary, flagged_rows).
    """
    df = df.copy()
    flagged_rows: list[dict] = []
    summary: dict[str, dict[str, int]] = {}

    for col, ontology in columns.items():
        if col not in df.columns:
            print(f"  WARNING: column '{col}' not found in input; skipping", file=sys.stderr)
            continue

        memo: dict[str, list[dict]] = {}
        ids, labels, similarities, flags, candidates_json = [], [], [], [], []
        n_matched = n_flagged = n_no_candidates = 0

        print(f"  Annotating '{col}' -> {ontology.upper()} ({df[col].nunique(dropna=True)} unique values)")

        for row_idx, raw_value in df[col].items():
            if raw_value is None or (isinstance(raw_value, float) and raw_value != raw_value):
                ids.append("")
                labels.append("")
                similarities.append(0.0)
                flags.append(True)
                candidates_json.append("[]")
                n_no_candidates += 1
                continue

            value = str(raw_value)
            memo_key = f"{ontology}:{value.strip().lower()}"
            if memo_key not in memo:
                result = client.search(value, ontology)
                memo[memo_key] = rank_candidates(value, result["docs"], top_n=3)

            candidates = memo[memo_key]
            top, is_flagged = best_candidate(candidates, min_similarity)

            ids.append(top["ontology_id"] if top else "")
            labels.append(top["label"] if top else "")
            similarities.append(top["string_similarity"] if top else 0.0)
            flags.append(is_flagged)
            candidates_json.append(json.dumps(candidates))

            if not candidates:
                n_no_candidates += 1
                reason = "No OLS4 candidates found"
            elif is_flagged:
                n_flagged += 1
                reason = f"Top string similarity {top['string_similarity']} < min {min_similarity}; review"
            else:
                n_matched += 1
                reason = None

            if is_flagged:
                flagged_rows.append(
                    {
                        "row_index": row_idx,
                        "column": col,
                        "value": value,
                        "top_label": top["label"] if top else "—",
                        "top_string_similarity": top["string_similarity"] if top else "—",
                        "reason": reason,
                    }
                )

        df[f"{col}_ontology_id"] = ids
        df[f"{col}_label"] = labels
        df[f"{col}_string_similarity"] = similarities
        df[f"{col}_needs_review"] = flags
        df[f"{col}_candidates"] = candidates_json

        summary[col] = {
            "n_matched": n_matched,
            "n_flagged": n_flagged,
            "n_no_candidates": n_no_candidates,
        }

    return df, summary, flagged_rows


def run_annotation(
    input_path: Path,
    output_dir: Path,
    columns_arg: str | None,
    min_similarity: float,
    cache_dir: Path = DEFAULT_CACHE_DIR,
    use_cache: bool = True,
    fixtures_path: Path | None = None,
    demo: bool = False,
) -> dict:
    """Run the full ontology-annotator pipeline. Returns the result envelope dict."""
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"  Loading input: {input_path}")
    df = load_table(input_path)
    print(f"  Loaded {len(df)} rows, {len(df.columns)} columns")

    if columns_arg:
        columns = parse_columns_arg(columns_arg)
    else:
        columns = infer_columns(df)
        if not columns:
            raise ValueError(
                "Could not auto-detect columns to annotate. Pass --columns "
                "col:ontology,col2:ontology2 (e.g. tissue:uberon,disease:mondo)."
            )
        print(f"  Auto-detected columns: {columns}")

    fixtures = None
    if fixtures_path is not None:
        print(f"  Offline mode: loading fixtures from {fixtures_path.name}")
        fixtures = load_fixtures(fixtures_path)

    client = OLS4Client(cache_dir=cache_dir, use_cache=use_cache, fixtures=fixtures)

    print(f"  Annotating {len(columns)} column(s) against OLS4 (min_similarity={min_similarity})...")
    annotated_df, summary, flagged_rows = annotate_table(df, columns, client, min_similarity)

    annotated_path = output_dir / "annotated.csv"
    annotated_df.to_csv(annotated_path, index=False)
    print(f"  Wrote {annotated_path}")

    print("  Writing report...")
    report_md = generate_markdown(
        annotated_columns=columns,
        row_summaries=summary,
        flagged_rows=flagged_rows,
        n_rows=len(df),
        min_similarity=min_similarity,
        catalog_path=CATALOG_PATH,
    )
    (output_dir / "report.md").write_text(report_md)

    total_matched = sum(s["n_matched"] for s in summary.values())
    total_flagged = sum(s["n_flagged"] for s in summary.values())
    total_no_candidates = sum(s["n_no_candidates"] for s in summary.values())

    print("  Writing result.json...")
    write_result_json(
        output_dir=output_dir,
        skill="ontology-annotator",
        version="0.1.0",
        summary={
            "input": str(input_path),
            "n_rows": len(df),
            "columns": columns,
            "min_similarity": min_similarity,
            "n_matched": total_matched,
            "n_flagged": total_flagged,
            "n_no_candidates": total_no_candidates,
        },
        data={"column_summary": summary, "flagged_rows": flagged_rows},
    )

    print("  Writing reproducibility bundle...")
    write_environment_yml(
        output_dir,
        env_name="clawbio-ontology-annotator",
        pip_deps=["pandas", "requests"],
        python_version="3.11",
    )
    args: list[str | ReproPath] = ["--demo"] if demo else [
        "--input",
        ReproPath(
            input_path,
            "repo_root" if input_path.is_relative_to(_PROJECT_ROOT) else "auto",
        ),
    ]
    args += ["--output", ReproPath(output_dir, "output_dir")]
    if columns_arg:
        args += ["--columns", columns_arg]
    args += ["--min-similarity", str(min_similarity)]
    write_portable_commands_sh(
        output_dir,
        ReproCommand(
            script_path=Path("skills/ontology-annotator/ontology_annotator.py"),
            args=args,
            comment="Reproduce this ontology-annotator run",
        ),
        repo_root=_PROJECT_ROOT,
    )
    write_checksums(
        [output_dir / "report.md", output_dir / "result.json", output_dir / "annotated.csv"],
        output_dir,
        anchor=output_dir,
    )

    print(f"\n  Matched: {total_matched}  Flagged: {total_flagged}  No candidates: {total_no_candidates}")
    print(f"  Report: {output_dir / 'report.md'}")
    print(f"  Full output: {output_dir}/")
    print(f"\n  {DISCLAIMER}")

    return {
        "annotated_path": str(annotated_path),
        "summary": summary,
        "flagged_rows": flagged_rows,
    }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Ontology Annotator — map metadata column values to UBERON/CL/MONDO/EFO IDs via EBI OLS4"
    )
    parser.add_argument("--input", help="Input CSV, TSV, or .h5ad file")
    parser.add_argument("--output", "-o", required=True, help="Output directory")
    parser.add_argument(
        "--columns",
        default=None,
        help="Comma-separated col:ontology pairs (e.g. tissue:uberon,disease:mondo). "
        "If omitted, auto-detects tissue/cell_type/disease/trait columns.",
    )
    parser.add_argument(
        "--min-similarity",
        type=float,
        default=DEFAULT_MIN_SIMILARITY,
        help=f"Minimum string similarity (label/synonym vs input, not biological confidence) below which a row is flagged for review (default: {DEFAULT_MIN_SIMILARITY})",
    )
    parser.add_argument("--demo", action="store_true", help="Run with bundled demo data, fully offline")
    parser.add_argument("--no-cache", action="store_true", help="Bypass local OLS4 cache")
    parser.add_argument("--cache-dir", default=str(DEFAULT_CACHE_DIR), help="Cache directory")

    args = parser.parse_args(argv)

    if not args.input and not args.demo:
        parser.print_help()
        print("\nError: provide --input or --demo")
        return 1

    output_dir = Path(args.output)

    print("Ontology Annotator")
    print("=" * 60)
    print()

    if args.demo:
        input_path = DEMO_INPUT_PATH
        fixtures_path = DEMO_FIXTURES_PATH
        columns_arg = args.columns  # allow --demo --columns override for testing, else auto-detect
        print(f"  Demo mode: {input_path.name} + recorded OLS4 fixtures (fully offline)")
        print()
    else:
        input_path = Path(args.input)
        if not input_path.exists():
            print(f"Error: input file not found: {input_path}")
            return 1
        fixtures_path = None
        columns_arg = args.columns

    try:
        run_annotation(
            input_path=input_path,
            output_dir=output_dir,
            columns_arg=columns_arg,
            min_similarity=args.min_similarity,
            cache_dir=Path(args.cache_dir),
            use_cache=not args.no_cache,
            fixtures_path=fixtures_path,
            demo=args.demo,
        )
    except ValueError as exc:
        print(f"Error: {exc}")
        return 1
    except RuntimeError as exc:
        print(f"Error: {exc}")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
