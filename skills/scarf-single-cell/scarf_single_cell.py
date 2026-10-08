#!/usr/bin/env python3
"""scarf-single-cell: a baseline out-of-core scRNA-seq pass with Scarf (NygenAnalytics/scarf).

Converts a raw-count H5AD into a Scarf Zarr store (or opens an existing store), runs the Scarf
RNA pipeline under an explicit memory budget, audits QC retention per sample, and writes
report.md, result.json, figures, tables, a raw-count H5AD handoff for scrna-orchestrator or
scrna-embedding, and a reproducibility bundle.

Usage:
    python skills/scarf-single-cell/scarf_single_cell.py --demo --output /tmp/scarf_demo
    python skills/scarf-single-cell/scarf_single_cell.py --input counts.h5ad --output out/ \\
        --sample-column donor_id
    python skills/scarf-single-cell/scarf_single_cell.py --input store.zarr --output out/ \\
        --label baseline_v2

Requires Python 3.12+ and scarf 1.0.0rc17 or newer: pip install "scarf[extra]>=1.0.0rc17".
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import os
import shlex
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

SKILL_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SKILL_DIR.parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

# Only the light shared helpers: Scarf's own environment need not carry ClawBio's full stack.
from clawbio.common.checksums import sha256_file  # noqa: E402
from clawbio.common.reproducibility import (  # noqa: E402
    write_checksums,
    write_commands_sh,
    write_environment_yml,
)

SKILL_NAME = "scarf-single-cell"
SKILL_VERSION = "0.3.0"
SCARF_REQUIREMENT = "scarf[extra]>=1.0.0rc17"
INSTALL_HINT = (
    f"scarf-single-cell needs Python 3.12+ and {SCARF_REQUIREMENT}. Install it with:\n"
    f'  pip install "{SCARF_REQUIREMENT}"\n'
    "A bare `pip install scarf[extra]` resolves to the old 0.32 series, which this skill "
    "does not support."
)
DEMO_LABEL = "clawbio_demo"
DEMO_SAMPLE_COLUMN = "donor_id"
DEMO_HOLDOUT_COLUMN = "author_cell_type"
HANDOFF_NAME = "run_cells.h5ad"
EXIT_MISSING_SCARF = 2
DISCLAIMER = (
    "ClawBio is a research and educational tool. It is not a medical device and does not "
    "provide clinical diagnoses. Consult a healthcare professional before making any medical "
    "decisions."
)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Baseline out-of-core scRNA-seq pass with Scarf (NygenAnalytics/scarf)."
    )
    parser.add_argument(
        "--input",
        type=Path,
        dest="input_path",
        help="Raw-count .h5ad (converted into OUTPUT/store.zarr) or an existing Scarf .zarr store",
    )
    parser.add_argument("--output", type=Path, help="Output directory")
    parser.add_argument("--demo", action="store_true", help="Run on the bundled synthetic dataset")
    parser.add_argument("--label", help="Pipeline run label; must be new for the store")
    parser.add_argument(
        "--n-mads", type=float, default=5.0, help="MAD bound for the pipeline QC filter (default 5)"
    )
    parser.add_argument(
        "--sample-column",
        help="Cell column with the sample or donor unit; frozen into the run and used for audits",
    )
    parser.add_argument(
        "--holdout-column",
        help="Author label column compared only after the analysis is finished (never used upstream)",
    )
    parser.add_argument(
        "--hvg-blacklist",
        help="Regex that replaces Scarf's default HVG blacklist (see references/gene-blacklists.md)",
    )
    parser.add_argument(
        "--no-cell-cycle",
        action="store_true",
        help="Skip cell-cycle scoring (needed when gene names are not human or mouse symbols)",
    )
    parser.add_argument(
        "--skip-handoff",
        action="store_true",
        help="Do not write handoff/run_cells.h5ad (it holds X in memory; skip for atlas-scale runs)",
    )
    parser.add_argument(
        "--mem-budget",
        default=os.environ.get("SCARF_MEM_BUDGET", "2G"),
        help="Scarf memory budget with a unit, e.g. 2G (default: $SCARF_MEM_BUDGET or 2G)",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=int(os.environ.get("SCARF_WORKERS", "2")),
        help="Scarf worker count (default: $SCARF_WORKERS or 2)",
    )
    args = parser.parse_args(argv)
    if not args.demo and args.input_path is None:
        parser.error("provide --input <file.h5ad|store.zarr> or --demo")
    if args.demo:
        args.label = args.label or DEMO_LABEL
        args.sample_column = args.sample_column or DEMO_SAMPLE_COLUMN
        args.holdout_column = args.holdout_column or DEMO_HOLDOUT_COLUMN
        args.no_cell_cycle = True  # synthetic gene names contain no cell-cycle genes
    args.label = args.label or "clawbio_baseline"
    if args.holdout_column and args.holdout_column == args.sample_column:
        parser.error("--holdout-column must differ from --sample-column")
    if args.output is None:
        stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
        mode = "demo" if args.demo else "run"
        args.output = Path(tempfile.gettempdir()) / "scarf_single_cell" / f"{mode}-{stamp}"
    return args


def import_scarf():
    """Configure resources from the environment, then import Scarf or exit with an install hint."""
    if sys.version_info < (3, 12):
        print(INSTALL_HINT, file=sys.stderr)
        sys.exit(EXIT_MISSING_SCARF)
    os.environ.setdefault("MPLBACKEND", "Agg")
    try:
        import scarf
    except ImportError:
        print(INSTALL_HINT, file=sys.stderr)
        sys.exit(EXIT_MISSING_SCARF)
    version = getattr(scarf, "__version__", "0")
    if version.split(".")[0].isdigit() and int(version.split(".")[0]) < 1:
        print(f"Found scarf {version}. " + INSTALL_HINT, file=sys.stderr)
        sys.exit(EXIT_MISSING_SCARF)
    scarf.configure_output(level="WARNING", progress=False)
    return scarf


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def refuse_existing(output: Path, paths: list[Path]) -> None:
    existing = [p for p in paths if p.exists()]
    if existing:
        names = ", ".join(str(p.relative_to(output)) for p in existing)
        print(
            f"Error: {output} already holds {names}. Choose a new --output directory; "
            "this skill does not overwrite earlier results.",
            file=sys.stderr,
        )
        sys.exit(1)


def prepare_store(scarf, args: argparse.Namespace, output: Path) -> tuple[Path, Path | None]:
    """Return (store, source_h5ad). H5AD input is converted into OUTPUT/store.zarr."""
    if args.demo:
        generator = load_module(
            "scarf_single_cell_demo_data", SKILL_DIR / "examples" / "make_demo_data.py"
        )
        source = generator.make_demo_h5ad(output / "demo_input.h5ad")
    else:
        source = args.input_path.resolve()
        if not source.exists():
            print(f"Error: input not found: {source}", file=sys.stderr)
            sys.exit(1)
        if source.is_dir() and source.suffix == ".zarr":
            return source, None
        if source.suffix != ".h5ad":
            print(
                "Error: --input must be a .h5ad file or a Scarf .zarr store. Convert 10x, MTX "
                "or Seurat inputs first (references/data-access.md).",
                file=sys.stderr,
            )
            sys.exit(1)
    store = output / "store.zarr"
    reader = scarf.H5adReader.from_inspect(scarf.inspect_h5ad(str(source)))
    scarf.H5adToZarr(reader, zarr_loc=str(store)).dump()
    return store, source


def inspect_store(store: Path, output: Path) -> dict:
    """Run scripts/inspect_store.py read-only and keep its text and JSON output."""
    profile = output / "tables" / "store_profile.json"
    done = subprocess.run(
        [sys.executable, str(SKILL_DIR / "scripts" / "inspect_store.py"), str(store),
         "--json", str(profile)],
        capture_output=True, text=True, env=os.environ.copy(),
    )
    (output / "tables" / "inspect_store.txt").write_text(done.stdout + done.stderr)
    if done.returncode != 0:
        print(done.stderr, file=sys.stderr)
        sys.exit(1)
    return json.loads(profile.read_text())


def retention_table(pd, total_mask, run_mask, groups):
    frame = pd.DataFrame({"group": groups, "active": total_mask, "kept": run_mask})
    frame = frame[frame["active"]]
    table = frame.groupby("group", sort=True).agg(n_cells=("kept", "size"), n_kept=("kept", "sum"))
    table["frac_kept"] = (table["n_kept"] / table["n_cells"]).round(4)
    return table.reset_index()


def write_handoff(ds, run, path: Path) -> None:
    """Write the run's cells as a raw-count H5AD that the ClawBio Scanpy/scVI skills load.

    `ds.to_anndata(run=run)` holds X in memory, so this is for runs that fit in RAM.
    Scarf stores small counts as narrow unsigned integers (uint8 on the demo); Scanpy's QC
    kernels reject those, so X is cast to float32 (values stay integer counts). Nullable
    string indexes and columns become plain object strings for the same reason.
    """
    import numpy as np
    import pandas as pd

    adata = ds.to_anndata(run=run)
    adata.X = adata.X.astype(np.float32)

    def plain(index):
        return pd.Index(np.asarray(index, dtype=object), dtype=object, name=index.name)

    adata.obs_names = plain(adata.obs_names)
    adata.var_names = plain(adata.var_names)
    for frame in (adata.obs, adata.var):
        for column in frame.columns:
            if isinstance(frame[column].dtype, pd.StringDtype):
                frame[column] = frame[column].astype(object)
    adata.write_h5ad(path)


def raw_count_guard(handoff: Path) -> str | None:
    """Load the handoff with the loader scrna-orchestrator and scrna-embedding use.

    Returns None when ClawBio's raw-count guard accepts the file, else the rejection reason.
    """
    import anndata as ad

    from clawbio.common.scrna_io import load_count_adata

    try:
        load_count_adata(handoff, h5ad_loader=ad.read_h5ad, expected_input="raw-count .h5ad")
    except ValueError as exc:
        return str(exc)
    return None


def md_table(frame) -> str:
    cols = list(frame.columns)
    lines = ["| " + " | ".join(map(str, cols)) + " |", "|" + "---|" * len(cols)]
    for row in frame.itertuples(index=False):
        lines.append("| " + " | ".join(str(v) for v in row) + " |")
    return "\n".join(lines)


def run_analysis(args: argparse.Namespace) -> dict:
    os.environ["SCARF_MEM_BUDGET"] = str(args.mem_budget)
    os.environ["SCARF_WORKERS"] = str(args.workers)
    scarf = import_scarf()
    import numpy as np
    import pandas as pd
    from sklearn.metrics import adjusted_rand_score

    output = args.output.resolve()
    refuse_existing(output, [output / "report.md", output / "result.json", output / "store.zarr"])
    for sub in ("figures", "tables", "handoff"):
        (output / sub).mkdir(parents=True, exist_ok=True)

    store, source = prepare_store(scarf, args, output)
    # A writable open prepares a fresh store; -1 keeps every cell active (golden rule 1).
    ds = scarf.DataStore(str(store), min_features_per_cell=-1)
    if args.sample_column and args.sample_column not in ds.cells.columns:
        print(f"Error: --sample-column {args.sample_column!r} is not a cell column", file=sys.stderr)
        sys.exit(1)
    if args.holdout_column and args.holdout_column not in ds.cells.columns:
        print(f"Error: --holdout-column {args.holdout_column!r} is not a cell column", file=sys.stderr)
        sys.exit(1)
    profile = inspect_store(store, output)

    params = {"hvg": {"blacklist": args.hvg_blacklist}} if args.hvg_blacklist else None
    try:
        run = ds.pipeline.run(
            label=args.label,
            filtering={"method": "mad", "n_mads": args.n_mads},
            snapshot_columns=[args.sample_column] if args.sample_column else [],
            cell_cycle=not args.no_cell_cycle,
            params=params,
        )
    except ValueError as exc:
        print(f"Error: {exc}. Labels are immutable: pass a new --label.", file=sys.stderr)
        sys.exit(1)

    (output / "run_report.md").write_text(run.report(format="markdown"))
    (output / "lineage.md").write_text(ds.lineage({"clusters": run["clusters"]}).to_markdown())
    refs = {name: run[name].to_dict() for name in run.keys()}
    handoff_record = {
        "store": str(store),
        "label": args.label,
        "run_id": run.run_id,
        "status": str(run.status),
        "refs": refs,
    }
    (output / "handoff.json").write_text(json.dumps(handoff_record, indent=2, default=str))

    total_mask = np.asarray(ds.cells.fetch_all("I"), dtype=bool)
    run_mask = np.asarray(run.cells.fetch_all("I"), dtype=bool)
    clusters_all = np.asarray(run.cells.fetch_all("clusters"))
    clusters = clusters_all[run_mask]
    sizes = pd.Series(clusters).value_counts().sort_index()
    sizes_table = sizes.rename_axis("cluster").reset_index(name="n_cells")
    sizes_table.to_csv(output / "tables" / "cluster_sizes.csv", index=False)

    markers = ds.get_markers(marker=run["markers"])
    markers.to_csv(output / "tables" / "markers.csv", index=False)
    top = (
        markers.groupby("group_id", sort=False).head(5)
        .groupby("group_id", sort=False)["feature_name"].agg(", ".join)
    )

    retention = composition = None
    if args.sample_column:
        groups = np.asarray(ds.cells.fetch_all(args.sample_column)).astype(str)
        retention = retention_table(pd, total_mask, run_mask, groups)
        retention.to_csv(output / "tables" / "qc_retention_by_sample.csv", index=False)
        composition = pd.crosstab(
            pd.Series(clusters, name="cluster"), pd.Series(groups[run_mask], name=args.sample_column)
        )
        composition.to_csv(output / "tables" / "cluster_by_sample.csv")

    umap = ds.plots.embedding(run=run, color_by="clusters", show=False)
    umap.save(output / "figures" / "umap_clusters.png", dpi=120)
    umap.close()

    handoff = output / "handoff" / HANDOFF_NAME
    guard = "skipped (--skip-handoff)"
    if not args.skip_handoff:
        write_handoff(ds, run, handoff)
        guard = raw_count_guard(handoff)

    holdout = None
    if args.holdout_column:
        # Read only now, after every analysis decision above (golden rule 14).
        labels = np.asarray(ds.cells.fetch_all(args.holdout_column)).astype(str)[run_mask]
        crosstab = pd.crosstab(
            pd.Series(clusters, name="cluster"), pd.Series(labels, name=args.holdout_column)
        )
        crosstab.to_csv(output / "tables" / "holdout_crosstab.csv")
        holdout = {
            "column": args.holdout_column,
            "adjusted_rand_index": round(float(adjusted_rand_score(labels, clusters)), 4),
        }

    n_kept = int(run_mask.sum())
    summary = {
        "label": args.label,
        "run_status": str(run.status),
        "cells_total": int(ds.cells.N),
        "cells_in_run": n_kept,
        "clusters": int(sizes.size),
        "handoff_raw_count_guard": (
            "accepted" if guard is None else "skipped" if args.skip_handoff else "rejected"
        ),
    }
    data = {
        "scarf_version": scarf.__version__,
        "input": str(source if source is not None else store),
        "store": str(store),
        "run_id": run.run_id,
        "resources": {"mem_budget": args.mem_budget, "workers": args.workers},
        "filtering": {"method": "mad", "n_mads": args.n_mads},
        "cell_cycle": not args.no_cell_cycle,
        "hvg_blacklist": args.hvg_blacklist or "scarf default",
        "cluster_sizes": {str(k): int(v) for k, v in sizes.items()},
        "qc_retention_by_sample": retention.to_dict(orient="records") if retention is not None else None,
        "matrix_hints": profile.get("hints", []),
        "handoff": {
            "path": None if args.skip_handoff else f"handoff/{HANDOFF_NAME}",
            "raw_count_guard_reason": guard,
        },
        "holdout": holdout,
    }
    return {
        "summary": summary, "data": data, "top_markers": top, "sizes": sizes_table,
        "retention": retention, "composition": composition, "source": source, "output": output,
    }


def write_report(result: dict, args: argparse.Namespace) -> None:
    output: Path = result["output"]
    summary, data = result["summary"], result["data"]
    lines = [
        "# Scarf single-cell baseline report",
        "",
        f"**Date**: {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}",
        f"**Skill**: {SKILL_NAME} {SKILL_VERSION}",
        f"**Scarf**: {data['scarf_version']}",
        f"**Run**: `{summary['label']}` ({summary['run_status']})",
        f"**Resources**: SCARF_MEM_BUDGET={args.mem_budget}, SCARF_WORKERS={args.workers}",
    ]
    if result["source"] is not None:
        source = Path(result["source"])
        lines.append(f"**Input**: `{source.name}` (sha256 `{sha256_file(source)}`)")
    else:
        lines.append(f"**Input**: Scarf store `{data['store']}`")
    lines += ["", "---", ""]
    if args.demo:
        lines += [
            "> Demo mode: a seeded synthetic dataset (2,000 cells, 1,500 genes, four planted",
            "> groups, two donors). Numbers below describe that toy matrix only.",
            "",
        ]
    lines += [
        "## Summary",
        "",
        f"- Cells: {summary['cells_in_run']:,} of {summary['cells_total']:,} kept by the "
        f"MAD {args.n_mads:g} filter",
        f"- Clusters (pipeline silhouette pick): {summary['clusters']}",
        f"- Handoff H5AD and ClawBio's raw-count guard: {summary['handoff_raw_count_guard']}",
        "",
        "## Matrix check (scripts/inspect_store.py)",
        "",
    ]
    lines += [f"- {hint}" for hint in data["matrix_hints"]] or [
        "- No sign of corrected or pre-filtered counts in these checks."
    ]
    if result["retention"] is not None:
        lines += ["", f"## QC retention by `{args.sample_column}`", "", md_table(result["retention"])]
    lines += ["", "## Cluster sizes", "", md_table(result["sizes"]), "", "## Top markers", ""]
    lines += [f"- Cluster {gid}: {genes}" for gid, genes in result["top_markers"].items()]
    if result["composition"] is not None:
        comp = result["composition"].reset_index()
        lines += ["", f"## Composition by `{args.sample_column}`", "", md_table(comp)]
    if args.skip_handoff:
        lines += ["", "## Handoff", "", "Skipped (`--skip-handoff`).", ""]
    else:
        lines += [
            "",
            "## Handoff",
            "",
            f"`handoff/{HANDOFF_NAME}` holds the run's cells with raw counts in `X`, Scarf `clusters`,",
            "Leiden and Paris labels in `obs`, and `X_umap`. Continue in Scanpy or scVI with:",
            "",
            "```bash",
            f"python clawbio.py run scrna --input {output}/handoff/{HANDOFF_NAME} --output <dir>",
            f"python clawbio.py run scrna-embedding --input {output}/handoff/{HANDOFF_NAME} --output <dir>",
            "```",
            "",
            "Those skills recompute QC, normalisation and clustering; Scarf labels stay in `obs`.",
            "",
        ]
    if data["holdout"]:
        lines += [
            f"## Held-out comparison: `{data['holdout']['column']}`",
            "",
            "Read only after the analysis above was finished; nothing upstream used it.",
            "",
            f"- Adjusted Rand index (clusters vs held-out labels): "
            f"{data['holdout']['adjusted_rand_index']}",
            "- Crosstab: `tables/holdout_crosstab.csv`",
            "",
        ]
    lines += [
        "## Caveats",
        "",
        "- This is a baseline pass. Audit QC removals, cluster stability and markers before",
        "  naming cell types (SKILL.md, Gotchas).",
        "- Cells are not replicates: condition claims need donor-level aggregation.",
        "- Run details: `run_report.md`; provenance: `lineage.md`, `handoff.json`.",
        "",
        "---",
        "",
        "## Disclaimer",
        "",
        f"*{DISCLAIMER}*",
        "",
    ]
    (output / "report.md").write_text("\n".join(lines))
    checksum = sha256_file(result["source"]) if result["source"] is not None else ""
    envelope = {
        "skill": SKILL_NAME,
        "version": SKILL_VERSION,
        "completed_at": datetime.now(timezone.utc).isoformat(),
        "input_checksum": f"sha256:{checksum}" if checksum else "",
        "summary": summary,
        "data": data,
    }
    (output / "result.json").write_text(json.dumps(envelope, indent=2, default=str))


def write_reproducibility_bundle(args: argparse.Namespace, output: Path) -> None:
    if args.demo:
        invocation = "--demo --output \"$OUTPUT_DIR\""
    else:
        invocation = (
            f"--input {shlex.quote(str(args.input_path))} --output \"$OUTPUT_DIR\" "
            f"--label {shlex.quote(args.label)}"
        )
        for flag, value in (
            ("--sample-column", args.sample_column),
            ("--holdout-column", args.holdout_column),
            ("--hvg-blacklist", args.hvg_blacklist),
        ):
            if value:
                invocation += f" {flag} {shlex.quote(value)}"
        if args.no_cell_cycle:
            invocation += " --no-cell-cycle"
        if args.skip_handoff:
            invocation += " --skip-handoff"
    invocation += f" --n-mads {args.n_mads:g}"
    write_commands_sh(
        output,
        "# scarf-single-cell - reproducibility\n"
        "set -euo pipefail\n"
        'OUTPUT_DIR="${OUTPUT_DIR:-./scarf_single_cell_rerun}"\n'
        "\n"
        "# 1. Recreate the environment (Python 3.12+)\n"
        'conda env create -f "$(dirname "$0")/environment.yml"\n'
        "conda activate scarf-single-cell\n"
        "\n"
        "# 2. Re-run the analysis from the ClawBio repository root, same budget\n"
        f"SCARF_MEM_BUDGET={shlex.quote(str(args.mem_budget))} SCARF_WORKERS={args.workers} \\\n"
        f"  python skills/scarf-single-cell/scarf_single_cell.py {invocation}\n"
        "\n"
        "# 3. Verify this bundle's original outputs (labels are relative to its output directory)\n"
        '( cd "$(dirname "$0")/.." && sha256sum -c reproducibility/checksums.sha256 )',
    )
    write_environment_yml(
        output,
        env_name="scarf-single-cell",
        python_version="3.12",
        pip_deps=[SCARF_REQUIREMENT],
    )
    tracked = [output / "report.md", output / "result.json", output / "run_report.md",
               output / "lineage.md", output / "handoff.json"]
    tracked += sorted((output / "tables").glob("*")) + sorted((output / "figures").glob("*"))
    tracked += sorted((output / "handoff").glob("*.h5ad"))
    write_checksums(tracked, output, anchor=output)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    result = run_analysis(args)
    write_report(result, args)
    write_reproducibility_bundle(args, result["output"])
    output = result["output"]
    print(f"Report written to {output / 'report.md'}")
    print(f"Results written to {output / 'result.json'}")
    if not args.skip_handoff:
        print(f"Handoff H5AD written to {output / 'handoff' / HANDOFF_NAME}")
    print(f"Reproducibility bundle written to {output / 'reproducibility'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
