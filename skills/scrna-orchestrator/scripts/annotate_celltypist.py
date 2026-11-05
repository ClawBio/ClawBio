#!/usr/bin/env python3
"""Step 6/7 of scrna_orchestrator: OPTIONAL local CellTypist cluster annotation.

Faithful port of `resolve_celltypist_model_path()`, `run_celltypist_annotation()`,
and `aggregate_cluster_annotations()` from the 1.x skill. Putative, model-based
labels aggregated to cluster-level majority vote with support/confidence.

Declared `optional: true` in pipeline.yaml: CellTypist needs a locally
installed human-immune model (runtime downloads are intentionally disabled,
matching the source skill's safety posture) and is not always relevant (e.g.
non-immune tissue, no local model available). `report.py` discovers whether
this step ran the same way `preprocess.py` discovers `detect_doublets`: by
checking for `cluster_annotations.csv` next to its mandatory inputs in the
shared `{out}` directory, not via a `requires` gate — see `report.py`.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import import_scanpy  # noqa: E402

DEFAULT_CELLTYPIST_MODEL = "Immune_All_Low"


def _cluster_sort_key(cluster: str):
    try:
        return (0, int(cluster))
    except ValueError:
        return (1, cluster)


def _normalize_celltypist_model_name(model_name: str) -> str:
    normalized = model_name.strip()
    if not normalized:
        return f"{DEFAULT_CELLTYPIST_MODEL}.pkl"
    if "/" not in normalized and "\\" not in normalized and not normalized.endswith(".pkl"):
        normalized = f"{normalized}.pkl"
    return normalized


def resolve_celltypist_model_path(celltypist, model_name: str) -> Path:
    """Resolve a CellTypist model to a local path without triggering downloads. Verbatim from 1.x skill."""
    normalized = _normalize_celltypist_model_name(model_name)
    if "/" in normalized or "\\" in normalized:
        path = Path(normalized).expanduser()
    else:
        path = Path(celltypist.models.models_path) / normalized

    if path.exists():
        return path

    raise RuntimeError(
        f"CellTypist model '{normalized}' was not found locally at {path}. Install it first with: "
        f"python -c \"import celltypist; celltypist.models.download_models(model='{normalized}')\". "
        "Runtime downloads are disabled for this skill."
    )


def aggregate_cluster_annotations(clusters, predicted_labels, conf_scores, model_name: str) -> pd.DataFrame:
    """Aggregate per-cell CellTypist predictions to cluster-level summaries. Verbatim from 1.x skill."""
    frame = pd.DataFrame(
        {
            "cluster": clusters.astype(str),
            "predicted_cell_type": predicted_labels.astype(str),
            "conf_score": conf_scores.astype(float),
        },
        index=clusters.index,
    )

    rows = []
    cluster_order = sorted(frame["cluster"].unique().tolist(), key=_cluster_sort_key)
    for cluster in cluster_order:
        group = frame.loc[frame["cluster"] == cluster]
        counts = group["predicted_cell_type"].value_counts()
        max_count = int(counts.max())
        winning_labels = sorted(counts[counts == max_count].index.astype(str).tolist())
        majority_label = winning_labels[0]
        majority_mask = group["predicted_cell_type"] == majority_label
        rows.append(
            {
                "cluster": cluster,
                "n_cells": int(group.shape[0]),
                "predicted_cell_type": majority_label,
                "support_fraction": round(float(majority_mask.mean()), 4),
                "mean_confidence": round(float(group.loc[majority_mask, "conf_score"].mean()), 4),
                "annotation_model": model_name,
            }
        )
    return pd.DataFrame(rows)


def run_celltypist_annotation(adata, model_name: str):
    """Run CellTypist and return cluster-level summaries. Verbatim from 1.x skill."""
    import celltypist  # type: ignore

    model_path = resolve_celltypist_model_path(celltypist, model_name)
    model = celltypist.models.Model.load(str(model_path))

    overlap = int(np.isin(np.asarray(adata.var_names, dtype=object), model.features).sum())
    if overlap == 0:
        raise RuntimeError(
            "CellTypist annotation requires human gene symbols overlapping the local model. "
            f"Found 0 overlapping genes with {model_path.name}."
        )

    adata = adata.copy()
    adata.var_names_make_unique()
    result = celltypist.annotate(adata, model=model, majority_voting=False)

    predicted = result.predicted_labels["predicted_labels"].reindex(adata.obs_names)
    conf_scores = result.probability_matrix.max(axis=1).reindex(adata.obs_names)
    annotations = aggregate_cluster_annotations(
        clusters=adata.obs["leiden"], predicted_labels=predicted, conf_scores=conf_scores, model_name=model_path.name
    )
    metadata = {
        "backend": "celltypist",
        "model": model_path.name,
        "model_path": str(model_path),
        "overlap_genes": overlap,
        "putative": True,
        "n_clusters_annotated": int(annotations.shape[0]),
    }
    return annotations, metadata


def main() -> None:
    parser = argparse.ArgumentParser(description="scrna_orchestrator: annotate_celltypist step (optional)")
    parser.add_argument("--markers-h5ad", required=True, help="markers.h5ad from the find_markers step")
    parser.add_argument("--annotation-model", default=DEFAULT_CELLTYPIST_MODEL)
    parser.add_argument("--out", required=True, help="output path for cluster_annotations.csv")
    parser.add_argument("--info-out", required=True, help="output path for annotation_info.json")
    args = parser.parse_args()

    sc = import_scanpy()
    adata = sc.read_h5ad(Path(args.markers_h5ad))

    try:
        annotations, info = run_celltypist_annotation(adata, model_name=args.annotation_model)
    except (RuntimeError, ImportError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    annotations.to_csv(out_path, index=False)

    info_path = Path(args.info_out)
    info_path.parent.mkdir(parents=True, exist_ok=True)
    info_path.write_text(json.dumps(info, indent=2))

    print(f"annotate_celltypist complete: {info['n_clusters_annotated']} clusters annotated with {info['model']}")
    print(f"  wrote {out_path}")
    print(f"  wrote {info_path}")


if __name__ == "__main__":
    main()
