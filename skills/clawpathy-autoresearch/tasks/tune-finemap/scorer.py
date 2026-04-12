"""Scorer for tune-finemap task.

Expects skill_output shape:
    {"L1": {"credible_set": [rsid, ...], "L": int, ...}, "L2": {...}, "L3": {...}}

Loss = 1 - mean Jaccard(predicted_set, truth_set) across loci.
Perfect reproduction -> 0.0. Worst -> 1.0.
"""
from __future__ import annotations

import json
from pathlib import Path


def _jaccard(a, b):
    sa, sb = set(a), set(b)
    if not sa and not sb:
        return 1.0
    union = sa | sb
    if not union:
        return 1.0
    return len(sa & sb) / len(union)


def score(skill_output, task_dir):
    task_dir = Path(task_dir)
    truth = json.loads((task_dir / "ground_truth.json").read_text())
    breakdown = {}
    jaccards = []
    for locus_id, truth_rsids in truth.items():
        pred = skill_output.get(locus_id) or {}
        pred_rsids = pred.get("credible_set") or []
        if not isinstance(pred_rsids, list):
            pred_rsids = []
        j = _jaccard(pred_rsids, truth_rsids)
        jaccards.append(j)
        breakdown[f"{locus_id}_jaccard"] = j
        breakdown[f"{locus_id}_pred_size"] = len(pred_rsids)
        breakdown[f"{locus_id}_truth_size"] = len(truth_rsids)
    mean_j = sum(jaccards) / len(jaccards) if jaccards else 0.0
    total = 1.0 - mean_j
    breakdown["mean_jaccard"] = mean_j
    return total, breakdown
