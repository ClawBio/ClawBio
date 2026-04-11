"""GWAS reproduction scorer for the demo workspace.

Score = weighted mean of six error components. Lower is better. 0 = perfect.
"""
from __future__ import annotations

MISSING_VARIANT_PENALTY = 1.0
DIRECTION_PENALTY = 1.0
WEIGHTS = {
    "p_value": 0.20,
    "or": 0.25,
    "freq": 0.10,
    "locus_count": 0.15,
    "missing": 0.20,
    "direction": 0.10,
}


def _score_paper(gt_paper: dict, agent_paper: dict) -> float:
    """Score a single paper's reproduction."""
    gt_variants = gt_paper.get("lead_variants", [])
    agent_variants = agent_paper.get("variants_found", [])
    gt_loci = gt_paper.get("total_loci", 0)
    agent_loci = agent_paper.get("total_loci_reported", 0)

    if not gt_variants:
        lce = abs(gt_loci - agent_loci) / max(gt_loci, 1)
        return lce * WEIGHTS["locus_count"]

    agent_by_rsid = {v["rsid"]: v for v in agent_variants}

    p_errors, or_errors, freq_errors, direction_errors = [], [], [], []
    missing_count = 0

    for gt_v in gt_variants:
        rsid = gt_v["rsid"]
        if rsid not in agent_by_rsid:
            missing_count += 1
            continue

        av = agent_by_rsid[rsid]

        gt_p = gt_v.get("neg_log10_p", 0)
        ag_p = av.get("neg_log10_p", 0)
        p_errors.append(abs(gt_p - ag_p) / gt_p if gt_p > 0 else abs(gt_p - ag_p))

        gt_or = gt_v.get("odds_ratio", 1.0)
        ag_or = av.get("odds_ratio", 1.0)
        or_errors.append(abs(gt_or - ag_or))

        gt_freq = gt_v.get("effect_allele_freq", 0.0)
        ag_freq = av.get("effect_allele_freq", 0.0)
        freq_errors.append(abs(gt_freq - ag_freq))

        gt_dir = gt_v.get("effect_direction", "")
        ag_dir = av.get("effect_direction", "")
        direction_errors.append(
            DIRECTION_PENALTY if gt_dir and ag_dir and gt_dir != ag_dir else 0.0
        )

    n = len(gt_variants)
    p_val_err = sum(p_errors) / len(p_errors) if p_errors else 0.0
    or_err = sum(or_errors) / len(or_errors) if or_errors else 0.0
    freq_err = sum(freq_errors) / len(freq_errors) if freq_errors else 0.0
    dir_err = sum(direction_errors) / n
    missing_err = missing_count * MISSING_VARIANT_PENALTY / n
    locus_err = abs(gt_loci - agent_loci) / gt_loci if gt_loci > 0 else float(agent_loci)

    return (
        WEIGHTS["p_value"] * p_val_err
        + WEIGHTS["or"] * or_err
        + WEIGHTS["freq"] * freq_err
        + WEIGHTS["locus_count"] * locus_err
        + WEIGHTS["missing"] * missing_err
        + WEIGHTS["direction"] * dir_err
    )


def score(skill_output: dict, ground_truth: dict) -> float:
    """Score the agent's full output against all papers.

    Args:
        skill_output: {"papers": {"paper_id": {"variants_found": [...], "total_loci_reported": int}}}
        ground_truth: {"targets": [{"paper_id": str, "lead_variants": [...], "total_loci": int}]}

    Returns:
        Mean error across all papers. Lower is better. 0 = perfect.
    """
    targets = ground_truth.get("targets", [])
    if not targets:
        return 0.0

    agent_papers = skill_output.get("papers", {})
    paper_errors = []

    for target in targets:
        pid = target["paper_id"]
        agent_paper = agent_papers.get(pid, {"variants_found": [], "total_loci_reported": 0})
        paper_errors.append(_score_paper(target, agent_paper))

    return sum(paper_errors) / len(paper_errors)
