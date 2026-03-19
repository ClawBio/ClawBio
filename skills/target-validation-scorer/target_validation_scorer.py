#!/usr/bin/env python3
"""target-validation-scorer: Evidence-grounded target validation for drug discovery."""

import argparse
import json
from pathlib import Path
from datetime import datetime, timezone

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

SKILL_DIR = Path(__file__).resolve().parent
DEMO_DATA = SKILL_DIR / "demo_input.json"

# ---------------------------------------------------------------------------
# Scoring logic
# ---------------------------------------------------------------------------


def score_disease_association(ev):
    if ev is None:
        return {"score": None, "confidence": "low", "source": "Open Targets"}
    raw = ev.get("overall_score", 0)
    if raw >= 0.7:
        s = 20
    elif raw >= 0.4:
        s = 10
    else:
        s = 0
    conf = "high" if s == 20 else "medium" if s == 10 else "low"
    return {"score": s, "confidence": conf, "source": "Open Targets", "detail": ev.get("detail", "")}


def score_druggability(ev):
    if ev is None:
        return {"score": None, "confidence": "low", "source": "ChEMBL + UniProt"}
    s = 0
    tc = ev.get("target_class", "").lower()
    if tc in ("kinase", "gpcr", "nuclear receptor", "protease", "ion channel"):
        s += 10
    elif tc:
        s += 5
    n_lig = ev.get("known_ligands", 0)
    if n_lig > 500:
        s += 6
    elif n_lig > 100:
        s += 5
    elif n_lig > 10:
        s += 3
    elif n_lig > 0:
        s += 1
    s = min(20, s)
    conf = "high" if s >= 14 else "medium" if s >= 8 else "low"
    return {"score": s, "confidence": conf, "source": "ChEMBL + UniProt", "detail": ev.get("detail", "")}


def score_chemical_matter(ev):
    if ev is None:
        return {"score": None, "confidence": "low", "source": "ChEMBL"}
    s = 0
    n = ev.get("bioactive_compounds", 0)
    if n > 500:
        s += 8
    elif n > 100:
        s += 6
    elif n > 10:
        s += 3
    elif n > 0:
        s += 1
    n_potent = ev.get("compounds_sub_100nM", 0)
    if n_potent > 50:
        s += 7
    elif n_potent > 10:
        s += 5
    elif n_potent > 0:
        s += 3
    s = min(20, s)
    conf = "high" if s >= 14 else "medium" if s >= 8 else "low"
    return {"score": s, "confidence": conf, "source": "ChEMBL", "detail": ev.get("detail", "")}


def score_clinical_precedent(ev):
    if ev is None:
        return {"score": None, "confidence": "low", "source": "ChEMBL + ClinicalTrials.gov"}
    phase = ev.get("max_phase", 0)
    if phase >= 1:
        s = 20
    elif ev.get("trials_completed", 0) > 0 or ev.get("trials_active", 0) > 0:
        s = 10
    else:
        s = 0
    conf = "high" if phase >= 2 else "medium" if phase >= 1 else "low"
    return {"score": s, "confidence": conf, "source": "ChEMBL + ClinicalTrials.gov", "detail": ev.get("detail", "")}


def score_structural_data(ev):
    if ev is None:
        return {"score": None, "confidence": "low", "source": "PDB + AlphaFold"}
    has_pdb = ev.get("pdb_structures", 0) > 0
    has_cocrystal = ev.get("co_crystal_ligands", False)
    good_res = (ev.get("best_resolution_A") or 99) < 2.5
    has_af = ev.get("alphafold_available", False)

    if has_pdb and has_cocrystal and good_res:
        s = 20
    elif has_pdb or has_af:
        s = 10
    else:
        s = 0
    conf = "high" if s == 20 else "medium" if s == 10 else "low"
    return {"score": s, "confidence": conf, "source": "PDB + AlphaFold", "detail": ev.get("detail", "")}


def compute_safety_penalty(signals):
    if not signals:
        return 0, []
    total = 0
    flags = []
    for sig in signals:
        p = sig.get("penalty", 0)
        total += p
        flags.append({
            "risk": sig.get("risk", "unknown"),
            "penalty": p,
            "tier": sig.get("tier", "T4"),
            "evidence": sig.get("evidence", ""),
            "mitigation": sig.get("mitigation", ""),
        })
    return total, flags


DECISION_TIERS = [
    (75, "GO", "Strong evidence across multiple dimensions"),
    (50, "CONDITIONAL_GO", "Proceed with explicit risk mitigation plan"),
    (25, "REVIEW", "Insufficient evidence; needs more data before commitment"),
    (0, "NO_GO", "Target lacks fundamental validation"),
]


def decide(adjusted_score):
    for threshold, decision, meaning in DECISION_TIERS:
        if adjusted_score >= threshold:
            return decision, meaning
    return "NO_GO", "Target lacks fundamental validation"


def build_rationale(target, disease, sub_scores, safety_flags, decision):
    parts = []
    strong = [k for k, v in sub_scores.items() if v.get("score") and v["score"] >= 14]
    weak = [k for k, v in sub_scores.items() if v.get("score") is not None and v["score"] < 8]
    missing = [k for k, v in sub_scores.items() if v.get("score") is None]

    if strong:
        parts.append(f"Strong evidence in: {', '.join(s.replace('_', ' ') for s in strong)}.")
    if weak:
        parts.append(f"Weak evidence in: {', '.join(w.replace('_', ' ') for w in weak)}.")
    if missing:
        parts.append(f"No data for: {', '.join(m.replace('_', ' ') for m in missing)}.")
    if safety_flags:
        risks = [f["risk"] for f in safety_flags if f["penalty"] < 0]
        if risks:
            parts.append(f"Safety concerns: {'; '.join(risks)}.")

    parts.append(f"Decision: {decision}" + (f" for {disease}." if disease else "."))

    if decision == "CONDITIONAL_GO":
        mits = [f["mitigation"] for f in safety_flags if f.get("mitigation")]
        if mits:
            parts.append(f"Required mitigation: {mits[0]}")

    return " ".join(parts)


# ---------------------------------------------------------------------------
# Main validation
# ---------------------------------------------------------------------------


def validate_target(data):
    target = data.get("target", "unknown")
    disease = data.get("disease")
    ev = data.get("evidence", {})

    sub_scores = {
        "disease_association": score_disease_association(ev.get("disease_association")),
        "druggability": score_druggability(ev.get("druggability")),
        "chemical_matter": score_chemical_matter(ev.get("chemical_matter")),
        "clinical_precedent": score_clinical_precedent(ev.get("clinical_precedent")),
        "structural_data": score_structural_data(ev.get("structural_data")),
    }

    raw_score = sum(v["score"] for v in sub_scores.values() if v["score"] is not None)
    penalty, safety_flags = compute_safety_penalty(ev.get("safety_signals", []))
    adjusted = max(0, raw_score + penalty)
    decision, meaning = decide(adjusted)
    rationale = build_rationale(target, disease, sub_scores, safety_flags, decision)

    evidence_trail = []
    for dim, obj in sub_scores.items():
        if obj["score"] is not None:
            # Prefer tier from input evidence, fallback to confidence mapping
            input_ev = ev.get(dim) or {}
            tier = input_ev.get("tier") or (
                "T1" if obj["confidence"] == "high" else
                "T2" if obj["confidence"] == "medium" else "T3"
            )
            evidence_trail.append({
                "type": dim,
                "source": obj["source"],
                "detail": obj.get("detail", ""),
                "tier": tier,
                "confidence": obj["confidence"],
            })

    return {
        "target": target,
        "disease": disease,
        "raw_score": raw_score,
        "safety_penalty": penalty,
        "adjusted_score": adjusted,
        "decision": decision,
        "decision_meaning": meaning,
        "sub_scores": sub_scores,
        "safety_flags": safety_flags,
        "evidence": evidence_trail,
        "rationale": rationale,
        "retrieved_at": datetime.now(timezone.utc).isoformat(),
    }


# ---------------------------------------------------------------------------
# Report generation (ClawBio convention: output/report.md)
# ---------------------------------------------------------------------------


def generate_report(result, output_dir: Path) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    report_path = output_dir / "report.md"

    r = result
    lines = [
        "# Target Validation Report",
        "",
        f"**Target**: {r['target']}",
    ]
    if r["disease"]:
        lines.append(f"**Disease**: {r['disease']}")
    lines += [
        f"**Generated**: {r['retrieved_at']}",
        "",
        "## Scoring summary",
        "",
        f"| Metric | Value |",
        f"|--------|-------|",
        f"| Raw score | {r['raw_score']} / 100 |",
        f"| Safety penalty | {r['safety_penalty']} |",
        f"| **Adjusted score** | **{r['adjusted_score']}** |",
        f"| **Decision** | **{r['decision']}** |",
        "",
        "## Sub-scores",
        "",
        "| Dimension | Score | Confidence | Source |",
        "|-----------|-------|------------|--------|",
    ]
    for dim, obj in r["sub_scores"].items():
        sc = str(obj["score"]) if obj["score"] is not None else "N/A"
        lines.append(f"| {dim.replace('_', ' ').title()} | {sc} / 20 | {obj['confidence']} | {obj['source']} |")

    if r["safety_flags"]:
        lines += [
            "",
            "## Safety flags",
            "",
        ]
        for f in r["safety_flags"]:
            lines.append(f"- **{f['risk']}** (penalty: {f['penalty']}, tier: {f['tier']})")
            if f.get("evidence"):
                lines.append(f"  - Evidence: {f['evidence']}")
            if f.get("mitigation"):
                lines.append(f"  - Mitigation: {f['mitigation']}")

    lines += [
        "",
        "## Rationale",
        "",
        r["rationale"],
        "",
        "---",
        "",
        "*ClawBio is a research and educational tool. It is not a medical device "
        "and does not provide clinical diagnoses. Consult a healthcare professional "
        "before making any medical decisions.*",
    ]

    report_path.write_text("\n".join(lines))
    return report_path


def plot_scoring_summary(result, output_dir: Path) -> Path:
    """Generate a horizontal bar chart of sub-scores with safety penalty."""
    fig_dir = output_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    dims = list(result["sub_scores"].keys())
    scores = [result["sub_scores"][d]["score"] or 0 for d in dims]
    labels = [d.replace("_", " ").title() for d in dims]

    # Colors: green for high, amber for medium, gray for low/null
    colors = []
    for d in dims:
        conf = result["sub_scores"][d]["confidence"]
        if conf == "high":
            colors.append("#1D9E75")
        elif conf == "medium":
            colors.append("#EF9F27")
        else:
            colors.append("#888780")

    fig, ax = plt.subplots(figsize=(9, 5))

    y_pos = np.arange(len(labels))
    bars = ax.barh(y_pos, scores, color=colors, height=0.5, edgecolor="white", linewidth=0.5)

    # Add score labels on bars
    for bar, score in zip(bars, scores):
        ax.text(bar.get_width() + 0.3, bar.get_y() + bar.get_height() / 2,
                f"{score}/20", va="center", fontsize=10, fontweight="medium")

    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=11)
    ax.set_xlim(0, 24)
    ax.set_xlabel("Score", fontsize=11)
    ax.invert_yaxis()

    # Add vertical lines for reference
    ax.axvline(x=14, color="#cccccc", linestyle="--", linewidth=0.8, label="Strong threshold (14)")
    ax.axvline(x=8, color="#eeeeee", linestyle="--", linewidth=0.8, label="Medium threshold (8)")

    # Title with decision
    penalty = result["safety_penalty"]
    ax.set_title(
        f"{result['target']} — Raw: {result['raw_score']}  |  "
        f"Safety: {penalty}  |  Adjusted: {result['adjusted_score']}  →  {result['decision']}",
        fontsize=12, fontweight="medium", pad=12
    )

    # Legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor="#1D9E75", label="High confidence"),
        Patch(facecolor="#EF9F27", label="Medium confidence"),
        Patch(facecolor="#888780", label="Low / no data"),
    ]
    if penalty < 0:
        legend_elements.append(Patch(facecolor="#E24B4A", label=f"Safety penalty ({penalty})"))
    ax.legend(handles=legend_elements, loc="lower right", fontsize=9, framealpha=0.9)

    plt.tight_layout()
    fig_path = fig_dir / "scoring_summary.png"
    fig.savefig(fig_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return fig_path


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(description="Target Validation Scorer")
    parser.add_argument("--input", type=Path, required=False, help="Input JSON file")
    parser.add_argument("--output", type=Path, default=Path("/tmp/target_validation_output"), help="Output directory")
    parser.add_argument("--demo", action="store_true", help="Run with built-in TGFBR1/IPF demo data")
    args = parser.parse_args()

    if args.demo:
        input_path = DEMO_DATA
    elif args.input:
        input_path = args.input
    else:
        parser.error("Provide --input <file> or use --demo")

    print(f"Reading: {input_path}")
    with open(input_path) as f:
        data = json.load(f)

    print(f"Validating target: {data.get('target', 'unknown')}")
    result = validate_target(data)

    # Write JSON report
    args.output.mkdir(parents=True, exist_ok=True)
    json_path = args.output / "validation_report.json"
    with open(json_path, "w") as f:
        json.dump(result, f, indent=2)

    # Write markdown report
    report_path = generate_report(result, args.output)

    # Generate figure
    fig_path = plot_scoring_summary(result, args.output)

    # Print summary to stdout
    print(f"\n{'='*60}")
    print(f"TARGET VALIDATION REPORT")
    print(f"{'='*60}")
    print(f"Target:          {result['target']}")
    if result["disease"]:
        print(f"Disease:         {result['disease']}")
    print(f"Raw score:       {result['raw_score']} / 100")
    print(f"Safety penalty:  {result['safety_penalty']}")
    print(f"Adjusted score:  {result['adjusted_score']}")
    print(f"Decision:        {result['decision']}")
    print(f"{'='*60}")
    print(f"\nSub-scores:")
    for dim, obj in result["sub_scores"].items():
        sc = str(obj["score"]) if obj["score"] is not None else "N/A"
        print(f"  {dim.replace('_', ' ').title():25s}  {sc:>4s} / 20  ({obj['confidence']})")
    if result["safety_flags"]:
        print(f"\nSafety flags:")
        for flag in result["safety_flags"]:
            print(f"  - {flag['risk']} (penalty: {flag['penalty']})")
    print(f"\nRationale:")
    print(f"  {result['rationale']}")
    print(f"\nReport:  {report_path}")
    print(f"JSON:    {json_path}")
    print(f"Figure:  {fig_path}")


if __name__ == "__main__":
    main()
