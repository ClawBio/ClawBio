"""
generate_report.py — Markdown report + matplotlib figures for NutriGx Advisor
"""

import os
import json
from datetime import datetime, timezone
from pathlib import Path


DOMAIN_LABELS = {
    "folate": "Folate / B-Vitamins",
    "vitamin_d": "Vitamin D",
    "omega3": "Omega-3 / LC-PUFA",
    "vitamin_a": "Vitamin A (Beta-carotene)",
    "vitamin_c": "Vitamin C",
    "vitamin_b6": "Vitamin B6",
    "fat_metabolism": "Fat Metabolism",
    "carbohydrate": "Carbohydrate Metabolism",
    "caffeine": "Caffeine Metabolism",
    "alcohol": "Alcohol Metabolism",
    "lactose": "Lactose Tolerance",
    "gluten": "Gluten Sensitivity",
    "antioxidant": "Antioxidant / Detox",
}

RECOMMENDATIONS = {
    "folate": {
        "Low": "Standard dietary folate (leafy greens, legumes) is sufficient.",
        "Moderate": "Consider increasing dietary folate. If homocysteine is elevated, 5-MTHF (methylfolate) supplement preferred over folic acid.",
        "Elevated": "Clinical assessment of homocysteine and B12 recommended. 5-MTHF + methylcobalamin supplementation often indicated. Avoid high-dose synthetic folic acid.",
    },
    "vitamin_d": {
        "Low": "Standard sun exposure and dietary vitamin D3 are likely adequate.",
        "Moderate": "Serum 25(OH)D testing recommended. Supplement D3 (1000–2000 IU/day) if <75 nmol/L.",
        "Elevated": "Genetic predisposition to low vitamin D binding efficiency. Test and maintain 25(OH)D >100 nmol/L. D3+K2 co-supplementation advisable.",
    },
    "omega3": {
        "Low": "Dietary ALA (flaxseed, walnuts) efficiently converted. Two oily fish servings/week sufficient.",
        "Moderate": "Reduced LC-PUFA synthesis capacity. Increase direct EPA/DHA sources (oily fish, algae oil). Consider omega-6:omega-3 ratio <4:1.",
        "Elevated": "Markedly reduced FADS1/2 or ELOVL2 activity. Direct EPA+DHA supplementation (1–3 g/day algae or fish oil) recommended. Minimise linoleic acid (LA) competition.",
    },
    "vitamin_a": {
        "Low": "Beta-carotene conversion is efficient. Plant-based vitamin A sources adequate.",
        "Moderate": "Reduced BCMO1 activity. Include preformed vitamin A (liver, eggs, dairy) alongside carotenoid-rich vegetables.",
        "Elevated": "Substantially impaired beta-carotene → retinol conversion. Preformed vitamin A (retinol) from animal sources or supplements essential.",
    },
    "vitamin_c": {
        "Low": "Standard dietary vitamin C (citrus, peppers, kiwi) sufficient.",
        "Moderate": "May benefit from higher dietary vitamin C intake (200–500 mg/day from food).",
        "Elevated": "Reduced renal reabsorption of vitamin C. Higher dietary intake or supplementation (500–1000 mg/day) may be warranted.",
    },
    "fat_metabolism": {
        "Low": "No specific fat quality concerns from genetics alone.",
        "Moderate": "Monitor LDL response to saturated fat intake. Mediterranean-style fat profile advisable.",
        "Elevated": "APOE ε4 or PPARG variants detected. Significant LDL-C elevation on high saturated fat diet likely. Prioritise MUFA, PUFA; limit SFA to <7% of energy.",
    },
    "carbohydrate": {
        "Low": "Standard carbohydrate intake appropriate.",
        "Moderate": "Glycaemic index awareness recommended. Favour whole-grain, low-GI carbohydrates.",
        "Elevated": "TCF7L2 / FTO variants suggest elevated metabolic risk on high-GI carbohydrate intake. Restrict refined carbohydrates; increase dietary fibre.",
    },
    "caffeine": {
        "Low": "Fast CYP1A2 metaboliser. Moderate caffeine (≤400 mg/day) poses low cardiovascular risk.",
        "Moderate": "Intermediate metaboliser. Caffeine ≤200 mg/day recommended.",
        "Elevated": "Slow metaboliser. Caffeine intake >200 mg/day associated with increased MI risk in carriers. Consider limiting to 1 cup coffee/day.",
    },
    "alcohol": {
        "Low": "Standard alcohol metabolism. General guidelines apply.",
        "Moderate": "Some acetaldehyde accumulation risk. Limit alcohol consumption.",
        "Elevated": "ADH1B / ALDH2 variant(s) detected. Significantly increased acetaldehyde toxicity risk. Strong recommendation to minimise or avoid alcohol.",
    },
    "lactose": {
        "Low": "Likely lactase-persistent. Can tolerate dairy normally.",
        "Moderate": "Possible partial lactase non-persistence. Monitor tolerance to high-lactose foods.",
        "Elevated": "Predicted lactase non-persistence. Limit fresh milk; fermented dairy (yoghurt, hard cheese) typically well tolerated. Calcium intake from non-dairy sources or supplements.",
    },
    "gluten": {
        "Low": "No HLA-DQ risk haplotypes detected. Gluten restriction not genetically indicated.",
        "Moderate": "HLA-DQ risk haplotype(s) present. Genetic predisposition to coeliac disease. Clinical testing (tTG-IgA) recommended if symptomatic.",
        "Elevated": "High-risk HLA-DQ2/DQ8 haplotype(s). Coeliac disease screening strongly recommended. Do not begin gluten-free diet before testing.",
    },
    "antioxidant": {
        "Low": "Standard antioxidant-rich diet (fruits, vegetables) sufficient.",
        "Moderate": "Some oxidative stress pathway variants. Increase dietary antioxidants; selenium and CoQ10 food sources beneficial.",
        "Elevated": "Multiple antioxidant enzyme variants detected. Consider assessed supplementation: selenium, CoQ10, and magnesium. Limit pro-oxidant exposures (smoking, excess alcohol).",
    },
}


def generate_report(snp_calls, risk_scores, snp_panel, output_dir, figures=True, input_file=""):
    """Generate Markdown report and optional figures. Returns path to report file."""
    output_dir = Path(output_dir)
    timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")

    lines = []

    # ── Header ────────────────────────────────────────────────────────────────
    lines += [
        "# NutriGx Personalised Nutrition Report",
        "",
        f"**Generated**: {timestamp}  ",
        f"**Tool**: ClawBio NutriGx Advisor v0.2.0  ",
        f"**Input**: `{Path(input_file).name}`  ",
        "",
        "> **Disclaimer**: This report is for research and educational purposes only. "
        "It does not constitute medical advice. Consult a registered dietitian or clinical "
        "geneticist before making significant dietary changes or starting supplements.",
        "",
    ]

    # ── Executive Summary ─────────────────────────────────────────────────────
    elevated = [(d, v) for d, v in risk_scores.items() if v["category"] == "Elevated"]
    moderate = [(d, v) for d, v in risk_scores.items() if v["category"] == "Moderate"]

    lines += ["## Executive Summary", ""]
    total_tested = sum(v["tested_snps"] for v in risk_scores.values())
    total_panel = len(snp_panel)
    lines.append(f"Analysed **{total_tested}** of **{total_panel}** panel SNPs from your genetic data.")
    lines.append("")

    if elevated:
        lines.append(f"**{len(elevated)} nutrient domain(s) at elevated genetic risk:**")
        for d, v in elevated:
            label = DOMAIN_LABELS.get(d, d)
            lines.append(f"- {label} (score: {v['score']}/10)")
        lines.append("")

    if moderate:
        lines.append(f"**{len(moderate)} domain(s) at moderate genetic risk:**")
        for d, v in moderate:
            label = DOMAIN_LABELS.get(d, d)
            lines.append(f"- {label} (score: {v['score']}/10)")
        lines.append("")

    lines += [
        "---",
        "",
        "## Nutrient Risk Score Overview",
        "",
        "| Nutrient Domain | Score (0–10) | Risk Category |",
        "|-----------------|:------------:|:-------------:|",
    ]
    for domain, data in sorted(risk_scores.items(), key=lambda x: -(x[1]["score"] or 0)):
        label = DOMAIN_LABELS.get(domain, domain)
        score = data["score"] if data["score"] is not None else "N/A"
        cat = data["category"]
        coverage = data.get("coverage", "")
        emoji = {"Low": "🟢", "Moderate": "🟡", "Elevated": "🔴", "Unknown": "⚪"}.get(cat, "")
        cat_display = f"{cat} ({coverage})" if coverage else cat
        lines.append(f"| {label} | {score} | {emoji} {cat_display} |")
    lines += ["", "---", ""]

    # ── Per-Domain Sections ───────────────────────────────────────────────────
    lines.append("## Detailed Findings by Nutrient Domain")
    lines.append("")

    for domain, data in sorted(risk_scores.items(), key=lambda x: -(x[1]["score"] or 0)):
        label = DOMAIN_LABELS.get(domain, domain)
        score = data["score"]
        cat = data["category"]
        if cat == "Unknown":
            rec = "Insufficient genetic data to assess this domain."
        else:
            rec = RECOMMENDATIONS.get(domain, {}).get(cat, "No specific recommendation available.")

        coverage = data.get("coverage", "")
        lines += [
            f"### {label}",
            "",
            f"**Risk Score**: {score}/10 — **{cat}**  ",
            f"**SNPs tested**: {data['tested_snps']} | **Not on chip**: {data['missing_snps']} | **Coverage**: {coverage}",
            "",
        ]

        if data["contributing_snps"]:
            lines += [
                "| Gene | rsID | Genotype | Risk Alleles | Effect |",
                "|------|------|----------|:------------:|--------|",
            ]
            for s in data["contributing_snps"]:
                effect = s["effect_direction"].replace("_", " ").title()
                lines.append(
                    f"| {s['gene']} | {s['rsid']} | `{s['genotype']}` "
                    f"| {s['risk_count']}/2 | {effect} |"
                )
            lines.append("")

        lines += [
            "**Recommendation**",
            "",
            f"> {rec}",
            "",
            "---",
            "",
        ]

    # ── Supplement Interactions ───────────────────────────────────────────────
    lines += [
        "## Supplement Interaction Notes",
        "",
        "| Supplement | Relevant Genes | Caution |",
        "|------------|---------------|---------|",
        "| Folic acid (synthetic) | MTHFR | Prefer 5-MTHF if MTHFR variant present |",
        "| Vitamin D3 | VDR, GC | K2 co-supplementation improves utilisation |",
        "| Fish oil / EPA+DHA | FADS1, FADS2, ELOVL2 | Higher dose needed if FADS variants present |",
        "| Vitamin A (retinol) | BCMO1 | Do not exceed UL (3000 μg RAE/day) |",
        "| Iron | HFE (not in panel) | Not assessed; request separate iron panel |",
        "| CoQ10 | NQO1 | Ubiquinol form preferred if NQO1 P187S homozygous |",
        "",
        "---",
        "",
    ]

    # ── Footer ────────────────────────────────────────────────────────────────
    lines += [
        "## Reproducibility",
        "",
        "This report was generated deterministically. See `commands.sh` and `environment.yml`",
        "in the output directory to reproduce this analysis on any machine.",
        "",
        "## References",
        "",
        "SNP-nutrient associations sourced from GWAS Catalog, ClinVar, and CPIC guidelines.",
        "Full citations available in `skills/nutrigx-advisor/SKILL.md`.",
        "",
    ]

    report_text = "\n".join(lines)
    report_path = output_dir / "nutrigx_report.md"
    report_path.write_text(report_text)

    if figures:
        _generate_figures(risk_scores, output_dir)

    return str(report_path)


def _generate_figures(risk_scores: dict, output_dir: Path):
    """Generate radar chart and heatmap."""
    try:
        import numpy as np
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("[NutriGx] matplotlib/numpy not available — skipping figures")
        return

    # ── Radar Chart ───────────────────────────────────────────────────────────
    scored = {d: v for d, v in risk_scores.items() if v["score"] is not None}
    labels = [DOMAIN_LABELS.get(d, d) for d in scored]
    values = [v["score"] for v in scored.values()]

    if len(labels) < 3:
        return

    N = len(labels)
    angles = [n / float(N) * 2 * np.pi for n in range(N)]
    angles += angles[:1]
    values_plot = values + values[:1]

    fig, ax = plt.subplots(figsize=(8, 8), subplot_kw=dict(polar=True))
    ax.set_theta_offset(np.pi / 2)
    ax.set_theta_direction(-1)
    plt.xticks(angles[:-1], labels, color="grey", size=8)
    ax.set_rlabel_position(0)
    plt.yticks([2, 4, 6, 8, 10], ["2", "4", "6", "8", "10"], color="grey", size=7)
    plt.ylim(0, 10)
    ax.plot(angles, values_plot, linewidth=2, linestyle="solid", color="#2196F3")
    ax.fill(angles, values_plot, alpha=0.25, color="#2196F3")

    # Colour zones
    for r, color, alpha in [(3.5, "green", 0.05), (6.5, "orange", 0.05), (10, "red", 0.05)]:
        ax.fill_between(angles, 0, r, alpha=alpha, color=color)

    plt.title("NutriGx Nutrient Risk Profile", size=14, fontweight="bold", pad=20)
    plt.tight_layout()
    fig.savefig(output_dir / "nutrigx_radar.png", dpi=150, bbox_inches="tight")
    plt.close()

    # ── Heatmap ───────────────────────────────────────────────────────────────
    try:
        import seaborn as sns
        import pandas as pd

        rows = []
        for domain, data in risk_scores.items():
            for snp in data["contributing_snps"]:
                rows.append({
                    "Gene": snp["gene"],
                    "Nutrient": DOMAIN_LABELS.get(domain, domain),
                    "Score": snp["raw_score"],
                })

        if rows:
            df = pd.DataFrame(rows)
            pivot = df.pivot_table(index="Gene", columns="Nutrient", values="Score", aggfunc="mean")
            pivot = pivot.fillna(0)

            fig, ax = plt.subplots(figsize=(max(8, len(pivot.columns)), max(6, len(pivot.index) * 0.6)))
            sns.heatmap(pivot, annot=True, fmt=".1f", cmap="RdYlGn_r",
                        vmin=0, vmax=1, linewidths=0.5, ax=ax,
                        cbar_kws={"label": "Risk Score (0=Ref, 0.5=Het, 1=Hom Risk)"})
            ax.set_title("Gene × Nutrient Risk Heatmap", fontsize=13, fontweight="bold")
            plt.tight_layout()
            fig.savefig(output_dir / "nutrigx_heatmap.png", dpi=150, bbox_inches="tight")
            plt.close()
    except ImportError:
        pass
