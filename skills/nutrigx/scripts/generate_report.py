#!/usr/bin/env python3
"""
generate_report.py — DAG step 4/4: Markdown report + figures + result.json.

Ported from ClawBio 1.x skills/nutrigx/generate_report.py (report/figure
generation logic unchanged) plus the result.json envelope that 1.x's
nutrigx.py main() wrote via clawbio.common.report.write_result_json — that
shared-library helper is inlined here as _write_result_json (dependency-free,
no import of the 1.x-only clawbio.common package). This step also folds in
the SHA-256 input checksum that 1.x computed via clawbio.common.checksums.

1.x's repro_bundle.py (commands.sh / environment.yml / checksums.txt /
provenance.json) was dropped entirely — this repo's own `clawbio bundle`
is the authoritative, run-level reproducibility record, so a per-skill
duplicate would be redundant (same reasoning as claw_metagenomics and
rnaseq_de's ports).

Deliberate deviation from 1.x for CLI-contract reasons: 1.x's
_generate_figures() silently skipped the radar chart when fewer than 3
domains had a score, and silently skipped the heatmap when there were no
contributing SNPs. Under `clawbio done --output`, a declared `emits` key
whose file doesn't exist is a hard failure — there is no "this step
produced fewer outputs than declared" escape hatch (see pipeline-schema.md).
So both guards were relaxed to always write a PNG: a real chart when there
is enough data, otherwise a small placeholder figure stating why. The
underlying plotting logic for the "enough data" case is untouched.

Usage:
    python3 generate_report.py --snp-calls snp_calls.json --risk-scores risk_scores.json \\
        --panel snp_panel.json --input-file genome.csv --out OUT_DIR
"""

import argparse
import hashlib
import json
import sys
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
        "Moderate": "Serum 25(OH)D testing recommended. Supplement D3 (1000-2000 IU/day) if <75 nmol/L.",
        "Elevated": "Genetic predisposition to low vitamin D binding efficiency. Test and maintain 25(OH)D >100 nmol/L. D3+K2 co-supplementation advisable.",
    },
    "omega3": {
        "Low": "Dietary ALA (flaxseed, walnuts) efficiently converted. Two oily fish servings/week sufficient.",
        "Moderate": "Reduced LC-PUFA synthesis capacity. Increase direct EPA/DHA sources (oily fish, algae oil). Consider omega-6:omega-3 ratio <4:1.",
        "Elevated": "Markedly reduced FADS1/2 or ELOVL2 activity. Direct EPA+DHA supplementation (1-3 g/day algae or fish oil) recommended. Minimise linoleic acid (LA) competition.",
    },
    "vitamin_a": {
        "Low": "Beta-carotene conversion is efficient. Plant-based vitamin A sources adequate.",
        "Moderate": "Reduced BCMO1 activity. Include preformed vitamin A (liver, eggs, dairy) alongside carotenoid-rich vegetables.",
        "Elevated": "Substantially impaired beta-carotene -> retinol conversion. Preformed vitamin A (retinol) from animal sources or supplements essential.",
    },
    "vitamin_c": {
        "Low": "Standard dietary vitamin C (citrus, peppers, kiwi) sufficient.",
        "Moderate": "May benefit from higher dietary vitamin C intake (200-500 mg/day from food).",
        "Elevated": "Reduced renal reabsorption of vitamin C. Higher dietary intake or supplementation (500-1000 mg/day) may be warranted.",
    },
    "fat_metabolism": {
        "Low": "No specific fat quality concerns from genetics alone.",
        "Moderate": "Monitor LDL response to saturated fat intake. Mediterranean-style fat profile advisable.",
        "Elevated": "APOE e4 or PPARG variants detected. Significant LDL-C elevation on high saturated fat diet likely. Prioritise MUFA, PUFA; limit SFA to <7% of energy.",
    },
    "carbohydrate": {
        "Low": "Standard carbohydrate intake appropriate.",
        "Moderate": "Glycaemic index awareness recommended. Favour whole-grain, low-GI carbohydrates.",
        "Elevated": "TCF7L2 / FTO variants suggest elevated metabolic risk on high-GI carbohydrate intake. Restrict refined carbohydrates; increase dietary fibre.",
    },
    "caffeine": {
        "Low": "Fast CYP1A2 metaboliser. Moderate caffeine (<=400 mg/day) poses low cardiovascular risk.",
        "Moderate": "Intermediate metaboliser. Caffeine <=200 mg/day recommended.",
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


def sha256_file(filepath: str) -> str:
    h = hashlib.sha256()
    with open(filepath, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def generate_report(snp_calls, risk_scores, snp_panel, output_dir, input_file=""):
    """Generate Markdown report. Returns path to report file."""
    output_dir = Path(output_dir)
    timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")

    lines = []

    # -- Header -------------------------------------------------------------
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

    # -- Executive Summary ----------------------------------------------------
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
        "| Nutrient Domain | Score (0-10) | Risk Category |",
        "|-----------------|:------------:|:-------------:|",
    ]
    for domain, data in sorted(risk_scores.items(), key=lambda x: -(x[1]["score"] or 0)):
        label = DOMAIN_LABELS.get(domain, domain)
        score = data["score"] if data["score"] is not None else "N/A"
        cat = data["category"]
        coverage = data.get("coverage", "")
        emoji = {"Low": "\U0001F7E2", "Moderate": "\U0001F7E1", "Elevated": "\U0001F534", "Unknown": "⚪"}.get(cat, "")
        cat_display = f"{cat} ({coverage})" if coverage else cat
        lines.append(f"| {label} | {score} | {emoji} {cat_display} |")
    lines += ["", "![NutriGx Nutrient Risk Profile](nutrigx_radar.png)", "", "![Gene x Nutrient Risk Heatmap](nutrigx_heatmap.png)", "", "---", ""]

    # -- Per-Domain Sections ----------------------------------------------------
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
            f"**Risk Score**: {score}/10 -- **{cat}**  ",
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

    # -- Supplement Interactions ------------------------------------------------
    lines += [
        "## Supplement Interaction Notes",
        "",
        "| Supplement | Relevant Genes | Caution |",
        "|------------|---------------|---------|",
        "| Folic acid (synthetic) | MTHFR | Prefer 5-MTHF if MTHFR variant present |",
        "| Vitamin D3 | VDR, GC | K2 co-supplementation improves utilisation |",
        "| Fish oil / EPA+DHA | FADS1, FADS2, ELOVL2 | Higher dose needed if FADS variants present |",
        "| Vitamin A (retinol) | BCMO1 | Do not exceed UL (3000 ug RAE/day) |",
        "| Iron | HFE (not in panel) | Not assessed; request separate iron panel |",
        "| CoQ10 | NQO1 | Ubiquinol form preferred if NQO1 P187S homozygous |",
        "",
        "---",
        "",
    ]

    # -- Footer -------------------------------------------------------------
    lines += [
        "## Reproducibility",
        "",
        "This report was generated deterministically as one step of a `clawbio` run. See",
        "`clawbio bundle`'s `commands.sh` / `checksums.sha256` / `environment.yml` for the",
        "full reproducibility record of this run.",
        "",
        "## References",
        "",
        "SNP-nutrient associations sourced from GWAS Catalog, ClinVar, and CPIC guidelines.",
        "Full citations available in `skills/nutrigx/SKILL.md`.",
        "",
    ]

    report_text = "\n".join(lines)
    report_path = output_dir / "nutrigx_report.md"
    report_path.write_text(report_text)

    return str(report_path)


def _placeholder_figure(path: Path, title: str, message: str):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.axis("off")
    ax.set_title(title, fontsize=13, fontweight="bold")
    ax.text(0.5, 0.5, message, ha="center", va="center", wrap=True, fontsize=11)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def generate_figures(risk_scores: dict, output_dir: Path):
    """Generate radar chart and heatmap. Always writes both files (see module
    docstring for why the 1.x early-return guards were relaxed to placeholders
    instead of silently skipping)."""
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    radar_path = output_dir / "nutrigx_radar.png"
    heatmap_path = output_dir / "nutrigx_heatmap.png"

    # -- Radar Chart ----------------------------------------------------------
    scored = {d: v for d, v in risk_scores.items() if v["score"] is not None}
    labels = [DOMAIN_LABELS.get(d, d) for d in scored]
    values = [v["score"] for v in scored.values()]

    if len(labels) < 3:
        _placeholder_figure(
            radar_path,
            "NutriGx Nutrient Risk Profile",
            f"Only {len(labels)} nutrient domain(s) had a score (need >= 3 for a radar chart).\n"
            "Increase panel coverage for a full risk profile.",
        )
    else:
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
        fig.savefig(radar_path, dpi=150, bbox_inches="tight")
        plt.close(fig)

    # -- Heatmap ----------------------------------------------------------------
    rows = []
    for domain, data in risk_scores.items():
        for snp in data["contributing_snps"]:
            rows.append({
                "Gene": snp["gene"],
                "Nutrient": DOMAIN_LABELS.get(domain, domain),
                "Score": snp["raw_score"],
            })

    if not rows:
        _placeholder_figure(
            heatmap_path,
            "Gene x Nutrient Risk Heatmap",
            "No contributing SNPs found in the input data.\n"
            "No gene x nutrient risk relationships to plot.",
        )
    else:
        import seaborn as sns
        import pandas as pd

        df = pd.DataFrame(rows)
        pivot = df.pivot_table(index="Gene", columns="Nutrient", values="Score", aggfunc="mean")
        pivot = pivot.fillna(0)

        fig, ax = plt.subplots(figsize=(max(8, len(pivot.columns)), max(6, len(pivot.index) * 0.6)))
        sns.heatmap(pivot, annot=True, fmt=".1f", cmap="RdYlGn_r",
                    vmin=0, vmax=1, linewidths=0.5, ax=ax,
                    cbar_kws={"label": "Risk Score (0=Ref, 0.5=Het, 1=Hom Risk)"})
        ax.set_title("Gene x Nutrient Risk Heatmap", fontsize=13, fontweight="bold")
        plt.tight_layout()
        fig.savefig(heatmap_path, dpi=150, bbox_inches="tight")
        plt.close(fig)

    return str(radar_path), str(heatmap_path)


def _write_result_json(output_dir, skill, version, summary, data, input_checksum=""):
    """Inlined equivalent of 1.x's clawbio.common.report.write_result_json —
    that module also pulled in a home-directory audit log writer
    (clawbio.common.audit) which has no place in a reproducible, sandboxed
    `clawbio` step, so only the result.json envelope itself was ported."""
    output_dir = Path(output_dir)
    envelope = {
        "skill": skill,
        "version": version,
        "completed_at": datetime.now(timezone.utc).isoformat(),
        "input_checksum": f"sha256:{input_checksum}" if input_checksum else "",
        "summary": summary,
        "data": data,
    }
    result_path = output_dir / "result.json"
    result_path.write_text(json.dumps(envelope, indent=2, default=str))
    return result_path


def main():
    ap = argparse.ArgumentParser(description="Render NutriGx Markdown report + figures + result.json")
    ap.add_argument("--snp-calls", required=True, help="Path to snp_calls.json from extract_genotypes.py")
    ap.add_argument("--risk-scores", required=True, help="Path to risk_scores.json from score_variants.py")
    ap.add_argument("--panel", required=True, help="Path to snp_panel.json")
    ap.add_argument("--input-file", required=True, help="Original genetic data file path (for report header + checksum)")
    ap.add_argument("--out", required=True, help="Output directory")
    args = ap.parse_args()

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    with open(args.snp_calls) as f:
        snp_calls = json.load(f)
    with open(args.risk_scores) as f:
        risk_scores = json.load(f)
    with open(args.panel) as f:
        snp_panel = json.load(f)

    print("[generate_report] Generating report ...")
    report_path = generate_report(
        snp_calls=snp_calls,
        risk_scores=risk_scores,
        snp_panel=snp_panel,
        output_dir=str(out_dir),
        input_file=args.input_file,
    )

    print("[generate_report] Generating figures ...")
    generate_figures(risk_scores, out_dir)

    print("[generate_report] Writing result.json ...")
    elevated_domains = [d for d, v in risk_scores.items() if v["category"] == "Elevated"]
    moderate_domains = [d for d, v in risk_scores.items() if v["category"] == "Moderate"]
    total_tested = sum(v["tested_snps"] for v in risk_scores.values())
    total_missing = sum(v["missing_snps"] for v in risk_scores.values())

    input_checksum = ""
    if Path(args.input_file).exists():
        input_checksum = sha256_file(args.input_file)

    _write_result_json(
        output_dir=out_dir,
        skill="nutrigx",
        version="0.2.0",
        summary={
            "panel_snps_tested": total_tested,
            "panel_snps_missing": total_missing,
            "panel_size": len(snp_panel),
            "elevated_domains": elevated_domains,
            "moderate_domains": moderate_domains,
            "domains_assessed": len(risk_scores),
        },
        data={
            "risk_scores": risk_scores,
            "snp_calls": snp_calls,
        },
        input_checksum=input_checksum,
    )

    print(f"[generate_report] Done. Report: {report_path}")
    print(f"[generate_report] Results in: {out_dir}/")


if __name__ == "__main__":
    main()
