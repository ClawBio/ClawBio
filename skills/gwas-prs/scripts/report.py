#!/usr/bin/env python3
"""Stage 4 of the gwas_prs DAG: render the Markdown report + result.json envelope.

Port of ClawBio 1.x `gwas-prs` (skills/gwas-prs/gwas_prs.py::generate_report).
Report structure and content are unchanged from the 1.x source; only the
I/O boundary changed (reads scripts/calculate_prs.py's prs_results.json
instead of holding results in-process).

No dependency on the `clawbio.common` shared library (that package does not
exist in this repo) — the small write_result_json/DISCLAIMER envelope it
provided is inlined below, same pattern as skills/claw_ancestry_pca and
skills/rnaseq_de.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

DISCLAIMER = (
    "ClawBio is a research and educational tool. It is not a medical device "
    "and does not provide clinical diagnoses. Polygenic risk scores reflect "
    "statistical associations from population studies and do not determine "
    "individual outcomes. Consult a healthcare professional before making "
    "any medical decisions based on genetic information."
)


def sha256_hex(filepath: str | Path, length: int = 16) -> str:
    h = hashlib.sha256()
    with open(filepath, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()[:length]


def write_result_json(output_path: Path, skill: str, version: str, summary: dict,
                       data: dict, input_checksum: str = "") -> None:
    envelope = {
        "skill": skill,
        "version": version,
        "completed_at": datetime.now(timezone.utc).isoformat(),
        "input_checksum": f"sha256:{input_checksum}" if input_checksum else "",
        "datasets": {},
        "summary": summary,
        "data": data,
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(envelope, indent=2, default=str))


def generate_report(results: list[dict], input_info: dict) -> str:
    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    lines = [
        "# GWAS Polygenic Risk Score Report",
        "",
        f"**Date**: {now}",
        f"**Tool**: ClawBio GWAS-PRS (2.0 port)",
        f"**Input file**: {input_info.get('filepath', 'N/A')}",
        f"**Format detected**: {input_info.get('format', 'N/A')}",
        f"**Total SNPs in file**: {input_info.get('total_snps', 'N/A'):,}",
        f"**Genome build**: GRCh37",
        "",
        "---",
        "",
        "## Summary",
        "",
        "| PGS ID | Trait | Raw Score | Percentile | Risk Category | Overlap |",
        "|--------|-------|-----------|------------|---------------|---------|",
    ]

    for r in results:
        pgs_id = r["pgs_id"]
        trait = r.get("trait", "Unknown")
        prs = r["prs"]
        pct_info = r["percentile_info"]
        raw = f"{prs['raw_score']:.4f}"
        pct = f"{pct_info['percentile']:.1f}%" if pct_info["percentile"] is not None else "N/A"
        risk = pct_info.get("risk_category") or "N/A"
        overlap = f"{prs['overlap_fraction'] * 100:.1f}%"
        lines.append(f"| {pgs_id} | {trait} | {raw} | {pct} | {risk} | {overlap} |")

    lines.append("")

    low_overlap = [r for r in results if r["prs"]["overlap_fraction"] < 0.70]
    if low_overlap:
        lines.append("### Overlap Warnings")
        lines.append("")
        lines.append(
            "> **Warning**: The following scores have < 70% SNP overlap with "
            "your genotype file. Results should be interpreted with caution."
        )
        lines.append("")
        for r in low_overlap:
            prs = r["prs"]
            lines.append(
                f"- **{r['pgs_id']}** ({r.get('trait', '')}): "
                f"{prs['variants_used']}/{prs['variants_total']} variants "
                f"({prs['overlap_fraction'] * 100:.1f}% overlap)"
            )
        lines.append("")

    lines.append("---")
    lines.append("")
    lines.append("## Score Details")
    lines.append("")

    for r in results:
        pgs_id = r["pgs_id"]
        trait = r.get("trait", "Unknown")
        prs = r["prs"]
        pct_info = r["percentile_info"]
        metadata = r.get("metadata", {})

        lines.append(f"### {pgs_id} — {trait}")
        lines.append("")
        lines.append("| Property | Value |")
        lines.append("|----------|-------|")
        lines.append(f"| **PGS ID** | {pgs_id} |")
        lines.append(f"| **Trait** | {trait} |")
        if metadata.get("publication"):
            lines.append(f"| **Publication** | {metadata['publication']} |")
        lines.append(f"| **Variants in score** | {prs['variants_total']} |")
        lines.append(f"| **Variants matched** | {prs['variants_used']} |")
        lines.append(f"| **Variants missing** | {prs['variants_missing']} |")
        lines.append(f"| **Overlap** | {prs['overlap_fraction'] * 100:.1f}% |")
        lines.append(f"| **Raw PRS** | {prs['raw_score']:.6f} |")
        if pct_info["percentile"] is not None:
            lines.append(f"| **Percentile** | {pct_info['percentile']:.1f}% |")
            lines.append(f"| **Risk category** | {pct_info['risk_category']} |")
            lines.append(f"| **Z-score** | {pct_info['z_score']:.3f} |")
            lines.append(f"| **Estimation method** | {pct_info['method']} |")
            if pct_info.get("reference_population"):
                lines.append(f"| **Reference population** | {pct_info['reference_population']} |")
        else:
            lines.append("| **Percentile** | N/A (no reference available) |")
        lines.append("")

        scored = [v for v in prs["per_variant"] if v["status"] == "scored"]
        if scored:
            scored_sorted = sorted(scored, key=lambda v: abs(v["contribution"]), reverse=True)
            show_count = min(len(scored_sorted), 20)
            lines.append(f"**Top {show_count} contributing variants** (of {len(scored_sorted)} scored):")
            lines.append("")
            lines.append("| rsID | Genotype | Effect Allele | Dosage | Weight | Contribution |")
            lines.append("|------|----------|---------------|--------|--------|--------------|")
            for v in scored_sorted[:show_count]:
                lines.append(
                    f"| {v['rsid']} | {v['genotype']} | {v['effect_allele']} | "
                    f"{v['dosage']} | {v['weight']:.4f} | {v['contribution']:.4f} |"
                )
            if len(scored_sorted) > show_count:
                lines.append(f"| ... | | | | | ({len(scored_sorted) - show_count} more variants) |")
            lines.append("")
        else:
            lines.append("*No variants from this score were found in your genotype file.*")
            lines.append("")

    lines.append("---")
    lines.append("")
    lines.append("## Methods")
    lines.append("")
    lines.append(
        "Polygenic Risk Scores (PRS) were calculated as the weighted sum of "
        "effect allele dosages across all available variants:"
    )
    lines.append("")
    lines.append("```")
    lines.append("PRS = sum(effect_weight_i * dosage_i)")
    lines.append("```")
    lines.append("")
    lines.append(
        "Where `dosage_i` is the number of copies (0, 1, or 2) of the effect "
        "allele at variant *i*, and `effect_weight_i` is the log odds ratio "
        "(or beta coefficient) from the original GWAS."
    )
    lines.append("")
    lines.append(
        "Percentile estimates use one of two methods: (1) curated EUR "
        "reference distributions (mean and SD from published cohorts), or "
        "(2) population-level allele frequency estimates from the PGS Catalog "
        "scoring file to compute expected score and variance under "
        "Hardy-Weinberg equilibrium."
    )
    lines.append("")
    lines.append("**Scoring files**: PGS Catalog (https://www.pgscatalog.org/)")
    lines.append("")
    lines.append("**Genome build**: GRCh37")
    lines.append("")

    lines.append("## Limitations")
    lines.append("")
    lines.append(
        "- PRS calculated from DTC genotype data typically covers only a "
        "subset of the variants in the full score, reducing predictive power."
    )
    lines.append(
        "- Percentile estimates assume a European (EUR) reference population "
        "unless otherwise noted. PRS accuracy and transferability vary "
        "significantly across ancestries."
    )
    lines.append(
        "- DTC genotype arrays may not include all variants in a PGS, "
        "especially rare variants or indels."
    )
    lines.append(
        "- Risk categories are simplified and do not account for family "
        "history, lifestyle, or environmental factors."
    )
    lines.append("")

    lines.append("## Disclaimer")
    lines.append("")
    lines.append(f"*{DISCLAIMER}*")
    lines.append("")

    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description="Render PRS Markdown report + result.json")
    parser.add_argument("--results", required=True, help="Path to prs_results.json (from calculate_prs.py)")
    parser.add_argument("--genotype-json", required=True, help="Path to genotype.json (for input checksum)")
    parser.add_argument("--out-report", required=True, help="Output path for report.md")
    parser.add_argument("--out-result-json", required=True, help="Output path for result.json")
    args = parser.parse_args()

    payload = json.loads(Path(args.results).read_text())
    results = payload["results"]
    input_info = payload["input_info"]

    report_text = generate_report(results, input_info)
    out_report = Path(args.out_report)
    out_report.parent.mkdir(parents=True, exist_ok=True)
    out_report.write_text(report_text)
    print(f"Report written to {out_report}")

    input_path = Path(input_info["filepath"])
    first = results[0] if results else {}
    first_pct = first.get("percentile_info", {})
    out_result_json = Path(args.out_result_json)
    write_result_json(
        output_path=out_result_json,
        skill="gwas_prs",
        version="0.2.0",
        summary={
            "scores_calculated": len(results),
            "trait": first.get("trait", ""),
            "pgs_id": first.get("pgs_id", ""),
            "raw_score": first.get("prs", {}).get("raw_score"),
            "percentile": first_pct.get("percentile"),
            "risk_category": first_pct.get("risk_category"),
            "overlap_fraction": first.get("prs", {}).get("overlap_fraction"),
        },
        data={
            "input_info": input_info,
            "results": [
                {
                    "pgs_id": r["pgs_id"],
                    "trait": r["trait"],
                    "raw_score": r["prs"]["raw_score"],
                    "variants_used": r["prs"]["variants_used"],
                    "variants_total": r["prs"]["variants_total"],
                    "overlap_fraction": r["prs"]["overlap_fraction"],
                    "percentile": r["percentile_info"]["percentile"],
                    "risk_category": r["percentile_info"]["risk_category"],
                    "z_score": r["percentile_info"]["z_score"],
                    "method": r["percentile_info"]["method"],
                    "reference_population": r["percentile_info"].get("reference_population"),
                }
                for r in results
            ],
        },
        input_checksum=sha256_hex(input_path) if input_path.exists() else "",
    )
    print(f"Result envelope written to {out_result_json}")

    print()
    print("=" * 60)
    print("PRS SUMMARY")
    print("=" * 60)
    for r in results:
        pct = r["percentile_info"]
        pct_str = (
            f"{pct['percentile']:.1f}% ({pct['risk_category']})"
            if pct["percentile"] is not None else "N/A"
        )
        print(f"  {r['pgs_id']:12s} {r['trait']:<30s} "
              f"Score: {r['prs']['raw_score']:.4f}  Percentile: {pct_str}")
    print()
    print(DISCLAIMER)


if __name__ == "__main__":
    main()
