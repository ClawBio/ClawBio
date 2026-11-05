#!/usr/bin/env python3
"""Profile Report — unified personal genomic profile report generator.

Port of ClawBio 1.x `profile-report` for the ClawBio 2.0 `clawbio` CLI.
Reads a PatientProfile JSON and synthesizes all skill results into a single
human-readable markdown document. No re-computation — presentation only.

CLI surface differs slightly from the 1.x script for `clawbio` DAG fit:
  - `--demo` is dropped; the CLI drives inputs explicitly via
    `clawbio start --input profile=<path>`, so the demo profile is just an
    ordinary bundled fixture passed like any other `--profile` value (see
    `demo_full_profile.json` (skill root), copied verbatim from the 1.x repo).
  - No dependency on the `clawbio.common` shared library (that package does
    not exist in this repo) — the small set of helpers it provided
    (DISCLAIMER, write_result_json) are inlined below, following the
    precedent in skills/claw_ancestry_pca/scripts/ancestry_pca.py.

The core synthesis logic (profile loading, section renderers, cross-domain
insight detection, report assembly) is unchanged from the 1.x source.

Usage:
    python profile_report.py --profile <profile.json> --output <dir>
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

# ---------------------------------------------------------------------------
# Inlined helpers (ported from ClawBio 1.x clawbio/common/report.py)
# ---------------------------------------------------------------------------

DISCLAIMER = (
    "ClawBio is a research and educational tool. It is not a medical device "
    "and does not provide clinical diagnoses. Consult a healthcare "
    "professional before making any medical decisions."
)


def write_result_json(
    output_dir: str | Path,
    skill: str,
    version: str,
    summary: dict,
    data: dict,
    input_checksum: str = "",
) -> Path:
    """Write the standardized result.json envelope alongside report.md."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    envelope = {
        "skill": skill,
        "version": version,
        "completed_at": datetime.now(timezone.utc).isoformat(),
        "input_checksum": f"sha256:{input_checksum}" if input_checksum else "",
        "datasets": {},
        "summary": summary,
        "data": data,
    }

    result_path = output_dir / "result.json"
    result_path.write_text(json.dumps(envelope, indent=2, default=str))
    return result_path


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

VERSION = "0.1.0"
ALL_SKILLS = ["pharmgx", "nutrigx", "prs", "compare"]

# Cross-domain gene mapping: genes that appear in multiple skill contexts
CROSS_DOMAIN_GENES = {
    "CYP1A2": {
        "pharmgx": "Metabolizes clozapine, fluvoxamine, caffeine-containing drugs",
        "nutrigx": "Caffeine metabolism — determines slow vs fast metabolizer status",
    },
    "CYP2D6": {
        "pharmgx": "Major drug metabolizer — codeine, tramadol, TCAs, SSRIs",
        "nutrigx": "May influence metabolism of dietary amines",
    },
    "APOE": {
        "nutrigx": "Omega-3 and fat metabolism (APOE e2/e3/e4 status)",
        "prs": "Risk factor for Alzheimer's disease and cardiovascular traits",
    },
    "FTO": {
        "nutrigx": "Carbohydrate metabolism and satiety signalling",
        "prs": "BMI and obesity polygenic risk",
    },
    "TCF7L2": {
        "nutrigx": "Carbohydrate metabolism and insulin signalling",
        "prs": "Type 2 diabetes polygenic risk",
    },
}


# ---------------------------------------------------------------------------
# Profile loader
# ---------------------------------------------------------------------------


def load_profile(profile_path: str | Path) -> dict:
    """Load and validate a PatientProfile JSON file.

    Returns the raw profile dict with keys: metadata, genotypes, ancestry, skill_results.
    Raises ValueError if the file is missing or invalid.
    """
    path = Path(profile_path)
    if not path.exists():
        raise ValueError(f"Profile not found: {path}")

    data = json.loads(path.read_text())

    # Validate minimum structure
    if "metadata" not in data:
        raise ValueError("Profile missing 'metadata' section")
    if "genotypes" not in data:
        raise ValueError("Profile missing 'genotypes' section")

    # Ensure skill_results exists
    if "skill_results" not in data:
        data["skill_results"] = {}

    return data


def get_completed_skills(profile: dict) -> list[str]:
    """Return list of skill names that have results in the profile."""
    results = profile.get("skill_results", {})
    return [s for s in ALL_SKILLS if s in results and results[s]]


# ---------------------------------------------------------------------------
# Section renderers
# ---------------------------------------------------------------------------


def render_executive_summary(profile: dict) -> str:
    """Generate bullet-point executive summary of all available results."""
    lines = ["## Executive Summary", ""]
    results = profile.get("skill_results", {})

    # Pharmacogenomics summary
    pharmgx = _get_nested(results, "pharmgx", "data", "summary")
    if pharmgx:
        genes = pharmgx.get("genes_profiled", 0)
        avoid = pharmgx.get("drugs_avoid", 0)
        caution = pharmgx.get("drugs_caution", 0)
        standard = pharmgx.get("drugs_standard", 0)
        lines.append(
            f"- **Pharmacogenomics**: {genes} genes profiled "
            f"· {avoid} drugs to avoid · {caution} use with caution · {standard} standard use"
        )
    else:
        lines.append("- **Pharmacogenomics**: Not yet assessed")

    # PRS summary
    prs_data = _get_nested(results, "prs", "data", "data")
    if prs_data and "scores" in prs_data:
        scores = prs_data["scores"]
        trait_summaries = []
        for score in scores:
            trait = score.get("trait", "Unknown")
            percentile = score.get("percentile")
            category = score.get("risk_category", "")
            if percentile is not None:
                trait_summaries.append(f"{category} {trait} risk ({percentile:.0f}th percentile)")
        if trait_summaries:
            lines.append(f"- **Disease Risk**: {' · '.join(trait_summaries[:3])}")
            if len(trait_summaries) > 3:
                lines.append(f"  _(+{len(trait_summaries) - 3} more traits assessed)_")
        else:
            lines.append("- **Disease Risk**: Scores computed but no summary available")
    else:
        lines.append("- **Disease Risk**: Not yet assessed")

    # NutriGx summary
    nutrigx = _get_nested(results, "nutrigx", "data", "summary")
    if nutrigx:
        domains = nutrigx.get("domains_assessed", 0)
        elevated = nutrigx.get("elevated_domains", [])
        moderate = nutrigx.get("moderate_domains", [])
        highlights = []
        if elevated:
            highlights.append(f"elevated: {', '.join(elevated)}")
        if moderate:
            highlights.append(f"moderate: {', '.join(moderate)}")
        detail = " · ".join(highlights) if highlights else "no elevated domains"
        lines.append(f"- **Nutrition**: {domains} domains assessed · {detail}")
    else:
        lines.append("- **Nutrition**: Not yet assessed")

    # Ancestry summary — handle both data formats
    compare_data = _get_nested(results, "compare", "data", "data")
    if compare_data:
        ancestry = compare_data.get("ancestry_estimation", {})
        if not ancestry:
            ancestry_block = compare_data.get("ancestry", {})
            if isinstance(ancestry_block, dict):
                ancestry = ancestry_block.get("continental", {})
        if ancestry:
            numeric_pops = {k: v for k, v in ancestry.items() if isinstance(v, (int, float))}
            if numeric_pops:
                top_pop = max(numeric_pops.items(), key=lambda x: x[1])
                lines.append(
                    f"- **Ancestry**: Top component: {top_pop[0]} ({top_pop[1]:.1%})"
                )
            else:
                lines.append("- **Ancestry**: Estimation available")
        else:
            ibs = compare_data.get("ibs_summary", {})
            ibs_score = compare_data.get("ibs_score")
            if ibs or ibs_score is not None:
                lines.append("- **Ancestry**: IBS comparison completed")
            else:
                lines.append("- **Ancestry**: Results available")
    else:
        lines.append("- **Ancestry**: Not yet assessed")

    lines.append("")
    return "\n".join(lines)


def find_cross_domain_insights(profile: dict) -> str:
    """Identify genes/variants that appear across multiple skill results."""
    lines = ["### Cross-Domain Insights", ""]
    results = profile.get("skill_results", {})
    completed = set(get_completed_skills(profile))
    insights = []

    # Check PGx gene profiles
    pgx_genes = set()
    gene_profiles = _get_nested(results, "pharmgx", "data", "data", "gene_profiles")
    if gene_profiles:
        pgx_genes = set(gene_profiles.keys())

    # Check NutriGx SNP genes
    nutrigx_genes = set()
    snp_calls = _get_nested(results, "nutrigx", "data", "data", "snp_calls")
    if snp_calls:
        for rsid, call in snp_calls.items():
            gene = call.get("gene", "")
            if gene and call.get("status") == "found":
                nutrigx_genes.add(gene)

    # Find cross-domain genes actually present in results
    for gene, domains in CROSS_DOMAIN_GENES.items():
        present_in = []
        details = []
        if gene in pgx_genes and "pharmgx" in completed:
            present_in.append("pharmgx")
            gp = gene_profiles.get(gene, {})
            phenotype = gp.get("phenotype", "")
            details.append(f"PGx: {phenotype}" if phenotype else f"PGx: {domains.get('pharmgx', '')}")
        if gene in nutrigx_genes and "nutrigx" in completed:
            present_in.append("nutrigx")
            details.append(f"Nutrition: {domains.get('nutrigx', '')}")
        if len(present_in) >= 2:
            insights.append(f"- **{gene}** appears in {' and '.join(present_in)}: {' | '.join(details)}")

    if insights:
        lines.extend(insights)
    else:
        lines.append("_No cross-domain gene overlaps detected in current results._")

    lines.append("")
    return "\n".join(lines)


def render_pharmgx_section(profile: dict) -> str:
    """Render the pharmacogenomics section from profile data."""
    results = profile.get("skill_results", {})
    pgx = _get_nested(results, "pharmgx", "data", "data")

    if not pgx:
        return _missing_section(
            "Pharmacogenomics",
            "Run `clawbio.py run pharmgx --profile <your_profile.json>` to generate pharmacogenomic results.",
        )

    lines = ["## Pharmacogenomics", ""]

    # Gene profile table
    gene_profiles = pgx.get("gene_profiles", {})
    if gene_profiles:
        lines.append("### Gene Profiles")
        lines.append("")
        lines.append("| Gene | Diplotype | Phenotype |")
        lines.append("|------|-----------|-----------|")
        for gene, info in sorted(gene_profiles.items()):
            diplotype = info.get("diplotype", "—")
            phenotype = info.get("phenotype", "—")
            lines.append(f"| {gene} | {diplotype} | {phenotype} |")
        lines.append("")

    # Drug recommendations by classification
    drug_recs = pgx.get("drug_recommendations", {})
    for classification in ["avoid", "caution", "standard"]:
        drugs = drug_recs.get(classification, [])
        if not drugs:
            continue

        header_map = {
            "avoid": "Drugs to Avoid",
            "caution": "Drugs Requiring Caution",
            "standard": "Standard Use Drugs",
        }
        emoji_map = {
            "avoid": "🔴",
            "caution": "🟡",
            "standard": "🟢",
        }

        lines.append(f"### {emoji_map[classification]} {header_map[classification]} ({len(drugs)})")
        lines.append("")
        lines.append("| Drug | Brand | Class | Gene | Recommendation |")
        lines.append("|------|-------|-------|------|----------------|")
        for drug in drugs:
            lines.append(
                f"| {drug.get('drug', '—')} "
                f"| {drug.get('brand', '—')} "
                f"| {drug.get('class', '—')} "
                f"| {drug.get('gene', '—')} "
                f"| {drug.get('recommendation', '—')} |"
            )
        lines.append("")

    return "\n".join(lines)


def render_prs_section(profile: dict) -> str:
    """Render the polygenic risk scores section from profile data."""
    results = profile.get("skill_results", {})
    prs = _get_nested(results, "prs", "data", "data")

    if not prs:
        return _missing_section(
            "Polygenic Risk Scores",
            "Run `clawbio.py run prs --profile <your_profile.json>` to calculate polygenic risk scores.",
        )

    lines = ["## Polygenic Risk Scores", ""]

    scores = prs.get("scores", [])
    if scores:
        lines.append("| PGS ID | Trait | Raw Score | Percentile | Risk Category |")
        lines.append("|--------|-------|-----------|------------|---------------|")
        for score in scores:
            pgs_id = score.get("pgs_id", "—")
            trait = score.get("trait", "—")
            raw = score.get("raw_score")
            raw_str = f"{raw:.4f}" if raw is not None else "—"
            percentile = score.get("percentile")
            pct_str = f"{percentile:.0f}th" if percentile is not None else "—"
            category = score.get("risk_category", "—")
            lines.append(f"| {pgs_id} | {trait} | {raw_str} | {pct_str} | {category} |")
        lines.append("")

        # Per-trait detail
        for score in scores:
            percentile = score.get("percentile")
            if percentile is not None and percentile >= 90:
                trait = score.get("trait", "Unknown")
                pgs_id = score.get("pgs_id", "")
                lines.append(
                    f"> **Elevated risk**: {trait} at {percentile:.0f}th percentile ({pgs_id}). "
                    f"This is a statistical association, not a diagnosis."
                )
                lines.append("")
    else:
        lines.append("_PRS results present but no individual scores found._")
        lines.append("")

    return "\n".join(lines)


def render_nutrigx_section(profile: dict) -> str:
    """Render the nutrigenomics section from profile data."""
    results = profile.get("skill_results", {})
    nutrigx = _get_nested(results, "nutrigx", "data", "data")

    if not nutrigx:
        return _missing_section(
            "Nutrigenomics",
            "Run `clawbio.py run nutrigx --profile <your_profile.json>` to generate nutrigenomic results.",
        )

    lines = ["## Nutrigenomics", ""]

    risk_scores = nutrigx.get("risk_scores", {})
    if risk_scores:
        lines.append("### Domain Risk Assessment")
        lines.append("")
        lines.append("| Domain | Score | Category | Coverage |")
        lines.append("|--------|-------|----------|----------|")
        for domain, info in sorted(risk_scores.items()):
            score = info.get("score")
            score_str = f"{score:.1f}" if score is not None else "—"
            category = info.get("category", "Unknown")
            coverage = info.get("coverage", "—")
            display_domain = domain.replace("_", " ").title()
            lines.append(f"| {display_domain} | {score_str} | {category} | {coverage} |")
        lines.append("")

        # Highlight domains with data
        unknown = [d for d, info in risk_scores.items() if info.get("category") == "Unknown"]
        if unknown:
            lines.append(
                f"_Note: {len(unknown)} of {len(risk_scores)} domains lack sufficient SNP coverage. "
                f"A broader genotyping panel would improve coverage._"
            )
            lines.append("")

    return "\n".join(lines)


def render_ancestry_section(profile: dict) -> str:
    """Render the ancestry & relatedness section from profile data."""
    results = profile.get("skill_results", {})
    compare = _get_nested(results, "compare", "data", "data")

    if not compare:
        return _missing_section(
            "Ancestry & Relatedness",
            "Run `clawbio.py run compare --profile <your_profile.json>` to generate ancestry results.",
        )

    lines = ["## Ancestry & Relatedness", ""]

    # IBS summary — handle both formats:
    # Format A (genome-compare v2): ibs_summary.{total_snps_compared, ibs2_count, ...}
    # Format B (genome-compare v1): {ibs_score, n_overlap, n_concordant}
    ibs = compare.get("ibs_summary", {})
    if ibs:
        lines.append("### Identity-by-State (IBS) Comparison")
        lines.append("")
        total = ibs.get("total_snps_compared", 0)
        shared = ibs.get("ibs2_count", 0)
        ibs2_pct = ibs.get("ibs2_proportion", 0)
        lines.append(f"- **SNPs compared**: {total:,}")
        lines.append(f"- **IBS2 (identical genotypes)**: {shared:,} ({ibs2_pct:.1%})")
        ibs1 = ibs.get("ibs1_count", 0)
        ibs1_pct = ibs.get("ibs1_proportion", 0)
        if ibs1:
            lines.append(f"- **IBS1 (one allele shared)**: {ibs1:,} ({ibs1_pct:.1%})")
        ibs0 = ibs.get("ibs0_count", 0)
        ibs0_pct = ibs.get("ibs0_proportion", 0)
        if ibs0:
            lines.append(f"- **IBS0 (no alleles shared)**: {ibs0:,} ({ibs0_pct:.1%})")
        lines.append("")
    elif "ibs_score" in compare:
        # Format B: simple IBS score
        lines.append("### Identity-by-State (IBS) Comparison")
        lines.append("")
        ibs_score = compare.get("ibs_score", 0)
        n_overlap = compare.get("n_overlap", 0)
        n_concordant = compare.get("n_concordant", 0)
        lines.append(f"- **IBS score**: {ibs_score:.4f}")
        lines.append(f"- **Overlapping SNPs**: {n_overlap:,}")
        lines.append(f"- **Concordant genotypes**: {n_concordant:,}")
        lines.append("")

    # Ancestry estimation — handle both formats:
    # Format A: ancestry_estimation.{population: proportion}
    # Format B: ancestry.continental.{AFR, EUR, EAS, SAS, AMR}
    ancestry = compare.get("ancestry_estimation", {})
    if not ancestry:
        # Try format B
        ancestry_block = compare.get("ancestry", {})
        if isinstance(ancestry_block, dict):
            ancestry = ancestry_block.get("continental", {})
    if ancestry:
        lines.append("### Ancestry Estimation")
        lines.append("")
        lines.append("| Population | Proportion |")
        lines.append("|------------|------------|")
        for pop, proportion in sorted(ancestry.items(), key=lambda x: -x[1]):
            if isinstance(proportion, (int, float)):
                lines.append(f"| {pop} | {proportion:.1%} |")
        lines.append("")

    # Reference comparison
    ref = compare.get("reference_comparison", {})
    if ref:
        ref_name = ref.get("reference_name", "Reference")
        lines.append(f"### Comparison to {ref_name}")
        lines.append("")
        for key, val in ref.items():
            if key != "reference_name":
                lines.append(f"- **{key.replace('_', ' ').title()}**: {val}")
        lines.append("")

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Report assembler
# ---------------------------------------------------------------------------


def generate_profile_report(profile: dict) -> str:
    """Combine all sections into a unified profile report."""
    meta = profile.get("metadata", {})
    patient_id = meta.get("patient_id", "Unknown")
    genotype_count = len(profile.get("genotypes", {}))
    completed = get_completed_skills(profile)
    now = datetime.now(timezone.utc).strftime("%Y-%m-%d")

    # Header
    lines = [
        "# Your Genomic Profile",
        "",
        f"**Patient ID**: {patient_id}  |  **Date**: {now}  |  **SNPs**: {genotype_count}",
        f"**Skills completed**: {', '.join(completed) if completed else 'none'} ({len(completed)} of {len(ALL_SKILLS)})",
        "",
        "---",
        "",
    ]

    # Executive summary
    lines.append(render_executive_summary(profile))

    # Cross-domain insights (only if 2+ skills completed)
    if len(completed) >= 2:
        lines.append(find_cross_domain_insights(profile))

    lines.append("---")
    lines.append("")

    # Skill sections
    lines.append(render_pharmgx_section(profile))
    lines.append(render_prs_section(profile))
    lines.append(render_nutrigx_section(profile))
    lines.append(render_ancestry_section(profile))

    # Methods and disclaimer
    lines.extend([
        "---",
        "",
        "## Methods",
        "",
        "This report was generated by the ClawBio Profile Report skill (v{version}). "
        "Results are synthesized from previously computed skill outputs stored in the "
        "PatientProfile JSON. No re-computation was performed.".format(version=VERSION),
        "",
        "**Skills used**:",
    ])
    for skill in ALL_SKILLS:
        status = "completed" if skill in completed else "not run"
        lines.append(f"- `{skill}`: {status}")
    lines.extend([
        "",
        "---",
        "",
        "## Disclaimer",
        "",
        f"*{DISCLAIMER}*",
        "",
    ])

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _get_nested(d: dict, *keys: str) -> Any:
    """Safely traverse nested dicts, returning None if any key is missing."""
    for key in keys:
        if not isinstance(d, dict):
            return None
        d = d.get(key)
        if d is None:
            return None
    return d


def _missing_section(title: str, instruction: str) -> str:
    """Generate a placeholder section for a skill that hasn't been run."""
    return f"""## {title}

_Not yet assessed._ {instruction}

"""


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Profile Report — unified personal genomic profile report",
    )
    parser.add_argument(
        "--profile", "-p",
        required=True,
        help="Path to PatientProfile JSON file",
    )
    parser.add_argument(
        "--output", "-o",
        required=True,
        help="Output directory for report files",
    )
    args = parser.parse_args()

    print(f"Loading profile: {args.profile}")
    try:
        profile = load_profile(args.profile)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

    # Generate report
    completed = get_completed_skills(profile)
    patient_id = profile.get("metadata", {}).get("patient_id", "Unknown")
    print(f"Patient: {patient_id} | Skills completed: {completed}")

    report_md = generate_profile_report(profile)

    # Write output
    out_dir = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)

    report_path = out_dir / "profile_report.md"
    report_path.write_text(report_md)
    print(f"Report written: {report_path}")

    # Write result.json
    write_result_json(
        output_dir=out_dir,
        skill="profile_report",
        version=VERSION,
        summary={
            "patient_id": patient_id,
            "skills_completed": completed,
            "skills_total": len(ALL_SKILLS),
            "genotype_count": len(profile.get("genotypes", {})),
        },
        data={
            "report_file": "profile_report.md",
            "skills_completed": completed,
            "skills_missing": [s for s in ALL_SKILLS if s not in completed],
        },
    )
    print(f"Result JSON written: {out_dir / 'result.json'}")
    print("Done.")


if __name__ == "__main__":
    main()
