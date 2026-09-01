#!/usr/bin/env python3
"""
Ancestry-Aware Disease Risk Profiler.

Infers genetic super-population ancestry from a 23andMe/AncestryDNA file,
then computes ancestry-stratified odds ratios and an exploratory Ancestry
Elevation Score (AES) comparing ancestry-specific GWAS effect sizes to
European reference estimates.

NOTE: This skill does NOT compute absolute lifetime risk percentages.
Applying individual-variant ORs on top of a population baseline prevalence
double-counts allele contributions already embedded in that baseline.
For validated absolute risk estimates, use ClawBio's gwas-prs skill with
an ancestry-appropriate PGS Catalog score.

Usage:
    python ancestry_risk_profiler.py --input <23andme_file> --output <dir>
    python ancestry_risk_profiler.py --input <23andme_file> --ancestry SAS --output <dir>
    python ancestry_risk_profiler.py --demo --output /tmp/ancestry_risk_demo
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

_PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from clawbio.common.parsers import parse_genetic_file, genotypes_to_simple
from clawbio.common.report import DISCLAIMER

SKILL_DIR = Path(__file__).resolve().parent
DATA_DIR = SKILL_DIR / "data"

SUPERPOPULATIONS = ("AFR", "AMR", "EAS", "EUR", "SAS")
SUPERPOP_LABELS = {
    "AFR": "African",
    "AMR": "Admixed American / Latino",
    "EAS": "East Asian",
    "EUR": "European",
    "SAS": "South Asian",
}

# Minimum AISNP panel hits required before ancestry inference is attempted.
# 30 markers is the lower bound validated in published minimal AISNP panels
# (Kosoy et al. 2009, Hum Genet) for reliable continental-level super-population
# assignment. Below this the Hardy-Weinberg log-likelihood gap between populations
# is too small to be trustworthy.
MIN_AISNP_COVERAGE = 30

# AES display thresholds — for colouring only; not validated clinical cutoffs.
_AES_DISPLAY_ELEVATED = 1.3
_AES_DISPLAY_REDUCED = 0.77

# Minimum top-match posterior probability required before disease risk scoring is performed.
# Below this threshold the ancestry assignment is too ambiguous (likely admixed or mixed-
# signal) to reliably gate disease ORs on a single population label.
LOW_POSTERIOR_ABSTAIN = 0.35


# ---------------------------------------------------------------------------
# Custom exceptions
# ---------------------------------------------------------------------------


class InsufficientCoverageError(ValueError):
    """Raised when too few AISNPs match to infer ancestry reliably."""


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------


@dataclass
class VariantHit:
    rsid: str
    genotype: str
    dosage: int
    or_ancestry: float
    or_eur: float
    delta_log_or: float
    note: str


@dataclass
class DiseaseRisk:
    disease: str
    combined_or: float       # product of ancestry-specific ORs across risk alleles carried
    or_eur_combined: float   # same calculation using EUR reference ORs — for direct comparison
    aes: float               # exploratory directional indicator; not a validated clinical score
    n_variants: int = 0      # number of variants contributing to combined_or
    variants_found: list = field(default_factory=list)
    aes_category: str = "neutral"
    variants_detail: list = field(default_factory=list)


# ---------------------------------------------------------------------------
# Ancestry inference
# ---------------------------------------------------------------------------


def load_aisnp_panel(path: Path) -> dict[str, dict]:
    """Load AISNP panel CSV → {rsid: {alt, AFR, AMR, EAS, EUR, SAS}}."""
    panel: dict[str, dict] = {}
    with open(path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            rsid = row["rsid"].strip()
            panel[rsid] = {
                "ref": row["ref"].strip(),
                "alt": row["alt"].strip(),
                "AFR": float(row["AFR"]),
                "AMR": float(row["AMR"]),
                "EAS": float(row["EAS"]),
                "EUR": float(row["EUR"]),
                "SAS": float(row["SAS"]),
            }
    return panel


def _genotype_dosage(genotype: str, allele: str) -> int:
    """Count copies of `allele` in a two-character genotype."""
    return genotype.count(allele)


def infer_ancestry(
    genotypes: dict[str, str],
    panel: dict[str, dict],
    ancestry_override: Optional[str] = None,
) -> dict:
    """Compute genetic super-population likelihood from AISNP panel.

    Raises InsufficientCoverageError if fewer than MIN_AISNP_COVERAGE panel
    markers are present and no ancestry_override is provided — a sparse genotype
    file cannot support reliable super-population assignment.

    Returns dict with inferred_ancestry, confidence, scores, overridden flag.
    """
    if ancestry_override:
        ancestry_override = ancestry_override.upper()
        if ancestry_override not in SUPERPOPULATIONS:
            raise ValueError(f"Unknown ancestry: {ancestry_override}. Choose from {SUPERPOPULATIONS}")
        return {
            "inferred_ancestry": ancestry_override,
            "confidence": "user-supplied",
            "scores": {p: 0.0 for p in SUPERPOPULATIONS},
            "posterior": {p: (1.0 if p == ancestry_override else 0.0) for p in SUPERPOPULATIONS},
            "aisnp_coverage": 0,
            "overridden": True,
        }

    log_likelihoods: dict[str, float] = {p: 0.0 for p in SUPERPOPULATIONS}
    hits = 0

    for rsid, info in panel.items():
        if rsid not in genotypes:
            continue
        genotype = genotypes[rsid]
        alt = info["alt"]
        dosage = _genotype_dosage(genotype, alt)
        hits += 1

        for pop in SUPERPOPULATIONS:
            p = info[pop]
            p = max(min(p, 0.9999), 0.0001)
            q = 1.0 - p
            if dosage == 0:
                prob = q * q
            elif dosage == 1:
                prob = 2 * p * q
            else:
                prob = p * p
            log_likelihoods[pop] += math.log(max(prob, 1e-15))

    if hits < MIN_AISNP_COVERAGE:
        raise InsufficientCoverageError(
            f"Only {hits} AISNP panel marker(s) matched (minimum required: {MIN_AISNP_COVERAGE}). "
            f"Ancestry cannot be reliably inferred from this file. "
            f"Re-run with --ancestry AFR|AMR|EAS|EUR|SAS to specify your documented genetic ancestry."
        )

    best = max(log_likelihoods, key=lambda p: log_likelihoods[p])
    sorted_ll = sorted(log_likelihoods.values(), reverse=True)
    gap = sorted_ll[0] - sorted_ll[1]

    if gap > 50:
        confidence = "high"
    elif gap > 15:
        confidence = "medium"
    else:
        confidence = "low"

    # Numerically-stable softmax over log-likelihoods → posterior probability per population.
    # A hard best-match label is kept for downstream association lookup, but the posterior
    # is reported alongside it so low-confidence or admixed cases are not over-called.
    max_ll = max(log_likelihoods.values())
    exp_lls = {p: math.exp(log_likelihoods[p] - max_ll) for p in SUPERPOPULATIONS}
    total_exp = sum(exp_lls.values())
    posterior = {p: round(exp_lls[p] / total_exp, 4) for p in SUPERPOPULATIONS}

    return {
        "inferred_ancestry": best,
        "confidence": confidence,
        "scores": log_likelihoods,
        "posterior": posterior,
        "aisnp_coverage": hits,
        "overridden": False,
    }


# ---------------------------------------------------------------------------
# Risk scoring
# ---------------------------------------------------------------------------


def load_associations(path: Path) -> dict:
    """Load ancestry_risk_associations.json."""
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def _aes_category(aes: float) -> str:
    if aes >= _AES_DISPLAY_ELEVATED:
        return "elevated by ancestry"
    if aes <= _AES_DISPLAY_REDUCED:
        return "reduced by ancestry"
    return "neutral"


def _handle_compound_recessive_groups(
    genotypes: dict[str, str],
    ancestry: str,
    associations_data: dict,
) -> list[DiseaseRisk]:
    """Score recessive compound loci (e.g. APOL1) that require two risk alleles.

    Per-allele log-additive OR is biologically wrong for these loci: APOL1
    nephropathy risk is recessive and requires two high-risk alleles (G1+G2,
    G1/G1, or G2/G2). This function counts total risk alleles across the group
    and applies the validated compound OR rather than multiplying per-allele ORs.
    """
    groups = associations_data.get("compound_recessive_groups", {})
    results: list[DiseaseRisk] = []

    for group_name, group in groups.items():
        if group.get("ancestry") != ancestry:
            continue

        disease = group["disease"]
        group_assocs = [
            a for a in associations_data["associations"]
            if a.get("model") == "recessive_compound"
            and a.get("compound_group") == group_name
            and a["ancestry"] == ancestry
        ]

        risk_alleles_by_rsid: dict[str, str] = {
            a["rsid"]: a["risk_allele"] for a in group_assocs
        }

        total_risk_alleles = 0
        carrier_details: list[dict] = []
        for rsid, risk_allele in risk_alleles_by_rsid.items():
            if rsid not in genotypes:
                continue
            dosage = _genotype_dosage(genotypes[rsid], risk_allele)
            if dosage > 0:
                total_risk_alleles += dosage
                carrier_details.append({
                    "rsid": rsid,
                    "genotype": genotypes[rsid],
                    "dosage": dosage,
                })

        if total_risk_alleles == 0:
            continue

        eur_or = group.get("eur_two_allele_or", 1.0)
        if total_risk_alleles == 1:
            compound_or = group.get("one_allele_or", 1.3)
            note = group.get("one_allele_note", "Single risk allele: carrier status")
        else:
            compound_or = group["two_allele_or"]
            note = group.get("note", f"Compound recessive: {total_risk_alleles} risk alleles")

        aes = compound_or / eur_or if eur_or > 0 else compound_or

        variant_details = [
            {
                "rsid": d["rsid"],
                "genotype": d["genotype"],
                "dosage": d["dosage"],
                "or_ancestry": round(compound_or, 3),
                "or_eur": round(eur_or, 3),
                "delta_log_or": round(
                    math.log(compound_or / eur_or) if eur_or > 0 else math.log(compound_or), 4
                ),
                "note": note,
            }
            for d in carrier_details
        ]

        results.append(DiseaseRisk(
            disease=disease,
            combined_or=round(compound_or, 3),
            or_eur_combined=round(eur_or, 3),
            aes=round(aes, 3),
            n_variants=total_risk_alleles,
            variants_found=[d["rsid"] for d in carrier_details],
            aes_category=_aes_category(aes),
            variants_detail=variant_details,
        ))

    return results


def compute_disease_risks(
    genotypes: dict[str, str],
    ancestry: str,
    associations_data: dict,
) -> list[DiseaseRisk]:
    """Compute DiseaseRisk objects for all diseases in the association panel.

    For each disease:
    - combined_or: log-additive OR across risk alleles carried (ancestry-specific)
    - or_eur_combined: same calculation using EUR reference ORs
    - aes: exp(sum of delta_log_or vs EUR) — exploratory ancestry elevation indicator
    - n_variants: count of variants contributing to combined_or

    Entries marked model="recessive_compound" (e.g. APOL1) are handled by
    _handle_compound_recessive_groups instead of the additive loop here.
    """
    ancestry = ancestry.upper()
    all_assocs = associations_data["associations"]

    by_disease: dict[str, list[dict]] = {}
    for assoc in all_assocs:
        if assoc["ancestry"] != ancestry:
            continue
        if assoc.get("model") == "recessive_compound":
            continue  # handled separately by _handle_compound_recessive_groups
        d = assoc["disease"]
        by_disease.setdefault(d, []).append(assoc)

    results: list[DiseaseRisk] = []

    for disease, assocs in by_disease.items():
        sum_log_or = 0.0
        sum_log_or_eur = 0.0
        sum_delta_log_or = 0.0
        variants_found: list[str] = []
        variant_details: list[dict] = []

        for assoc in assocs:
            rsid = assoc["rsid"]
            risk_allele = assoc["risk_allele"]
            or_ancestry = assoc["or"]
            or_eur = assoc.get("eur_or", or_ancestry)

            if rsid not in genotypes:
                continue
            genotype = genotypes[rsid]
            dosage = _genotype_dosage(genotype, risk_allele)
            if dosage == 0:
                continue

            log_or = math.log(or_ancestry) * dosage
            log_or_eur = math.log(or_eur) * dosage
            delta_log_or = (math.log(or_ancestry) - math.log(or_eur)) * dosage

            sum_log_or += log_or
            sum_log_or_eur += log_or_eur
            sum_delta_log_or += delta_log_or
            variants_found.append(rsid)
            variant_details.append({
                "rsid": rsid,
                "genotype": genotype,
                "dosage": dosage,
                "or_ancestry": round(or_ancestry, 3),
                "or_eur": round(or_eur, 3),
                "delta_log_or": round(delta_log_or, 4),
                "note": assoc.get("note", ""),
            })

        if not variants_found:
            continue

        combined_or = math.exp(sum_log_or)
        or_eur_combined = math.exp(sum_log_or_eur)
        aes = math.exp(sum_delta_log_or)

        results.append(DiseaseRisk(
            disease=disease,
            combined_or=round(combined_or, 3),
            or_eur_combined=round(or_eur_combined, 3),
            aes=round(aes, 3),
            n_variants=len(variants_found),
            variants_found=variants_found,
            aes_category=_aes_category(aes),
            variants_detail=variant_details,
        ))

    results.extend(_handle_compound_recessive_groups(genotypes, ancestry, associations_data))
    results.sort(key=lambda r: r.aes, reverse=True)
    return results


def _get_risks_with_confidence_gating(
    genotypes: dict[str, str],
    ancestry_result: dict,
    associations: dict,
) -> tuple[list[DiseaseRisk], Optional[str]]:
    """Apply soft posterior gating before disease risk scoring.

    Three outcomes:
    - confidence != "low" or user override: score normally, no note.
    - confidence == "low", top posterior >= LOW_POSTERIOR_ABSTAIN: score with a caveat note.
    - confidence == "low", top posterior < LOW_POSTERIOR_ABSTAIN: abstain (return empty list
      and an explanatory note). The ancestry assignment is too ambiguous to hard-gate ORs on.
    """
    ancestry = ancestry_result["inferred_ancestry"]
    confidence = ancestry_result["confidence"]
    overridden = ancestry_result.get("overridden", False)

    if overridden or confidence != "low":
        return compute_disease_risks(genotypes, ancestry, associations), None

    posterior = ancestry_result.get("posterior", {})
    top_posterior = max(posterior.values()) if posterior else 1.0

    if top_posterior < LOW_POSTERIOR_ABSTAIN:
        note = (
            f"Disease risk scoring was not performed. "
            f"Ancestry assignment confidence is too low "
            f"(top posterior for {ancestry}: {top_posterior:.1%}). "
            f"This typically indicates admixed or ambiguous ancestry signal. "
            f"Re-run with `--ancestry {{AFR|AMR|EAS|EUR|SAS}}` to specify your "
            f"documented genetic ancestry and enable disease risk scoring."
        )
        return [], note

    note = (
        f"Ancestry assignment has low confidence "
        f"(top posterior for {ancestry}: {top_posterior:.1%}). "
        f"Disease risk results below use this tentative assignment. "
        f"Consider re-running with `--ancestry {ancestry}` (or another population) "
        f"if you know your documented genetic ancestry."
    )
    return compute_disease_risks(genotypes, ancestry, associations), note


# ---------------------------------------------------------------------------
# Report generation
# ---------------------------------------------------------------------------

_AES_EMOJI = {
    "elevated by ancestry": "🔴",
    "neutral": "🟡",
    "reduced by ancestry": "🟢",
}


def _bar(value: float, max_val: float, width: int = 20) -> str:
    filled = int(round(value / max_val * width)) if max_val > 0 else 0
    return "█" * filled + "░" * (width - filled)


def generate_report(
    risks: list[DiseaseRisk],
    ancestry_result: dict,
    output_dir: Path,
    scoring_note: Optional[str] = None,
) -> None:
    """Write ancestry_risk_report.md, ancestry_risk_result.json, and figures/."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    figs_dir = output_dir / "figures"
    figs_dir.mkdir(exist_ok=True)

    ancestry = ancestry_result["inferred_ancestry"]
    ancestry_label = SUPERPOP_LABELS.get(ancestry, ancestry)
    confidence = ancestry_result["confidence"]
    coverage = ancestry_result.get("aisnp_coverage", "n/a")
    overridden = ancestry_result.get("overridden", False)
    timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")

    lines: list[str] = []

    lines += [
        "# Ancestry-Aware Disease Risk Profile",
        "",
        f"**Generated**: {timestamp}  ",
        f"**Skill**: ancestry-risk-profiler v1.3.0  ",
        "",
    ]

    # Section 1: Inferred genetic super-population ancestry
    override_note = " *(user-supplied)*" if overridden else f" *(confidence: {confidence}, AISNPs matched: {coverage})*"
    lines += [
        "---",
        "## 1. Inferred Genetic Super-Population Ancestry",
        "",
        "> **Note**: This is an inference of genetic super-population based on allele frequencies at",
        "> ~80 ancestry-informative SNPs. It is **not** self-reported ethnicity, cultural identity,",
        "> or nationality. Super-population labels (AFR, EAS, EUR, SAS, AMR) are analytical",
        "> categories from the 1000 Genomes Project — not ethnic identifiers.",
        "",
        f"| Field | Value |",
        f"|---|---|",
        f"| Genetic super-population | **{ancestry}** — {ancestry_label}{override_note} |",
        f"| Confidence | {confidence} |",
        f"| AISNPs matched | {coverage} |",
        "",
    ]

    if confidence == "low" and not overridden:
        posterior = ancestry_result.get("posterior", {})
        sorted_pops = sorted(posterior.items(), key=lambda x: x[1], reverse=True)
        lines += [
            "> ⚠️ **Low confidence**: The log-likelihood gap between the top two super-populations",
            "> is small, which can indicate admixed ancestry or insufficient marker overlap.",
            "> Consider re-running with `--ancestry` to specify your documented genetic ancestry.",
            ">",
            "> **Estimated posterior probability across super-populations** (risk scoring uses the top match):",
            ">",
            "> | Population | Posterior Probability |",
            "> |---|---|",
        ]
        for pop, prob in sorted_pops:
            label = SUPERPOP_LABELS.get(pop, pop)
            lines.append(f"> | {pop} — {label} | {prob * 100:.1f}% |")
        lines.append("")

    # Section 2: Disease risk table (OR-based, no absolute lifetime risk %)
    lines += [
        "---",
        "## 2. Ancestry-Stratified Disease Risk Summary",
        "",
        "> **What these numbers mean:**",
        "> - **Ancestry OR**: combined odds ratio using published GWAS effect sizes for your inferred",
        ">   super-population, across the risk variants you carry (log-additive model).",
        "> - **EUR ref OR**: the same calculation using European reference effect sizes for those",
        ">   same variants — allows direct comparison.",
        "> - **AES** (Ancestry Elevation Score): exp(Σ[log OR_ancestry − log OR_EUR]) — an",
        ">   **exploratory directional indicator** showing where ancestry-specific effect sizes",
        ">   diverge from European estimates. AES > 1 = ancestry amplifies the signal; AES < 1 =",
        ">   ancestry attenuates it. **AES has not been externally validated and is not a clinical",
        ">   risk score.** Use it as a signal to explore further, not a probability.",
        ">",
        "> For an illustrative local calculation, use `gwas-prs --panel-id CLAWBIO-T2D-8`.",
        "> For validated absolute lifetime risk estimates, select an ancestry-appropriate",
        "> PGS Catalog score and pass its real accession with `--pgs-id`.",
        "",
    ]

    if scoring_note:
        lines += [
            f"> ⚠️ **Ancestry scoring note**: {scoring_note}",
            "",
        ]

    if not risks:
        lines += [
            "> No disease risk results available. See the ancestry scoring note above for details.",
            "",
        ]
    else:
        lines += [
            "| Disease | AES (exploratory) | Direction | N variants | Per-allele ORs (see detail) |",
            "|---|---|---|---|---|",
        ]
        for r in risks:
            emoji = _AES_EMOJI.get(r.aes_category, "🟡")
            or_summary = ", ".join(
                f"{v['rsid']} {v['or_ancestry']:.2f}x"
                for v in r.variants_detail
            )
            lines.append(
                f"| {r.disease} "
                f"| {r.aes:.2f} "
                f"| {emoji} {r.aes_category.title()} "
                f"| {r.n_variants} "
                f"| {or_summary} |"
            )
        lines += [
            "",
            "> **Note on per-allele ORs**: values above are individual published GWAS effect sizes,",
            "> not a combined disease risk estimate. Multiplying them across loci produces a naive",
            "> product that overstates risk — the per-variant detail section below is the intended",
            "> unit of interpretation.",
            "",
        ]

    # Section 3: AES visual
    lines += [
        "---",
        "## 3. Ancestry Elevation Score — Top Diseases (Exploratory)",
        "",
        "> AES is a directional indicator, not a validated clinical metric.",
        "",
        "```",
        f"{'Disease':<35} {'AES':>5}   {'Bar':<22} Direction",
        "-" * 75,
    ]
    max_aes = max((r.aes for r in risks), default=2.0)
    for r in risks[:10]:
        bar = _bar(r.aes, max(max_aes, 2.0))
        lines.append(
            f"{r.disease:<35} {r.aes:>5.2f}   {bar:<22} {r.aes_category}"
        )
    lines += ["```", ""]

    # Section 4: Variant detail
    lines += [
        "---",
        "## 4. Variant Detail",
        "",
        "| rsID | Genotype | Disease | OR (your ancestry) | OR (EUR ref) | Delta log-OR | Note |",
        "|---|---|---|---|---|---|---|",
    ]
    for r in risks:
        for v in r.variants_detail:
            lines.append(
                f"| {v['rsid']} | {v['genotype']} | {r.disease} "
                f"| {v['or_ancestry']:.2f} | {v['or_eur']:.2f} "
                f"| {v['delta_log_or']:+.3f} | {v['note'][:80]} |"
            )
    lines.append("")

    # Section 5: Methodology and disclaimer
    lines += [
        "---",
        "## Methodology",
        "",
        f"**Ancestry inference**: Hardy-Weinberg log-likelihood at ~80 ancestry-informative SNPs",
        f"(AISNPs) across 5 super-populations (AFR, AMR, EAS, EUR, SAS). Minimum {MIN_AISNP_COVERAGE} matched",
        "markers required before inference is attempted.",
        "",
        "**Ancestry Elevation Score (AES)**: exp(Σ[log OR_ancestry_i − log OR_EUR_i] × dosage_i)",
        "summed across risk variants for each disease. This is an **exploratory metric** with no",
        "published external validation. It shows divergence from European GWAS estimates — it does",
        "not represent a calibrated disease probability.",
        "",
        "**Why no absolute lifetime risk %?** Applying individual-variant ORs on top of a",
        "population baseline prevalence double-counts the allele contribution already embedded in",
        "that baseline. For validated absolute risk, use `gwas-prs` with a PGS Catalog score",
        "calibrated to your ancestry.",
        "",
        "**SNP independence assumed**: variants are treated as independent (log-additive model);",
        "linkage disequilibrium is not modelled.",
        "",
        "---",
        "## Disclaimer",
        "",
        f"*{DISCLAIMER}*",
        "",
    ]

    report_text = "\n".join(lines)
    report_path = output_dir / "ancestry_risk_report.md"
    report_path.write_text(report_text, encoding="utf-8")

    # JSON result
    result = {
        "inferred_ancestry": ancestry,
        "ancestry_label": ancestry_label,
        "confidence": confidence,
        "posterior": ancestry_result.get("posterior", {}),
        "aisnp_coverage": coverage,
        "overridden": overridden,
        "scoring_note": scoring_note,
        "risks": [
            {
                "disease": r.disease,
                "combined_or": r.combined_or,
                "or_eur_combined": r.or_eur_combined,
                "aes": r.aes,
                "n_variants": r.n_variants,
                "aes_category": r.aes_category,
                "variants_found": r.variants_found,
                "variants_detail": r.variants_detail,
            }
            for r in risks
        ],
        "timestamp": timestamp,
        "skill_version": "1.3.0",
    }
    json_path = output_dir / "ancestry_risk_result.json"
    json_path.write_text(json.dumps(result, indent=2), encoding="utf-8")

    _try_write_aes_chart(risks, figs_dir)


def _try_write_aes_chart(risks: list[DiseaseRisk], figs_dir: Path) -> None:
    """Write a horizontal AES bar chart; silently skip if matplotlib unavailable."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return

    top = risks[:10]
    if not top:
        return

    names = [r.disease for r in reversed(top)]
    aes_vals = [r.aes for r in reversed(top)]
    colors = []
    for r in reversed(top):
        if r.aes_category == "elevated by ancestry":
            colors.append("#d62728")
        elif r.aes_category == "reduced by ancestry":
            colors.append("#2ca02c")
        else:
            colors.append("#ff7f0e")

    fig, ax = plt.subplots(figsize=(10, max(4, len(top) * 0.55)))
    bars = ax.barh(names, aes_vals, color=colors, edgecolor="white", height=0.6)
    ax.axvline(x=1.0, color="black", linestyle="--", linewidth=1.0, label="AES = 1 (EUR parity)")
    ax.axvline(x=_AES_DISPLAY_ELEVATED, color="#d62728", linestyle=":", linewidth=0.8)
    ax.axvline(x=_AES_DISPLAY_REDUCED, color="#2ca02c", linestyle=":", linewidth=0.8)

    for bar, val in zip(bars, aes_vals):
        ax.text(val + 0.02, bar.get_y() + bar.get_height() / 2,
                f"{val:.2f}", va="center", fontsize=8)

    ax.set_xlabel("Ancestry Elevation Score (AES) — exploratory")
    ax.set_title("Ancestry-Stratified Disease Signal\n(AES > 1 = ancestry-specific OR exceeds EUR reference; not a validated risk score)")
    ax.legend(fontsize=8)
    plt.tight_layout()
    fig.savefig(figs_dir / "aes_chart.png", dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Demo
# ---------------------------------------------------------------------------


def run_demo(output_dir: Path) -> None:
    """Run with built-in South Asian demo patient."""
    output_dir = Path(output_dir)
    demo_file = DATA_DIR / "demo_patient_south_asian.txt"
    panel = load_aisnp_panel(DATA_DIR / "aisnp_panel.csv")
    associations = load_associations(DATA_DIR / "ancestry_risk_associations.json")

    genotypes = genotypes_to_simple(parse_genetic_file(demo_file))
    ancestry_result = infer_ancestry(genotypes, panel)
    risks, scoring_note = _get_risks_with_confidence_gating(genotypes, ancestry_result, associations)
    generate_report(risks, ancestry_result, output_dir, scoring_note=scoring_note)

    print(f"\n[ancestry-risk-profiler] Demo complete → {output_dir}/ancestry_risk_report.md")
    print(f"Inferred genetic super-population: {ancestry_result['inferred_ancestry']} (confidence: {ancestry_result['confidence']})")
    if scoring_note:
        print(f"\n  Note: {scoring_note}")
    if risks:
        print(f"\nTop diseases by ancestry elevation (exploratory AES):")
        for r in risks[:5]:
            print(f"  {r.disease:<35} Ancestry OR={r.combined_or:.2f}x  EUR ref OR={r.or_eur_combined:.2f}x  AES={r.aes:.2f}  [{r.aes_category}]")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Ancestry-aware disease risk profiler",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--input", help="23andMe or AncestryDNA raw file")
    p.add_argument("--ancestry", help="Override ancestry inference (AFR/AMR/EAS/EUR/SAS)")
    p.add_argument("--output", default="./ancestry_risk_output", help="Output directory")
    p.add_argument("--demo", action="store_true", help="Run with built-in South Asian demo patient")
    return p


def main(argv: list[str] | None = None) -> None:
    args = _build_parser().parse_args(argv)
    output_dir = Path(args.output)

    if args.demo:
        run_demo(output_dir)
        return

    if not args.input:
        print("ERROR: Provide --input <file> or --demo", file=sys.stderr)
        sys.exit(1)

    input_path = Path(args.input)
    if not input_path.exists():
        print(f"ERROR: File not found: {input_path}", file=sys.stderr)
        sys.exit(1)

    panel = load_aisnp_panel(DATA_DIR / "aisnp_panel.csv")
    associations = load_associations(DATA_DIR / "ancestry_risk_associations.json")

    print(f"[ancestry-risk-profiler] Parsing {input_path.name} …")
    genotypes = genotypes_to_simple(parse_genetic_file(input_path))
    print(f"  {len(genotypes):,} SNPs loaded")

    try:
        ancestry_result = infer_ancestry(genotypes, panel, ancestry_override=args.ancestry)
    except InsufficientCoverageError as e:
        print(f"\nERROR: {e}", file=sys.stderr)
        sys.exit(1)

    ancestry = ancestry_result["inferred_ancestry"]
    confidence = ancestry_result["confidence"]
    print(f"  Genetic super-population: {ancestry} (confidence: {confidence})")

    risks, scoring_note = _get_risks_with_confidence_gating(genotypes, ancestry_result, associations)

    if scoring_note:
        print(f"\n  NOTE: {scoring_note}")

    generate_report(risks, ancestry_result, output_dir, scoring_note=scoring_note)

    print(f"\n[ancestry-risk-profiler] Done → {output_dir}/ancestry_risk_report.md")
    if risks:
        print(f"\nTop 5 diseases by ancestry elevation (exploratory AES):")
        for r in risks[:5]:
            print(f"  {r.disease:<35} Ancestry OR={r.combined_or:.2f}x  EUR ref OR={r.or_eur_combined:.2f}x  AES={r.aes:.2f}  [{r.aes_category}]")


if __name__ == "__main__":
    main()
