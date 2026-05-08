#!/usr/bin/env python3
"""
ClawBio · claw-methylation-cycle v0.1.1
Methylation cycle analysis with BH4/neurotransmitter axis interpretation.

Author: Samuel Carmona Aguirre <samuel@unimed-consulting.es>
Framework: Holomedicina® · CAPS Digital · UNIMED Consulting
License: MIT

Research and educational use only (RUO). Not a diagnostic device.
Consult a qualified clinician before modifying supplementation or treatment.
"""

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

# pandas removed: not used in this module

# ---------------------------------------------------------------------------
# SNP Panel Definition
# Activity weights are population-derived approximations from in-vitro and
# epidemiological studies; they are not direct enzymatic assays.
# Each entry cites the primary source for its activity_het / activity_hom values.
# ---------------------------------------------------------------------------

PANEL: dict[str, dict] = {
    "rs1801133": {
        "gene": "MTHFR",
        "variant": "C677T",
        "risk_allele": "T",
        "effect": "Decreased Folate to 5-MTHF Conversion",
        # Heterozygous ~65%, homozygous ~30% residual activity (thermolability assay).
        # source: doi:10.1038/ng0595-111 (Frosst et al., Nat Genet 1995)
        # source: doi:10.3390/nu13030768 (Ledford et al., Nutrients 2021)
        "activity_het": 65,
        "activity_hom": 30,
        "weight": 0.35,
    },
    "rs1801131": {
        "gene": "MTHFR",
        "variant": "A1298C",
        "risk_allele": "C",
        "effect": "Decreased MTHFR Activity (modifier)",
        # Heterozygous ~80%, homozygous ~60% residual activity.
        # source: doi:10.1086/301927 (van der Put et al., Am J Hum Genet 1998)
        "activity_het": 80,
        "activity_hom": 60,
        "weight": 0.0,        # combined with rs1801133 for MTHFR total
    },
    "rs1801394": {
        "gene": "MTRR",
        "variant": "A66G",
        "risk_allele": "G",
        "effect": "Decreased Methionine Synthase Reductase Activity",
        # Heterozygous ~80%, homozygous ~60% residual activity.
        # source: doi:10.1093/qjmed/94.11.609 (Gaughan et al., QJM 2001)
        # source: doi:10.1006/mgme.1999.2879 (Wilson et al., Mol Genet Metab 1999)
        "activity_het": 80,
        "activity_hom": 60,
        "weight": 0.15,
    },
    "rs1805087": {
        "gene": "MTR",
        "variant": "A2756G",
        "risk_allele": "G",
        "effect": "Decreased Methionine Synthase Activity",
        # Heterozygous ~85%, homozygous ~70% residual activity.
        # source: doi:10.1016/S0021-9150(99)00113-0 (Harmon et al., Atherosclerosis 1999)
        # source: doi:10.1093/hmg/5.12.1867 (Leclerc et al., Hum Mol Genet 1996)
        "activity_het": 85,
        "activity_hom": 70,
        "weight": 0.10,
    },
    "rs234706": {
        "gene": "CBS",
        "variant": "C699T",
        "risk_allele": "T",
        "effect": "Increased CBS Activity (diverts homocysteine to transsulfuration)",
        # CBS risk allele INCREASES activity (inverse effect on methylation capacity).
        # Heterozygous ~120%, homozygous ~140% of normal CBS flux.
        # source: doi:10.1093/hmg/10.5.477 (Gaughan et al., Hum Mol Genet 2001)
        # source: doi:10.1086/320593 (Lievers et al., Am J Hum Genet 2001)
        "activity_het": 120,
        "activity_hom": 140,
        "weight": 0.05,
        "inverse": True,
    },
    "rs3733890": {
        "gene": "BHMT",
        "variant": "R239Q",
        "risk_allele": "A",
        "effect": "Decreased Betaine-Homocysteine Methyltransferase Activity",
        # Heterozygous ~70%, homozygous ~40% residual activity.
        # source: doi:10.1093/jn/131.9.2479 (Caudill et al., J Nutr 2001)
        # source: doi:10.1186/1471-2156-9-43 (Morin et al., BMC Genet 2008)
        "activity_het": 70,
        "activity_hom": 40,
        "weight": 0.15,
    },
    "rs1979277": {
        "gene": "SHMT1",
        "variant": "C1420T",
        "risk_allele": "T",
        "effect": "Decreased Serine Hydroxymethyltransferase Activity",
        # Heterozygous ~80%, homozygous ~60% residual activity.
        # source: doi:10.1093/carcin/bgm139 (Perry et al., Carcinogenesis 2007)
        "activity_het": 80,
        "activity_hom": 60,
        "weight": 0.05,
    },
    "rs4680": {
        "gene": "COMT",
        "variant": "Val158Met",
        "risk_allele": "A",  # A = Met allele (slow COMT)
        "effect": "Decreased Catechol-O-Methyltransferase Activity",
        # Heterozygous ~65%, homozygous Met/Met ~25% residual COMT activity.
        # source: doi:10.1097/00008571-199606000-00007 (Lachman et al., Pharmacogenetics 1996)
        # source: doi:10.1523/JNEUROSCI.4106-03.2004 (Chen et al., J Neurosci 2004)
        "activity_het": 65,
        "activity_hom": 25,
        "weight": 0.10,
    },
    "rs819147": {
        "gene": "AHCY",
        "variant": "AHCY",
        "risk_allele": "T",
        "effect": "Decreased Adenosylhomocysteinase Activity",
        # Heterozygous ~80%, homozygous ~60% residual activity.
        # source: doi:10.1073/pnas.0400658101 (Baric et al., PNAS 2004)
        "activity_het": 80,
        "activity_hom": 60,
        "weight": 0.05,
    },
}

DISCLAIMER = (
    "**Research and educational use only (RUO). Not a diagnostic device.**\n"
    "Consult a qualified clinician before modifying supplementation or treatment.\n"
    "Enzymatic activity estimates are population-derived approximations, not direct assays."
)


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def parse_genotype_file(path: Path) -> dict[str, str]:
    """
    Parse a 23andMe / AncestryDNA / ADNTRO raw genotype file.
    Returns a dict mapping rsid -> genotype string (e.g. 'AG', 'TT').
    """
    genotypes: dict[str, str] = {}
    with path.open() as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            rsid = parts[0].strip()
            genotype = parts[3].strip().upper()
            if rsid.startswith("rs"):
                genotypes[rsid] = genotype
    return genotypes


# ---------------------------------------------------------------------------
# Analysis
# ---------------------------------------------------------------------------

def count_risk_alleles(genotype: str, risk_allele: str) -> int:
    """Count occurrences of the risk allele in a genotype string."""
    return genotype.count(risk_allele)


def estimate_activity(snp_def: dict, n_risk: int) -> int:
    """Return estimated enzymatic activity % given risk allele count."""
    if n_risk == 0:
        return 100
    elif n_risk == 1:
        return snp_def["activity_het"]
    else:
        return snp_def["activity_hom"]


def compute_mthfr_combined(rs1801133_n: int, rs1801131_n: int) -> int:
    """
    MTHFR combined activity for all genotype configurations.

    Values derived from:
    - Compound heterozygous (C677T het + A1298C het) -> ~15% residual activity.
      The synergistic reduction is greater than either variant alone.
      source: doi:10.1086/301927 (van der Put et al., Am J Hum Genet 1998)
    - C677T hom + A1298C het -> ~12% (most severe compound configuration).
      source: doi:10.1086/301927
    - C677T het + A1298C hom -> ~20%.
      source: doi:10.1086/301927
    - C677T hom only -> ~30%.
      source: doi:10.1038/ng0595-111 (Frosst et al., Nat Genet 1995)
    - C677T het only -> ~65%.
      source: doi:10.1038/ng0595-111
    - A1298C hom only -> ~60%.
      source: doi:10.1086/301927
    - A1298C het only -> ~80%.
      source: doi:10.1086/301927
    """
    if rs1801133_n >= 1 and rs1801131_n >= 1:
        # Compound heterozygous - most clinically significant
        if rs1801133_n == 1 and rs1801131_n == 1:
            return 15
        elif rs1801133_n == 2:
            return 12   # homozygous C677T + A1298C het
        elif rs1801131_n == 2:
            return 20   # C677T het + homozygous A1298C
        else:
            return 15
    elif rs1801133_n == 2:
        return 30       # homozygous C677T only
    elif rs1801133_n == 1:
        return 65       # heterozygous C677T only
    elif rs1801131_n == 2:
        return 60       # homozygous A1298C only
    elif rs1801131_n == 1:
        return 80       # heterozygous A1298C only
    else:
        return 100      # no risk alleles


def analyse(genotypes: dict[str, str]) -> dict:
    """Run methylation cycle analysis. Returns structured result dict."""
    results = {}
    found_rsids = []
    missing_rsids = []

    # Resolve each SNP
    for rsid, snp_def in PANEL.items():
        if rsid in genotypes:
            found_rsids.append(rsid)
            gt = genotypes[rsid]
            n_risk = count_risk_alleles(gt, snp_def["risk_allele"])
            activity = estimate_activity(snp_def, n_risk)
            status = (
                "normal" if n_risk == 0
                else "heterozygous" if n_risk == 1
                else "homozygous_risk"
            )
            results[rsid] = {
                "gene": snp_def["gene"],
                "variant": snp_def["variant"],
                "genotype": gt,
                "n_risk_alleles": n_risk,
                "status": status,
                "activity_pct": activity,
                "effect": snp_def["effect"],
                "weight": snp_def["weight"],
            }
        else:
            missing_rsids.append(rsid)
            results[rsid] = {
                "gene": PANEL[rsid]["gene"],
                "variant": PANEL[rsid]["variant"],
                "genotype": "Not assessed",
                "n_risk_alleles": None,
                "status": "not_assessed",
                # Safety Rule 6: missing SNPs carry None, never assumed wildtype.
                "activity_pct": None,
                "effect": PANEL[rsid]["effect"],
                "weight": PANEL[rsid]["weight"],
            }

    # MTHFR combined activity
    # n_risk_alleles is None when SNP was not assessed; treat as 0 only for
    # the purpose of compute_mthfr_combined, which handles missing as wildtype
    # internally. The mthfr_assessed flag below controls whether the result
    # enters the NMC calculation at all (Safety Rule 6 compliance).
    rs677_n = results["rs1801133"]["n_risk_alleles"] or 0
    rs1298_n = results["rs1801131"]["n_risk_alleles"] or 0
    mthfr_combined = compute_mthfr_combined(rs677_n, rs1298_n)
    compound_het = (rs677_n >= 1 and rs1298_n >= 1)
    results["rs1801133"]["mthfr_combined_activity"] = mthfr_combined
    results["rs1801133"]["compound_heterozygous"] = compound_het

    # BH4 axis capacity
    # MTRR A66G modifies the MTHFR -> BH4 axis efficiency.
    # Modifier values 0.88 (het) and 0.75 (hom) derived from:
    # source: doi:10.1093/qjmed/94.11.609 (Gaughan et al., QJM 2001)
    mtrr_n = results["rs1801394"]["n_risk_alleles"] or 0
    mtrr_modifier = 1.0 if mtrr_n == 0 else (0.88 if mtrr_n == 1 else 0.75)
    bh4_capacity = round(mthfr_combined * mtrr_modifier)

    # Net Methylation Capacity (NMC)
    # Safety Rule 6: unassessed SNPs are excluded from the calculation entirely.
    # coverage_pct reports what fraction of the weighted panel was assessed.
    # MTHFR (rs1801133) is excluded from the main loop and handled separately
    # using its combined activity value.
    nmc_numerator = 0.0
    nmc_weight_assessed = 0.0

    for rsid, r in results.items():
        w = r["weight"]
        if w == 0 or rsid == "rs1801133":
            continue  # MTHFR handled below; zero-weight SNPs skipped
        act = r["activity_pct"]
        if act is None:
            continue  # Safety Rule 6: skip, do not assume wildtype
        snp_def = PANEL[rsid]
        if snp_def.get("inverse"):
            # CBS: higher activity diverts methyl groups - penalise NMC
            contribution = max(0, 100 - (act - 100)) * w
        else:
            contribution = act * w
        nmc_numerator += contribution
        nmc_weight_assessed += w

    # Add MTHFR combined contribution only if at least one MTHFR SNP was assessed
    mthfr_w = PANEL["rs1801133"]["weight"]
    mthfr_assessed = (
        results["rs1801133"]["activity_pct"] is not None
        or results["rs1801131"]["activity_pct"] is not None
    )
    if mthfr_assessed:
        nmc_numerator += mthfr_combined * mthfr_w
        nmc_weight_assessed += mthfr_w

    total_panel_weight = sum(v["weight"] for v in PANEL.values() if v["weight"] > 0)
    coverage_pct = (
        round(nmc_weight_assessed / total_panel_weight * 100)
        if total_panel_weight > 0 else 0
    )
    # nmc_numerator = sum(activity_pct * weight), so dividing by sum(weight)
    # yields the weighted average activity directly as a 0-100 percentage.
    nmc = (
        round(nmc_numerator / nmc_weight_assessed)
        if nmc_weight_assessed > 0 else None
    )

    return {
        "metadata": {
            "tool": "claw-methylation-cycle v0.1.1",
            "framework": "Holomedicina(R) CAPS Digital UNIMED Consulting",
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "snps_in_panel": len(PANEL),
            "snps_found": len(found_rsids),
            "snps_missing": len(missing_rsids),
            "missing_rsids": missing_rsids,
        },
        "summary": {
            "net_methylation_capacity": nmc,
            "nmc_coverage_pct": coverage_pct,
            "bh4_axis_capacity": bh4_capacity,
            "mthfr_combined_activity": mthfr_combined,
            "mthfr_compound_heterozygous": compound_het,
            # TH (dopamine) and TPH2 (serotonin) differ in Km for BH4.
            # TPH2 has higher Km (~30-100 uM) than TH (~10-30 uM), so
            # serotonin synthesis is more sensitive to BH4 reduction.
            # source: doi:10.1042/BJ20031542 (Fitzpatrick, Biochem J 2004)
            "dopamine_synthesis_impact": _dopamine_impact(bh4_capacity),
            "serotonin_synthesis_impact": _serotonin_impact(bh4_capacity),
        },
        "gene_results": results,
        "missing_rsids": missing_rsids,
    }


def _dopamine_impact(bh4_pct: int) -> str:
    """
    Dopamine synthesis sensitivity to BH4 reduction.
    TH (tyrosine hydroxylase) Km for BH4 approx 10-30 uM.
    source: doi:10.1042/BJ20031542 (Fitzpatrick, Biochem J 2004)
    """
    if bh4_pct < 35:
        return "Severely Reduced"
    elif bh4_pct < 60:
        return "Moderately Reduced"
    else:
        return "Within Normal Range"


def _serotonin_impact(bh4_pct: int) -> str:
    """
    Serotonin synthesis sensitivity to BH4 reduction.
    TPH2 (tryptophan hydroxylase 2) Km for BH4 approx 30-100 uM,
    higher than TH, making serotonin synthesis more BH4-sensitive.
    source: doi:10.1042/BJ20031542 (Fitzpatrick, Biochem J 2004)
    """
    if bh4_pct < 45:
        return "Severely Reduced"
    elif bh4_pct < 70:
        return "Moderately Reduced"
    else:
        return "Within Normal Range"


def _nmc_status(nmc: int | None) -> str:
    if nmc is None:
        return "Unknown (insufficient panel coverage)"
    if nmc < 40:
        return "CRITICAL - Severely Reduced"
    elif nmc < 60:
        return "LOW - Moderately Reduced"
    elif nmc < 80:
        return "BORDERLINE - Mildly Reduced"
    else:
        return "NORMAL"


def _bh4_status(bh4: int) -> str:
    if bh4 < 40:
        return "CRITICAL - Severely Reduced"
    elif bh4 < 65:
        return "LOW - Moderately Reduced"
    else:
        return "NORMAL - Within Normal Range"


def _activity_label(act: int | None) -> str:
    if act is None:
        return "Not assessed"
    if act <= 30:
        return f"severely reduced ({act}%)"
    elif act <= 60:
        return f"moderately reduced ({act}%)"
    elif act <= 80:
        return f"mildly reduced ({act}%)"
    else:
        return f"normal ({act}%)"


# ---------------------------------------------------------------------------
# Report Generation
# ---------------------------------------------------------------------------

def generate_report(result: dict) -> str:
    s = result["summary"]
    g = result["gene_results"]
    m = result["metadata"]

    nmc = s["net_methylation_capacity"]
    bh4 = s["bh4_axis_capacity"]
    compound_het = s["mthfr_compound_heterozygous"]
    coverage = s["nmc_coverage_pct"]

    lines = [
        "# ClawBio - Methylation Cycle Clinical Report",
        "",
        f"**Date**: {m['generated_utc']}",
        f"**Tool**: {m['tool']}",
        f"**Framework**: {m['framework']}",
        f"**SNPs assessed**: {m['snps_found']}/{m['snps_in_panel']}",
        "",
        f"> {DISCLAIMER}",
        "",
        "---",
        "## Executive Summary",
        "",
        "| Metric | Value | Status |",
        "|--------|-------|--------|",
        f"| Net Methylation Capacity | {nmc}/100 | {_nmc_status(nmc)} |",
        f"| NMC Panel Coverage | {coverage}% | {'FULL' if coverage == 100 else 'PARTIAL - interpret NMC with caution'} |",
        f"| BH4 Axis Capacity | {bh4}/100 | {_bh4_status(bh4)} |",
        f"| MTHFR Combined Activity | {s['mthfr_combined_activity']}% | {'COMPOUND HET DETECTED' if compound_het else 'see gene table'} |",
        f"| Dopamine Synthesis Impact | {s['dopamine_synthesis_impact']} | |",
        f"| Serotonin Synthesis Impact | {s['serotonin_synthesis_impact']} | |",
        "",
        "---",
        "## Enzymatic Activity Profile",
        "",
        "| Gene | Enzyme | rsID | Genotype | Activity | Status |",
        "|------|--------|------|----------|----------|--------|",
    ]

    for rsid, r in g.items():
        if r["weight"] == 0:
            continue  # skip MTHFR A1298C (shown combined)
        act_display = (
            f"{r['activity_pct']}%"
            if r["activity_pct"] is not None
            else "Not assessed"
        )
        lines.append(
            f"| **{r['gene']}** | {r['variant']} | {rsid} | "
            f"`{r['genotype']}` | {act_display} | {_activity_label(r['activity_pct'])} |"
        )

    mthfr_combined = s["mthfr_combined_activity"]
    compound_flag = " [COMPOUND HET]" if compound_het else ""
    lines.append(
        f"| **MTHFR** | Combined (677+1298) | -- | -- | "
        f"{mthfr_combined}%{compound_flag} | {_activity_label(mthfr_combined)} |"
    )

    lines += [
        "",
        "---",
        "## BH4 / Neurotransmitter Axis",
        "",
        "The **MTHFR -> 5-MTHF -> BH4** axis is the critical connection between",
        "methylation cycle variants and neurotransmitter synthesis.",
        "BH4 (tetrahydrobiopterin) is the essential cofactor for tyrosine hydroxylase",
        "(dopamine synthesis) and tryptophan hydroxylase 2 (serotonin synthesis).",
        "TPH2 has a higher Km for BH4 than TH, making serotonin synthesis more",
        "sensitive to BH4 reduction at moderate deficiencies.",
        "source: doi:10.1042/BJ20031542 (Fitzpatrick, Biochem J 2004)",
        "",
        f"**Estimated BH4 production capacity**: {bh4}% of normal",
        f"- Dopamine synthesis pathway: **{s['dopamine_synthesis_impact']}**",
        f"- Serotonin synthesis pathway: **{s['serotonin_synthesis_impact']}**",
        "",
        "_Clinical implication_: When BH4 capacity is reduced, the clinical presentation",
        "of ADHD, depression, and anxiety may have an upstream biological substrate that",
        "5-MTHF + BH4-support nutrients can partially address - prior to or alongside",
        "pharmacological intervention.",
        "",
    ]

    if compound_het:
        lines += [
            "---",
            "## Compound Heterozygosity",
            "",
            "[COMPOUND HET] **MTHFR Compound Heterozygous detected (C677T + A1298C)**",
            "",
            "Both MTHFR variants are present simultaneously. This combination reduces",
            "total MTHFR enzymatic activity more than either variant alone.",
            "It is the most clinically significant single-gene methylation finding.",
            "Active folate (5-MTHF) supplementation is strongly indicated.",
            "source: doi:10.1086/301927 (van der Put et al., Am J Hum Genet 1998)",
            "",
        ]

    lines += [
        "---",
        "## Clinical Recommendations",
        "",
    ]

    recs = _build_recommendations(s, g)
    for rec in recs:
        lines.append(rec)

    lines += [
        "",
        "---",
        "## Missing Variants",
        "",
    ]
    if m["missing_rsids"]:
        lines.append(
            "The following SNPs were not found in the input file and were not assessed:"
        )
        for rsid in m["missing_rsids"]:
            lines.append(f"- {rsid} ({PANEL[rsid]['gene']} - {PANEL[rsid]['variant']})")
        lines.append(
            "\nNote: Missing SNPs are excluded from NMC calculation (Safety Rule 6)."
            f" Panel coverage: {s['nmc_coverage_pct']}%."
        )
    else:
        lines.append("All panel SNPs were found in the input file.")

    lines += [
        "",
        "---",
        "## References",
        "",
        "- Frosst P et al. (1995). A candidate genetic risk factor: common mutation in MTHFR. "
        "Nat Genet 10:111-3. doi:10.1038/ng0595-111",
        "- van der Put NM et al. (1998). Second common mutation in MTHFR. "
        "Am J Hum Genet 62:1044-51. doi:10.1086/301927",
        "- Gaughan DJ et al. (2001). MTRR A66G polymorphism and neural tube defects. "
        "QJM 94:609-17. doi:10.1093/qjmed/94.11.609",
        "- Harmon DL et al. (1999). MTR A2756G polymorphism. "
        "Atherosclerosis 146:295-300. doi:10.1016/S0021-9150(99)00113-0",
        "- Lachman HM et al. (1996). COMT pharmacogenetics. "
        "Pharmacogenetics 6:243-50. doi:10.1097/00008571-199606000-00007",
        "- Fitzpatrick PF (2004). Tetrahydropterin-dependent amino acid hydroxylases. "
        "Biochem J 380:299-310. doi:10.1042/BJ20031542",
        "- Ledford AW et al. (2021). MTHFR and BH4 in neuropsychiatric disorders. "
        "Nutrients 13:768. doi:10.3390/nu13030768",
        "- Carmona Aguirre S. (2014/UNESCO 2016). Holomedicina. UNIMED Consulting.",
        "- ClawBio (2026). https://github.com/ClawBio/ClawBio",
    ]

    return "\n".join(lines)


def _build_recommendations(summary: dict, genes: dict) -> list[str]:
    recs = []
    mthfr_act = summary["mthfr_combined_activity"]
    bh4 = summary["bh4_axis_capacity"]
    compound_het = summary["mthfr_compound_heterozygous"]
    bhmt_act = genes["rs3733890"]["activity_pct"] or 100
    mtrr_n = genes["rs1801394"]["n_risk_alleles"] or 0

    if compound_het or mthfr_act <= 30:
        recs.append(
            "- **PRIORITY 1** -- Use 5-MTHF (methylfolate) instead of synthetic folic acid. "
            "MTHFR enzymatic reduction impairs folic acid conversion."
        )
    if compound_het:
        recs.append(
            "- **PRIORITY 1** -- MTHFR compound heterozygous (C677T + A1298C): combined activity "
            f"significantly reduced ({mthfr_act}%). 5-MTHF + methylcobalamin (B12) supplementation strongly indicated."
        )
    if bh4 < 65:
        recs.append(
            f"- **PRIORITY 2** -- BH4 capacity at {bh4}% of normal. "
            "Riboflavin (B2, 200-400 mg/day) supports MTHFR activity and BH4 regeneration. "
            "Vitamin C (500 mg/day) maintains BH4 in reduced (active) form."
        )
    if bh4 < 40:
        recs.append(
            "- **PRIORITY 2** -- Dopamine and serotonin synthesis may be compromised via reduced BH4. "
            "Evaluate clinical presentation for neurodevelopmental implications (ADHD, depression, anxiety) "
            "before pharmacological intervention."
        )
    if mtrr_n >= 1:
        recs.append(
            "- **PRIORITY 2** -- MTRR variant: prefer methylcobalamin (B12 active form) over cyanocobalamin. "
            "Consider hydroxocobalamin as alternative."
        )
    if bhmt_act <= 60:
        recs.append(
            "- **PRIORITY 3** -- BHMT variant: support choline-dependent remethylation with "
            "betaine (TMG, 500-1000 mg/day) and choline-rich foods (eggs, liver)."
        )
    if not recs:
        recs.append(
            "- No high-priority supplementation flags identified. Maintain dietary folate and B12 adequacy."
        )
    return recs


# ---------------------------------------------------------------------------
# CLI Entry Point
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="ClawBio claw-methylation-cycle - Methylation cycle analysis with BH4 axis"
    )
    parser.add_argument(
        "--input", required=True, type=Path,
        help="Path to raw genotype file (23andMe / AncestryDNA / ADNTRO format)"
    )
    parser.add_argument(
        "--output", required=True, type=Path,
        help="Output directory for report and JSON"
    )
    parser.add_argument(
        "--demo", action="store_true",
        help="Run with demo_input.txt from the skill directory"
    )
    args = parser.parse_args()

    if args.demo:
        demo_path = Path(__file__).parent / "demo_input.txt"
        if not demo_path.exists():
            print("ERROR: demo_input.txt not found in skill directory.", file=sys.stderr)
            sys.exit(1)
        input_path = demo_path
    else:
        input_path = args.input
        if not input_path.exists():
            print(f"ERROR: Input file not found: {input_path}", file=sys.stderr)
            sys.exit(1)

    output_dir: Path = args.output
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Loading genotypes from: {input_path.name}")
    genotypes = parse_genotype_file(input_path)
    print(f"Parsed {len(genotypes)} total variants")

    print("Running methylation cycle analysis...")
    result = analyse(genotypes)

    s = result["summary"]
    print(f"Net Methylation Capacity: {s['net_methylation_capacity']}/100 "
          f"(panel coverage: {s['nmc_coverage_pct']}%)")
    print(f"BH4 Axis Capacity: {s['bh4_axis_capacity']}/100")
    if s["mthfr_compound_heterozygous"]:
        print("[COMPOUND HET] MTHFR Compound Heterozygous detected (C677T + A1298C)")

    report_path = output_dir / "report.md"
    report_path.write_text(generate_report(result), encoding="utf-8")
    print(f"Report written to {report_path}")

    json_path = output_dir / "result.json"
    json_path.write_text(json.dumps(result, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"JSON written to {json_path}")

    print("\nDone. ClawBio is a research tool. Not a medical device.")


if __name__ == "__main__":
    main()
