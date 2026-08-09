"""Prioritisation engine for clinical-variant-prioritizer.

Implements the pathogenicity-screening logic of Corpas et al. 2021
(Front Genet 12:535123): screen a genotype set against curated clinical-gene
panels (OMIM-morbid, ACMG-SF, Hereditary-Cancer), then rank carried variants by
ClinVar significance, gnomAD frequency, inheritance model and zygosity.

The skill does not re-query ClinVar/gnomAD per call; the curated panel ships the
catalogued annotations so the logic is deterministic and offline-reproducible.
"""
from __future__ import annotations

import json
from pathlib import Path

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}

# Lower rank = higher clinical priority (shown first).
CATEGORY_RANK = {
    "affected": 0,
    "actionable": 1,
    "carrier": 2,
    "uncertain": 3,
    "benign": 4,
    "reference": 5,
}

_FINDING_FIELDS = (
    "id", "rsid", "gene", "hgvs_c", "hgvs_p", "consequence",
    "clinvar_significance", "clinvar_review", "clinvar_id", "gnomad_af", "vep",
    "condition", "inheritance", "panels", "source",
)


def load_panel(path: str | Path) -> list[dict]:
    with open(path) as fh:
        return json.load(fh)


def call_zygosity(genotype: str, ref: str, alt: str) -> str | None:
    """Return 'ref' | 'het' | 'hom' for a diploid SNV call, or None if the
    genotype does not resolve to the ref/alt alleles on either strand."""
    g = (genotype or "").strip().upper()
    if len(g) < 2:
        return None
    bases = [g[0], g[1]]
    if set(bases) <= {ref, alt}:
        n_alt = bases.count(alt)
    else:
        flipped = [COMPLEMENT.get(b, b) for b in bases]
        if set(flipped) <= {ref, alt}:
            n_alt = flipped.count(alt)
        else:
            return None
    return {0: "ref", 1: "het", 2: "hom"}[n_alt]


def classify(entry: dict, zygosity: str) -> tuple[str, str]:
    """Map a carried variant to a clinical category + plain-language rationale."""
    if zygosity == "ref":
        return "reference", "Screened; reference genotype, variant allele not carried."

    sig = (entry.get("clinvar_significance") or "").lower()
    inh = (entry.get("inheritance") or "").lower()
    sev = (entry.get("condition_severity") or "").lower()
    is_pathogenic = "pathogenic" in sig and "conflicting" not in sig
    is_benign = "benign" in sig and "pathogenic" not in sig
    is_uncertain = ("uncertain" in sig) or ("conflicting" in sig)

    if is_benign:
        return "benign", "Carried, but ClinVar classifies the allele benign; no action."

    if is_uncertain:
        if inh == "ar" and zygosity == "het" and sev in ("benign_trait", "low_penetrance"):
            return ("carrier",
                    "Heterozygous carrier of a recessive, largely benign-spectrum trait; "
                    "conflicting ClinVar status, common in population. Not personally actionable.")
        return ("uncertain",
                "Variant of uncertain/conflicting significance; insufficient evidence to act, "
                "flagged for transparency and future reclassification.")

    if is_pathogenic:
        if inh == "ar":
            if zygosity == "hom":
                return ("affected",
                        "Homozygous for a recessive pathogenic allele; clinically relevant, confirm with a clinical assay.")
            return ("carrier",
                    "Heterozygous carrier of a recessive pathogenic allele; reproductive-risk relevant, not personally actionable.")
        if inh in ("ad", "risk_factor", "dominant"):
            return ("actionable",
                    "Pathogenic allele carried in a dominant/risk gene; clinically actionable, confirm with a clinical assay.")
        return ("actionable", "Pathogenic allele carried; clinical follow-up warranted.")

    return "uncertain", "Carried; classification context incomplete."


def screen(genotypes: dict[str, str], panel: list[dict]) -> tuple[list[dict], dict]:
    """Screen genotypes against the panel and return (ranked findings, summary)."""
    counts = {c: 0 for c in ("actionable", "affected", "carrier", "uncertain", "benign", "reference")}
    not_tested = 0
    findings: list[dict] = []

    for entry in panel:
        key, rsid = entry["id"], entry.get("rsid")
        geno = genotypes.get(key)
        if geno is None and rsid:
            geno = genotypes.get(rsid)
        if geno is None:
            not_tested += 1
            continue

        zyg = call_zygosity(geno, entry["ref_allele"], entry["alt_allele"])
        if zyg is None:
            not_tested += 1
            continue

        category, rationale = classify(entry, zyg)
        counts[category] += 1
        finding = {k: entry.get(k) for k in _FINDING_FIELDS}
        finding.update({
            "genotype": (geno or "").strip().upper(),
            "zygosity": zyg,
            "category": category,
            "priority": CATEGORY_RANK[category],
            "rationale": rationale,
        })
        findings.append(finding)

    findings.sort(key=lambda f: (f["priority"], f["gene"]))

    loci_carried = sum(counts[c] for c in ("actionable", "affected", "carrier", "uncertain", "benign"))
    summary = {
        "panels": sorted({p for e in panel for p in e.get("panels", [])}),
        "panel_size": len(panel),
        "loci_tested": len(findings),
        "loci_carried": loci_carried,
        "reference": counts["reference"],
        "not_tested": not_tested,
        "actionable": counts["actionable"],
        "affected": counts["affected"],
        "carriers": counts["carrier"],
        "uncertain": counts["uncertain"],
        "benign": counts["benign"],
    }
    return findings, summary


def build_headline(summary: dict) -> str:
    if summary["loci_tested"] == 0:
        return "No catalogued panel loci were present in the input genotypes."
    flagged = []
    if summary["actionable"] or summary["affected"]:
        flagged.append(f"{summary['actionable'] + summary['affected']} actionable/affected")
    if summary["carriers"]:
        flagged.append(f"{summary['carriers']} carrier")
    if summary["uncertain"]:
        flagged.append(f"{summary['uncertain']} uncertain")
    lead = ("No actionable pathogenic variant carried"
            if (summary["actionable"] == 0 and summary["affected"] == 0)
            else "Actionable finding present")
    detail = ", ".join(flagged) if flagged else "nothing flagged"
    return (f"{lead}. Across {summary['loci_tested']} screened loci: "
            f"{detail}; {summary['reference']} reference.")
