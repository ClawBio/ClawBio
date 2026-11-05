#!/usr/bin/env python3
"""Step 4 of pharmgx_reporter: render the markdown report + result.json
from the outputs of `parse_and_profile.py`, `lookup_drugs.py`, and
`enrich_clinpgx.py`.

Ported from ClawBio 1.x `pharmgx_reporter.py`'s `generate_report()`
(DATA QUALITY WARNING, Panel Limitations, Gene Profiles, Detected Variants,
Complete Drug Recommendations, Disclaimer, Methods, References sections —
all unchanged) plus its `main()`'s `write_result_json()` call. Two changes
from 1.x:

  - No `generate_html_report()` — see pgx_common.py's module docstring for
    why the HTML report (and the ClinPGx `source_recs` guideline-table
    parsing that fed part of it) isn't ported.
  - A new "ClinPGx Evidence" section, populated from `enrich_clinpgx.py`'s
    output when non-empty. 1.x computed this evidence data but only ever
    rendered it into the (here, omitted) HTML report — folding it into the
    markdown report instead of dropping it keeps the enrichment step's
    output meaningfully consumed by the DAG rather than an orphaned emit.
    When the enrichment step ran offline (empty `enrichment` dict — see
    enrich_clinpgx.py), this section prints a one-line note instead of an
    empty table.

Usage:
    python report.py --profile /tmp/out/profile.json \\
        --drug-recommendations /tmp/out/drug_recommendations.json \\
        --enrichment /tmp/out/clinpgx_enrichment.json \\
        --out /tmp/out
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import pgx_common as pgx

ICON = {
    "standard": "STANDARD",
    "caution": "CAUTION",
    "avoid": "AVOID",
    "indeterminate": "INDETERMINATE - INSUFFICIENT DATA",
}

# Per-gene limitation warnings for genes with known CNV/SV/repeat issues that
# DTC SNP arrays cannot detect. Must appear in the report body for disclosure
# compliance, not just as a stderr warning.
_GENE_LIMITATIONS = {
    "CYP2D6": (
        "CYP2D6: Copy number variation (gene deletion CYP2D6*5, duplication CYP2D6*1xN/*2xN, "
        "hybrid alleles CYP2D6*13/*36) cannot be detected from DTC genotyping data. "
        "Phenotype assignment may be incomplete. Clinical-grade CYP2D6 testing includes CNV analysis."
    ),
    "UGT1A1": (
        "UGT1A1: The UGT1A1*28 allele (rs8175347) is a TA-repeat polymorphism that cannot be "
        "reliably genotyped by SNP arrays. If rs8175347 is present in the input, the reported "
        "genotype is a proxy SNP call and may not reflect true TA repeat count."
    ),
    "CYP3A5": (
        "CYP3A5: The CYP3A5*7 allele (rs41303343) is an insertion/deletion variant. "
        "DTC genotyping platforms use proxy SNP calls that may not accurately detect this allele."
    ),
}

_OUT_OF_PANEL = [
    ("HLA-B", "HLA-B*57:01", "Not assessed", "Indeterminate (not in panel)"),
    ("MT-RNR1", "MT-RNR1 (mitochondrial)", "Not assessed", "Indeterminate (not in panel)"),
    ("G6PD", "Glucose-6-Phosphate Dehydrogenase", "Not assessed", "Indeterminate (not in panel)"),
    ("HLA-A", "HLA-A*31:01", "Not assessed", "Indeterminate (not in panel)"),
]


def generate_report(profile_data, drug_results, enrichment) -> str:
    fmt = profile_data["format"]
    total_snps = profile_data["total_snps"]
    pgx_snps = profile_data["pgx_snps"]
    profiles = profile_data["profiles"]
    checksum = profile_data["input_checksum"]
    fname = Path(profile_data["input_path"]).name
    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")

    lines = []
    lines.append("# ClawBio PharmGx Report")
    lines.append("")
    lines.append(f"**Date**: {now}")
    lines.append(f"**Input**: `{fname}`")
    lines.append(f"**Format detected**: {fmt}")
    lines.append(f"**Checksum (SHA-256)**: `{checksum}`")
    lines.append(f"**Total SNPs in file**: {total_snps}")
    lines.append(f"**Pharmacogenomic SNPs found**: {len(pgx_snps)}/{len(pgx.PGX_SNPS)}")
    lines.append(f"**Genes profiled**: {len(profiles)}")
    lines.append(f"**Drugs assessed**: {sum(len(v) for v in drug_results.values())}")
    lines.append("")
    lines.append("---")
    lines.append("")

    # Data quality warning
    not_tested = [g for g, p in profiles.items() if p["diplotype"] == "NOT_TESTED"]
    unknown_pheno = [g for g, p in profiles.items()
                     if "unknown" in p["phenotype"].lower() or "indeterminate" in p["phenotype"].lower()]

    lines.append("## DATA QUALITY WARNING")
    lines.append("")
    if not_tested:
        lines.append(f"**{len(not_tested)} gene(s) could not be assessed** because the "
                     "relevant SNPs were not found in the input file: "
                     f"{', '.join(not_tested)}")
        lines.append("")
        lines.append("Drugs depending on these genes are marked INDETERMINATE below. "
                     "Do not assume normal metabolism for untested genes.")
        lines.append("")
    if unknown_pheno:
        unmapped = [g for g in unknown_pheno if g not in not_tested]
        if unmapped:
            lines.append(f"**{len(unmapped)} gene(s) have unmapped diplotypes**: "
                         f"{', '.join(unmapped)}. These diplotypes could not be matched "
                         "to a known phenotype. Clinical pharmacogenomic testing is recommended.")
            lines.append("")

    gene_warnings = []
    for gene_name in sorted(profiles.keys()):
        if gene_name in _GENE_LIMITATIONS:
            gene_warnings.append(_GENE_LIMITATIONS[gene_name])
    if gene_warnings:
        lines.append("**Gene-Specific Limitations:**")
        lines.append("")
        for warning in gene_warnings:
            lines.append(f"- {warning}")
        lines.append("")

    lines.append("---")
    lines.append("")

    # Panel limitations disclosure
    lines.append("## Panel Limitations")
    lines.append("")
    lines.append("This report uses **SNP-based genotyping only** from a panel of "
                 f"{len(pgx.PGX_SNPS)} pharmacogenomic variants. The following cannot be detected:")
    lines.append("")
    lines.append("- **Copy number variants (CNVs)**: Gene deletions (e.g. CYP2D6\\*5) "
                 "and duplications (e.g. CYP2D6\\*1xN, \\*2xN)")
    lines.append("- **Structural variants**: CYP2D6-CYP2D7 hybrid alleles (e.g. \\*13, \\*36)")
    lines.append("- **Repeat polymorphisms**: UGT1A1\\*28 (TA7 repeat in rs8175347)")
    lines.append("- **HLA typing**: HLA-B\\*57:01 (abacavir hypersensitivity)")
    lines.append("- **Mitochondrial variants**: MT-RNR1 m.1555A>G (aminoglycoside ototoxicity)")
    lines.append("- **G6PD deficiency**: G6PD A- and Mediterranean variants")
    lines.append("")
    lines.append("For CYP2D6, a result of 'Normal Metabolizer' does NOT rule out "
                 "gene deletions or duplications that alter metabolizer status. "
                 "Clinical-grade CYP2D6 testing includes CNV analysis.")
    lines.append("")
    lines.append("## Genes Not Assessed by This Panel")
    lines.append("")
    lines.append("The following clinically relevant pharmacogenomic genes are **not included** "
                 "in this panel. Consider targeted testing if indicated:")
    lines.append("")
    lines.append("| Gene | Clinical Relevance | CPIC Guideline |")
    lines.append("|------|-------------------|----------------|")
    lines.append("| HLA-B\\*57:01 | Abacavir hypersensitivity | CPIC Abacavir (2014, updated 2020) |")
    lines.append("| HLA-B\\*58:01 | Allopurinol hypersensitivity (SJS/TEN) | CPIC Allopurinol (Hershfield 2013, updated 2015) |")
    lines.append("| G6PD | Rasburicase, chloroquine contraindication | CPIC G6PD (2014) |")
    lines.append("| MT-RNR1 | Aminoglycoside ototoxicity | FDA label |")
    lines.append("| HLA-A\\*31:01 | Carbamazepine hypersensitivity | CPIC Carbamazepine (2017) |")
    lines.append("| CYP2B6 | Efavirenz metabolism | CPIC Efavirenz (2019) |")
    lines.append("")
    lines.append("---")
    lines.append("")

    # Summary counts
    n_std = len(drug_results["standard"])
    n_cau = len(drug_results["caution"])
    n_avo = len(drug_results["avoid"])
    n_ind = len(drug_results.get("indeterminate", []))
    lines.append("## Drug Response Summary")
    lines.append("")
    lines.append("| Category | Count |")
    lines.append("|----------|-------|")
    lines.append(f"| Standard dosing | {n_std} |")
    lines.append(f"| Use with caution | {n_cau} |")
    lines.append(f"| Avoid / use alternative | {n_avo} |")
    if n_ind > 0:
        lines.append(f"| Insufficient data | {n_ind} |")
    lines.append("")

    if n_avo > 0 or n_cau > 0:
        lines.append("### Actionable Alerts")
        lines.append("")
        if n_avo > 0:
            lines.append("**AVOID / USE ALTERNATIVE:**")
            lines.append("")
            for d in drug_results["avoid"]:
                lines.append(f"- **{d['drug']}** ({d['brand']}) [{d['gene']}]")
            lines.append("")
        if n_cau > 0:
            lines.append("**USE WITH CAUTION:**")
            lines.append("")
            for d in drug_results["caution"]:
                lines.append(f"- **{d['drug']}** ({d['brand']}) [{d['gene']}]")
            lines.append("")

    lines.append("---")
    lines.append("")

    # Gene profiles
    lines.append("## Gene Profiles")
    lines.append("")
    lines.append("| Gene | Full Name | Diplotype | Phenotype |")
    lines.append("|------|-----------|-----------|-----------|")
    for gene in pgx.GENE_DEFS:
        if gene in profiles:
            p = profiles[gene]
            lines.append(f"| {gene} | {pgx.GENE_DEFS[gene]['name']} | {p['diplotype']} | {p['phenotype']} |")
    for gene_sym, full_name, dip, pheno in _OUT_OF_PANEL:
        lines.append(f"| {gene_sym} | {full_name} | {dip} | {pheno} |")
    lines.append("")

    # Detected variants
    lines.append("## Detected Variants")
    lines.append("")
    lines.append("| rsID | Gene | Star Allele | Genotype | Effect |")
    lines.append("|------|------|-------------|----------|--------|")
    for rsid, info in sorted(pgx_snps.items(), key=lambda x: x[1]["gene"]):
        lines.append(f"| {rsid} | {info['gene']} | {info['allele']} | {info['genotype']} | {info['effect']} |")
    lines.append("")

    # Full drug table
    lines.append("---")
    lines.append("")
    lines.append("## Complete Drug Recommendations")
    lines.append("")
    lines.append("| Drug | Brand | Class | Gene | Status |")
    lines.append("|------|-------|-------|------|--------|")
    for cat in ["avoid", "caution", "indeterminate", "standard"]:
        for d in sorted(drug_results.get(cat, []), key=lambda x: x["drug"]):
            status = ICON.get(d["classification"], d["classification"].upper())
            note_suffix = f" — {d['note']}" if d.get("note") else ""
            lines.append(f"| {d['drug']} | {d['brand']} | {d['class']} | {d['gene']} | {status}{note_suffix} |")
    lines.append("")

    # ClinPGx evidence (2.0 addition — see module docstring)
    lines.append("---")
    lines.append("")
    lines.append("## ClinPGx Evidence")
    lines.append("")
    if enrichment:
        lines.append("Evidence level and guideline sources from the live ClinPGx API, for "
                     "drugs flagged AVOID or USE WITH CAUTION above:")
        lines.append("")
        lines.append("| Drug | Evidence Level | CPIC-Verified | Sources |")
        lines.append("|------|-----------------|----------------|---------|")
        flagged = {d["drug"].lower(): d["drug"] for cat in ("avoid", "caution") for d in drug_results.get(cat, [])}
        any_row = False
        for drug_key, drug_name in sorted(flagged.items()):
            entry = enrichment.get(drug_key)
            if not entry:
                continue
            any_row = True
            level = entry.get("evidence_level") or "N/A"
            verified = "Yes" if entry.get("verified") else "No"
            sources = ", ".join(entry.get("sources", [])) or "—"
            lines.append(f"| {drug_name} | {level} | {verified} | {sources} |")
        if not any_row:
            lines.append("| _(no ClinPGx evidence matched the flagged drugs)_ | | | |")
        lines.append("")
    else:
        lines.append("_ClinPGx enrichment was unavailable for this run (offline, or no genes "
                     "required lookup). Drug recommendations above are CPIC classification "
                     "summaries only — see `enrich_clinpgx.py`'s step output for details._")
        lines.append("")

    # Disclaimer
    lines.append("---")
    lines.append("")
    lines.append("## Disclaimer")
    lines.append("")
    lines.append("This report is for **research and educational purposes only**. "
                 "It is NOT a diagnostic device and should NOT be used to make medication decisions "
                 "without consulting a qualified healthcare professional.")
    lines.append("")
    lines.append("Pharmacogenomic recommendations are based on CPIC guidelines (cpicpgx.org). "
                 "DTC genetic tests have limitations: they may not detect all relevant variants, "
                 "and results should be confirmed by clinical-grade testing before clinical use.")
    lines.append("")
    lines.append(f"*{pgx.DISCLAIMER}*")
    lines.append("")

    # Methods
    lines.append("## Methods")
    lines.append("")
    lines.append("- **Tool**: ClawBio PharmGx Reporter (2.0 port)")
    lines.append(f"- **SNP panel**: {len(pgx.PGX_SNPS)} pharmacogenomic variants across {len(pgx.GENE_DEFS)} genes")
    lines.append("- **Star allele calling**: Simplified DTC-compatible algorithm (single-SNP per allele)")
    lines.append("- **Phenotype assignment**: CPIC-based diplotype-to-phenotype mapping")
    lines.append(f"- **Drug guidelines**: {len(pgx.GUIDELINES)} drugs from CPIC (cpicpgx.org), simplified for DTC context")
    lines.append("")

    # Reproducibility
    lines.append("## Reproducibility")
    lines.append("")
    lines.append("```bash")
    lines.append("clawbio start skills/pharmgx_reporter --input patient_file=" + fname)
    lines.append("clawbio next / open / done   # loop until bundle")
    lines.append("clawbio bundle")
    lines.append("```")
    lines.append("")
    lines.append(f"**Input checksum**: `{checksum}`")
    lines.append("")

    # References
    lines.append("## References")
    lines.append("")
    lines.append("- Corpas, M. ClawBio. https://github.com/ClawBio/ClawBio")
    lines.append("- CPIC. Clinical Pharmacogenetics Implementation Consortium. https://cpicpgx.org/")
    lines.append("- Caudle, K.E. et al. (2014). Standardizing terms for clinical pharmacogenetic test results. Genet Med, 16(9), 655-663.")
    lines.append("- PharmGKB / ClinPGx. https://www.pharmgkb.org/")
    lines.append("")

    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description="Render the pharmgx_reporter markdown report + result.json")
    parser.add_argument("--profile", required=True, help="Path to profile.json from parse_and_profile.py")
    parser.add_argument("--drug-recommendations", required=True, help="Path to drug_recommendations.json from lookup_drugs.py")
    parser.add_argument("--enrichment", required=True, help="Path to clinpgx_enrichment.json from enrich_clinpgx.py")
    parser.add_argument("--out", required=True, help="Output directory")
    args = parser.parse_args()

    profile_data = json.loads(Path(args.profile).read_text())
    drug_payload = json.loads(Path(args.drug_recommendations).read_text())
    drug_results = drug_payload["drug_recommendations"]
    enrich_payload = json.loads(Path(args.enrichment).read_text())
    enrichment = enrich_payload.get("enrichment", {})

    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)

    report = generate_report(profile_data, drug_results, enrichment)
    report_path = outdir / "report.md"
    report_path.write_text(report)

    n_std = len(drug_results["standard"])
    n_cau = len(drug_results["caution"])
    n_avo = len(drug_results["avoid"])
    n_ind = len(drug_results.get("indeterminate", []))
    total_assessed = n_std + n_cau + n_avo + n_ind

    result_json_path = pgx.write_result_json(
        output_dir=outdir,
        skill="pharmgx_reporter",
        version="2.0.0",
        summary={
            "total_snps_in_file": profile_data["total_snps"],
            "pgx_snps_found": profile_data["pgx_snps_found"],
            "pgx_snps_total": profile_data["pgx_snps_total"],
            "genes_profiled": len(profile_data["profiles"]),
            "drugs_assessed": total_assessed,
            "drugs_standard": n_std,
            "drugs_caution": n_cau,
            "drugs_avoid": n_avo,
            "drugs_indeterminate": n_ind,
            "clinpgx_enriched": len(enrichment),
            "ref_genome": profile_data["ref_genome"],
        },
        data={
            "gene_profiles": profile_data["profiles"],
            "drug_recommendations": drug_results,
            "clinpgx_enrichment": enrichment,
        },
        input_checksum=profile_data["input_checksum"],
    )

    print(f"Report saved: {report_path}")
    print(f"Result JSON:  {result_json_path}")
    print("Done.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
