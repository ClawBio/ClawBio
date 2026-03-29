#!/usr/bin/env python3
"""Clinical Variant Reporter — ACMG/AMP 2015 variant classification from VCF.

Classifies germline variants according to the ACMG/AMP 28-criteria evidence
framework and generates clinical-grade interpretation reports with evidence
audit trails and ACMG SF v3.2 secondary findings screening.

Usage:
    python clinical_variant_reporter.py --demo --output /tmp/acmg_demo
    python clinical_variant_reporter.py --input patient.vcf --output report_dir
    python clinical_variant_reporter.py --input patient.vcf --genes "BRCA1,BRCA2" --output report_dir
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent

sys.path.insert(0, str(SCRIPT_DIR))
from acmg_engine import (
    ACMG_SF_V32_GENES,
    ClassifiedVariant,
    VariantEvidence,
    classify_variant,
)

DISCLAIMER = (
    "ClawBio is a research and educational tool. It is not a medical device "
    "and does not provide clinical diagnoses. Consult a healthcare "
    "professional before making any medical decisions."
)

VEP_REST_URL = "https://rest.ensembl.org/vep/homo_sapiens/region"
VEP_BATCH_SIZE = 200
VEP_RATE_LIMIT_SECONDS = 0.07  # ~15 requests/second


# ---------------------------------------------------------------------------
# VCF parsing (lightweight, no pysam dependency)
# ---------------------------------------------------------------------------
@dataclass
class VcfRecord:
    chrom: str
    pos: int
    id: str
    ref: str
    alt: str
    qual: str
    filt: str
    info: dict[str, str]
    genotype: str | None = None


def parse_vcf(path: Path) -> list[VcfRecord]:
    """Parse a VCF file into a list of VcfRecord objects (one per ALT allele)."""
    records: list[VcfRecord] = []
    opener = open
    filepath = str(path)

    if filepath.endswith(".gz"):
        import gzip
        opener = gzip.open

    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 8:
                continue

            info: dict[str, str] = {}
            for item in fields[7].split(";"):
                if "=" in item:
                    k, v = item.split("=", 1)
                    info[k] = v
                else:
                    info[item] = "true"

            genotype = None
            if len(fields) >= 10:
                genotype = fields[9].split(":")[0]

            for alt_allele in fields[4].split(","):
                records.append(VcfRecord(
                    chrom=fields[0],
                    pos=int(fields[1]),
                    id=fields[2],
                    ref=fields[3],
                    alt=alt_allele.strip(),
                    qual=fields[5],
                    filt=fields[6],
                    info=info,
                    genotype=genotype,
                ))

    return records


# ---------------------------------------------------------------------------
# Evidence collection — demo mode (pre-cached, offline)
# ---------------------------------------------------------------------------
def load_demo_evidence_cache() -> dict[str, dict]:
    """Load pre-cached evidence for demo variants."""
    cache_path = SCRIPT_DIR / "example_data" / "demo_evidence_cache.json"
    with open(cache_path) as fh:
        return json.load(fh)


def build_evidence_from_cache(record: VcfRecord, cache: dict[str, dict]) -> VariantEvidence:
    """Build VariantEvidence from the pre-cached demo evidence."""
    key = f"{record.chrom}:{record.pos}:{record.ref}:{record.alt}"
    cached = cache.get(key, {})

    return VariantEvidence(
        chrom=record.chrom,
        pos=record.pos,
        ref=record.ref,
        alt=record.alt,
        rsid=record.id if record.id != "." else "",
        gene=cached.get("gene", record.info.get("GENE", "")),
        consequence=cached.get("consequence", ""),
        impact=cached.get("impact", ""),
        hgvsc=cached.get("hgvsc", ""),
        hgvsp=cached.get("hgvsp", ""),
        transcript=cached.get("transcript", ""),
        clinvar_significance=cached.get("clinvar_significance", ""),
        clinvar_review_stars=cached.get("clinvar_review_stars", 0),
        gnomad_af=cached.get("gnomad_af"),
        gnomad_af_popmax=cached.get("gnomad_af_popmax"),
        cadd_phred=cached.get("cadd_phred"),
        sift_prediction=cached.get("sift_prediction", ""),
        polyphen_prediction=cached.get("polyphen_prediction", ""),
        spliceai_max_delta=cached.get("spliceai_max_delta"),
        is_lof=cached.get("is_lof", False),
        is_missense=cached.get("is_missense", False),
        is_synonymous=cached.get("is_synonymous", False),
        is_inframe_indel=cached.get("is_inframe_indel", False),
    )


# ---------------------------------------------------------------------------
# Evidence collection — live mode (VEP REST API)
# ---------------------------------------------------------------------------
def _vcf_to_vep_region(record: VcfRecord) -> str:
    """Convert a VCF record to VEP region format: 'chr start end alleles strand'."""
    chrom = record.chrom.replace("chr", "")
    ref, alt = record.ref, record.alt

    if len(ref) == 1 and len(alt) == 1:
        return f"{chrom} {record.pos} {record.pos} {ref}/{alt} 1"

    if len(ref) > len(alt):
        start = record.pos + 1
        end = record.pos + len(ref) - 1
        deleted = ref[1:]
        return f"{chrom} {start} {end} {deleted}/- 1"

    if len(alt) > len(ref):
        start = record.pos + 1
        end = record.pos
        inserted = alt[1:]
        return f"{chrom} {start} {end} -/{inserted} 1"

    return f"{chrom} {record.pos} {record.pos + len(ref) - 1} {ref}/{alt} 1"


def _extract_evidence_from_vep(vep_result: dict, record: VcfRecord) -> VariantEvidence:
    """Extract VariantEvidence fields from a single VEP REST response entry."""
    most_severe = vep_result.get("most_severe_consequence", "")

    gene = ""
    consequence = most_severe
    impact = ""
    hgvsc = ""
    hgvsp = ""
    transcript_id = ""

    for tc in vep_result.get("transcript_consequences", []):
        if tc.get("consequence_terms") and most_severe in tc["consequence_terms"]:
            gene = tc.get("gene_symbol", gene)
            impact = tc.get("impact", impact)
            hgvsc = tc.get("hgvsc", hgvsc)
            hgvsp = tc.get("hgvsp", hgvsp)
            transcript_id = tc.get("transcript_id", transcript_id)
            break

    if not gene:
        for tc in vep_result.get("transcript_consequences", []):
            gene = tc.get("gene_symbol", "")
            impact = tc.get("impact", "")
            if gene:
                break

    clinvar_sig = ""
    clinvar_stars = 0
    gnomad_af = None
    gnomad_af_popmax = None

    for cv in vep_result.get("colocated_variants", []):
        if "clinvar" in cv.get("var_synonyms", {}) or cv.get("clin_sig"):
            sigs = cv.get("clin_sig", "")
            if sigs and not clinvar_sig:
                clinvar_sig = sigs
                clinvar_stars = cv.get("clin_sig_allele", {}).get("review_status_stars", 0)

        freq_data = cv.get("frequencies", {})
        if freq_data:
            for allele_freq in freq_data.values():
                af_val = allele_freq.get("gnomade", allele_freq.get("gnomad", None))
                if af_val is not None and (gnomad_af is None or af_val > gnomad_af):
                    gnomad_af = af_val

    sift = ""
    polyphen = ""
    cadd = None
    for tc in vep_result.get("transcript_consequences", []):
        if tc.get("sift_prediction") and not sift:
            sift = tc["sift_prediction"]
        if tc.get("polyphen_prediction") and not polyphen:
            polyphen = tc["polyphen_prediction"]
        if tc.get("cadd_phred") and cadd is None:
            cadd = tc["cadd_phred"]

    from acmg_engine import (
        INFRAME_CONSEQUENCES,
        LOF_CONSEQUENCES,
        MISSENSE_CONSEQUENCES,
        SYNONYMOUS_CONSEQUENCES,
    )

    return VariantEvidence(
        chrom=record.chrom,
        pos=record.pos,
        ref=record.ref,
        alt=record.alt,
        rsid=record.id if record.id != "." else vep_result.get("id", ""),
        gene=gene or record.info.get("GENE", ""),
        consequence=consequence,
        impact=impact,
        hgvsc=hgvsc,
        hgvsp=hgvsp,
        transcript=transcript_id,
        clinvar_significance=clinvar_sig,
        clinvar_review_stars=clinvar_stars,
        gnomad_af=gnomad_af,
        gnomad_af_popmax=gnomad_af_popmax,
        cadd_phred=cadd,
        sift_prediction=sift,
        polyphen_prediction=polyphen,
        is_lof=consequence in LOF_CONSEQUENCES,
        is_missense=consequence in MISSENSE_CONSEQUENCES,
        is_synonymous=consequence in SYNONYMOUS_CONSEQUENCES,
        is_inframe_indel=consequence in INFRAME_CONSEQUENCES,
    )


def annotate_variants_vep(
    records: list[VcfRecord],
    assembly: str = "GRCh38",
) -> list[VariantEvidence]:
    """Annotate variants via Ensembl VEP REST API and return evidence objects."""
    try:
        import requests
    except ImportError:
        print("ERROR: 'requests' package required for live mode. Install: pip install requests", file=sys.stderr)
        sys.exit(1)

    evidence_list: list[VariantEvidence] = []

    for batch_start in range(0, len(records), VEP_BATCH_SIZE):
        batch = records[batch_start:batch_start + VEP_BATCH_SIZE]
        regions = [_vcf_to_vep_region(r) for r in batch]

        try:
            resp = requests.post(
                VEP_REST_URL,
                json={"variants": regions},
                headers={"Content-Type": "application/json", "Accept": "application/json"},
                params={"assembly": assembly},
                timeout=60,
            )
            resp.raise_for_status()
            vep_results = resp.json()
        except Exception as exc:
            print(f"WARNING: VEP batch failed ({exc}). Marking {len(batch)} variants as unannotated.", file=sys.stderr)
            for rec in batch:
                evidence_list.append(VariantEvidence(
                    chrom=rec.chrom, pos=rec.pos, ref=rec.ref, alt=rec.alt,
                    gene=rec.info.get("GENE", ""),
                ))
            continue

        result_map: dict[str, dict] = {}
        for vr in vep_results:
            loc = vr.get("input", "")
            result_map[loc] = vr

        for rec, region in zip(batch, regions):
            vep_result = result_map.get(region, {})
            if vep_result:
                evidence_list.append(_extract_evidence_from_vep(vep_result, rec))
            else:
                evidence_list.append(VariantEvidence(
                    chrom=rec.chrom, pos=rec.pos, ref=rec.ref, alt=rec.alt,
                    gene=rec.info.get("GENE", ""),
                ))

        time.sleep(VEP_RATE_LIMIT_SECONDS)

    return evidence_list


# ---------------------------------------------------------------------------
# Classification pipeline
# ---------------------------------------------------------------------------
def run_classification(
    records: list[VcfRecord],
    demo: bool = False,
    gene_filter: set[str] | None = None,
    assembly: str = "GRCh38",
) -> list[ClassifiedVariant]:
    """Run the full ACMG classification pipeline on VCF records."""
    if demo:
        cache = load_demo_evidence_cache()
        evidence_list = [build_evidence_from_cache(r, cache) for r in records]
    else:
        evidence_list = annotate_variants_vep(records, assembly=assembly)

    if gene_filter:
        evidence_list = [e for e in evidence_list if e.gene in gene_filter]

    return [classify_variant(ev) for ev in evidence_list]


# ---------------------------------------------------------------------------
# Report generation
# ---------------------------------------------------------------------------
CLASS_ORDER = ["Pathogenic", "Likely Pathogenic", "Uncertain Significance", "Likely Benign", "Benign"]
CLASS_SHORT = {"Pathogenic": "P", "Likely Pathogenic": "LP", "Uncertain Significance": "VUS", "Likely Benign": "LB", "Benign": "B"}


def generate_report(
    classified: list[ClassifiedVariant],
    output_dir: Path,
    demo: bool = False,
    assembly: str = "GRCh38",
    input_path: str = "demo",
) -> None:
    """Generate all output files: report.md, result.json, tables/, figures/, reproducibility/."""
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "tables").mkdir(exist_ok=True)
    (output_dir / "figures").mkdir(exist_ok=True)
    (output_dir / "reproducibility").mkdir(exist_ok=True)

    timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")

    counts = {c: 0 for c in CLASS_ORDER}
    for cv in classified:
        counts[cv.classification] = counts.get(cv.classification, 0) + 1

    sf_variants = [cv for cv in classified if cv.is_secondary_finding]

    _write_markdown_report(classified, counts, sf_variants, output_dir, timestamp, demo, assembly, input_path)
    _write_classification_table(classified, output_dir)
    _write_secondary_findings_table(sf_variants, output_dir)
    _write_result_json(classified, counts, sf_variants, output_dir, timestamp, demo, assembly)
    _write_classification_figure(counts, output_dir)
    _write_reproducibility(output_dir, demo, assembly, input_path, timestamp)


def _write_markdown_report(
    classified: list[ClassifiedVariant],
    counts: dict[str, int],
    sf_variants: list[ClassifiedVariant],
    output_dir: Path,
    timestamp: str,
    demo: bool,
    assembly: str,
    input_path: str,
) -> None:
    lines: list[str] = []
    lines.append("# Clinical Variant Report — ACMG/AMP Classification")
    lines.append("")
    lines.append(f"**Generated**: {timestamp}")
    lines.append(f"**Input**: {input_path}")
    lines.append(f"**Assembly**: {assembly}")
    lines.append(f"**Total variants classified**: {len(classified)}")
    lines.append(f"**Mode**: {'Demo (pre-cached evidence)' if demo else 'Live (VEP REST API)'}")
    lines.append("")

    lines.append("## Classification Summary")
    lines.append("")
    lines.append("| Classification | Count |")
    lines.append("|----------------|-------|")
    for cls in CLASS_ORDER:
        lines.append(f"| {cls} | {counts.get(cls, 0)} |")
    lines.append("")

    actionable = [cv for cv in classified if cv.classification in ("Pathogenic", "Likely Pathogenic")]
    if actionable:
        lines.append("## Actionable Variants (Pathogenic / Likely Pathogenic)")
        lines.append("")
        for cv in actionable:
            ev = cv.evidence
            sf_tag = " **[SF]**" if cv.is_secondary_finding else ""
            lines.append(f"### {ev.gene} — {ev.hgvsp or ev.hgvsc or f'{ev.ref}>{ev.alt}'}{sf_tag}")
            lines.append("")
            lines.append(f"- **Classification**: {cv.classification}")
            lines.append(f"- **Position**: {ev.chrom}:{ev.pos}")
            lines.append(f"- **rsID**: {ev.rsid or 'N/A'}")
            lines.append(f"- **Consequence**: {ev.consequence}")
            lines.append(f"- **ClinVar**: {ev.clinvar_significance or 'N/A'} (stars: {ev.clinvar_review_stars})")
            lines.append(f"- **gnomAD AF**: {ev.gnomad_af if ev.gnomad_af is not None else 'N/A'}")
            lines.append(f"- **Evidence codes**: {cv.evidence_summary}")
            lines.append("")
            lines.append("| Criterion | Triggered | Strength | Direction | Source |")
            lines.append("|-----------|-----------|----------|-----------|--------|")
            for c in cv.criteria:
                tick = "Yes" if c.triggered else "No"
                lines.append(f"| {c.code} | {tick} | {c.strength} | {c.direction} | {c.source} |")
            lines.append("")

    vus_list = [cv for cv in classified if cv.classification == "Uncertain Significance"]
    if vus_list:
        lines.append("## Variants of Uncertain Significance (VUS)")
        lines.append("")
        lines.append("| Gene | Variant | Position | gnomAD AF | Evidence Codes | SF Gene |")
        lines.append("|------|---------|----------|-----------|----------------|---------|")
        for cv in vus_list:
            ev = cv.evidence
            af = f"{ev.gnomad_af:.6f}" if ev.gnomad_af is not None else "N/A"
            sf = "Yes" if cv.is_secondary_finding else "No"
            lines.append(f"| {ev.gene} | {ev.hgvsp or ev.hgvsc or f'{ev.ref}>{ev.alt}'} | {ev.chrom}:{ev.pos} | {af} | {cv.evidence_summary} | {sf} |")
        lines.append("")

    benign_list = [cv for cv in classified if cv.classification in ("Likely Benign", "Benign")]
    if benign_list:
        lines.append("## Benign / Likely Benign Variants")
        lines.append("")
        lines.append("| Gene | Variant | Classification | gnomAD AF | Evidence Codes |")
        lines.append("|------|---------|----------------|-----------|----------------|")
        for cv in benign_list:
            ev = cv.evidence
            af = f"{ev.gnomad_af:.4f}" if ev.gnomad_af is not None else "N/A"
            lines.append(f"| {ev.gene} | {ev.hgvsp or ev.hgvsc or f'{ev.ref}>{ev.alt}'} | {cv.classification} | {af} | {cv.evidence_summary} |")
        lines.append("")

    if sf_variants:
        lines.append("## ACMG SF v3.2 Secondary Findings Screening")
        lines.append("")
        lines.append(f"**{len(sf_variants)}** variant(s) found in ACMG SF v3.2 genes ({len(ACMG_SF_V32_GENES)} genes screened).")
        lines.append("")
        lines.append("| Gene | Variant | Classification | Evidence Codes |")
        lines.append("|------|---------|----------------|----------------|")
        for cv in sf_variants:
            ev = cv.evidence
            lines.append(f"| {ev.gene} | {ev.hgvsp or ev.hgvsc or f'{ev.ref}>{ev.alt}'} | {cv.classification} | {cv.evidence_summary} |")
        lines.append("")

    lines.append("## Methodology")
    lines.append("")
    lines.append("Variants were classified according to the ACMG/AMP 2015 standards and guidelines ")
    lines.append("(Richards et al., *Genet Med* 2015; PMID 25741868). Evidence was collected from Ensembl VEP ")
    lines.append("(consequence annotation, ClinVar, gnomAD colocated frequencies, SIFT, PolyPhen). ")
    lines.append("The ACMG combining rules were applied to assign one of five classifications: ")
    lines.append("Pathogenic, Likely Pathogenic, Uncertain Significance, Likely Benign, or Benign. ")
    lines.append("Secondary findings were screened against ACMG SF v3.2 (Miller et al., 2023; 81 genes).")
    lines.append("")
    lines.append("### Criteria Not Automatically Assessed")
    lines.append("")
    lines.append("The following ACMG criteria require additional data (family studies, functional assays, etc.) ")
    lines.append("and were not evaluated in this automated run:")
    lines.append("")
    lines.append("- **PS2/PM6**: De novo status (requires parental samples)")
    lines.append("- **PS3/BS3**: Functional studies (requires experimental data)")
    lines.append("- **PS4**: Case-control prevalence (requires cohort data)")
    lines.append("- **PM3**: In trans with pathogenic variant (requires phased data)")
    lines.append("- **PP1/BS4**: Family segregation (requires pedigree)")
    lines.append("- **PP2/BP1**: Gene-level missense constraint (planned)")
    lines.append("- **PP4**: Phenotype specificity (requires HPO terms)")
    lines.append("- **BP2/BP3/BP5/BS2**: Require additional contextual data")
    lines.append("")

    lines.append("## Limitations")
    lines.append("")
    lines.append("- Not all 28 ACMG/AMP criteria can be evaluated automatically; manual review is recommended for actionable variants")
    lines.append("- In silico predictor scores may not be available for all variants")
    lines.append("- ClinVar assertions reflect submitter interpretations and may change over time")
    lines.append("- gnomAD does not include all populations equally; AF may underestimate prevalence in underrepresented groups")
    lines.append("")

    lines.append("---")
    lines.append("")
    lines.append(f"*{DISCLAIMER}*")
    lines.append("")

    (output_dir / "report.md").write_text("\n".join(lines), encoding="utf-8")


def _write_classification_table(classified: list[ClassifiedVariant], output_dir: Path) -> None:
    table_path = output_dir / "tables" / "acmg_classifications.tsv"
    with open(table_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow([
            "chrom", "pos", "ref", "alt", "rsid", "gene", "consequence",
            "clinvar_significance", "clinvar_stars", "gnomad_af",
            "cadd_phred", "acmg_classification", "acmg_short",
            "evidence_codes", "is_secondary_finding",
        ])
        for cv in classified:
            ev = cv.evidence
            writer.writerow([
                ev.chrom, ev.pos, ev.ref, ev.alt, ev.rsid, ev.gene,
                ev.consequence, ev.clinvar_significance, ev.clinvar_review_stars,
                ev.gnomad_af if ev.gnomad_af is not None else "",
                ev.cadd_phred if ev.cadd_phred is not None else "",
                cv.classification, CLASS_SHORT.get(cv.classification, "?"),
                cv.evidence_summary, cv.is_secondary_finding,
            ])


def _write_secondary_findings_table(sf_variants: list[ClassifiedVariant], output_dir: Path) -> None:
    table_path = output_dir / "tables" / "secondary_findings.tsv"
    with open(table_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow([
            "gene", "chrom", "pos", "ref", "alt", "consequence",
            "acmg_classification", "evidence_codes",
        ])
        for cv in sf_variants:
            ev = cv.evidence
            writer.writerow([
                ev.gene, ev.chrom, ev.pos, ev.ref, ev.alt,
                ev.consequence, cv.classification, cv.evidence_summary,
            ])


def _write_result_json(
    classified: list[ClassifiedVariant],
    counts: dict[str, int],
    sf_variants: list[ClassifiedVariant],
    output_dir: Path,
    timestamp: str,
    demo: bool,
    assembly: str,
) -> None:
    result = {
        "tool": "ClawBio Clinical Variant Reporter",
        "version": "0.1.0",
        "framework": "ACMG/AMP 2015 (Richards et al., PMID 25741868)",
        "sf_list": "ACMG SF v3.2 (Miller et al., 2023)",
        "assembly": assembly,
        "timestamp": timestamp,
        "mode": "demo" if demo else "live",
        "total_variants": len(classified),
        "classification_counts": counts,
        "secondary_findings_count": len(sf_variants),
        "variants": [
            {
                "chrom": cv.evidence.chrom,
                "pos": cv.evidence.pos,
                "ref": cv.evidence.ref,
                "alt": cv.evidence.alt,
                "rsid": cv.evidence.rsid,
                "gene": cv.evidence.gene,
                "consequence": cv.evidence.consequence,
                "classification": cv.classification,
                "is_secondary_finding": cv.is_secondary_finding,
                "triggered_criteria": cv.triggered_codes,
                "evidence_summary": cv.evidence_summary,
            }
            for cv in classified
        ],
        "disclaimer": DISCLAIMER,
    }
    (output_dir / "result.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False), encoding="utf-8",
    )


def _write_classification_figure(counts: dict[str, int], output_dir: Path) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        labels = list(CLASS_SHORT.values())
        values = [counts.get(c, 0) for c in CLASS_ORDER]
        colours = ["#d32f2f", "#ff9800", "#9e9e9e", "#4caf50", "#2196f3"]

        fig, ax = plt.subplots(figsize=(8, 4))
        bars = ax.bar(labels, values, color=colours, edgecolor="white", linewidth=0.8)
        ax.set_ylabel("Variant count")
        ax.set_title("ACMG Classification Summary")
        ax.set_ylim(0, max(values) + 2 if values else 5)

        for bar, val in zip(bars, values):
            if val > 0:
                ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.3,
                        str(val), ha="center", va="bottom", fontweight="bold")

        fig.tight_layout()
        fig.savefig(output_dir / "figures" / "classification_summary.png", dpi=150)
        plt.close(fig)
    except ImportError:
        pass


def _write_reproducibility(
    output_dir: Path,
    demo: bool,
    assembly: str,
    input_path: str,
    timestamp: str,
) -> None:
    cmd = f"python {Path(__file__).name}"
    if demo:
        cmd += " --demo"
    else:
        cmd += f" --input {input_path}"
    cmd += f" --output {output_dir} --assembly {assembly}"
    (output_dir / "reproducibility" / "commands.sh").write_text(
        f"#!/usr/bin/env bash\n# Reproducibility command — generated {timestamp}\n{cmd}\n",
        encoding="utf-8",
    )

    db_versions = {
        "acmg_framework": "Richards et al. 2015 (PMID 25741868)",
        "sf_list": "ACMG SF v3.2 (Miller et al. 2023)",
        "sf_gene_count": len(ACMG_SF_V32_GENES),
        "annotation_backend": "demo_cache" if demo else "Ensembl VEP REST (GRCh38)",
        "generated": timestamp,
    }
    (output_dir / "reproducibility" / "database_versions.json").write_text(
        json.dumps(db_versions, indent=2), encoding="utf-8",
    )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main() -> None:
    parser = argparse.ArgumentParser(
        description="Clinical Variant Reporter — ACMG/AMP 2015 classification",
    )
    parser.add_argument("--input", type=str, help="Input VCF/BCF file path")
    parser.add_argument("--output", type=str, help="Output directory for reports")
    parser.add_argument("--demo", action="store_true", help="Run with built-in GIAB demo panel")
    parser.add_argument("--genes", type=str, help="Comma-separated gene list to filter (e.g. BRCA1,BRCA2)")
    parser.add_argument("--assembly", type=str, default="GRCh38", choices=["GRCh37", "GRCh38"])

    args = parser.parse_args()

    if not args.demo and not args.input:
        parser.error("Provide --input <vcf> or --demo")
    if not args.output:
        parser.error("Provide --output <directory>")

    output_dir = Path(args.output)
    if output_dir.exists() and any(output_dir.iterdir()):
        print(f"WARNING: Output directory '{output_dir}' is not empty — files may be overwritten.", file=sys.stderr)

    if args.demo:
        vcf_path = SCRIPT_DIR / "example_data" / "giab_acmg_panel.vcf"
        input_label = str(vcf_path)
    else:
        vcf_path = Path(args.input)
        input_label = str(vcf_path)

    if not vcf_path.exists():
        print(f"ERROR: Input file not found: {vcf_path}", file=sys.stderr)
        sys.exit(1)

    print(f"[CVR] Parsing VCF: {vcf_path}")
    records = parse_vcf(vcf_path)
    print(f"[CVR] Found {len(records)} variant record(s)")

    gene_filter = None
    if args.genes:
        gene_filter = {g.strip() for g in args.genes.split(",")}
        print(f"[CVR] Filtering to genes: {', '.join(sorted(gene_filter))}")

    print(f"[CVR] Running ACMG classification ({'demo mode' if args.demo else 'live mode'})...")
    classified = run_classification(
        records, demo=args.demo, gene_filter=gene_filter, assembly=args.assembly,
    )

    print(f"[CVR] Generating report in: {output_dir}")
    generate_report(classified, output_dir, demo=args.demo, assembly=args.assembly, input_path=input_label)

    counts = {}
    for cv in classified:
        counts[cv.classification] = counts.get(cv.classification, 0) + 1
    sf_count = sum(1 for cv in classified if cv.is_secondary_finding)

    print(f"[CVR] Classification complete:")
    for cls in CLASS_ORDER:
        print(f"       {cls}: {counts.get(cls, 0)}")
    print(f"       Secondary findings (SF v3.2): {sf_count}")
    print(f"[CVR] Report written to: {output_dir / 'report.md'}")


if __name__ == "__main__":
    main()
