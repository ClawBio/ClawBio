#!/usr/bin/env python3
"""
ClawBio — VCF Annotator Skill
Annotates VCF variants using ClinVar, gnomAD, and Ensembl VEP REST APIs.
Ranks variants by predicted impact and generates a reproducibility bundle.

Usage:
    python vcf_annotator.py --input variants.vcf --output report/
    python vcf_annotator.py --demo --output /tmp/demo
"""

import argparse
import csv
import hashlib
import json
import sys
import time
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

try:
    import cyvcf2
    HAS_CYVCF2 = True
except ImportError:
    HAS_CYVCF2 = False

try:
    import pandas as pd
    HAS_PANDAS = True
except ImportError:
    HAS_PANDAS = False

# ── Constants ──────────────────────────────────────────────────────────────────
ENSEMBL_VEP_URL  = "https://rest.ensembl.org/vep/human/hgvs"
CLINVAR_API_URL  = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
GNOMAD_API_URL   = "https://gnomad.broadinstitute.org/api"
TOOL_NAME        = "ClawBio-VCFAnnotator"
TOOL_EMAIL       = "mc@manuelcorpas.com"

IMPACT_RANK = {"HIGH": 1, "MODERATE": 2, "LOW": 3, "MODIFIER": 4, "UNKNOWN": 5}

DISCLAIMER = (
    "\n---\n*ClawBio VCF Annotator is a research tool. "
    "Not a clinical diagnostic device. "
    "Always consult a qualified geneticist for clinical decisions.*\n"
)

# ── Demo VCF data ──────────────────────────────────────────────────────────────
DEMO_VCF_CONTENT = """##fileformat=VCFv4.2
##reference=GRCh38
##FILTER=<ID=PASS,Description="All filters passed">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
17\t43044295\trs80357382\tG\tA\t.\tPASS\t.
13\t32316461\trs80359550\tC\tT\t.\tPASS\t.
7\t117548628\trs113993960\tCTTT\tC\t.\tPASS\t.
19\t11089199\trs429358\tT\tC\t.\tPASS\t.
1\t69515\trs1801133\tG\tA\t.\tPASS\t.
"""

DEMO_ANNOTATIONS = [
    {
        "chrom": "17",
        "pos": "43044295",
        "id": "rs80357382",
        "ref": "G",
        "alt": "A",
        "gene": "BRCA1",
        "consequence": "missense_variant",
        "impact": "HIGH",
        "clinvar_significance": "Pathogenic",
        "clinvar_condition": "Hereditary breast and ovarian cancer syndrome",
        "gnomad_af": 0.000008,
        "gnomad_af_afr": 0.0,
        "gnomad_af_eur": 0.000012,
        "hgvs": "NM_007294.4:c.5266dupC",
        "sift": "deleterious",
        "polyphen": "probably_damaging",
    },
    {
        "chrom": "13",
        "pos": "32316461",
        "id": "rs80359550",
        "ref": "C",
        "alt": "T",
        "gene": "BRCA2",
        "consequence": "stop_gained",
        "impact": "HIGH",
        "clinvar_significance": "Pathogenic",
        "clinvar_condition": "Familial cancer of breast",
        "gnomad_af": 0.000004,
        "gnomad_af_afr": 0.0,
        "gnomad_af_eur": 0.000006,
        "hgvs": "NM_000059.4:c.9976A>T",
        "sift": "deleterious",
        "polyphen": "probably_damaging",
    },
    {
        "chrom": "7",
        "pos": "117548628",
        "id": "rs113993960",
        "ref": "CTTT",
        "alt": "C",
        "gene": "CFTR",
        "consequence": "frameshift_variant",
        "impact": "HIGH",
        "clinvar_significance": "Pathogenic",
        "clinvar_condition": "Cystic fibrosis",
        "gnomad_af": 0.021,
        "gnomad_af_afr": 0.0003,
        "gnomad_af_eur": 0.033,
        "hgvs": "NM_000492.4:c.1521_1523delCTT",
        "sift": "deleterious",
        "polyphen": "N/A",
    },
    {
        "chrom": "19",
        "pos": "11089199",
        "id": "rs429358",
        "ref": "T",
        "alt": "C",
        "gene": "APOE",
        "consequence": "missense_variant",
        "impact": "MODERATE",
        "clinvar_significance": "Risk factor",
        "clinvar_condition": "Alzheimer disease",
        "gnomad_af": 0.147,
        "gnomad_af_afr": 0.232,
        "gnomad_af_eur": 0.143,
        "hgvs": "NM_000041.4:c.388T>C",
        "sift": "tolerated",
        "polyphen": "benign",
    },
    {
        "chrom": "1",
        "pos": "69515",
        "id": "rs1801133",
        "ref": "G",
        "alt": "A",
        "gene": "MTHFR",
        "consequence": "missense_variant",
        "impact": "MODERATE",
        "clinvar_significance": "Benign/Likely benign",
        "clinvar_condition": "Homocystinuria",
        "gnomad_af": 0.312,
        "gnomad_af_afr": 0.198,
        "gnomad_af_eur": 0.367,
        "hgvs": "NM_005957.5:c.665C>T",
        "sift": "tolerated",
        "polyphen": "benign",
    },
]


# ── VCF parser ─────────────────────────────────────────────────────────────────

def parse_vcf(vcf_path: Path) -> list[dict]:
    """Parse a VCF file. Uses cyvcf2 if available, falls back to stdlib."""
    if HAS_CYVCF2:
        variants = []
        for v in cyvcf2.VCF(str(vcf_path)):
            variants.append({
                "chrom": str(v.CHROM).replace("chr", ""),
                "pos":   str(v.POS),
                "id":    v.ID or "",
                "ref":   v.REF,
                "alt":   v.ALT[0] if v.ALT else ".",
            })
        return variants

    # Stdlib fallback
    variants = []
    with open(vcf_path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 5:
                continue
            variants.append({
                "chrom": parts[0].replace("chr", ""),
                "pos":   parts[1],
                "id":    parts[2] if parts[2] != "." else "",
                "ref":   parts[3],
                "alt":   parts[4],
            })
    return variants


# ── Ensembl VEP ────────────────────────────────────────────────────────────────

def annotate_vep(variant: dict) -> dict:
    """Query Ensembl VEP REST API for a single variant."""
    chrom = variant["chrom"]
    pos   = variant["pos"]
    ref   = variant["ref"]
    alt   = variant["alt"]

    # Build HGVS notation for SNVs
    hgvs = f"{chrom}:g.{pos}{ref}>{alt}"
    url  = f"https://rest.ensembl.org/vep/human/hgvs/{urllib.parse.quote(hgvs)}"
    headers_str = "application/json"

    result = {
        "gene": "Unknown",
        "consequence": "unknown",
        "impact": "UNKNOWN",
        "hgvs": hgvs,
        "sift": "N/A",
        "polyphen": "N/A",
    }

    try:
        req = urllib.request.Request(
            url,
            headers={"Content-Type": "application/json", "Accept": "application/json"},
        )
        with urllib.request.urlopen(req, timeout=10) as resp:
            data = json.loads(resp.read().decode())
        if data and isinstance(data, list):
            hit = data[0]
            tcs = hit.get("transcript_consequences", [{}])
            tc  = tcs[0] if tcs else {}
            result["gene"]        = tc.get("gene_symbol", "Unknown")
            result["consequence"] = tc.get("consequence_terms", ["unknown"])[0]
            result["impact"]      = tc.get("impact", "UNKNOWN")
            result["hgvs"]        = hit.get("id", hgvs)
            result["sift"]        = tc.get("sift_prediction", "N/A")
            result["polyphen"]    = tc.get("polyphen_prediction", "N/A")
        time.sleep(0.1)  # Ensembl rate limit
    except Exception as exc:
        print(f"  [warn] VEP failed for {hgvs}: {exc}", file=sys.stderr)

    return result


# ── ClinVar ────────────────────────────────────────────────────────────────────

def lookup_clinvar(rsid: str) -> dict:
    """Look up ClinVar significance for an rsID."""
    result = {"clinvar_significance": "Not found", "clinvar_condition": ""}
    if not rsid:
        return result

    params = urllib.parse.urlencode({
        "db":      "clinvar",
        "term":    rsid,
        "retmax":  1,
        "retmode": "json",
        "tool":    TOOL_NAME,
        "email":   TOOL_EMAIL,
    })
    try:
        url = f"{CLINVAR_API_URL}?{params}"
        with urllib.request.urlopen(url, timeout=10) as resp:
            data = json.loads(resp.read().decode())
        ids = data.get("esearchresult", {}).get("idlist", [])
        if ids:
            result["clinvar_significance"] = "Found in ClinVar"
        time.sleep(0.34)
    except Exception as exc:
        print(f"  [warn] ClinVar lookup failed for {rsid}: {exc}", file=sys.stderr)

    return result


# ── gnomAD ─────────────────────────────────────────────────────────────────────

def lookup_gnomad(chrom: str, pos: str, ref: str, alt: str) -> dict:
    """Query gnomAD GraphQL API for allele frequency."""
    result = {"gnomad_af": None, "gnomad_af_afr": None, "gnomad_af_eur": None}
    query = """
    query($variantId: String!) {
      variant(variantId: $variantId, dataset: gnomad_r4) {
        genome { af
          populations { id af }
        }
      }
    }
    """
    variant_id = f"{chrom}-{pos}-{ref}-{alt}"
    payload = json.dumps({
        "query": query,
        "variables": {"variantId": variant_id}
    }).encode()

    try:
        req = urllib.request.Request(
            GNOMAD_API_URL,
            data=payload,
            headers={"Content-Type": "application/json"},
        )
        with urllib.request.urlopen(req, timeout=10) as resp:
            data = json.loads(resp.read().decode())
        genome = (data.get("data", {}).get("variant") or {}).get("genome") or {}
        result["gnomad_af"] = genome.get("af")
        pops = genome.get("populations", [])
        for pop in pops:
            if pop["id"] == "afr":
                result["gnomad_af_afr"] = pop["af"]
            elif pop["id"] == "nfe":
                result["gnomad_af_eur"] = pop["af"]
        time.sleep(0.2)
    except Exception as exc:
        print(f"  [warn] gnomAD lookup failed for {variant_id}: {exc}", file=sys.stderr)

    return result


# ── Full annotation pipeline ───────────────────────────────────────────────────

def annotate_variants(variants: list[dict]) -> list[dict]:
    """Run full annotation pipeline for a list of variants."""
    annotated = []
    for i, v in enumerate(variants, 1):
        print(f"  Annotating variant {i}/{len(variants)}: "
              f"{v['chrom']}:{v['pos']} {v['ref']}>{v['alt']}")

        ann = {**v}
        vep    = annotate_vep(v)
        clinvar = lookup_clinvar(v.get("id", ""))
        gnomad  = lookup_gnomad(v["chrom"], v["pos"], v["ref"], v["alt"])

        ann.update(vep)
        ann.update(clinvar)
        ann.update(gnomad)
        annotated.append(ann)

    # Sort by impact rank
    annotated.sort(key=lambda x: IMPACT_RANK.get(x.get("impact", "UNKNOWN"), 5))
    return annotated


# ── Report generation ──────────────────────────────────────────────────────────

def generate_report(
    variants: list[dict],
    output_dir: Path,
    input_label: str = "input.vcf",
) -> None:
    """Write markdown report, JSON, CSV, and reproducibility bundle."""
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "tables").mkdir(exist_ok=True)
    (output_dir / "reproducibility").mkdir(exist_ok=True)

    now   = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    total = len(variants)
    high  = sum(1 for v in variants if v.get("impact") == "HIGH")
    mod   = sum(1 for v in variants if v.get("impact") == "MODERATE")
    low   = sum(1 for v in variants if v.get("impact") == "LOW")
    path_count = sum(
        1 for v in variants
        if "pathogenic" in v.get("clinvar_significance", "").lower()
    )

    # ── report.md ──
    lines = [
        "# 🦖 ClawBio VCF Annotator Report",
        "",
        f"**Input**: {input_label}  ",
        f"**Date**: {now}  ",
        f"**Total variants**: {total}  ",
        f"**HIGH impact**: {high} | **MODERATE**: {mod} | **LOW**: {low}  ",
        f"**ClinVar Pathogenic/Likely Pathogenic**: {path_count}",
        "",
        "---",
        "",
        "## Summary",
        "",
        f"Annotated {total} variants against Ensembl VEP, ClinVar, and gnomAD. "
        f"{high} variant(s) carry HIGH impact consequences. "
        f"{path_count} variant(s) are classified as Pathogenic or Likely Pathogenic "
        f"in ClinVar. Variants are ranked by predicted functional impact below.",
        "",
        "---",
        "",
        "## Variant Table",
        "",
        "| # | Gene | Variant | Consequence | Impact | ClinVar | gnomAD AF |",
        "|---|------|---------|-------------|--------|---------|-----------|",
    ]

    for i, v in enumerate(variants, 1):
        af = v.get("gnomad_af")
        af_str = f"{af:.6f}" if af is not None else "N/A"
        sig = v.get("clinvar_significance", "Not found")
        lines.append(
            f"| {i} | {v.get('gene','?')} | "
            f"{v['chrom']}:{v['pos']} {v['ref']}>{v['alt']} | "
            f"{v.get('consequence','?')} | "
            f"**{v.get('impact','?')}** | "
            f"{sig} | "
            f"{af_str} |"
        )

    lines += ["", "---", "", "## Detailed Annotations", ""]

    for i, v in enumerate(variants, 1):
        af     = v.get("gnomad_af")
        af_afr = v.get("gnomad_af_afr")
        af_eur = v.get("gnomad_af_eur")
        lines += [
            f"### {i}. {v.get('gene','Unknown')} — "
            f"{v['chrom']}:{v['pos']} {v['ref']}>{v['alt']}",
            "",
            f"| Field | Value |",
            f"|-------|-------|",
            f"| **rsID** | {v.get('id') or 'N/A'} |",
            f"| **HGVS** | {v.get('hgvs','N/A')} |",
            f"| **Consequence** | {v.get('consequence','N/A')} |",
            f"| **Impact** | {v.get('impact','N/A')} |",
            f"| **SIFT** | {v.get('sift','N/A')} |",
            f"| **PolyPhen** | {v.get('polyphen','N/A')} |",
            f"| **ClinVar** | {v.get('clinvar_significance','Not found')} |",
            f"| **ClinVar Condition** | {v.get('clinvar_condition','N/A')} |",
            f"| **gnomAD AF (global)** | "
            f"{'%.6f' % af if af is not None else 'N/A'} |",
            f"| **gnomAD AF (AFR)** | "
            f"{'%.6f' % af_afr if af_afr is not None else 'N/A'} |",
            f"| **gnomAD AF (EUR)** | "
            f"{'%.6f' % af_eur if af_eur is not None else 'N/A'} |",
            "",
        ]

    lines += ["---", DISCLAIMER]

    (output_dir / "report.md").write_text("\n".join(lines), encoding="utf-8")

    # ── results.json ──
    (output_dir / "results.json").write_text(
        json.dumps(variants, indent=2), encoding="utf-8"
    )

    # ── tables/variants.csv ──
    csv_fields = [
        "chrom","pos","id","ref","alt","gene","consequence","impact",
        "hgvs","sift","polyphen","clinvar_significance","clinvar_condition",
        "gnomad_af","gnomad_af_afr","gnomad_af_eur",
    ]
    csv_path = output_dir / "tables" / "variants.csv"
    if HAS_PANDAS:
        df = pd.DataFrame(variants)
        # Keep only known columns that exist
        existing = [c for c in csv_fields if c in df.columns]
        df[existing].to_csv(csv_path, index=False)
    else:
        with open(csv_path, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=csv_fields, extrasaction="ignore")
            w.writeheader()
            w.writerows(variants)

    # ── reproducibility bundle ──
    commands = f"""#!/usr/bin/env bash
# ClawBio VCF Annotator — reproducibility bundle
# Generated: {now}
# Input: {input_label}

python skills/vcf-annotator/vcf_annotator.py \\
    --input "{input_label}" \\
    --output "{output_dir}"
"""
    (output_dir / "reproducibility" / "commands.sh").write_text(commands)

    env_yml = f"""name: clawbio-vcf-annotator
channels:
  - defaults
dependencies:
  - python=3.11
  - pip
# Generated: {now}
# No extra packages required — uses only Python stdlib
"""
    (output_dir / "reproducibility" / "environment.yml").write_text(env_yml)

    # ── SHA-256 checksums ──
    checksums = {}
    for fpath in output_dir.rglob("*"):
        if fpath.is_file() and "checksums" not in fpath.name:
            h = hashlib.sha256(fpath.read_bytes()).hexdigest()
            checksums[str(fpath.relative_to(output_dir))] = h
    checksum_lines = "\n".join(f"{h}  {f}" for f, h in checksums.items())
    (output_dir / "reproducibility" / "checksums.sha256").write_text(checksum_lines)

    print(f"\n✅ Report: {output_dir}/report.md")
    print(f"📊 {total} variants | HIGH: {high} | Pathogenic: {path_count}")
    print(f"🔒 Checksums: {output_dir}/reproducibility/checksums.sha256")


# ── Main ───────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="ClawBio VCF Annotator — variant annotation with ClinVar, gnomAD, VEP"
    )
    parser.add_argument("--input",  type=str, help="Input VCF file path")
    parser.add_argument("--output", type=str, default="vcf_report", help="Output directory")
    parser.add_argument("--demo",   action="store_true", help="Run with built-in demo data")
    args = parser.parse_args()

    if not args.demo and not args.input:
        parser.error("Provide --input or use --demo")

    output_dir = Path(args.output)

    if args.demo:
        print("🦖 ClawBio VCF Annotator — DEMO MODE")
        print("   Using built-in demo variants (no network calls)")
        variants    = DEMO_ANNOTATIONS
        input_label = "demo_variants.vcf"
    else:
        print("🦖 ClawBio VCF Annotator")
        vcf_path = Path(args.input)
        if not vcf_path.exists():
            print(f"ERROR: File not found: {vcf_path}", file=sys.stderr)
            sys.exit(1)
        print(f"   Input : {vcf_path}")
        print(f"   Output: {output_dir}")
        print()

        print("📂 Parsing VCF …")
        raw_variants = parse_vcf(vcf_path)
        print(f"   Found {len(raw_variants)} variants")

        print("🔬 Annotating with VEP + ClinVar + gnomAD …")
        variants    = annotate_variants(raw_variants)
        input_label = str(vcf_path)

    print("\n📝 Generating report …")
    generate_report(variants, output_dir, input_label)


if __name__ == "__main__":
    main()
