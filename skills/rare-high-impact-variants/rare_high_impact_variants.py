#!/usr/bin/env python3
"""Rare High Impact Variants - count rare, high-impact loss-of-function variants carried in a VCF.

A carried variant is counted as "rare high-impact" when:
  - the genotype carries the ALT allele (heterozygous or homozygous), AND
  - the molecular consequence is loss-of-function (nonsense / stop-gained,
    frameshift, splice donor/acceptor, start-lost, stop-lost), AND
  - the population allele frequency is below the rarity threshold (default 1 per
    cent), or the variant is absent from the reference frequency panels.

The input VCF must be annotated with molecular consequence (ClinVar `MC`, or a
VEP/SnpEff consequence) and population frequency (`AF_TGP`, `AF_EXAC`, `AF_ESP`,
or `gnomAD_AF`). Scope: this counts high-impact variants that are catalogued and
annotated; genome-wide novel loss-of-function calling needs a consequence
predictor (VEP / SnpEff / bcftools csq) plus gnomAD and is out of scope for v0.
"""

import argparse
import csv
import gzip
import json
import shlex
import sys
from pathlib import Path

SKILL_DIR = Path(__file__).resolve().parent
_PROJECT_ROOT = SKILL_DIR.parents[1]
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from clawbio.common.reproducibility import (  # noqa: E402
    write_checksums,
    write_commands_sh,
    write_environment_yml,
)
VERSION = "0.1.0"
DISCLAIMER = ("ClawBio is a research and educational tool. It is not a medical device "
             "and does not provide clinical diagnoses. Consult a healthcare professional "
             "before making any medical decisions.")

# Loss-of-function / high-impact consequence substrings (matched case-insensitively
# against the consequence string, covering ClinVar MC and VEP/SnpEff naming).
HIGH_IMPACT_TERMS = (
    "nonsense", "stop_gained", "frameshift", "splice_donor", "splice_acceptor",
    "start_lost", "initiator_codon", "stop_lost",
)
AF_KEYS = ("AF_TGP", "AF_EXAC", "AF_ESP", "gnomAD_AF", "AF_grpmax", "gnomad_af")
DEFAULT_MAX_AF = 0.01
ULTRA_RARE_AF = 0.001


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", type=Path, dest="input_file", help="Annotated VCF (.vcf or .vcf.gz)")
    parser.add_argument("--output", type=Path, help="Output directory")
    parser.add_argument("--max-af", type=float, default=DEFAULT_MAX_AF,
                        help=f"Rarity threshold on population AF (default {DEFAULT_MAX_AF})")
    parser.add_argument("--demo", action="store_true", help="Run with synthetic demo data")
    return parser.parse_args()


def _open(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def _parse_info(info: str) -> dict:
    out = {}
    if info == ".":
        return out
    for field in info.split(";"):
        if "=" in field:
            k, v = field.split("=", 1)
            out[k] = v
        else:
            out[field] = True
    return out


def _carries_alt(gt: str) -> bool:
    alleles = gt.replace("|", "/").split("/")
    return any(a not in ("0", ".", "") for a in alleles)


def _zygosity(gt: str) -> str:
    a = gt.replace("|", "/").split("/")
    if len(a) == 2 and a[0] == a[1] and a[0] not in ("0", "."):
        return "hom"
    return "het"


def _pop_af(info: dict):
    afs = []
    for k in AF_KEYS:
        v = info.get(k)
        if v not in (None, ".", "", True):
            try:
                afs.append(float(v))
            except ValueError:
                pass
    return max(afs) if afs else None


def _is_high_impact(consequence: str) -> bool:
    c = consequence.lower()
    return any(term in c for term in HIGH_IMPACT_TERMS)


def validate_input(input_path: Path) -> dict:
    """Parse an annotated VCF. Returns {'variants': [...], 'source': str}."""
    if not input_path.exists():
        print(f"Error: input file not found: {input_path}", file=sys.stderr)
        sys.exit(1)
    variants = []
    sample_idx = 9  # first sample column in a standard VCF
    saw_header = False
    with _open(input_path) as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if not line:
                continue
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                saw_header = True
                continue
            cols = line.split("\t")
            if len(cols) < 8:
                continue  # not a valid VCF data line
            info = _parse_info(cols[7])
            gt = "1/1"  # default: assume carried if no genotype column present
            if len(cols) > sample_idx:
                gt = cols[sample_idx].split(":")[0]
            consequence = info.get("MC", "") or info.get("Consequence", "") or info.get("ANN", "")
            variants.append({
                "chrom": cols[0], "pos": cols[1], "ref": cols[3], "alt": cols[4],
                "gt": gt, "consequence": consequence,
                "gene": (info.get("GENEINFO", "") or info.get("SYMBOL", "")).split(":")[0].split("|")[0],
                "clnsig": info.get("CLNSIG", ""),
                "pop_af": _pop_af(info),
                "rsid": cols[2] if cols[2] not in (".", "") else "",
            })
    return {"variants": variants, "source": str(input_path), "saw_header": saw_header}


def run_analysis(data: dict, max_af: float = DEFAULT_MAX_AF) -> dict:
    """Count carried rare high-impact (loss-of-function) variants.

    A variant is only called "rare" when it has a DOCUMENTED population frequency
    below the threshold. Variants with no frequency in the source are reported in
    a separate "frequency unknown" bucket: absence of a frequency annotation is
    not evidence of rarity (many are common loss-of-function polymorphisms, often
    ClinVar-benign), so they are not counted toward the headline number.
    """
    variants = data.get("variants", [])
    carried = 0
    high_impact_carried = 0
    findings = []            # documented-rare, the headline set
    by_consequence = {}
    by_rarity = {"ultra_rare": 0, "rare": 0}
    freq_unknown = []        # high-impact, no frequency data: cannot confirm rare
    common = 0               # high-impact but common (documented AF >= threshold)

    for v in variants:
        if not _carries_alt(v["gt"]):
            continue
        carried += 1
        if not _is_high_impact(v["consequence"]):
            continue
        high_impact_carried += 1
        af = v["pop_af"]
        cons = next((t for t in HIGH_IMPACT_TERMS if t in v["consequence"].lower()), "high_impact")
        rec = {
            "gene": v["gene"], "rsid": v["rsid"],
            "locus": f"{v['chrom']}:{v['pos']} {v['ref']}>{v['alt']}",
            "consequence": cons, "zygosity": _zygosity(v["gt"]),
            "population_af": af, "clnsig": v["clnsig"],
        }
        if af is None:
            rec["rarity"] = "frequency_unknown"
            freq_unknown.append(rec)
            continue
        if af >= max_af:
            common += 1
            continue
        rec["rarity"] = "ultra_rare" if af < ULTRA_RARE_AF else "rare"
        by_rarity[rec["rarity"]] += 1
        by_consequence[cons] = by_consequence.get(cons, 0) + 1
        findings.append(rec)

    findings.sort(key=lambda f: (f["population_af"], f["gene"]))
    return {
        "skill": "rare-high-impact-variants",
        "version": VERSION,
        "source": data.get("source", "unknown"),
        "rarity_threshold": max_af,
        "variants_processed": len(variants),
        "carried_variants": carried,
        "high_impact_carried": high_impact_carried,
        "rare_high_impact_count": len(findings),
        "high_impact_common": common,
        "high_impact_frequency_unknown": len(freq_unknown),
        "by_rarity": by_rarity,
        "by_consequence": by_consequence,
        "findings": findings,
        "frequency_unknown_genes": sorted({r["gene"] for r in freq_unknown if r["gene"]}),
        "scope_note": ("Counts high-impact (loss-of-function) variants annotated with molecular "
                       "consequence and population frequency in the input VCF. 'Rare' requires a "
                       "documented frequency below the threshold; variants with no frequency are "
                       "reported separately and not called rare. Genome-wide novel LoF calling "
                       "(VEP / SnpEff / bcftools csq) and a complete frequency reference (gnomAD) "
                       "are out of scope for v0.1.0."),
        "disclaimer": DISCLAIMER,
    }


def write_report(result: dict, output_dir: Path) -> None:
    """Write report.md and result.json to output_dir."""
    output_dir.mkdir(parents=True, exist_ok=True)
    with open(output_dir / "result.json", "w") as f:
        json.dump(result, f, indent=2)

    n = result["rare_high_impact_count"]
    r = result["by_rarity"]
    unknown = result["high_impact_frequency_unknown"]
    report = [
        "# Rare High-Impact Variants Report",
        "",
        f"**Input**: {result.get('source', 'unknown')}",
        f"**Rarity threshold**: population AF < {result['rarity_threshold']}",
        "",
        f"## {n} rare high-impact variant{'s' if n != 1 else ''} carried",
        "",
        f"Of {result['carried_variants']} carried, annotated variants, "
        f"{result['high_impact_carried']} are high-impact (loss-of-function). Of those:",
        "",
        f"- **{n} rare** with documented population frequency below {result['rarity_threshold']} "
        f"(ultra-rare AF < {ULTRA_RARE_AF}: {r['ultra_rare']}; rare: {r['rare']})",
        f"- {result['high_impact_common']} common (documented AF at or above the threshold)",
        f"- {unknown} with no population-frequency data, so they cannot be confirmed rare "
        f"(absence of a frequency is not evidence of rarity; many are common LoF polymorphisms)",
        "",
    ]
    if result["findings"]:
        report += ["## Variants", "",
                   "| Gene | Locus | Consequence | Zygosity | Population AF | ClinVar |",
                   "|------|-------|-------------|----------|---------------|---------|"]
        for f_ in result["findings"]:
            af = "absent" if f_["population_af"] is None else f"{f_['population_af']:.2g}"
            report.append(f"| {f_['gene'] or '-'} | {f_['locus']} | {f_['consequence']} | "
                          f"{f_['zygosity']} | {af} | {f_['clnsig'] or '-'} |")
        report.append("")
    report += ["## Scope", "", result["scope_note"], "", f"*{result['disclaimer']}*", ""]
    with open(output_dir / "report.md", "w") as f:
        f.write("\n".join(report))

    # tables/results.csv (ClawBio output convention)
    tables = output_dir / "tables"
    tables.mkdir(exist_ok=True)
    with open(tables / "results.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["gene", "locus", "consequence", "zygosity", "population_af", "rarity", "clnsig"])
        for x in result["findings"]:
            w.writerow([x["gene"], x["locus"], x["consequence"], x["zygosity"],
                        "" if x["population_af"] is None else x["population_af"],
                        x["rarity"], x["clnsig"]])

    print(f"Report written to {output_dir / 'report.md'}")
    print(f"Results written to {output_dir / 'result.json'}")


def write_reproducibility(output_dir: Path, args: argparse.Namespace) -> None:
    """Write reproducibility/ (commands.sh, environment.yml, checksums.sha256)."""
    cmd = ["python", "skills/rare-high-impact-variants/rare_high_impact_variants.py"]
    if args.demo:
        cmd.append("--demo")
    else:
        cmd += ["--input", str(args.input_file)]
    cmd += ["--output", str(output_dir), "--max-af", str(args.max_af)]
    write_commands_sh(output_dir, shlex.join(cmd))
    env_path = write_environment_yml(output_dir, "clawbio-rare-high-impact-variants", [])
    outputs = [output_dir / "report.md", output_dir / "result.json",
               output_dir / "tables" / "results.csv", env_path]
    write_checksums(outputs, output_dir, anchor=output_dir)


def run_demo(output_dir: Path, max_af: float = DEFAULT_MAX_AF) -> None:
    demo_input = SKILL_DIR / "demo_input.txt"
    if not demo_input.exists():
        print("Error: demo data not found", file=sys.stderr)
        sys.exit(1)
    data = validate_input(demo_input)
    write_report(run_analysis(data, max_af), output_dir)


def main():
    args = parse_args()
    if args.demo:
        output = args.output or Path("/tmp") / "rare_high_impact_variants" / "demo"
        run_demo(output, args.max_af)
    elif args.input_file:
        data = validate_input(args.input_file)
        result = run_analysis(data, args.max_af)
        output = args.output or args.input_file.parent / "output"
        write_report(result, output)
    else:
        print("Error: provide --input <file> or --demo", file=sys.stderr)
        sys.exit(1)
    write_reproducibility(output, args)


if __name__ == "__main__":
    main()
