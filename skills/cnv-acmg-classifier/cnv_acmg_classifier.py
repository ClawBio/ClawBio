#!/usr/bin/env python3
"""CNV ACMG Classifier.

Deterministic copy-number variant (CNV / structural variant) classification
using the ClinGen / ACMG 2019 technical standard point framework
(Riggs et al., Genet Med 2020; PMID 31690835).

The skill scores the computationally-derivable evidence categories from CNV
coordinates, a dosage-sensitivity map, and a protein-coding gene model, then
maps the summed score to the five-tier ACMG classification. Curated case-level
evidence (Section 4) and inheritance (Section 5) are supplied per-CNV by the
analyst; they cannot be auto-mined and are never fabricated.

Local-first, no network, stdlib-only.
"""

from __future__ import annotations

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

SCRIPT_REL = "skills/cnv-acmg-classifier/cnv_acmg_classifier.py"
DISCLAIMER = (
    "ClawBio is a research and educational tool. It is not a medical device "
    "and does not provide clinical diagnoses. Consult a healthcare professional "
    "before making any medical decisions."
)

REQUIRED_CNV_COLUMNS = {"cnv_id", "chrom", "start", "end", "type"}
OPTIONAL_CNV_COLUMNS = {"inheritance", "case_evidence_points"}

# ClinGen / ACMG 2020 classification thresholds (loss and gain share the table).
# Symmetric with the pathogenic side: Pathogenic >= +0.99 / Benign <= -0.99.
PATHOGENIC = 0.99
LIKELY_PATHOGENIC = 0.90
LIKELY_BENIGN = -0.90
BENIGN = -0.99

# Section 5 inheritance points (analyst-confirmed).
DE_NOVO_POINTS = 0.45        # 5A: observed de novo (both parents tested)
INHERITED_POINTS = -0.30     # 5B: inherited from a clinically unaffected parent

CASE_EVIDENCE_CAP = 0.90     # Section 4 clamp (per-CNV aggregate curator evidence)


# --------------------------------------------------------------------------- #
# Parsing helpers
# --------------------------------------------------------------------------- #
def normalize_chrom(value: str) -> str:
    return str(value).strip().lower().removeprefix("chr")


def normalize_type(value: str, chrom: str | None = None) -> str:
    """Map raw SV/CNV type tokens to 'loss' or 'gain'.

    When ``chrom`` is supplied, copy-number CN1 (single copy) on a sex
    chromosome is rejected as ambiguous, because CN1 on chrX/chrY is the normal
    hemizygous state in males and cannot be called a deletion without sex
    information. Use an explicit DEL/DUP in that case.
    """
    token = str(value).strip().upper().lstrip("<").rstrip(">")
    if chrom is not None and token == "CN1" and normalize_chrom(chrom) in {"x", "y"}:
        raise ValueError(
            "CN1 on a sex chromosome is ambiguous (normal hemizygous state in males); "
            "specify DEL or DUP explicitly"
        )
    if token in {"DEL", "LOSS", "DELETION"} or token in {"CN0", "CN1"}:
        return "loss"
    if token in {"DUP", "GAIN", "DUPLICATION"} or (token.startswith("CN") and token[2:].isdigit() and int(token[2:]) >= 3):
        return "gain"
    raise ValueError(f"Unsupported CNV type {value!r}; expected a deletion (DEL/LOSS) or duplication (DUP/GAIN)")


def _open(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open(newline="", encoding="utf-8")


def _parse_vcf(path: Path) -> list[dict]:
    rows: list[dict] = []
    with _open(path) as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            chrom, pos, vid, _ref, alt, _qual, _filt, info = fields[:8]
            info_map = {}
            for item in info.split(";"):
                if "=" in item:
                    k, v = item.split("=", 1)
                    info_map[k] = v
            svtype = info_map.get("SVTYPE") or alt
            end = info_map.get("END")
            if end is None:
                raise ValueError(f"VCF record {vid or chrom + ':' + pos} has no END in INFO; cannot size the CNV")
            rows.append({
                "cnv_id": vid if vid not in (".", "") else f"{chrom}:{pos}",
                "chrom": chrom,
                "start": int(pos),
                "end": int(end),
                "type": svtype,
                "inheritance": "unknown",
                "case_evidence_points": 0.0,
            })
    if not rows:
        raise ValueError("VCF contains no structural-variant records")
    return rows


def _parse_table(path: Path) -> list[dict]:
    with _open(path) as handle:
        sample = handle.read(4096)
        handle.seek(0)
        delimiter = "\t" if ("\t" in sample and sample.count("\t") >= sample.count(",")) else ","
        reader = csv.DictReader(handle, delimiter=delimiter)
        cols = set(reader.fieldnames or [])
        missing = REQUIRED_CNV_COLUMNS - cols
        if missing:
            raise ValueError(f"Input is missing required columns: {sorted(missing)}")
        rows: list[dict] = []
        for raw in reader:
            try:
                row = {
                    "cnv_id": raw["cnv_id"],
                    "chrom": raw["chrom"],
                    "start": int(float(raw["start"])),
                    "end": int(float(raw["end"])),
                    "type": raw["type"],
                    "inheritance": (raw.get("inheritance") or "unknown").strip() or "unknown",
                    "case_evidence_points": float(raw.get("case_evidence_points") or 0.0),
                }
            except (TypeError, ValueError) as exc:
                raise ValueError(f"Malformed CNV row {raw.get('cnv_id', '?')}: {exc}") from exc
            if row["end"] < row["start"]:
                raise ValueError(f"CNV {row['cnv_id']} has end < start")
            rows.append(row)
    if not rows:
        raise ValueError("Input contains no CNV records")
    return rows


def load_cnvs(path: Path) -> list[dict]:
    path = Path(path)
    name = path.name.lower()
    if name.endswith(".vcf") or name.endswith(".vcf.gz"):
        return _parse_vcf(path)
    return _parse_table(path)


def load_dosage_map(path: Path) -> list[dict]:
    rows: list[dict] = []
    with _open(Path(path)) as handle:
        for raw in csv.DictReader(handle):
            rows.append({
                "chrom": raw["chrom"],
                "start": int(raw["start"]),
                "end": int(raw["end"]),
                "name": raw["name"],
                "hi_score": int(float(raw.get("hi_score") or 0)),
                "ts_score": int(float(raw.get("ts_score") or 0)),
                "benign": str(raw.get("benign", "")).strip().lower() in {"1", "true", "yes"},
                "element_type": (raw.get("element_type") or "gene").strip().lower(),
                "strand": (raw.get("strand") or "+").strip() or "+",
                "cds_start": int(float(raw["cds_start"])) if (raw.get("cds_start") or "").strip() else None,
                "cds_end": int(float(raw["cds_end"])) if (raw.get("cds_end") or "").strip() else None,
            })
    return rows


def load_gene_model(path: Path) -> list[dict]:
    rows: list[dict] = []
    with _open(Path(path)) as handle:
        for raw in csv.DictReader(handle):
            rows.append({
                "chrom": raw["chrom"],
                "start": int(raw["start"]),
                "end": int(raw["end"]),
                "gene": raw["gene"],
            })
    return rows


# --------------------------------------------------------------------------- #
# Geometry
# --------------------------------------------------------------------------- #
def _overlaps(a_start: int, a_end: int, b_start: int, b_end: int) -> bool:
    return a_start <= b_end and a_end >= b_start


def _contains(outer_start: int, outer_end: int, inner_start: int, inner_end: int) -> bool:
    return outer_start <= inner_start and outer_end >= inner_end


# --------------------------------------------------------------------------- #
# Classification
# --------------------------------------------------------------------------- #
def classify_score(score: float) -> str:
    if score >= PATHOGENIC:
        return "Pathogenic"
    if score >= LIKELY_PATHOGENIC:
        return "Likely pathogenic"
    if score > LIKELY_BENIGN:
        return "Variant of uncertain significance"
    if score > BENIGN:
        return "Likely benign"
    return "Benign"


def _partial_gene_call(start: int, end: int, gene: dict) -> tuple[str, float, str]:
    """ClinGen Section 2 sub-call for a copy-number LOSS that partially overlaps an
    established haploinsufficient *gene*, derived from breakpoint geometry.

    Uses the gene's strand and coding boundaries to distinguish:
      2C-1 (+0.90) 5' end deleted, coding sequence involved -> predicted LoF
      2C-2 (0.00)  5' end deleted, 5'UTR only (no coding)
      2D-4 (+0.90) 3' end deleted, coding exon(s) of an established HI gene involved
      2D-1 (0.00)  3' end deleted, 3'UTR only
      2E   (+0.90 / 0.00) intragenic (both breakpoints inside the gene), coding vs not
    No free-text override is consulted; the call is a function of coordinates only.
    """
    gs, ge = gene["start"], gene["end"]
    strand = gene.get("strand", "+")
    cds_s = gene.get("cds_start") if gene.get("cds_start") is not None else gs
    cds_e = gene.get("cds_end") if gene.get("cds_end") is not None else ge
    five_prime = gs if strand == "+" else ge
    three_prime = ge if strand == "+" else gs
    covers_5p = start <= five_prime <= end
    covers_3p = start <= three_prime <= end
    cds_involved = start <= cds_e and end >= cds_s
    name = gene["name"]
    if covers_5p and not covers_3p:
        if cds_involved:
            return "2C-1", 0.90, f"5' end of established HI gene {name} deleted, coding sequence involved (predicted LoF)"
        return "2C-2", 0.0, f"5' end of {name} deleted, 5'UTR only (no coding sequence) — uncertain"
    if covers_3p and not covers_5p:
        if cds_involved:
            return "2D-4", 0.90, f"3' end of established HI gene {name} deleted, coding exon(s) involved (predicted LoF)"
        return "2D-1", 0.0, f"3' end of {name} deleted, 3'UTR only — uncertain"
    # both breakpoints inside the gene -> intragenic
    if cds_involved:
        return "2E", 0.90, f"Intragenic deletion within {name} coding sequence (predicted LoF)"
    return "2E", 0.0, f"Intragenic deletion in {name}, no coding sequence involved — uncertain"


def _gene_count_points(cnv_type: str, n_genes: int) -> tuple[str, float]:
    if cnv_type == "loss":
        if n_genes >= 35:
            return "3C", 0.90
        if n_genes >= 25:
            return "3B", 0.45
        return "3A", 0.0
    # gain thresholds are more permissive
    if n_genes >= 50:
        return "3C", 0.90
    if n_genes >= 35:
        return "3B", 0.45
    return "3A", 0.0


def score_cnv(cnv: dict, dosage_map: list[dict], gene_model: list[dict]) -> dict:
    cnv_type = normalize_type(cnv["type"], cnv.get("chrom"))
    chrom = normalize_chrom(cnv["chrom"])
    start, end = int(cnv["start"]), int(cnv["end"])

    genes = [g["gene"] for g in gene_model
             if normalize_chrom(g["chrom"]) == chrom and _overlaps(start, end, g["start"], g["end"])]
    dosage_hits = [d for d in dosage_map
                   if normalize_chrom(d["chrom"]) == chrom and _overlaps(start, end, d["start"], d["end"])]

    evidence: list[dict] = []

    def add(code: str, points: float, detail: str) -> None:
        evidence.append({"code": code, "points": round(points, 3), "detail": detail})

    # --- Section 1: genomic content ---
    has_content = bool(genes) or bool(dosage_hits)
    if not has_content:
        add("1B", -0.60, "CNV contains no protein-coding genes or known functional elements")
        total = -0.60
        return _result(cnv, cnv_type, genes, evidence, total)
    add("1A", 0.0, "CNV contains protein-coding or known functionally important elements")

    # --- Section 2: overlap with established dosage-sensitive / benign regions ---
    # Per Riggs 2020 the sections are ADDITIVE: Section 2 contributes one option's
    # points and the run continues. There is no early return for 2A or 2F — Sections
    # 4 (case) and 5 (inheritance) are always summed, so e.g. a complete 2A deletion
    # inherited from an unaffected parent is 1.00 + (-0.30) = 0.70 = VUS, exactly the
    # ClinGen worked example.
    relevant_key = "hi_score" if cnv_type == "loss" else "ts_score"
    benign_container = next((d for d in dosage_hits if d["benign"] and _contains(d["start"], d["end"], start, end)), None)
    established = [d for d in dosage_hits if not d["benign"] and d[relevant_key] >= 3]
    full = next((d for d in established if _contains(start, end, d["start"], d["end"])), None)
    section2_definitive = False  # only a full established/benign call excludes Section 3
    if benign_container is not None:
        add("2F", -1.00, f"CNV completely contained within established benign region {benign_container['name']}")
        section2_definitive = True
    elif full is not None:
        add("2A", 1.00, f"Complete overlap of established dosage-sensitive {full.get('element_type', 'gene')} {full['name']}")
        section2_definitive = True
    elif established:
        # Partial overlap of an established element.
        gene_entries = [d for d in established if d.get("element_type", "gene") == "gene"]
        if cnv_type == "loss" and gene_entries:
            code2, pts2, detail2 = _partial_gene_call(start, end, gene_entries[0])
            add(code2, pts2, detail2)
        else:
            # region-level partial overlap, or any gain partial: not scored positively
            add("2B", 0.0, f"Partial overlap of established dosage region {established[0]['name']}; critical gene/dose effect not established (uncertain)")

    # --- Section 3: number of protein-coding genes ---
    # ClinGen note: gene-count evaluates content NOT already covered by an established
    # dosage-sensitive gene/region. It is omitted only when a *complete* established
    # call (2A / 2F) was made; a 0-scoring partial (2B / 2C-2 / 2D-1 / 2E) does not
    # suppress it, so a deletion that merely clips an established region but spans
    # many genes is still scored on size.
    if not section2_definitive:
        code3, pts3 = _gene_count_points(cnv_type, len(genes))
        add(code3, pts3, f"{len(genes)} protein-coding gene(s) wholly or partially included")

    # --- Section 4: curator case/literature evidence (clamped, always summed) ---
    case_pts = float(cnv.get("case_evidence_points") or 0.0)
    if case_pts:
        clamped = max(-CASE_EVIDENCE_CAP, min(CASE_EVIDENCE_CAP, case_pts))
        add("4", clamped, "Analyst-supplied case-level / literature evidence")

    # --- Section 5: inheritance (always summed) ---
    inh = str(cnv.get("inheritance", "unknown")).strip().lower()
    if inh in {"de_novo", "denovo", "de novo"}:
        add("5A", DE_NOVO_POINTS, "Observed de novo (both parents tested)")
    elif inh in {"inherited", "maternal", "paternal", "inherited_unaffected"}:
        add("5B", INHERITED_POINTS, "Inherited from a clinically unaffected parent")

    total = sum(e["points"] for e in evidence)
    return _result(cnv, cnv_type, genes, evidence, total)


def _result(cnv: dict, cnv_type: str, genes: list[str], evidence: list[dict], total: float) -> dict:
    total = round(total, 2)
    return {
        "cnv_id": cnv["cnv_id"],
        "chrom": cnv["chrom"],
        "start": int(cnv["start"]),
        "end": int(cnv["end"]),
        "cnv_type": cnv_type,
        "gene_count": len(genes),
        "genes": genes,
        "evidence": evidence,
        "total_score": total,
        "classification": classify_score(total),
    }


def classify_all(cnvs: list[dict], dosage_map: list[dict], gene_model: list[dict]) -> dict:
    classifications = [score_cnv(c, dosage_map, gene_model) for c in cnvs]
    tier_counts: dict[str, int] = {}
    for c in classifications:
        tier_counts[c["classification"]] = tier_counts.get(c["classification"], 0) + 1
    return {
        "skill": "cnv-acmg-classifier",
        "framework": "ClinGen/ACMG 2019 (Riggs et al. 2020)",
        "summary": {"cnv_count": len(classifications), "tier_counts": tier_counts},
        "classifications": classifications,
        "disclaimer": DISCLAIMER,
    }


# --------------------------------------------------------------------------- #
# Output
# --------------------------------------------------------------------------- #
def write_outputs(result: dict, input_path: Path, output_dir: Path, command: list[str], demo: bool) -> None:
    output_dir = Path(output_dir)
    if output_dir.exists() and any(output_dir.iterdir()):
        print(f"WARNING: output directory exists; files may be overwritten: {output_dir}", file=sys.stderr)
    (output_dir / "tables").mkdir(parents=True, exist_ok=True)

    with (output_dir / "tables" / "cnv_classifications.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=[
            "cnv_id", "chrom", "start", "end", "cnv_type", "gene_count",
            "total_score", "classification", "evidence_codes"])
        writer.writeheader()
        for c in result["classifications"]:
            writer.writerow({
                "cnv_id": c["cnv_id"], "chrom": c["chrom"], "start": c["start"], "end": c["end"],
                "cnv_type": c["cnv_type"], "gene_count": c["gene_count"],
                "total_score": c["total_score"], "classification": c["classification"],
                "evidence_codes": ";".join(e["code"] for e in c["evidence"]),
            })

    (output_dir / "result.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")

    lines = [
        "# CNV ACMG Classification Report",
        "",
        f"**Input**: `{input_path}`",
        f"**Mode**: {'Synthetic demo data' if demo else 'User-provided local data'}",
        f"**Framework**: {result['framework']}",
        f"**CNVs classified**: {result['summary']['cnv_count']}",
        "",
        "| CNV | Region | Type | Genes | Score | Classification | Evidence |",
        "|---|---|---|---:|---:|---|---|",
    ]
    for c in result["classifications"]:
        region = f"{c['chrom']}:{c['start']:,}-{c['end']:,}"
        codes = ", ".join(e["code"] for e in c["evidence"])
        lines.append(
            f"| {c['cnv_id']} | {region} | {c['cnv_type']} | {c['gene_count']} | "
            f"{c['total_score']:.2f} | {c['classification']} | {codes} |"
        )
    lines += [
        "",
        "## Tier counts",
        "",
    ]
    for tier, n in result["summary"]["tier_counts"].items():
        lines.append(f"- {tier}: {n}")
    lines += [
        "",
        "## Interpretation",
        "",
        "Scores follow the ClinGen/ACMG copy-number point framework: Pathogenic ≥ 0.99, "
        "Likely pathogenic 0.90–0.98, VUS −0.89 to 0.89, Likely benign −0.90 to −0.98, "
        "Benign ≤ −0.99. Section 3 gene-count points apply only when no established "
        "dosage gene/region is overlapped. Sections 4 (case/literature) and 5 "
        "(inheritance) reflect analyst-supplied evidence and are never auto-generated.",
        "",
        DISCLAIMER,
        "",
    ]
    (output_dir / "report.md").write_text("\n".join(lines), encoding="utf-8")

    # --- reproducibility bundle (ClawBio contract: commands.sh + environment.yml + checksums) ---
    write_commands_sh(output_dir, shlex.join(command))
    # stdlib-only skill: no pip deps
    env_path = write_environment_yml(
        output_dir, "clawbio-cnv-acmg-classifier", [], channels=["conda-forge", "nodefaults"])
    write_checksums(
        [output_dir / "report.md", output_dir / "result.json",
         output_dir / "tables" / "cnv_classifications.csv", env_path],
        output_dir, anchor=output_dir)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="CNV ACMG Classifier (ClinGen/ACMG 2019 point framework)")
    parser.add_argument("--input", type=Path, help="CNV calls: VCF (SVTYPE/END) or CSV/TSV")
    parser.add_argument("--output", type=Path, default=Path("cnv_acmg_out"))
    parser.add_argument("--dosage-map", type=Path, help="Dosage-sensitivity map CSV (overrides bundled curated map)")
    parser.add_argument("--gene-model", type=Path, help="Protein-coding gene model CSV (overrides bundled demo model)")
    parser.add_argument("--demo", action="store_true")
    args = parser.parse_args(argv)

    input_path = SKILL_DIR / "demo_cnv_calls.csv" if args.demo else args.input
    if input_path is None:
        parser.error("--input is required unless --demo is used")

    dosage_path = args.dosage_map or (SKILL_DIR / "data" / "curated_dosage_map.csv")
    gene_path = args.gene_model or (SKILL_DIR / "data" / "curated_gene_model.csv")

    try:
        cnvs = load_cnvs(input_path)
        dosage_map = load_dosage_map(dosage_path)
        gene_model = load_gene_model(gene_path)
        result = classify_all(cnvs, dosage_map, gene_model)
    except (ValueError, FileNotFoundError, KeyError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    command = ["python", SCRIPT_REL]
    if args.demo:
        command.append("--demo")
    else:
        command += ["--input", str(args.input)]
    if args.dosage_map:
        command += ["--dosage-map", str(args.dosage_map)]
    if args.gene_model:
        command += ["--gene-model", str(args.gene_model)]
    command += ["--output", str(args.output)]
    write_outputs(result, input_path, args.output, command, args.demo)
    print(f"CNV ACMG Classifier wrote {args.output / 'report.md'} ({result['summary']['cnv_count']} CNVs)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
