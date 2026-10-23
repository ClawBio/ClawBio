#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

_PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from clawbio.common.reproducibility import (  # noqa: E402
    ReproCommand,
    ReproPath,
    write_checksums,
    write_environment_yml,
    write_portable_commands_sh,
)

SKILL_DIR = Path(__file__).resolve().parent

# Layer weights: contact evidence (TAD/loop) is weighted higher because it is
# harder to obtain spuriously than a single accessibility or histone peak call.
LAYER_WEIGHTS = {"accessibility": 1.0, "histone": 1.0, "contact": 2.0}
HIGH_SCORE_THRESHOLD = 0.75
HIGH_MIN_DISTINCT_LAYERS = 2

# Promoter-proximal TSS window (standard convention, e.g. HOMER/ChIPseeker
# defaults in this range): upstream is more generous than downstream because
# most regulatory elements cluster just upstream of the TSS.
TSS_UPSTREAM_BP = 1000
TSS_DOWNSTREAM_BP = 250

Interval = tuple[str, int, int, str]
LoopAnchor = tuple[str, int, int, str, int, int, str]


def parse_bed(path: Path) -> list[Interval]:
    intervals: list[Interval] = []
    with path.open(encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, start=1):
            line = line.strip()
            if not line or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                raise ValueError(f"{path}:{line_no}: BED record needs >=3 tab-separated fields")
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            name = fields[3] if len(fields) > 3 else f"{chrom}:{start}-{end}"
            intervals.append((chrom, start, end, name))
    return intervals


def parse_gene_targets(path: Path) -> list[Interval]:
    """Parse gene-targets BED (chrom, start, end, gene_symbol, strand) into
    strand-aware TSS windows: [TSS-1000bp, TSS+250bp).

    The window, not the full gene body, IS the gene's region from here on —
    nearest_gene() and contact_genes() just see it as an ordinary Interval and
    need no changes. This replaced an earlier "nearest-TSS" label that was
    actually measuring distance to the nearer edge of the whole gene body
    (wrong for any minus-strand gene, since BED carries no strand of its own
    without this column).
    """
    genes: list[Interval] = []
    with path.open(encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, start=1):
            line = line.strip()
            if not line or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split("\t")
            if len(fields) < 5:
                raise ValueError(
                    f"{path}:{line_no}: gene-targets BED needs >=5 tab-separated "
                    f"fields (chrom, start, end, gene_symbol, strand)"
                )
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            name, strand = fields[3], fields[4]
            if strand == "+":
                tss = start
                win_start, win_end = tss - TSS_UPSTREAM_BP, tss + TSS_DOWNSTREAM_BP
            elif strand == "-":
                tss = end
                win_start, win_end = tss - TSS_DOWNSTREAM_BP, tss + TSS_UPSTREAM_BP
            else:
                raise ValueError(f"{path}:{line_no}: strand must be '+' or '-', got {strand!r}")
            genes.append((chrom, max(win_start, 0), win_end, name))
    return genes


def parse_bedpe(path: Path) -> list[LoopAnchor]:
    loops: list[LoopAnchor] = []
    with path.open(encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, start=1):
            line = line.strip()
            if not line or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split("\t")
            if len(fields) < 6:
                raise ValueError(f"{path}:{line_no}: BEDPE record needs >=6 tab-separated fields")
            chrom1, start1, end1 = fields[0], int(fields[1]), int(fields[2])
            chrom2, start2, end2 = fields[3], int(fields[4]), int(fields[5])
            name = fields[6] if len(fields) > 6 else f"loop_{line_no}"
            loops.append((chrom1, start1, end1, chrom2, start2, end2, name))
    return loops


def _overlaps(chrom: str, start: int, end: int, chrom2: str, start2: int, end2: int) -> bool:
    return chrom == chrom2 and start < end2 and start2 < end


def chrom_spans(intervals: list[Interval]) -> dict[str, tuple[int, int]]:
    spans: dict[str, tuple[int, int]] = {}
    for chrom, start, end, _name in intervals:
        lo, hi = spans.get(chrom, (start, end))
        spans[chrom] = (min(lo, start), max(hi, end))
    return spans


def loop_chrom_spans(loops: list[LoopAnchor]) -> dict[str, tuple[int, int]]:
    spans: dict[str, tuple[int, int]] = {}
    for chrom1, start1, end1, chrom2, start2, end2, _name in loops:
        for chrom, start, end in ((chrom1, start1, end1), (chrom2, start2, end2)):
            lo, hi = spans.get(chrom, (start, end))
            spans[chrom] = (min(lo, start), max(hi, end))
    return spans


def within_span(chrom: str, start: int, end: int, spans: dict[str, tuple[int, int]]) -> bool:
    if chrom not in spans:
        return False
    lo, hi = spans[chrom]
    return start < hi and end > lo


def overlapping_names(chrom: str, start: int, end: int, intervals: list[Interval]) -> list[str]:
    return [name for (c, s, e, name) in intervals if _overlaps(chrom, start, end, c, s, e)]


def nearest_gene(chrom: str, start: int, end: int, genes: list[Interval]) -> tuple[str | None, int | None]:
    best_name, best_dist = None, None
    midpoint = (start + end) // 2
    for g_chrom, g_start, g_end, g_name in genes:
        if g_chrom != chrom:
            continue
        if g_start <= midpoint < g_end:
            dist = 0
        else:
            dist = min(abs(midpoint - g_start), abs(midpoint - g_end))
        if best_dist is None or dist < best_dist:
            best_name, best_dist = g_name, dist
    return best_name, best_dist


def locus_in_any_loop(chrom: str, start: int, end: int, loops: list[LoopAnchor]) -> bool:
    """Pure interval overlap: does this locus sit in ANY loop anchor at all?

    Deliberately independent of whether a gene is known at the other anchor —
    "is there contact evidence here" and "do we know the gene it reaches" are
    two separate questions (see contact_genes below). Conflating them let a
    sparse --gene-targets list silently hide real loop evidence (PR #678).
    """
    for chrom1, start1, end1, chrom2, start2, end2, _name in loops:
        if _overlaps(chrom, start, end, chrom1, start1, end1) or _overlaps(
            chrom, start, end, chrom2, start2, end2
        ):
            return True
    return False


def contact_genes(
    chrom: str, start: int, end: int, loops: list[LoopAnchor], genes: list[Interval]
) -> list[tuple[str, str]]:
    """Genes contacted via a loop whose OTHER anchor overlaps this locus."""
    hits: list[tuple[str, str]] = []
    for chrom1, start1, end1, chrom2, start2, end2, loop_name in loops:
        other_anchor = None
        if _overlaps(chrom, start, end, chrom1, start1, end1):
            other_anchor = (chrom2, start2, end2)
        elif _overlaps(chrom, start, end, chrom2, start2, end2):
            other_anchor = (chrom1, start1, end1)
        if other_anchor is None:
            continue
        o_chrom, o_start, o_end = other_anchor
        for g_chrom, g_start, g_end, g_name in genes:
            if _overlaps(o_chrom, o_start, o_end, g_chrom, g_start, g_end):
                hits.append((g_name, loop_name))
    return hits


def score_locus(
    chrom: str,
    start: int,
    end: int,
    name: str,
    layers: dict[str, list[Interval]],
    layer_spans: dict[str, dict[str, tuple[int, int]]],
    loops: list[LoopAnchor],
    loop_span: dict[str, tuple[int, int]],
    genes: list[Interval],
) -> dict:
    checked: list[str] = []
    supported: list[str] = []
    not_assayed: list[str] = []

    for layer_name, intervals in layers.items():
        if not within_span(chrom, start, end, layer_spans[layer_name]):
            not_assayed.append(layer_name)
            continue
        checked.append(layer_name)
        if overlapping_names(chrom, start, end, intervals):
            supported.append(layer_name)

    if loops:
        if within_span(chrom, start, end, loop_span):
            checked.append("contact")
            if locus_in_any_loop(chrom, start, end, loops):
                supported.append("contact")
        else:
            not_assayed.append("contact")

    checked_weight = sum(LAYER_WEIGHTS[layer] for layer in checked)
    supported_weight = sum(LAYER_WEIGHTS[layer] for layer in supported)

    if not checked:
        tier = "not_assayed"
        score = None
    else:
        score = round(supported_weight / checked_weight, 3)
        if score >= HIGH_SCORE_THRESHOLD and len(supported) >= HIGH_MIN_DISTINCT_LAYERS:
            tier = "high"
        elif supported:
            tier = "partial"
        else:
            tier = "unsupported"

    gene_nearest, dist_nearest = nearest_gene(chrom, start, end, genes)
    contacts = contact_genes(chrom, start, end, loops, genes) if "contact" in checked else []
    gene_contact = contacts[0][0] if contacts else None
    loop_name = contacts[0][1] if contacts else None
    disagreement = bool(gene_contact and gene_nearest and gene_contact != gene_nearest)

    return {
        "locus_id": name,
        "chrom": chrom,
        "start": start,
        "end": end,
        "layers_checked": sorted(checked),
        "layers_supported": sorted(supported),
        "layers_not_assayed": sorted(not_assayed),
        "score": score,
        "tier": tier,
        "target_gene_nearest": gene_nearest,
        "target_gene_nearest_distance_bp": dist_nearest,
        "target_gene_contact": gene_contact,
        "contact_loop_id": loop_name,
        "nearest_vs_contact_disagreement": disagreement,
    }


def analyze(
    loci: list[Interval],
    accessibility: list[Interval],
    histone: list[Interval],
    loops: list[LoopAnchor],
    genes: list[Interval],
) -> dict:
    layers = {}
    if accessibility:
        layers["accessibility"] = accessibility
    if histone:
        layers["histone"] = histone
    layer_spans = {name: chrom_spans(intervals) for name, intervals in layers.items()}
    loop_span = loop_chrom_spans(loops)

    results = [
        score_locus(chrom, start, end, name, layers, layer_spans, loops, loop_span, genes)
        for (chrom, start, end, name) in loci
    ]

    tier_counts: dict[str, int] = {}
    for r in results:
        tier_counts[r["tier"]] = tier_counts.get(r["tier"], 0) + 1

    return {
        "skill": "multiomics-locus-integrator",
        "summary": {
            "locus_count": len(results),
            "layers_supplied": sorted(list(layers.keys()) + (["contact"] if loops else [])),
            "tier_counts": tier_counts,
            "disagreement_count": sum(1 for r in results if r["nearest_vs_contact_disagreement"]),
        },
        "loci": results,
        "note": (
            "Locus scores are a deterministic weighted overlap heuristic over the "
            "evidence layers actually supplied, not a biological confidence probability. "
            "A 'not_assayed' layer means the supplied track did not cover this region, "
            "which is distinct from an 'unsupported' locus where every supplied layer "
            "was checked and found no overlap."
        ),
    }


def _write_evidence_table(path: Path, loci: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "locus_id", "chrom", "start", "end", "score", "tier",
        "layers_checked", "layers_supported", "layers_not_assayed",
        "target_gene_nearest", "target_gene_nearest_distance_bp",
        "target_gene_contact", "contact_loop_id", "nearest_vs_contact_disagreement",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in loci:
            flat = dict(row)
            for key in ("layers_checked", "layers_supported", "layers_not_assayed"):
                flat[key] = ",".join(row[key]) if row[key] else ""
            writer.writerow({field: flat.get(field, "") for field in fieldnames})


def repro_command(loci_path: Path, output_dir: Path, demo: bool, extra_inputs: dict[str, Path | None]) -> ReproCommand:
    def rp(path: Path) -> ReproPath:
        return ReproPath(path, "repo_root" if path.is_relative_to(_PROJECT_ROOT) else "auto")

    args: list[str | ReproPath] = ["--demo"] if demo else ["--loci", rp(loci_path)]
    if not demo:
        for flag, path in extra_inputs.items():
            if path is not None:
                args += [flag, rp(path)]
    args += ["--output", ReproPath(output_dir, "output_dir")]
    return ReproCommand(
        script_path=Path("skills/multiomics-locus-integrator/multiomics_locus_integrator.py"),
        args=args,
        comment="Reproduce this multiomics-locus-integrator run",
    )


def write_repro_bundle(
    loci_path: Path, output_dir: Path, demo: bool, extra_inputs: dict[str, Path | None]
) -> None:
    write_environment_yml(
        output_dir,
        env_name="clawbio-multiomics-locus-integrator",
        pip_deps=[],
        python_version="3.11",
    )
    write_portable_commands_sh(
        output_dir,
        repro_command(loci_path, output_dir, demo, extra_inputs),
        repo_root=_PROJECT_ROOT,
    )
    write_checksums(
        [
            output_dir / "report.md",
            output_dir / "result.json",
            output_dir / "tables" / "evidence_table.tsv",
        ],
        output_dir,
        anchor=output_dir,
    )


def write_outputs(
    result: dict, loci_path: Path, output_dir: Path, demo: bool, extra_inputs: dict[str, Path | None]
) -> None:
    if output_dir.exists() and any(output_dir.iterdir()):
        print(f"WARNING: output directory already exists and files may be overwritten: {output_dir}", file=sys.stderr)
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "tables").mkdir(exist_ok=True)
    (output_dir / "reproducibility").mkdir(exist_ok=True)

    _write_evidence_table(output_dir / "tables" / "evidence_table.tsv", result["loci"])
    (output_dir / "result.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")

    rows = [
        "# Multiomics Locus Integrator Report",
        "",
        f"**Input loci**: `{loci_path}`",
        f"**Mode**: {'Synthetic demo data' if demo else 'User-provided local data'}",
        f"**Loci evaluated**: {result['summary']['locus_count']}",
        f"**Layers supplied**: {', '.join(result['summary']['layers_supplied']) or 'none'}",
        f"**Tier counts**: {result['summary']['tier_counts']}",
        f"**Nearest-gene vs. contact-gene disagreements**: {result['summary']['disagreement_count']}",
        "",
        "| Locus | Coordinates | Score | Tier | Checked | Supported | Not assayed | Nearest gene | Contact gene | Disagreement |",
        "|---|---|---:|---|---|---|---|---|---|---|",
    ]
    for r in result["loci"]:
        score_cell = "" if r["score"] is None else r["score"]
        rows.append(
            f"| {r['locus_id']} | {r['chrom']}:{r['start']}-{r['end']} | {score_cell} | {r['tier']} | "
            f"{','.join(r['layers_checked']) or '-'} | {','.join(r['layers_supported']) or '-'} | "
            f"{','.join(r['layers_not_assayed']) or '-'} | {r['target_gene_nearest'] or '-'} | "
            f"{r['target_gene_contact'] or '-'} | {'YES' if r['nearest_vs_contact_disagreement'] else ''} |"
        )
    rows.extend([
        "",
        "## Interpretation",
        "",
        result["note"],
        "",
        "ClawBio is a research tool. Locus tiers are a starting point for manual review, "
        "not a substitute for experimental validation (e.g. CRISPRi/luciferase for enhancer "
        "activity, 4C/Capture-C for contact confirmation).",
        "",
    ])
    (output_dir / "report.md").write_text("\n".join(rows), encoding="utf-8")
    write_repro_bundle(loci_path, output_dir, demo, extra_inputs)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Multiomics Locus Integrator")
    parser.add_argument("--loci", type=Path)
    parser.add_argument("--accessibility-track", type=Path)
    parser.add_argument("--histone-track", type=Path)
    parser.add_argument("--contact-evidence", type=Path)
    parser.add_argument("--gene-targets", type=Path)
    parser.add_argument("--output", type=Path, default=Path("multiomics_locus_integrator_out"))
    parser.add_argument("--demo", action="store_true")
    args = parser.parse_args(argv)

    if args.demo:
        loci_path = SKILL_DIR / "demo_loci.bed"
        accessibility_path = SKILL_DIR / "demo_atac.bed"
        histone_path = SKILL_DIR / "demo_h3k27ac.bed"
        contact_path = SKILL_DIR / "demo_loops.bedpe"
        gene_path = SKILL_DIR / "demo_targets.bed"
    else:
        loci_path = args.loci
        accessibility_path = args.accessibility_track
        histone_path = args.histone_track
        contact_path = args.contact_evidence
        gene_path = args.gene_targets

    if loci_path is None:
        parser.error("--loci is required unless --demo is used")
    if gene_path is None:
        parser.error("--gene-targets is required unless --demo is used")

    try:
        loci = parse_bed(loci_path)
        genes = parse_gene_targets(gene_path)
        accessibility = parse_bed(accessibility_path) if accessibility_path else []
        histone = parse_bed(histone_path) if histone_path else []
        loops = parse_bedpe(contact_path) if contact_path else []
    except (ValueError, FileNotFoundError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    if not loci:
        print("ERROR: --loci contains no records", file=sys.stderr)
        return 2
    if not any([accessibility, histone, loops]):
        print("ERROR: at least one evidence layer (accessibility, histone, or contact) is required", file=sys.stderr)
        return 2

    result = analyze(loci, accessibility, histone, loops, genes)
    extra_inputs = {
        "--accessibility-track": accessibility_path,
        "--histone-track": histone_path,
        "--contact-evidence": contact_path,
        "--gene-targets": gene_path,
    }
    write_outputs(result, loci_path, args.output, args.demo, extra_inputs)
    print(f"Multiomics Locus Integrator wrote {args.output / 'report.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
