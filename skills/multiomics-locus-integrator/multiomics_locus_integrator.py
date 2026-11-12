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

# A sparse --gene-targets file can leave the only candidate gene
# megabases away. 1Mb is a generous but finite
# cis-regulatory distance (most enhancer-promoter contacts and TAD
# sizes fall well inside it); past it, the nearest gene is flagged as
# distant rather than silently reported as if it were meaningful. The
# gene and its real distance are still returned -- never hidden.
NEAREST_GENE_FAR_THRESHOLD_BP = 1_000_000

# Promoter-proximal TSS window: a deliberately chosen, declared
# convention (not a specific tool's default -- HOMER uses roughly
# -1000/+100bp and ChIPseeker defaults to about +-3kb, neither of
# which matches these numbers). Upstream is more generous than
# downstream because most regulatory elements cluster just upstream
# of the TSS.
TSS_UPSTREAM_BP = 1000
TSS_DOWNSTREAM_BP = 250

# Named out of the report-building list (CodeQL: implicit string
# concatenation inside a list literal looks like a missing comma, even
# though this one is intentional -- one long sentence split across source
# lines for readability). Naming it here removes the ambiguous pattern
# entirely rather than just suppressing the warning.
_RESEARCH_TOOL_NOTE = (
    "ClawBio is a research tool. Locus tiers are a starting point for manual review, "
    "not a substitute for experimental validation (e.g. CRISPRi/luciferase for enhancer "
    "activity, 4C/Capture-C for contact confirmation)."
)

Interval = tuple[str, int, int, str]
# Genes carry their strand-aware TSS point alongside the regulatory
# window, so ranking/distance can use the real TSS directly instead of
# window-edge distance -- the window alone could occasionally rank the
# wrong gene as nearest, not just misreport the distance number.
GeneWindow = tuple[str, int, int, str, int]
LoopAnchor = tuple[str, int, int, str, int, int, str]


def _field_count_error(path: Path, line_no: int, line: str, fields: list[str], minimum: int, context: str) -> ValueError:
    """Build a field-count error, naming a likely wrong-delimiter cause
    directly rather than just reporting the count.
    """
    if len(fields) < minimum and len(line.split()) >= minimum:
        return ValueError(
            f"{path}:{line_no}: {context} needs >={minimum} tab-separated "
            f"fields, but got {len(fields)} -- this line looks "
            f"space-separated instead of tab-separated. BED/BEDPE files "
            f"must use literal tab characters between fields."
        )
    return ValueError(f"{path}:{line_no}: {context} needs >={minimum} tab-separated fields")


def parse_bed(path: Path) -> list[Interval]:
    intervals: list[Interval] = []
    with path.open(encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, start=1):
            line = line.strip()
            if not line or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                raise _field_count_error(path, line_no, line, fields, 3, "BED record")
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            name = fields[3] if len(fields) > 3 else f"{chrom}:{start}-{end}"
            intervals.append((chrom, start, end, name))
    return intervals


def parse_gene_targets(path: Path) -> list[GeneWindow]:
    """Parse gene-targets BED into strand-aware TSS windows:
    [TSS-1000bp, TSS+250bp), carrying the TSS point itself alongside
    the window.

    Accepts two column layouts: this skill's
    own 5-column form (chrom, start, end, gene_symbol, strand), or a
    standard BED6 file (chrom, start, end, name, score, strand) so
    users do not have to convert a real BED6 file first -- column 6 is
    used as strand whenever it is a valid '+'/'-', falling back to
    column 5 otherwise.

    The window is still what contact_genes() overlaps against (a locus
    reaching anywhere in the regulatory window counts); the TSS point is
    what nearest_gene() ranks by, so ranking and the reported distance
    both reflect genuine proximity to the TSS, not to whichever window
    edge happens to be nearer.
    """
    genes: list[GeneWindow] = []
    with path.open(encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, start=1):
            line = line.strip()
            if not line or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split("\t")
            if len(fields) < 5:
                raise _field_count_error(
                    path, line_no, line, fields, 5,
                    "gene-targets BED (chrom, start, end, gene_symbol, strand)",
                )
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            name = fields[3]
            # BED6 carries strand in column 6 (fields[5]), with a score
            # in column 5 (fields[4]) that this skill does not use.
            # Prefer it when present and valid; fall back to this
            # skill's own column-5 convention otherwise.
            # Deterministic, not a fallback chain: a 6+-field line is
            # BED6-shaped, so column 6 IS strand -- never cross-check
            # column 5 (that is the score column in real BED6) if
            # column 6 turns out invalid. A 5-field line is this
            # skill's own format, where column 5 is strand. An
            # earlier version fell back to column 5 whenever column 6
            # was not a recognised strand, so a BED6 "." (unstranded)
            # line raised an error naming the score value instead of
            # saying which column was actually wrong.
            if len(fields) >= 6:
                strand, strand_col = fields[5], 6
            else:
                strand, strand_col = fields[4], 5
            if strand == "+":
                tss = start
                win_start, win_end = tss - TSS_UPSTREAM_BP, tss + TSS_DOWNSTREAM_BP
            elif strand == "-":
                # BED end is exclusive (half-open [start, end)), so the
                # feature's actual last base -- the TSS on the minus
                # strand -- is end-1, not end.
                tss = end - 1
                win_start, win_end = tss - TSS_DOWNSTREAM_BP, tss + TSS_UPSTREAM_BP
            else:
                raise ValueError(f"{path}:{line_no}: strand must be '+' or '-', got {strand!r} in column {strand_col}")
            genes.append((chrom, max(win_start, 0), win_end, name, tss))
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
                raise _field_count_error(path, line_no, line, fields, 6, "BEDPE record")
            chrom1, start1, end1 = fields[0], int(fields[1]), int(fields[2])
            chrom2, start2, end2 = fields[3], int(fields[4]), int(fields[5])
            name = fields[6] if len(fields) > 6 else f"loop_{line_no}"
            loops.append((chrom1, start1, end1, chrom2, start2, end2, name))
    return loops


def _normalize_chrom(chrom: str) -> str:
    """Normalize a 'chr' naming-convention difference for comparison only --
    'chr5' and '5' are the same chromosome under this convention, so every
    cross-file chromosome comparison in this skill goes through this first.
    A locus on '5' against a track on 'chr5' used to silently compare
    unequal, which reads as a false not_assayed/no-overlap rather than the
    real match. Every per-locus output field still reports the chromosome
    exactly as given in its own source file; only the comparison is
    normalized, never the data.

    This does NOT and cannot detect a genome-build mismatch (e.g. hg19 vs
    hg38 coordinates) -- that requires reference data this skill does not
    have, and a build mismatch can share identical 'chr5'-vs-'chr5' naming
    while meaning two different sequences. Verify all inputs share one
    assembly before trusting a result; this skill only guards against the
    naming-convention difference, not the build difference.
    """
    return chrom[3:] if len(chrom) > 3 and chrom[:3].lower() == "chr" else chrom


def _overlaps(chrom: str, start: int, end: int, chrom2: str, start2: int, end2: int) -> bool:
    return _normalize_chrom(chrom) == _normalize_chrom(chrom2) and start < end2 and start2 < end


def chrom_spans(intervals: list[Interval]) -> dict[str, tuple[int, int]]:
    spans: dict[str, tuple[int, int]] = {}
    for chrom, start, end, _name in intervals:
        lo, hi = spans.get(chrom, (start, end))
        spans[chrom] = (min(lo, start), max(hi, end))
    return spans


def layer_coverage_blocks(intervals: list[Interval]) -> list[tuple[str, int, int]]:
    """Merged contiguous blocks per chromosome, from only touching or
    overlapping intervals -- never a chromosome-wide min/max. Used ONLY
    for an explicit --coverage-bed override, where the file's whole
    purpose is to state the real assayed footprint (e.g. a capture
    panel's target regions): there, a gap between two entries genuinely
    means "never measured," unlike the default peak-track span (see
    chrom_spans above), where a gap between two peaks is a true negative
    on a whole-genome assay, not missing data.
    """
    by_chrom: dict[str, list[tuple[int, int]]] = {}
    for chrom, start, end, _name in intervals:
        by_chrom.setdefault(chrom, []).append((start, end))

    blocks: list[tuple[str, int, int]] = []
    for chrom, spans in by_chrom.items():
        spans.sort()
        cur_start, cur_end = spans[0]
        for start, end in spans[1:]:
            if start <= cur_end:
                cur_end = max(cur_end, end)
            else:
                blocks.append((chrom, cur_start, cur_end))
                cur_start, cur_end = start, end
        blocks.append((chrom, cur_start, cur_end))
    return blocks


def loop_chrom_spans(loops: list[LoopAnchor]) -> dict[str, tuple[int, int]]:
    """Per-chromosome min/max across every loop's anchors, for the
    informational "## Layer coverage" summary only -- NOT for deciding
    not_assayed (see loop_spans_per_loop below for why merging across
    unrelated loops is wrong for that).
    """
    spans: dict[str, tuple[int, int]] = {}
    for chrom1, start1, end1, chrom2, start2, end2, _name in loops:
        for chrom, start, end in ((chrom1, start1, end1), (chrom2, start2, end2)):
            lo, hi = spans.get(chrom, (start, end))
            spans[chrom] = (min(lo, start), max(hi, end))
    return spans


def loop_spans_per_loop(loops: list[LoopAnchor]) -> list[tuple[str, int, int]]:
    """One (chrom, lo, hi) span PER LOOP, covering just that loop's own
    two anchors -- never merged across different loops. A chromosome-wide
    merge would let a single distant loop widen the "assayed" region for
    a locus nowhere near any actual loop; checking
    against each loop's own span instead keeps not_assayed as tight as
    the real evidence.
    """
    spans: list[tuple[str, int, int]] = []
    for chrom1, start1, end1, chrom2, start2, end2, _name in loops:
        if _normalize_chrom(chrom1) == _normalize_chrom(chrom2):
            spans.append((chrom1, min(start1, start2), max(end1, end2)))
        else:
            spans.append((chrom1, start1, end1))
            spans.append((chrom2, start2, end2))
    return spans


def within_span(chrom: str, start: int, end: int, spans: dict[str, tuple[int, int]]) -> bool:
    key = _normalize_chrom(chrom)
    return any(
        _normalize_chrom(c) == key and start < hi and end > lo
        for c, (lo, hi) in spans.items()
    )


def within_any_span(chrom: str, start: int, end: int, spans: list[tuple[str, int, int]]) -> bool:
    key = _normalize_chrom(chrom)
    return any(_normalize_chrom(c) == key and start < hi and end > lo for c, lo, hi in spans)


def overlapping_names(chrom: str, start: int, end: int, intervals: list[Interval]) -> list[str]:
    return [name for (c, s, e, name) in intervals if _overlaps(chrom, start, end, c, s, e)]


def nearest_gene(chrom: str, start: int, end: int, genes: list[GeneWindow]) -> tuple[str | None, int | None]:
    best_name, best_dist = None, None
    midpoint = (start + end) // 2
    key = _normalize_chrom(chrom)
    for g_chrom, _g_start, _g_end, g_name, g_tss in genes:
        if _normalize_chrom(g_chrom) != key:
            continue
        dist = abs(midpoint - g_tss)
        if best_dist is None or dist < best_dist:
            best_name, best_dist = g_name, dist
    return best_name, best_dist


def locus_in_any_loop(chrom: str, start: int, end: int, loops: list[LoopAnchor]) -> bool:
    """Pure interval overlap: does this locus sit in ANY loop anchor at all?

    Deliberately independent of whether a gene is known at the other anchor —
    "is there contact evidence here" and "do we know the gene it reaches" are
    two separate questions (see contact_genes below). Conflating them let a
    sparse --gene-targets list silently hide real loop evidence.
    """
    for chrom1, start1, end1, chrom2, start2, end2, _name in loops:
        if _overlaps(chrom, start, end, chrom1, start1, end1) or _overlaps(
            chrom, start, end, chrom2, start2, end2
        ):
            return True
    return False


def contact_genes(
    chrom: str, start: int, end: int, loops: list[LoopAnchor], genes: list[GeneWindow]
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
        for g_chrom, g_start, g_end, g_name, _g_tss in genes:
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
    layer_coverage_overrides: dict[str, list[tuple[str, int, int]]],
    loops: list[LoopAnchor],
    loop_span: list[tuple[str, int, int]],
    contact_coverage_override: list[tuple[str, int, int]] | None,
    genes: list[GeneWindow],
) -> dict:
    checked: list[str] = []
    supported: list[str] = []
    not_assayed: list[str] = []

    for layer_name, intervals in layers.items():
        if layer_name in layer_coverage_overrides:
            assayed = within_any_span(chrom, start, end, layer_coverage_overrides[layer_name])
        else:
            assayed = within_span(chrom, start, end, layer_spans[layer_name])
        if not assayed:
            not_assayed.append(layer_name)
            continue
        checked.append(layer_name)
        if overlapping_names(chrom, start, end, intervals):
            supported.append(layer_name)

    if loops:
        # An explicit --contact-coverage-bed is opt-in, same as the
        # accessibility/histone coverage overrides: when supplied, it
        # states the real assayed footprint of the 3D-contact assay (e.g.
        # a Capture-C/HiChIP capture region), so a gap within it is a
        # true negative, not missing data -- unlike the default per-loop
        # span, which only ever covers each loop's own two anchors.
        effective_contact_span = (
            contact_coverage_override if contact_coverage_override is not None else loop_span
        )
        if within_any_span(chrom, start, end, effective_contact_span):
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
    # Report every contacted gene, not just the first: a
    # locus can have more than one loop, or a loop can reach more than one
    # gene at its far anchor, and dropping any of them could silently hide
    # a real nearest-vs-contact disagreement.
    contact_gene_names = sorted({name for name, _loop in contacts})
    contact_loop_ids = sorted({loop for _name, loop in contacts})
    gene_contact = ",".join(contact_gene_names) if contact_gene_names else None
    loop_name = ",".join(contact_loop_ids) if contact_loop_ids else None
    disagreement = bool(
        gene_nearest and any(name != gene_nearest for name in contact_gene_names)
    )

    nearest_is_distant = bool(
        dist_nearest is not None and dist_nearest > NEAREST_GENE_FAR_THRESHOLD_BP
    )

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
        "target_gene_nearest_is_distant": nearest_is_distant,
        "target_gene_contact": gene_contact,
        "contact_loop_id": loop_name,
        "nearest_vs_contact_disagreement": disagreement,
    }


def analyze(
    loci: list[Interval],
    accessibility: list[Interval],
    histone: list[Interval],
    loops: list[LoopAnchor],
    genes: list[GeneWindow],
    accessibility_coverage: list[Interval] | None = None,
    histone_coverage: list[Interval] | None = None,
    contact_coverage: list[Interval] | None = None,
) -> dict:
    layers = {}
    if accessibility:
        layers["accessibility"] = accessibility
    if histone:
        layers["histone"] = histone

    # An explicit --coverage-bed is opt-in, per layer, and never changes
    # the default chromosome-wide-span behaviour for a layer that does
    # not supply one.
    coverage_by_layer = {"accessibility": accessibility_coverage, "histone": histone_coverage}
    layer_coverage_overrides = {
        name: layer_coverage_blocks(coverage)
        for name, coverage in coverage_by_layer.items()
        if coverage
    }
    layer_spans = {
        name: chrom_spans(intervals)
        for name, intervals in layers.items()
        if name not in layer_coverage_overrides
    }
    loop_span = loop_spans_per_loop(loops)
    contact_coverage_override = layer_coverage_blocks(contact_coverage) if contact_coverage else None

    results = [
        score_locus(
            chrom, start, end, name, layers, layer_spans, layer_coverage_overrides,
            loops, loop_span, contact_coverage_override, genes,
        )
        for (chrom, start, end, name) in loci
    ]

    tier_counts: dict[str, int] = {}
    for r in results:
        tier_counts[r["tier"]] = tier_counts.get(r["tier"], 0) + 1

    # A layer's observed span can be silently stretched by one
    # stray distant interval, masking a
    # real coverage gap elsewhere. Surfacing the actual span per
    # layer lets a reader spot that directly. A layer with an explicit
    # --coverage-bed override is shown as its merged blocks instead,
    # since that is what was actually used to decide not_assayed for it.
    layer_coverage: dict[str, str] = {
        layer: ",".join(f"{chrom}:{lo}-{hi}" for chrom, (lo, hi) in sorted(spans.items()))
        for layer, spans in layer_spans.items()
    }
    layer_coverage.update({
        layer: ",".join(f"{chrom}:{lo}-{hi}" for chrom, lo, hi in sorted(blocks))
        for layer, blocks in layer_coverage_overrides.items()
    })
    if contact_coverage_override is not None:
        layer_coverage["contact"] = ",".join(
            f"{chrom}:{lo}-{hi}" for chrom, lo, hi in sorted(contact_coverage_override)
        )
    else:
        loop_display_span = loop_chrom_spans(loops)
        if loop_display_span:
            layer_coverage["contact"] = ",".join(
                f"{chrom}:{lo}-{hi}" for chrom, (lo, hi) in sorted(loop_display_span.items())
            )

    return {
        "skill": "multiomics-locus-integrator",
        "summary": {
            "locus_count": len(results),
            "layers_supplied": sorted(list(layers.keys()) + (["contact"] if loops else [])),
            "layer_coverage": layer_coverage,
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
        "target_gene_nearest_is_distant",
        "target_gene_contact", "contact_loop_id", "nearest_vs_contact_disagreement",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        # No leading comment line: the
        # distance column is now genuine TSS-point distance, so there is
        # no caveat left to document, and a plain csv.reader -- unlike
        # pandas.read_csv(comment="#") -- does not skip "#" lines, which
        # is a real risk for a chaining partner reading this table.
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
    ]
    coverage = result["summary"].get("layer_coverage") or {}
    if coverage:
        rows.append("## Layer coverage (observed span per supplied track)")
        rows.append("")
        rows.append(
            "A layer's `not_assayed` loci are only as trustworthy as this span -- "
            "one stray distant interval in a track inflates it and can mask a real "
            "coverage gap elsewhere."
        )
        rows.append("")
        for layer_name in sorted(coverage):
            rows.append(f"- **{layer_name}**: {coverage[layer_name]}")
        rows.append("")
    rows.append("| Locus | Coordinates | Score | Tier | Checked | Supported | Not assayed | Nearest gene | Distant | Contact gene | Disagreement |")
    rows.append("|---|---|---:|---|---|---|---|---|---|---|---|")
    for r in result["loci"]:
        score_cell = "" if r["score"] is None else r["score"]
        rows.append(
            f"| {r['locus_id']} | {r['chrom']}:{r['start']}-{r['end']} | {score_cell} | {r['tier']} | "
            f"{','.join(r['layers_checked']) or '-'} | {','.join(r['layers_supported']) or '-'} | "
            f"{','.join(r['layers_not_assayed']) or '-'} | {r['target_gene_nearest'] or '-'} | "
            f"{'YES' if r['target_gene_nearest_is_distant'] else ''} | "
            f"{r['target_gene_contact'] or '-'} | {'YES' if r['nearest_vs_contact_disagreement'] else ''} |"
        )
    rows.extend([
        "",
        "## Interpretation",
        "",
        result["note"],
        "",
        _RESEARCH_TOOL_NOTE,
        "",
    ])
    (output_dir / "report.md").write_text("\n".join(rows), encoding="utf-8")
    write_repro_bundle(loci_path, output_dir, demo, extra_inputs)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Multiomics Locus Integrator")
    parser.add_argument("--loci", type=Path)
    parser.add_argument("--accessibility-track", type=Path)
    parser.add_argument("--histone-track", type=Path)
    parser.add_argument(
        "--accessibility-coverage-bed", type=Path,
        help="Real assayed footprint for --accessibility-track (e.g. a capture "
             "panel's target regions). When supplied, not_assayed for this layer "
             "is decided from this file's own merged intervals instead of the "
             "track's chromosome-wide span -- use this for a targeted/capture "
             "track, where a gap really does mean 'never measured.'",
    )
    parser.add_argument(
        "--histone-coverage-bed", type=Path,
        help="Same as --accessibility-coverage-bed, for --histone-track.",
    )
    parser.add_argument("--contact-evidence", type=Path)
    parser.add_argument(
        "--contact-coverage-bed", type=Path,
        help="Real assayed footprint for --contact-evidence (e.g. a Capture-C/HiChIP "
             "capture region). When supplied, not_assayed for the contact layer is "
             "decided from this file's own merged intervals instead of each loop's "
             "own anchor span.",
    )
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

    accessibility_coverage_path = None if args.demo else args.accessibility_coverage_bed
    histone_coverage_path = None if args.demo else args.histone_coverage_bed
    contact_coverage_path = None if args.demo else args.contact_coverage_bed

    try:
        loci = parse_bed(loci_path)
        genes = parse_gene_targets(gene_path)
        accessibility = parse_bed(accessibility_path) if accessibility_path else []
        histone = parse_bed(histone_path) if histone_path else []
        loops = parse_bedpe(contact_path) if contact_path else []
        accessibility_coverage = parse_bed(accessibility_coverage_path) if accessibility_coverage_path else []
        histone_coverage = parse_bed(histone_coverage_path) if histone_coverage_path else []
        contact_coverage = parse_bed(contact_coverage_path) if contact_coverage_path else []
    except (ValueError, FileNotFoundError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    if not loci:
        print("ERROR: --loci contains no records", file=sys.stderr)
        return 2
    if not any([accessibility, histone, loops]):
        print("ERROR: at least one evidence layer (accessibility, histone, or contact) is required", file=sys.stderr)
        return 2

    result = analyze(
        loci, accessibility, histone, loops, genes,
        accessibility_coverage=accessibility_coverage,
        histone_coverage=histone_coverage,
        contact_coverage=contact_coverage,
    )
    extra_inputs = {
        "--accessibility-track": accessibility_path,
        "--histone-track": histone_path,
        "--accessibility-coverage-bed": accessibility_coverage_path,
        "--histone-coverage-bed": histone_coverage_path,
        "--contact-evidence": contact_path,
        "--contact-coverage-bed": contact_coverage_path,
        "--gene-targets": gene_path,
    }
    write_outputs(result, loci_path, args.output, args.demo, extra_inputs)
    print(f"Multiomics Locus Integrator wrote {args.output / 'report.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
