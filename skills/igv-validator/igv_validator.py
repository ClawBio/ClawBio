#!/usr/bin/env python3
"""
ClawBio — IGV Validator Skill
Read-level validation of somatic variant calls: counts tumor/normal read support in the
BAMs, flags artifact patterns, takes one tumor-over-normal IGV screenshot per variant and
writes a report. Numbers always come from the BAM, never from the screenshots.

Usage:
    python igv_validator.py --vcf calls.vcf.gz --tumor tumor.bam --normal normal.bam \
        --reference ref.fa --genes KRAS,BRAF --output report/
    python igv_validator.py --vcf calls.vcf.gz --tumor t.bam --normal n.bam \
        --variants list.tsv --no-igv --output report/
    python igv_validator.py --demo --output /tmp/igv_validator_demo
"""
from __future__ import annotations

import argparse
import base64
import csv
import hashlib
import html
import json
import os
import re
import shlex
import shutil
import socket
import subprocess
import sys
import tempfile
import time
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path

SKILL_DIR = Path(__file__).resolve().parent
_PROJECT_ROOT = SKILL_DIR.parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from clawbio.common.report import DISCLAIMER, write_result_json  # noqa: E402
from clawbio.common.reproducibility import (  # noqa: E402
    write_checksums,
    write_commands_sh,
    write_environment_yml,
)

try:
    import pysam
except ImportError:  # pragma: no cover - reported cleanly in main()
    pysam = None

VERSION = "0.1.0"
DEMO_DIR = SKILL_DIR / "demo"

# ── Thresholds (documented in SKILL.md; do not change without updating it) ─────
MIN_MAPQ = 20            # mapping quality floor
MIN_BQ = 20              # base quality floor (SNVs)
INDEL_WINDOW = 5         # an indel of the same length within ±5 bp counts
MIN_CLIP = 5             # soft clip length for indel clip rescue
SPLIT_WINDOW = 500       # SA tag must land within 500 bp of the partner breakpoint
PAIR_WINDOW = 1000       # discordant mate within 1 kb of the partner breakpoint
MIN_SV_LEN = 1000        # symbolic <DEL>/<DUP>/<INV> shorter than this are not evaluable
MIN_INV_LEN = 50         # inversions of any size are evaluated from same-strand read pairs
MIN_SUPPORT = 3          # fewer supporting reads -> insufficient
MIN_VAF = 5.0            # percent
STRAND_MIN = 4           # single-strand check only with at least this many reads
READ_END_BP = 10         # alt base within 10 bp of a read end
LOW_DEPTH = 10           # tumor or normal depth below this is flagged
GERMLINE_OTHER = 20.0    # percent of normal reads with a third allele at an SNV site
CALLER_VAF_DIFF = 10.0   # percentage points between caller-reported and BAM VAF
HIGH_DEPTH_X = 2.5       # depth above this multiple of the run's median suggests a repeat / mismapping
MAX_VARIANTS = 50

FLAG_TEXT = {
    "normal_support": "supporting reads in the normal",
    "low_support": f"fewer than {MIN_SUPPORT} supporting reads",
    "low_vaf": f"under {MIN_VAF:g}% of reads",
    "strand_bias": "all supporting reads on one strand",
    "read_end": f"alt bases mostly within {READ_END_BP} bp of read ends",
    "low_depth": f"under {LOW_DEPTH} reads in tumor or normal",
    "germline_site": f"the normal carries another allele here (>= {GERMLINE_OTHER:g}% of reads)",
    "caller_disagrees": f"the caller's reported counts (VCF AD) differ from the BAM by over {CALLER_VAF_DIFF:g} VAF points",
    "no_coverage": "no reads at this position in either BAM (region not in the BAM, or wrong BAM/contig)",
    "high_depth": f"depth over {HIGH_DEPTH_X:g}x the median of this run (possible repeat or mismapped reads)",
}


class InputError(Exception):
    """A problem with the user's inputs, reported without a traceback."""


@dataclass
class Variant:
    id: str
    chrom: str
    pos: int                 # 1-based VCF POS (anchor base for indels)
    ref: str
    alt: str
    kind: str                # snv | insertion | deletion | bnd | del | dup | inv
    gene: str = ""
    filter: str = "PASS"
    chrom2: str | None = None
    pos2: int | None = None
    genes: set = field(default_factory=set, repr=False)
    calls: dict = field(default_factory=dict, repr=False)  # VCF sample -> (alt reads, depth) the caller reported

    @property
    def is_sv(self) -> bool:
        return self.kind in ("bnd", "del", "dup", "inv")

    @property
    def label(self) -> str:
        if self.kind == "snv":
            return f"{self.chrom}:{self.pos:,} {self.ref}>{self.alt}"
        if self.kind in ("insertion", "deletion"):
            n = abs(len(self.alt) - len(self.ref))
            return f"{self.chrom}:{self.pos:,} {n} bp {self.kind}"
        if self.kind == "bnd":
            return f"{self.chrom}:{self.pos:,} <-> {self.chrom2}:{self.pos2:,} breakend"
        return f"{self.chrom}:{self.pos:,}-{self.pos2:,} {self.kind.upper()}"


# ── Loading and selecting variants ─────────────────────────────────────────────

GENE_KEYS = ("GENE", "GENEINFO", "SYMBOL",  # then ANNOVAR fields:
             "Gene.refGene", "Gene.refGeneWithVer", "Gene.ensGene", "Gene.knownGene")


def _info_genes(rec, csq_symbol: int | None, ann_gene: int | None) -> set:
    genes = set()
    for key in GENE_KEYS:
        if key in rec.info:
            val = rec.info[key]
            for g in (val if isinstance(val, tuple) else (val,)):
                # ANNOVAR joins genes with ';' written as the escape \x3b; others use | or ,
                for x in re.split(r"\\x3b|\\x2c|[;|,]", str(g)):
                    if x:
                        genes.add(x.split(":")[0])
    for key, idx in (("CSQ", csq_symbol), ("ANN", ann_gene)):
        if idx is not None and key in rec.info:
            for entry in rec.info[key]:
                parts = str(entry).split("|")
                if len(parts) > idx and parts[idx]:
                    genes.add(parts[idx])
    return {g.upper() for g in genes if g and g != "."}


def _format_index(header, key: str, name: str) -> int | None:
    if key not in header.info:
        return None
    m = re.search(r"Format:\s*([^\"]+)", header.info[key].description or "")
    if not m:
        return None
    cols = [c.strip() for c in m.group(1).replace("'", "").split("|")]
    return cols.index(name) if name in cols else None


def _raw_info(rec, key: str) -> str | None:
    """An INFO value read from the record text: pysam hides END when it is below POS (translocations)."""
    for item in str(rec).split("\t")[7].split(";"):
        k, _, v = item.partition("=")
        if k == key:
            return v
    return None


def _classify(rec, alt: str) -> tuple[str | None, str | None, int | None, str]:
    """(kind, chrom2, pos2, reason-if-not-evaluable)."""
    ref = rec.ref
    if alt in ("*", ".", None):
        return None, None, None, "no alternate allele"
    m = re.search(r"[\[\]]([^:\[\]]+):(\d+)[\[\]]", alt)
    if m:
        return "bnd", m.group(1), int(m.group(2)), ""
    if alt.startswith("<"):
        svtype = alt.strip("<>").split(":")[0].upper()
        if svtype in ("TRA", "BND", "CTX"):  # symbolic translocation (SURVIVOR, Delly): partner in CHR2/END
            c2, p2 = rec.info.get("CHR2"), _raw_info(rec, "END")
            if c2 and p2 and p2.isdigit():
                return "bnd", str(c2), int(p2), ""
            return None, None, None, f"symbolic {alt} without CHR2/END for the partner breakpoint"
        if svtype == "INV" and rec.stop - rec.pos >= MIN_INV_LEN:  # orientation, not distance, marks it
            return "inv", rec.chrom, rec.stop, ""
        if svtype in ("DEL", "DUP") and rec.stop - rec.pos >= MIN_SV_LEN:
            return svtype.lower(), rec.chrom, rec.stop, ""
        return None, None, None, (f"symbolic {alt} (only translocations and <DEL>/<DUP>/<INV> of at least "
                                  f"{MIN_SV_LEN} bp are evaluated)")
    if len(ref) == len(alt) == 1:
        return "snv", None, None, ""
    if len(ref) == len(alt):
        return None, None, None, "multi-base substitution (MNV) not supported"
    if ref[0] != alt[0]:
        return None, None, None, "complex indel without a shared anchor base"
    return ("insertion" if len(alt) > len(ref) else "deletion"), None, None, ""


def _short_id(chrom, pos, ref, alt, kind) -> str:
    """Readable ID for records without one: chr2:209694772:del28, chr19:2214506:G>A, chr9:22007648:bnd."""
    if kind == "snv":
        return f"{chrom}:{pos}:{ref}>{alt}"
    if kind in ("insertion", "deletion"):
        return f"{chrom}:{pos}:{kind[:3]}{abs(len(alt) - len(ref))}"
    return f"{chrom}:{pos}:{kind or 'other'}"


def _caller_counts(rec, alt_index: int) -> dict:
    """Per-sample (alt, depth) as the caller reported them: FORMAT AD, else AF x DP. Missing -> absent."""
    out = {}
    for name, smp in rec.samples.items():
        try:
            ad = smp.get("AD")
        except (KeyError, ValueError):
            ad = None
        if ad and len(ad) > alt_index + 1 and ad[alt_index + 1] is not None:
            depth = sum(x for x in ad if x is not None)
            out[name] = (ad[alt_index + 1], depth)
            continue
        try:
            af, dp = smp.get("AF"), smp.get("DP")
        except (KeyError, ValueError):
            af = dp = None
        af = af[alt_index] if isinstance(af, tuple) and len(af) > alt_index else af
        if isinstance(af, (int, float)) and isinstance(dp, int) and dp:
            out[name] = (round(af * dp), dp)
    return out


TOOL_SIGNS = [  # (substring in the VCF header, tool name); first match wins
    ("ID=Mutect2", "Mutect2"), ("source=Mutect2", "Mutect2"), ("source=FilterMutectCalls", "Mutect2"),
    ("source=SURVIVOR", "SURVIVOR"), ("##SURVIVOR", "SURVIVOR"), ("DRAGENCommandLine", "DRAGEN"),
    ("configManta", "Manta"), ("GenerateSVCandidates", "Manta"), ("purpleVersion", "PURPLE"),
    ("source=strelka", "Strelka"), ("source=DELLY", "Delly"), ("svaba", "SvABA"),
    ("HaplotypeCaller", "HaplotypeCaller"), ("source=freeBayes", "FreeBayes"), ("GRIDSS", "GRIDSS"),
    ("source=ESVEE", "ESVEE"), ("esveeVersion", "ESVEE")]


def vcf_tool(vcf: Path) -> str:
    """The tool that made the calls, from the VCF header; the file name when the header does not say."""
    pysam.set_verbosity(0)
    try:
        with pysam.VariantFile(str(vcf)) as vf:
            header = str(vf.header)
    except (OSError, ValueError):
        return Path(vcf).name
    low = header.lower()
    for sign, tool in TOOL_SIGNS:
        if sign.lower() in low:
            return tool + (" (ANNOVAR-annotated)" if "annovar" in low else "")
    return Path(vcf).name


def vcf_samples(vcf: Path) -> list[str]:
    pysam.set_verbosity(0)
    with pysam.VariantFile(str(vcf)) as vf:
        return list(vf.header.samples)


CHROM_COLS, POS_COLS = ("chrom", "chr", "#chrom", "contig", "chromosome"), ("pos", "start", "position")


def _annovar_to_vcf(chrom, pos, ref, alt):
    """ANNOVAR writes indels without the VCF anchor base: a deletion starts at the first deleted base
    (alt '-'), an insertion sits after `start` (ref '-'). Return the VCF-style (chrom, POS, ref, alt);
    ref/alt become None for indels, which are then matched on position alone."""
    if alt == "-":
        return chrom, pos - 1, None, None
    if ref == "-":
        return chrom, pos, None, None
    return chrom, pos, ref, alt


GENE_COLS = ("gene", "genes", "symbol", "gene_name", "gene.refgene", "gene.refgenewithver", "gene_symbol")


def _read_variant_list(path: Path) -> list[tuple]:
    """chrom pos [ref alt] rows, or any table with a header naming chrom/chr and pos/start (and ref/alt),
    such as an ANNOVAR-derived table; ANNOVAR's '-' indel notation is converted to VCF positions."""
    rows, cols = [], None
    for line in Path(path).read_text().splitlines():
        f = line.rstrip("\n").split("\t") if "\t" in line else line.split()
        f = [x.strip() for x in f]
        if not f or not f[0] or (f[0].startswith("#") and f[0].lower() != "#chrom"):
            continue
        low = [x.lower() for x in f]
        if cols is None and not (len(f) > 1 and f[1].isdigit()):
            ci = next((low.index(c) for c in CHROM_COLS if c in low), None)
            pi = next((low.index(c) for c in POS_COLS if c in low), None)
            if ci is not None and pi is not None:  # a header line: use its column names
                gi = next((low.index(c) for c in GENE_COLS if c in low), None)
                cols = (ci, pi, low.index("ref") if "ref" in low else None, low.index("alt") if "alt" in low else None,
                        gi)
            continue
        ci, pi, ri, ai, gi = cols or (0, 1, 2, 3, None)
        if len(f) <= max(ci, pi) or not f[pi].isdigit():
            continue
        ref = f[ri] if ri is not None and len(f) > ri else None
        alt = f[ai] if ai is not None and len(f) > ai else None
        gene = f[gi] if gi is not None and len(f) > gi and f[gi] not in ("", ".", "-") else None
        rows.append(_annovar_to_vcf(f[ci], int(f[pi]), ref, alt) + (gene,))
    if not rows:
        raise InputError(f"--variants file {path} has no rows (expected: chrom pos [ref alt])")
    return rows


def _read_regions(path: Path) -> dict:
    """BED (chrom, 0-based start, end[, name]) -> {chrom: [(start1, end, name)]} with 1-based inclusive starts."""
    regions = {}
    for line in Path(path).read_text().splitlines():
        f = line.split("\t") if "\t" in line else line.split()
        if not f or f[0].startswith(("#", "track", "browser")) or len(f) < 3:
            continue
        try:
            start, end = int(f[1]) + 1, int(f[2])
        except ValueError:
            continue
        regions.setdefault(f[0], []).append((start, end, f[3] if len(f) > 3 else f"{f[0]}:{start}-{end}"))
    if not regions:
        raise InputError(f"--regions file {path} has no BED rows (chrom start end [name])")
    return regions


def _in_regions(regions: dict, chrom, pos) -> set:
    return {name for a, b, name in regions.get(chrom, ()) if a <= pos <= b} if chrom and pos else set()


def _overlapping_regions(regions: dict, chrom, start, end) -> set:
    return {name for a, b, name in regions.get(chrom, ()) if a <= end and start <= b} if chrom else set()


def load_variants(vcf: Path, genes: set | None = None, variants_tsv: Path | None = None,
                  pass_only: bool = False, regions: Path | None = None,
                  allow_empty: bool = False) -> tuple[list[Variant], list[dict]]:
    """Parse and select variants. Returns (variants to validate, records not evaluated).

    Streams the VCF and keeps only selected records, so whole-genome VCFs of several GB are fine.
    Selectors combine with AND: --variants (positions), --regions (BED; region names become gene
    labels, and an SV matches if either breakpoint is inside), --genes (VCF gene names or region names).
    """
    vcf = Path(vcf)
    if not vcf.exists():
        raise InputError(f"VCF not found: {vcf}")
    if vcf.stat().st_size == 0:
        raise InputError(f"VCF is empty: {vcf}")
    pysam.set_verbosity(0)
    wanted = _read_variant_list(variants_tsv) if variants_tsv else None
    wanted_pos = {(w[0], w[1]) for w in wanted} if wanted else set()
    bed = _read_regions(regions) if regions else None
    genes = {g.upper() for g in genes} if genes else None

    variants, skipped, seen_bnd, found = [], [], set(), set()
    n_records, any_gene_names, genes_seen = 0, False, set()
    try:
        vf = pysam.VariantFile(str(vcf))
        csq_i = _format_index(vf.header, "CSQ", "SYMBOL")
        ann_i = _format_index(vf.header, "ANN", "Gene_Name")
        for rec in vf:
            n_records += 1
            if wanted is not None and (rec.chrom, rec.pos) not in wanted_pos and not any(
                    c in wanted_pos for c in _record_partners(rec)):
                continue  # fast path: not a listed position
            rec_genes = _info_genes(rec, csq_i, ann_i)
            any_gene_names |= bool(rec_genes)
            filt = ";".join(rec.filter.keys()) or "PASS"
            for i, alt in enumerate(rec.alts or ()):
                kind, c2, p2, why = _classify(rec, alt)
                labels = set(rec_genes)
                if bed is not None:
                    hit = _in_regions(bed, rec.chrom, rec.pos) | _in_regions(bed, c2, p2)
                    if kind in ("del", "dup", "inv") and c2 == rec.chrom:
                        # the span counts, as in AnnotSV: a deletion removing a whole gene has both
                        # breakpoints outside it. Translocations still need a breakpoint in the region.
                        hit |= _overlapping_regions(bed, rec.chrom, min(rec.pos, p2), max(rec.pos, p2))
                    if not hit:
                        continue  # outside every region
                    labels |= {x.upper() for x in hit}
                if genes is not None:
                    genes_seen |= labels & genes
                    if not labels & genes:
                        continue
                vid = rec.id or _short_id(rec.chrom, rec.pos, rec.ref, alt, kind)
                if len(rec.alts) > 1:
                    vid = f"{vid}_{i + 1}"
                if wanted is not None:
                    hits = [w for w in wanted if (rec.chrom == w[0] and rec.pos == w[1] and
                            (w[2] is None or (w[2] == rec.ref and w[3] == alt))) or (c2 == w[0] and p2 == w[1])]
                    if not hits:
                        continue
                    found.update(hits)
                    for w in hits:  # gene names from the list (e.g. ANNOVAR's column) label unannotated VCFs
                        if w[4]:
                            labels |= {g.upper() for g in re.split(r"[;,]", w[4]) if g}
                if kind is None:
                    skipped.append({"id": vid, "chrom": rec.chrom, "pos": rec.pos, "reason": why})
                    continue
                if kind == "bnd":  # mate breakends describe one event: keep the first
                    key = frozenset({(rec.chrom, rec.pos), (c2, p2)})
                    if key in seen_bnd:
                        continue
                    seen_bnd.add(key)
                calls = {} if kind in ("bnd", "del", "dup", "inv") else _caller_counts(rec, i)
                variants.append(Variant(vid, rec.chrom, rec.pos, rec.ref, alt, kind, ",".join(sorted(labels)),
                                        filt, c2, p2, labels, calls))
    except (OSError, ValueError) as e:
        raise InputError(f"could not read {vcf} as a VCF ({e})") from None
    if not n_records:
        raise InputError(f"no variant records in {vcf}")

    if wanted is not None:
        for w in wanted:
            if w not in found:
                skipped.append({"id": f"{w[0]}:{w[1]}", "chrom": w[0], "pos": w[1], "reason": "listed but not in the VCF"})
    if genes is not None:
        if not any_gene_names and bed is None:
            raise InputError("--genes needs gene names in the VCF (VEP CSQ, SnpEff ANN, ANNOVAR Gene.refGene or "
                             "INFO/GENE); use --regions with a BED of gene coordinates, or --variants")
        for g in sorted(genes - genes_seen):
            skipped.append({"id": g, "chrom": "", "pos": "", "reason": "gene requested but no variants in the VCF"})
    if pass_only:
        variants = [v for v in variants if v.filter == "PASS"]
    if not variants:
        if allow_empty:  # e.g. no SV touches the selected genes: a result, reported as such
            return [], skipped
        if genes is not None and not genes_seen:
            raise InputError(f"no variants in the VCF for gene(s): {', '.join(sorted(genes))}")
        raise InputError("no variants left to validate after selection")
    return variants, skipped


def _record_partners(rec) -> list[tuple]:
    """Partner breakpoint(s) of an SV record, for matching a listed position to either end."""
    out = []
    for alt in rec.alts or ():
        m = re.search(r"[\[\]]([^:\[\]]+):(\d+)[\[\]]", alt)
        if m:
            out.append((m.group(1), int(m.group(2))))
        elif alt.startswith("<") and "CHR2" in rec.info and (_raw_info(rec, "END") or "").isdigit():
            out.append((str(rec.info["CHR2"]), int(_raw_info(rec, "END"))))
        elif alt.startswith("<"):
            out.append((rec.chrom, rec.stop))
    return out


# ── Counting read support ──────────────────────────────────────────────────────

def _ok(r) -> bool:
    return not (r.is_unmapped or r.is_duplicate or r.is_secondary or r.is_supplementary
                or r.is_qcfail) and r.mapping_quality >= MIN_MAPQ


def _open_bam(path, fasta=None):
    return pysam.AlignmentFile(str(path), reference_filename=str(fasta) if fasta else None)


def _count_snv(bam, v: Variant) -> dict:
    depth, hits, other = 0, [], {}
    for col in bam.pileup(v.chrom, v.pos - 1, v.pos, truncate=True, min_base_quality=MIN_BQ,
                          stepper="nofilter", ignore_orphans=False, max_depth=1_000_000):
        for pr in col.pileups:  # overlapping mates: pysam zeroes one base quality, so a molecule counts once
            r = pr.alignment
            if not _ok(r) or pr.is_del or pr.is_refskip or pr.query_position is None:
                continue
            depth += 1
            qp = pr.query_position
            base = r.query_sequence[qp]
            if base == v.alt:
                hits.append((r.is_reverse, min(qp, r.query_length - 1 - qp), r.mapping_quality))
            elif base != v.ref:
                other[base] = other.get(base, 0) + 1
    return {"depth": depth, "hits": hits, "other": max(other.values(), default=0)}


def _mismatches(a: str, b: str) -> int:
    return sum(x != y for x, y in zip(a.upper(), b.upper())) + abs(len(a) - len(b))


def _clip_supports(r, v: Variant, fa) -> bool:
    """Soft-clipped exactly at the indel edge, with clipped bases matching ALT better than REF."""
    if fa is None or not r.cigartuples:
        return False
    pos, size = v.pos, len(v.alt) - len(v.ref)  # pos: 0-based index of the first base after the anchor
    L = abs(size)
    ins = v.alt[1:] if size > 0 else ""
    op0, n0 = r.cigartuples[0]
    op1, n1 = r.cigartuples[-1]
    q = r.query_sequence
    if op1 == 4 and n1 >= MIN_CLIP and r.reference_end == pos:          # right clip at the anchor
        k = min(n1, 20)
        clip = q[-n1:][:k]
        alt_seq = (ins + fa.fetch(v.chrom, pos + (L if size < 0 else 0), pos + (L if size < 0 else 0) + k))[:k]
        ref_seq = fa.fetch(v.chrom, pos, pos + k)
    elif op0 == 4 and n0 >= MIN_CLIP and r.reference_start == pos + (L if size < 0 else 0):  # left clip after it
        k = min(n0, 20)
        clip = q[:n0][-k:]
        alt_seq = (fa.fetch(v.chrom, max(0, pos - k), pos) + ins)[-k:]
        ref_seq = fa.fetch(v.chrom, max(0, pos + (L if size < 0 else 0) - k), pos + (L if size < 0 else 0))
    else:
        return False
    m_alt, m_ref = _mismatches(clip, alt_seq), _mismatches(clip, ref_seq)
    return m_alt <= max(1, int(0.15 * k)) and m_alt < m_ref


def _count_indel(bam, v: Variant, fa) -> dict:
    size = len(v.alt) - len(v.ref)  # negative = deletion
    pos = v.pos                      # 0-based coordinate of the first inserted/deleted base
    depth, hits = set(), {}
    # count molecules, not reads: overlapping mates would otherwise count one DNA fragment twice
    for r in bam.fetch(v.chrom, max(0, pos - 1), pos + abs(size) + 1):
        if not _ok(r):
            continue
        spans = r.reference_start <= pos - 1 and r.reference_end >= pos + 1
        rp, qp, hit = r.reference_start, 0, None
        for op, n in r.cigartuples:
            if op in (0, 7, 8):
                rp += n; qp += n
            elif op == 1:
                if size > 0 and n == size and abs(rp - pos) <= INDEL_WINDOW:
                    hit = qp; break
                qp += n
            elif op == 2:
                if size < 0 and n == -size and abs(rp - pos) <= INDEL_WINDOW:
                    hit = qp; break
                rp += n
            elif op == 4:
                qp += n
        if hit is None and _clip_supports(r, v, fa):
            clip = r.cigartuples[0][1] if r.cigartuples[0][0] == 4 and r.reference_start > pos - 1 else r.cigartuples[-1][1]
            hit = r.query_length - clip if r.reference_end == pos else clip
        if spans or hit is not None:
            depth.add(r.query_name)
        if hit is not None:
            hits.setdefault(r.query_name, (r.is_reverse, min(hit, r.query_length - hit), r.mapping_quality))
    return {"depth": len(depth), "hits": list(hits.values())}


def _sv_side(bam, chrom, pos, chrom2, pos2, inversion: bool = False):
    """Split reads and discordant pairs seen from one breakpoint, as sets of read names.

    Discordant pairs: mate near the partner breakpoint and, on the same chromosome, either both mates
    on the same strand (inversions, any size) or an improper pair with a long insert (DEL/DUP >= 1 kb).
    """
    split, disc, depth = set(), set(), 0
    same_chrom = chrom2 == chrom
    for r in bam.fetch(chrom, max(0, pos - SPLIT_WINDOW), pos + SPLIT_WINDOW):
        if not _ok(r):
            continue
        if r.reference_start <= pos <= r.reference_end:
            depth += 1
        if r.has_tag("SA"):
            for sa in r.get_tag("SA").strip(";").split(";"):
                ch, sp = sa.split(",")[:2]
                if ch == chrom2 and abs(int(sp) - pos2) < SPLIT_WINDOW:
                    split.add(r.query_name)
        if r.is_paired and not r.mate_is_unmapped and r.next_reference_name == chrom2 \
                and abs(r.next_reference_start - pos2) < PAIR_WINDOW:
            if not same_chrom:
                disc.add(r.query_name)
            elif inversion:
                if r.is_reverse == r.mate_is_reverse:
                    disc.add(r.query_name)
            elif not r.is_proper_pair and abs(r.template_length) >= MIN_SV_LEN // 2:
                disc.add(r.query_name)
    return split, disc, depth


def _count_sv(bam, v: Variant) -> dict:
    # Look from both breakpoints and merge by read name: each end sees a different subset of the
    # split reads, and a molecule seen from both ends must still count once.
    inv = v.kind == "inv"
    s1, d1, depth = _sv_side(bam, v.chrom, v.pos, v.chrom2, v.pos2, inv)
    s2, d2, _ = _sv_side(bam, v.chrom2, v.pos2, v.chrom, v.pos, inv)
    split, disc = s1 | s2, (d1 | d2) - (s1 | s2)
    return {"depth": depth, "split": len(split), "discordant": len(disc), "union": len(split | disc)}


def count_variant(v: Variant, bam_path, fasta=None) -> dict:
    """Tumor-or-normal support for one variant, from the BAM."""
    bam = _open_bam(bam_path, fasta)
    fa = pysam.FastaFile(str(fasta)) if fasta else None
    try:
        if v.is_sv:
            c = _count_sv(bam, v)
            return {"depth": c["depth"], "alt": c["union"], "vaf_pct": _pct(c["union"], c["depth"]),
                    "split": c["split"], "discordant": c["discordant"],
                    "alt_fwd": None, "alt_rev": None, "alt_near_read_end": None, "alt_mean_mapq": None,
                    "other_allele_pct": None}
        c = _count_snv(bam, v) if v.kind == "snv" else _count_indel(bam, v, fa)
        h = c["hits"]
        return {"depth": c["depth"], "alt": len(h), "vaf_pct": _pct(len(h), c["depth"]),
                "alt_fwd": sum(not x[0] for x in h), "alt_rev": sum(x[0] for x in h),
                "alt_near_read_end": sum(x[1] < READ_END_BP for x in h),
                "alt_mean_mapq": round(sum(x[2] for x in h) / len(h)) if h else None,
                "other_allele_pct": _pct(c["other"], c["depth"]) if "other" in c else None,
                "split": None, "discordant": None}
    finally:
        bam.close()
        if fa:
            fa.close()


def _pct(n, d):
    return round(100 * n / d, 1) if d else None


def flag_variant(v: Variant, t: dict, n: dict | None, caller: dict | None = None) -> tuple[list[str], str]:
    """Rule-based flags and status. `n` is None in tumor-only mode (normal-based checks are skipped)."""
    if t["depth"] == 0 and (n is None or n["depth"] == 0):  # nothing to judge: not "weak evidence"
        return ["no_coverage"], "no_coverage"
    flags = []
    if n is not None and n["alt"] > 0:
        flags.append("normal_support")
    if t["alt"] < MIN_SUPPORT:
        flags.append("low_support")
    if t["depth"] and 100 * t["alt"] / t["depth"] < MIN_VAF:
        flags.append("low_vaf")
    if not v.is_sv and t["alt"] >= STRAND_MIN and min(t["alt_fwd"], t["alt_rev"]) == 0:
        flags.append("strand_bias")
    if not v.is_sv and t["alt"] and t["alt_near_read_end"] / t["alt"] > 0.5:
        flags.append("read_end")
    if t["depth"] < LOW_DEPTH or (n is not None and n["depth"] < LOW_DEPTH):
        flags.append("low_depth")
    if n is not None and (n.get("other_allele_pct") or 0) >= GERMLINE_OTHER:
        flags.append("germline_site")
    if caller_disagrees(t, caller):
        flags.append("caller_disagrees")
    status = "insufficient" if t["alt"] < MIN_SUPPORT else ("flagged" if flags else "supported")
    return flags, status


def caller_disagrees(t: dict, caller: dict | None) -> bool:
    """The caller's own counts (VCF AD) and the BAM counts tell different stories."""
    if not caller or not caller.get("depth"):
        return False
    ours = 100 * t["alt"] / t["depth"] if t["depth"] else 0.0
    return abs(ours - caller["vaf_pct"]) > CALLER_VAF_DIFF or (caller["alt"] >= MIN_SUPPORT and t["alt"] == 0)


def flag_depth_outliers(results: dict) -> None:
    """Flag variants whose tumor or normal depth is far above the run's median (needs >= 3 variants)."""
    if len(results) < 3:
        return
    for role in ("tumor", "normal"):
        if any(r[role] is None for r in results.values()):
            continue  # tumor-only run
        depths = sorted(r[role]["depth"] for r in results.values())
        median = depths[len(depths) // 2]
        for r in results.values():
            if median and r[role]["depth"] > HIGH_DEPTH_X * median and "high_depth" not in r["flags"]:
                r["flags"].append("high_depth")
                if r["status"] == "supported":
                    r["status"] = "flagged"


# ── IGV screenshots ────────────────────────────────────────────────────────────

MAC_APP_DIRS = [Path("/Applications")]


def find_igv(explicit: str | None) -> tuple[list[str] | None, Path | None, str]:
    """(launcher argv prefix, working dir, reason if unavailable).

    Looks for an IGV .app (macOS), then an `igv.sh` or `igv` command on PATH (Linux, HPC modules).
    `explicit` may be a path or a command name.
    """
    if explicit:
        p = Path(explicit)
        cands = [p if p.exists() or not shutil.which(explicit) else Path(shutil.which(explicit))]
    else:
        cands = [a for d in MAC_APP_DIRS for a in sorted(d.glob("IGV*.app"), reverse=True)]
        cands += [Path(shutil.which(c)) for c in ("igv.sh", "igv") if shutil.which(c)]
    for c in cands:
        if c.suffix == ".app" and c.is_dir():
            java = sorted((c / "Contents").glob("jdk*/bin/java"))
            if java:
                return ([str(java[-1]), "--module-path=Java/lib", "-Xmx4g", "@Java/igv.args", "-Xdock:name=IGV",
                         "--module=org.igv/org.broad.igv.ui.Main"], c / "Contents", "")
        elif c.is_file():
            return [str(c)], None, ""
    return None, None, (f"IGV not found at {explicit}" if explicit else
                        "IGV desktop not found: load it first (e.g. `module load igv` on an HPC), install it, "
                        "or pass --igv-path")


def igv_location(explicit: str | None) -> str | None:
    """The IGV this run would use (.app bundle or command), recorded so --summarize can reuse it."""
    argv, cwd, _ = find_igv(explicit)
    if argv is None:
        return None
    return str(cwd.parent) if cwd else str(Path(argv[0]).resolve())


def overview_launcher(explicit: str | None, plan: list[dict]) -> str | None:
    """--igv-path if given, else the first IGV recorded by the runs that still exists (None: search as usual)."""
    if explicit:
        return explicit
    return next((p["igv"] for p in plan if p.get("igv") and Path(p["igv"]).exists()), None)


def launcher_failed(returncode: int | None) -> bool:
    """A launcher that exited with 0 may have started IGV in the background (HPC module wrappers)."""
    return returncode is not None and returncode != 0


def _free_port() -> int:
    with socket.socket() as s:
        s.bind(("127.0.0.1", 0))
        return s.getsockname()[1]


class IGVSession:
    """An IGV of our own: free port, temporary settings folder, closed on exit.

    Never talks to an IGV the user already has open, and never touches their IGV settings.
    """

    def __init__(self, argv, cwd, genome: Path, timeout: int = 300, log_copy: Path | None = None):
        self.port = _free_port()
        self.can_query_genome = True  # IGV <= 2.16 has no currentGenomePath; then the log is used instead
        self.genome = Path(genome)
        self.home = Path(tempfile.mkdtemp(prefix="igv_validator_"))
        # IGV opens its command port before it has loaded a genome, then loads its default one.
        # Make ours the default (so exactly one genome is ever loaded) and wait until IGV reports
        # it as current; sending anything earlier races the load and can leave IGV broken.
        (self.home / "prefs.properties").write_text(f"DEFAULT_GENOME_KEY={self.genome}\n")
        self.log = open(self.home / "igv.log", "w")
        self.proc = subprocess.Popen([*argv, "--port", str(self.port), "--igvDirectory", str(self.home),
                                      "-g", str(genome)], cwd=cwd, stdout=self.log, stderr=subprocess.STDOUT)
        self.f = None
        deadline = time.time() + timeout
        while time.time() < deadline:
            if launcher_failed(self.proc.poll()):
                self._fail(f"IGV exited during start-up", log_copy)
            try:
                if self.f is None:
                    self._connect()
                reply = self._raw("currentGenomePath")
                if reply and "COMMAND" in reply.upper():  # 'UNKOWN COMMAND' (sic) from IGV <= 2.16
                    self.can_query_genome = False
                if self.can_query_genome and reply and Path(reply) == self.genome:
                    return
                if not self.can_query_genome and str(self.genome) in self.genomes_loaded():
                    time.sleep(2)  # let the load finish; every image is re-checked for reads anyway
                    return
            except OSError:
                self.f = None
            time.sleep(0.5)
        self._fail(f"IGV did not load {self.genome.name} within {timeout} s (raise --igv-timeout on slow "
                   "file systems)", log_copy)

    def _fail(self, msg: str, log_copy: Path | None):
        """Keep IGV's own log before cleaning up, then raise."""
        if log_copy:
            try:
                self.log.flush()
                log_copy.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy(self.home / "igv.log", log_copy)
                msg += f" (IGV log: {log_copy})"
            except OSError:
                pass
        self.close()
        raise RuntimeError(msg)

    def _connect(self):
        s = socket.create_connection(("127.0.0.1", self.port), timeout=5)
        s.settimeout(60)  # a command that hangs longer than this means IGV is stuck
        self.f = s.makefile("rw")

    def _raw(self, cmd: str) -> str | None:
        """Send one command; None when IGV dropped the connection or never answered."""
        try:
            self.f.write(cmd + "\n"); self.f.flush()
            line = self.f.readline()
        except OSError:
            return None
        return line.strip() if line else None

    def current_genome(self) -> Path | None:
        if not self.can_query_genome:  # IGV <= 2.16: the final log check guards against a genome switch
            return self.genome
        reply = self._raw("currentGenomePath")
        return Path(reply) if reply else None

    def send(self, cmd: str, retry: bool = True) -> str:
        reply = self._raw(cmd)
        if reply is None and not launcher_failed(self.proc.poll()):
            # IGV occasionally drops the command connection: reconnect, and resend when the command
            # is safe to repeat (every command here is, except ones that may have crashed IGV's handler)
            time.sleep(1)
            self._connect()
            reply = self._raw(cmd) if retry else None
        if reply != "OK":
            raise RuntimeError(f"IGV command failed: {cmd} -> {reply or 'no reply'}")
        return reply

    def genomes_loaded(self) -> set[str]:
        """Every genome IGV has loaded so far, from its log."""
        self.log.flush()
        text = (self.home / "igv.log").read_text(errors="replace")
        return set(re.findall(r"Loading genome:\s*(\S+)", text))

    def close(self):
        if getattr(self, "_closed", False):
            return
        self._closed = True
        try:
            if self.f:
                self.f.write("exit\n"); self.f.flush()
        except OSError:
            pass
        try:
            self.proc.wait(timeout=15)
        except subprocess.TimeoutExpired:
            self.proc.kill()
        # a background-starting wrapper has already exited: wait for IGV itself to release its port
        for _ in range(15):
            try:
                socket.create_connection(("127.0.0.1", self.port), timeout=1).close()
                time.sleep(1)
            except OSError:
                break
        self.log.close()
        shutil.rmtree(self.home, ignore_errors=True)


def _loci(v: Variant) -> list[tuple[str, str, str]]:
    """(tag, goto locus, sort position) for each image of a variant."""
    if v.is_sv:
        return [("bp1", f"{v.chrom}:{v.pos - 40}-{v.pos + 40}", f"{v.chrom}:{v.pos}"),
                ("bp2", f"{v.chrom2}:{v.pos2 - 40}-{v.pos2 + 40}", f"{v.chrom2}:{v.pos2}")]
    span = max(len(v.ref), len(v.alt))
    start, end = v.pos - 20, v.pos + span + 20
    # sort on the anchor base: inside a deletion the deleted reads have no base and IGV's sort crashes
    return [("", f"{v.chrom}:{start}-{end}", f"{v.chrom}:{v.pos}")]


def _png_name(v: Variant, tag: str) -> str:
    safe = re.sub(r"[^A-Za-z0-9._-]+", "_", f"{v.id}_{v.gene.split(',')[0] or v.chrom}")
    return f"{safe}{'_' + tag if tag else ''}.png"


ANNOT_PAD = 50_000   # bp of annotation kept around each screenshot window


def _other_chr_name(c: str) -> str:
    return c[3:] if c.startswith("chr") else "chr" + c


def annotation_subset(path, windows: list[tuple[str, int, int]], out_stem: Path) -> Path | None:
    """The lines of a GTF/GFF/BED gene file (optionally .gz) that overlap the screenshot windows, with chromosome
    names rewritten to the BAM's style (chr2 vs 2). Returns the small file IGV loads, or None if nothing overlaps."""
    import gzip
    path = Path(path)
    name = path.name[:-3] if path.name.endswith(".gz") else path.name
    ext = next((e for e in (".gff3", ".gff", ".gtf", ".bed") if name.lower().endswith(e)), None)
    if ext is None:
        raise InputError(f"--annotation must be a .gtf, .gff/.gff3 or .bed file (optionally .gz): {path}")
    bed = ext == ".bed"
    wins: dict[str, list[tuple[int, int]]] = {}
    for c, a, b in windows:
        for key in (c, _other_chr_name(c)):
            wins.setdefault(key, []).append((max(0, a - ANNOT_PAD), b + ANNOT_PAD))
    bam_name = {k: c for c, _, _ in windows for k in (c, _other_chr_name(c))}
    keep = []
    with (gzip.open(path, "rt") if path.name.endswith(".gz") else open(path)) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < (3 if bed else 5) or f[0] not in wins:
                continue
            if not bed and f[2] == "gene":   # transcripts and exons draw the gene; a gene line only repeats it
                continue
            try:
                a, b = (int(f[1]), int(f[2])) if bed else (int(f[3]), int(f[4]))
            except ValueError:
                continue
            if any(a <= wb and b >= wa for wa, wb in wins[f[0]]):
                f[0] = bam_name[f[0]]
                keep.append("\t".join(f))
    if not keep:
        return None
    out = Path(str(out_stem) + ext)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(keep) + "\n")
    return out


def annotation_cmds(sub: Path | None) -> list[str]:
    return [f"load {Path(sub).resolve()} name=genes", "expand genes"] if sub else []


def igv_commands(v: Variant, tumor: Path, normal: Path | None, snapdir: Path, image: int = 0,
                 tumor_name: str = "tumor", normal_name: str | None = "normal",
                 annotation: Path | None = None) -> list[str]:
    """IGV batch commands for one image; a tumor-only run (normal None) loads a single track."""
    tag, locus, sortpos = _loci(v)[image]
    return ["new", f"snapshotDirectory {snapdir}",
            "preference SAM.DOWNSAMPLE_READS false", "preference SAM.ALLELE_THRESHOLD 0.02",
            f"preference SAM.SHOW_SOFT_CLIPPED {'false' if v.kind == 'snv' else 'true'}",
            "maxPanelHeight 900",
            f"load {Path(tumor).resolve()} name={tumor_name}",
            *([f"load {Path(normal).resolve()} name={normal_name or 'normal'}"] if normal else []),
            *annotation_cmds(annotation),
            f"goto {locus}", f"sort base {sortpos}", "expand", f"snapshot {_png_name(v, tag)}"]


def image_rendered(png: Path) -> bool:
    """False when IGV saved before drawing, or drew no reads.

    Checks that the ruler/header at the top is drawn, and that read-grey pixels fill part of the
    track area (a whole-genome or "zoom in to see alignments" view has none).
    """
    from PIL import Image
    im = Image.open(png).convert("L")
    w, h = im.size
    ruler = im.crop((int(w * 0.15), 0, int(w * 0.9), min(100, h))).histogram()
    tracks = im.crop((int(w * 0.15), min(150, h - 1), int(w * 0.9), h)).histogram()
    area = max(1, (int(w * 0.9) - int(w * 0.15)) * max(1, h - 150))
    return sum(ruler[:120]) >= 500 and sum(tracks[160:215]) / area >= 0.02


def _locus_window(locus: str) -> tuple[str, int, int]:
    c, rng = locus.rsplit(":", 1)
    a, b = rng.replace(",", "").split("-")
    return c, int(a), int(b)


def take_screenshots(variants, results, tumor, normal, genome, snapdir, igv_path, names,
                     log_copy: Path | None = None, timeout: int = 300, annotation=None) -> tuple[int, str]:
    argv, cwd, why = find_igv(igv_path)
    if argv is None:
        return 0, why
    if sys.platform.startswith("linux") and not __import__("os").environ.get("DISPLAY"):
        return 0, "no display available (run with xvfb-run, or use an HPC desktop session)"
    snapdir.mkdir(parents=True, exist_ok=True)
    try:
        igv = IGVSession(argv, cwd, genome, timeout, log_copy)
    except RuntimeError as e:
        return 0, str(e)
    taken, incomplete, unsorted, note = 0, 0, 0, ""
    sub = annotation_subset(annotation, [_locus_window(l) for v in variants for _, l, _ in _loci(v)],
                            snapdir / "_genes") if annotation else None
    try:
        for v in variants:
            for i in range(len(_loci(v))):
                if igv.current_genome() != Path(genome):
                    raise RuntimeError(f"IGV is no longer on {Path(genome).name}")
                cmds = igv_commands(v, tumor, normal, snapdir, i, *names, annotation=sub)
                goto, sort, expand, snap = cmds[-4:]
                for c in cmds[:-3]:          # everything up to and including goto
                    igv.send(c)
                # IGV loads reads in the background after goto; sorting before they arrive crashes it.
                # Probe with throwaway snapshots until reads are drawn (about 10 s at most).
                probe = snapdir / "_probe.png"
                for _ in range(20):
                    igv.send(goto)  # re-sent: harmless, and covers a genome that finished loading late
                    igv.send(f"snapshot {probe.name}")
                    if probe.exists() and image_rendered(probe):
                        break
                    time.sleep(0.5)
                probe.unlink(missing_ok=True)
                try:  # sorting only groups the alt reads on top; never let it stop the run
                    igv.send(sort, retry=False)
                except RuntimeError:
                    unsorted += 1
                igv.send(expand)
                png = snapdir / snap.split(" ", 1)[1]
                for attempt in range(3):   # short pause before the snapshot, then check it rendered
                    igv.send(goto)         # goto again forces a full redraw
                    time.sleep(0.5 * (attempt + 1))
                    igv.send(snap)
                    if png.exists() and image_rendered(png):
                        break
                else:
                    incomplete += 1
                results[v.id]["screenshots"].append(png)
                taken += 1
        other = {g for g in igv.genomes_loaded() if Path(g) != Path(genome)}
        if other:  # IGV switched away from the BAMs' reference: the images show the wrong genome
            for r in results.values():
                r["screenshots"] = []
            taken, note = 0, f"IGV loaded a different genome ({', '.join(sorted(other))}), so the images were discarded"
        else:
            notes = [f"{incomplete} image(s) may be incomplete (no reads drawn after 3 tries)"] if incomplete else []
            notes += [f"{unsorted} image(s) left unsorted (IGV could not sort there)"] if unsorted else []
            note = "; ".join(notes)
    except RuntimeError as e:
        note = f"stopped early: {e}"
    finally:
        if note and log_copy:  # keep IGV's own log whenever something went wrong
            log_copy.parent.mkdir(parents=True, exist_ok=True)
            igv.log.flush()
            shutil.copy(igv.home / "igv.log", log_copy)
            note += f" (IGV log: {log_copy.name})"
        igv.close()
    return taken, note


def _font(size, bold=False):
    from PIL import ImageFont
    for p in (["/System/Library/Fonts/Helvetica.ttc"] if sys.platform == "darwin" else []) + [
            "/usr/share/fonts/truetype/dejavu/DejaVuSans%s.ttf" % ("-Bold" if bold else ""),
            "/usr/share/fonts/dejavu/DejaVuSans%s.ttf" % ("-Bold" if bold else "")]:
        if Path(p).exists():
            return ImageFont.truetype(p, size, index=1 if bold and p.endswith(".ttc") else 0)
    return ImageFont.load_default(size)


def caption(png: Path, out: Path, lines: list[str]):
    from PIL import Image, ImageDraw
    shot = Image.open(png).convert("RGB")
    fonts = [_font(24, True)] + [_font(19)] * (len(lines) - 1)
    pad, gap = 16, 8
    h = pad * 2 + sum(f.size for f in fonts) + gap * (len(lines) - 1)
    img = Image.new("RGB", (shot.width, h + shot.height), "white")
    d = ImageDraw.Draw(img)
    y = pad
    for text, f in zip(lines, fonts):
        d.text((pad, y), text, fill="black", font=f); y += f.size + gap
    d.line([(0, h - 1), (shot.width, h - 1)], fill=(180, 180, 180), width=2)
    img.paste(shot, (0, h))
    img.save(out)


# ── Reports ────────────────────────────────────────────────────────────────────

def _support_text(c: dict, sv: bool) -> str:
    if sv:
        return f"{c['alt']} ({c['split']} split, {c['discordant']} pairs) / {c['depth']}"
    if not c["depth"]:
        return f"{c['alt']}/0 (no reads)"
    return f"{c['alt']}/{c['depth']} ({c['vaf_pct']}%)"


def _strand_text(c: dict) -> str:
    return f", {c['alt_fwd']} forward / {c['alt_rev']} reverse" if c.get("alt_fwd") is not None and c["alt"] else ""


def _reason_summary(skipped: list[dict]) -> str:
    counts = {}
    for s in skipped:
        counts[s["reason"]] = counts.get(s["reason"], 0) + 1
    return "; ".join(f"{n:,} {r}" for r, n in sorted(counts.items(), key=lambda x: -x[1])[:4]) + \
        ("; ..." if len(counts) > 4 else "")


def _caller_text(c: dict | None) -> str:
    return f"{c['alt']}/{c['depth']} ({c['vaf_pct']}%)" if c else "-"


def write_reports(out: Path, rows: list[dict], skipped: list[dict], meta: dict) -> list[Path]:
    paired = meta["mode"] == "tumor_normal"
    tables = out / "tables"; tables.mkdir(parents=True, exist_ok=True)
    tsv = tables / "support_counts.tsv"
    cols = ["id", "gene", "kind", "variant", "filter", "sample", "depth", "alt", "vaf_pct", "alt_fwd", "alt_rev",
            "alt_near_read_end", "alt_mean_mapq", "split", "discordant", "caller_alt", "caller_depth", "flags", "status"]
    with open(tsv, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t"); w.writerow(cols)
        for r in rows:
            for role in ("tumor", "normal"):
                c = r[role]
                if c is None:
                    continue  # tumor-only run
                cc = r["caller"] if role == "tumor" and r["caller"] else {}
                w.writerow([r["id"], r["gene"], r["kind"], r["label"], r["filter"], meta[f"{role}_name"]] +
                           ["" if c[k] is None else c[k] for k in cols[6:15]] +
                           [cc.get("alt", ""), cc.get("depth", ""), ";".join(r["flags"]), r["status"]])

    counts = {s: sum(r["status"] == s for r in rows) for s in ("supported", "flagged", "insufficient", "no_coverage")}
    shots = meta["screenshots"]
    shot_line = (f"IGV screenshots: {shots['taken']} image(s) in `figures/igv/`." if shots["taken"] and not shots["note"]
                 else f"IGV screenshots skipped: {shots['note']}." if not shots["taken"]
                 else f"IGV screenshots: {shots['taken']} image(s); {shots['note']}.")
    md = [f"# IGV Validator Report", "",
          f"**Input**: `{meta['vcf']}` ({len(rows)} variant(s) checked) · tumor: {meta['tumor_name']} · "
          + (f"normal: {meta['normal_name']}" if paired else "**tumor-only** (no normal)"),
          f"**Date**: {meta['date']}  ", f"**Skill**: igv-validator {VERSION}", "",
          f"**{counts['supported']} supported, {counts['flagged']} flagged, {counts['insufficient']} insufficient**"
          + (f", {counts['no_coverage']} no coverage" if counts["no_coverage"] else "")
          + (f"; {len(skipped):,} not evaluated ({_reason_summary(skipped)})" if skipped else "") + ".", ""]
    if not rows:
        md += ["**No calls in the selected regions or genes.** The VCF has no call touching them "
               "(for SVs: no breakpoint inside, and no deletion/duplication/inversion spanning them).", ""]
    if not paired:
        md += ["> **Tumor-only mode:** without a matched normal the skill cannot tell a somatic mutation from an "
               "inherited (germline) variant, and the normal-based checks (support in the normal, germline site, "
               "swapped samples) are skipped.", ""]
    if meta.get("caller_note"):
        md += [f"> Caller counts: {meta['caller_note']}.", ""]
    if meta.get("swap_warning"):
        md += [f"> **Warning:** {meta['swap_warning']}", ""]
    if counts["no_coverage"]:
        md += [f"> **Warning:** {counts['no_coverage']} variant(s) have no reads in either BAM. Check that the BAMs "
               "cover these regions (a sliced BAM, a different sample, or a VCF from another reference build).", ""]
    md += [
          shot_line, "",
          "Read counts come from the BAMs (MAPQ >= 20, base quality >= 20, duplicate/secondary/supplementary "
          "reads removed, each DNA molecule counted once), never from the screenshots. **Status** is a rule-based "
          "summary of the flags, not a verdict on whether the variant is real.", "",
          "| ID | Gene | Variant | Tumor support | " + ("Normal support | " if paired else "")
          + "Caller reported (VCF) | Flags | Status |",
          "|---|---|---|---|" + ("---|" if paired else "") + "---|---|---|"]
    for r in rows:
        md.append(f"| {r['id']} | {r['gene'] or '-'} | {r['label']} | {_support_text(r['tumor'], r['sv'])} | "
                  + (f"{_support_text(r['normal'], r['sv'])} | " if paired else "")
                  + f"{_caller_text(r['caller'])} | {'; '.join(r['flags']) or 'none'} | **{r['status']}** |")
    if skipped:
        md += ["", "## Not evaluated", "", "| Record | Reason |", "|---|---|"]
        md += [f"| {s['id']} | {s['reason']} |" for s in skipped[:30]]
        if len(skipped) > 30:
            md += ["", f"...and {len(skipped) - 30} more; the full list is in `result.json` (`not_evaluated`)."]
    md += ["", "## Flags", ""] + [f"- `{k}`: {t}" for k, t in FLAG_TEXT.items()]
    md += ["", "## Variants", ""]
    for r in rows:
        md += [f"### {r['id']} {r['gene']}: {r['status']}", "",
               f"- {meta['tumor_name']}: {_support_text(r['tumor'], r['sv'])}{_strand_text(r['tumor'])}"]
        if paired:
            md.append(f"- {meta['normal_name']}: {_support_text(r['normal'], r['sv'])}{_strand_text(r['normal'])}")
        md += [f"- Caller reported (VCF): {_caller_text(r['caller'])}",
               f"- Flags: {', '.join(FLAG_TEXT[f] for f in r['flags']) or 'none'}", ""]
        md += [f"![{r['id']} IGV]({p})" for p in r["figures"]] + [""]
    md += ["---", "", "## Disclaimer", "", f"*{DISCLAIMER}*", ""]
    (out / "report.md").write_text("\n".join(md))

    e = html.escape

    def img(p):
        return "data:image/png;base64," + base64.b64encode((out / p).read_bytes()).decode()

    body = "".join(
        f"<tr><td>{e(r['id'])}</td><td><b>{e(r['gene'] or '-')}</b></td><td>{e(r['label'])}</td>"
        f"<td class=n>{e(_support_text(r['tumor'], r['sv']))}</td>"
        + (f"<td class=n>{e(_support_text(r['normal'], r['sv']))}</td>" if paired else "")
        + f"<td class=n>{e(_caller_text(r['caller']))}</td>"
        f"<td>{e('; '.join(r['flags']) or 'none')}</td><td><span class='b {r['status']}'>{e(r['status'])}</span></td></tr>"
        for r in rows)
    cards = "".join(
        f"<section class=card><h3>{e(r['id'])} {e(r['gene'])} <span class=m>{e(r['label'])}</span> "
        f"<span class='b {r['status']}'>{e(r['status'])}</span></h3>"
        f"<p>{e(meta['tumor_name'])}: {e(_support_text(r['tumor'], r['sv']) + _strand_text(r['tumor']))}<br>"
        + (f"{e(meta['normal_name'])}: {e(_support_text(r['normal'], r['sv']) + _strand_text(r['normal']))}<br>"
           if paired else "")
        + f"Caller reported (VCF): {e(_caller_text(r['caller']))}<br>"
        f"Flags: {e(', '.join(FLAG_TEXT[f] for f in r['flags']) or 'none')}</p>"
        + "".join(f"<details><summary>IGV screenshot</summary><img src='{img(p)}' alt='{e(r['id'])}'></details>"
                  for p in r["figures"]) + "</section>" for r in rows)
    page = f"""<!doctype html><html lang=en><head><meta charset=utf-8><meta name=viewport content="width=device-width,initial-scale=1">
<title>IGV Validator Report</title><style>
:root{{--bg:#fbfaf8;--fg:#1d1d1b;--m:#6b6a66;--card:#fff;--line:#e4e1db;--g:#1f7a4d;--gb:#e3f3ea;--y:#8a6200;--yb:#f7edd2;--r:#a23b2a;--rb:#f8e4df}}
@media (prefers-color-scheme:dark){{:root{{--bg:#181817;--fg:#ecebe7;--m:#a09e98;--card:#222220;--line:#3a3935;--g:#7fd3a5;--gb:#1e3a2b;--y:#e8c46c;--yb:#3d3219;--r:#f0a292;--rb:#45251f}}}}
body{{margin:0;background:var(--bg);color:var(--fg);font:15px/1.5 -apple-system,system-ui,sans-serif}}main{{max-width:1100px;margin:0 auto;padding:28px 16px}}
.m{{color:var(--m);font-weight:400;font-size:14px}}.w{{overflow-x:auto;border:1px solid var(--line);border-radius:10px;background:var(--card)}}
table{{border-collapse:collapse;width:100%;font-size:13.5px}}th,td{{padding:7px 10px;border-bottom:1px solid var(--line);text-align:left}}td.n{{white-space:nowrap}}
.b{{padding:2px 8px;border-radius:999px;font-size:12px;font-weight:600}}.supported{{color:var(--g);background:var(--gb)}}.flagged{{color:var(--y);background:var(--yb)}}.insufficient,.no_coverage{{color:var(--r);background:var(--rb)}}
.card{{background:var(--card);border:1px solid var(--line);border-radius:12px;padding:14px 16px;margin:12px 0}}.card h3{{margin:0 0 6px;font-size:16px}}
img{{max-width:100%;margin-top:8px;border:1px solid var(--line)}}.disc{{margin-top:32px;padding:12px;border:1px solid var(--line);border-radius:8px;color:var(--m)}}
</style></head><body><main><h1>IGV Validator Report</h1>
<p class=m>{e(meta['vcf'])} · tumor {e(meta['tumor_name'])} · {e('normal ' + meta['normal_name'] if paired else 'tumor-only (no normal)')} · {e(meta['date'])}</p>
<p><b>{counts['supported']} supported, {counts['flagged']} flagged, {counts['insufficient']} insufficient</b>. {e(shot_line.replace('`', ''))}</p>
<p class=m>Counts come from the BAMs, never from the screenshots. Status is a rule-based summary of the flags, not a verdict.</p>
<div class=w><table><tr><th>ID</th><th>Gene</th><th>Variant</th><th>Tumor</th>{'<th>Normal</th>' if paired else ''}<th>Caller (VCF)</th><th>Flags</th><th>Status</th></tr>{body}</table></div>
{cards}<p class=disc>{e(DISCLAIMER)}</p></main></body></html>"""
    (out / "report.html").write_text(page)
    return [out / "report.md", out / "report.html", tsv]


# ── Copy-number mode ───────────────────────────────────────────────────────────
# A CNV is checked by depth, not by reads carrying a change: for each gene, read depth inside the gene
# is compared with the flanking regions and with the GATK segment call covering it.

CNV_FLANK_MIN, CNV_FLANK_MAX = 10_000, 1_000_000   # flank on each side: 3x the gene, within these bounds
CNV_BIN_MIN = 1000                                  # bp; genes are split into about 50 bins
CNV_DEEP_LOG2 = -2.0                                # a segment below this is a deep (homozygous) deletion
CNV_IGV_MAX = 100_000                               # IGV screenshot only for windows up to this size
CNV_LOG2_TOL = 0.5                                  # depth log2 within this of the segment log2 = agreement
CNV_FOCAL_MAX = 3_000_000                           # gain/loss segments up to this size move the local flanks
CNV_BASELINE_WINDOWS, CNV_BASELINE_BIN = 300, 10_000  # sample-wide depth: median of random 10 kb windows
CNV_FLAG_TEXT = {
    "cnv_disagrees": f"the depth log2 (gene vs the sample-wide depth) does not reproduce the segment log2 "
                     f"(off by more than {CNV_LOG2_TOL:g}, or deep on one side only)",
    "no_segment": "no segment covers the gene (it falls between segments); depth is still measured",
    "multiple_segments": "more than one segment overlaps the gene (a copy-number breakpoint inside it)",
    "low_mapq_remaining": "the few reads left in a deletion are mostly MAPQ < 20 without alternative hits "
                          "(in PDX tumors often mouse DNA)",
    "ambiguous_mapping": "most reads here have MAPQ < 20 with alternative hits (XA), e.g. on _alt haplotype contigs "
                         "(regions with alternate haplotype copies): the reads are present (IGV shows them), so a MAPQ-filtering caller can "
                         "report a false loss",
    "low_depth": "flank depth under 5x; the ratio is unreliable",
    "cn_step": "the depth changes sharply inside the gene (a copy-number breakpoint, often at an SV); the gene "
               "average hides it",
    "cnv_not_visible": "a gain/loss was called, but the read depth in the gene is within the normal range",
}


def read_cnv_segments(path: Path, sample: str | None = None) -> tuple[list[dict], str]:
    """Segments from an igv-validator/cnv_extracted TSV (sample contig start end log2 cn_call loh ...) or a
    GATK called.seg (CONTIG START END ... MEAN_LOG2_COPY_RATIO CALL). Returns (segments, sample used)."""
    lines = [l for l in Path(path).read_text().splitlines() if l.strip() and not l.startswith("@")]
    if not lines:
        raise InputError(f"--cnv file {path} is empty")
    head = [h.strip().lower() for h in lines[0].split("\t")]
    col = {h: i for i, h in enumerate(head)}
    chrom_c = next((col[c] for c in ("contig", "chrom", "chr", "chromosome") if c in col), None)
    log2_c = next((col[c] for c in ("log2", "mean_log2_copy_ratio", "log2_copy_ratio") if c in col), None)
    if chrom_c is None or "start" not in col or "end" not in col or log2_c is None:
        raise InputError(f"--cnv file {path} needs contig/start/end/log2 columns (got: {', '.join(head)})")
    call_c = col.get("cn_call", col.get("call"))
    rows = [l.split("\t") for l in lines[1:]]
    if "sample" in col:
        names = sorted({r[col["sample"]] for r in rows})
        if sample is None:
            if len(names) != 1:
                raise InputError(f"--cnv file has several samples ({', '.join(names)}); pass --cnv-sample")
            sample = names[0]
        elif sample not in names:
            raise InputError(f"--cnv-sample {sample} is not in {path} (samples: {', '.join(names)})")
        rows = [r for r in rows if r[col["sample"]] == sample]
    calls = {"+": "amp", "-": "del", "0": "neutral"}
    segs = []
    for r in rows:
        try:
            log2 = float(r[log2_c])
        except ValueError:
            log2 = float("nan")
        call = r[call_c].strip() if call_c is not None and call_c < len(r) else ""
        loh = r[col["loh"]].strip() if "loh" in col and col["loh"] < len(r) else ""
        segs.append({"contig": r[chrom_c], "start": int(r[col["start"]]), "end": int(r[col["end"]]),
                     "log2": log2, "cn_call": calls.get(call, call or "na"), "loh": loh in ("1", "true", "True")})
    return segs, sample or ""


def _depth_bins(bam, chrom, start, end, binsize):
    """Per-bin read starts (MAPQ >= 20 and all mapping qualities), duplicates/secondary/supplementary excluded."""
    n = max(1, (end - start + 1) // binsize)  # complete bins only: a partial last bin would look like a dip
    good, allq, lens = [0] * n, [0] * n, []
    for r in bam.fetch(chrom, max(0, start - 1), end):
        if r.is_unmapped or r.is_duplicate or r.is_secondary or r.is_supplementary or r.is_qcfail:
            continue
        i = (r.reference_start + 1 - start) // binsize
        if not 0 <= i < n:
            continue
        allq[i] += 1
        if r.mapping_quality >= MIN_MAPQ:
            good[i] += 1
            if len(lens) < 2000:
                lens.append(r.query_length or 100)
    return good, allq, (sum(lens) / len(lens)) if lens else 100.0


def cnv_log2_agrees(seg_log2: float, obs_log2: float) -> bool:
    """Does the log2 measured from depth reproduce the segment's log2?"""
    if seg_log2 <= CNV_DEEP_LOG2:            # a deep deletion was called: the reads must be nearly gone too
        return obs_log2 <= -1.5
    if obs_log2 <= CNV_DEEP_LOG2 and seg_log2 > -1:   # the reads are nearly gone but no loss was called
        return False
    # both a clear loss or both a clear gain (depth_state's words): the same change, even if the sizes differ
    if seg_log2 <= -0.5 and obs_log2 <= -0.5 or seg_log2 >= 0.4 and obs_log2 >= 0.4:
        return True
    return abs(obs_log2 - seg_log2) <= CNV_LOG2_TOL


STEP_MIN_LOG2 = 0.58     # one copy on two (2 -> 3); a smaller jump inside a gene is not reported as a step
STEP_MIN_SPREAD = 2.0   # a step must exceed twice the two sides' median absolute deviations
STEP_MIN_FRAC = 0.2      # each side of a step covers at least this share of the gene


def depth_step(vals: list[float], start: int, binsize: int) -> dict | None:
    """The strongest depth step inside a gene: {pos, left, right, log2_change}, or None.

    Medians on each side, and the two sides' middle halves (25th-75th percentile) must not overlap, so noise
    and a one-bin mappability dip are not steps."""
    import math
    import statistics as st
    n = len(vals)
    kmin = max(2, math.ceil(n * STEP_MIN_FRAC))
    best = None
    for k in range(kmin, n - kmin + 1):
        left, right = vals[:k], vals[k:]
        lm, rm = st.median(left), st.median(right)
        change = math.log2((rm + 0.5) / (lm + 0.5))
        q = lambda xs: (sorted(xs)[len(xs) // 4], sorted(xs)[(3 * len(xs)) // 4])
        (l25, l75), (r25, r75) = q(left), q(right)
        if not (r25 > l75 or l25 > r75):
            continue
        sse = sum((x - lm) ** 2 for x in left) + sum((x - rm) ** 2 for x in right)
        if best is None or sse < best[4]:   # the split that best separates two flat levels
            best = (change, k, lm, rm, sse)
    if best is None or abs(best[0]) < STEP_MIN_LOG2:
        return None
    change, k, lm, rm, _ = best
    mad = lambda xs, m: st.median(abs(x - m) for x in xs)
    if abs(rm - lm) < STEP_MIN_SPREAD * (mad(vals[:k], lm) + mad(vals[k:], rm) + 1):
        return None   # the jump is not large against the scatter on each side (noisy samples, spikes)
    return {"pos": start + k * binsize, "left": round(lm, 1), "right": round(rm, 1), "log2_change": round(change, 2)}


def segment_parts(vals: list[float], gstart: int, binsize: int, segs: list[dict], baseline: float) -> list[dict]:
    """For a gene covered by several segments: the depth log2 of the part of the gene inside each one."""
    import math
    parts = []
    for sg in sorted(segs, key=lambda x: x["start"]):
        idx = [i for i in range(len(vals)) if sg["start"] <= gstart + i * binsize + binsize // 2 <= sg["end"]]
        if not idx or not baseline:
            continue
        d = sum(vals[i] for i in idx) / len(idx)
        dl = round(max(math.log2(d / baseline), -10.0), 2) if d > 0 else -10.0
        parts.append({"start": max(sg["start"], gstart), "end": min(sg["end"], gstart + len(vals) * binsize - 1),
                      "cn_call": sg["cn_call"], "seg_log2": sg["log2"], "depth": round(d, 1), "depth_log2": dl,
                      "agrees": math.isnan(sg["log2"]) or cnv_log2_agrees(sg["log2"], dl)})
    return parts


def cnv_visible(cn_call: str, depth_log2: float | None) -> bool:
    """Does the read depth show the called direction (loss <= -0.5, gain >= 0.4, as depth_state words it)?"""
    if depth_log2 is None or cn_call not in ("del", "amp"):
        return True
    state = depth_state(depth_log2)
    return state in ("loss", "deep loss") if cn_call == "del" else state == "gain"


def cnv_flank_windows(gstart: int, gend: int, segs: list[dict], contig_len: int,
                      all_segments: list[dict] | None = None) -> list[tuple[int, int]]:
    """Local flank windows, placed outside a focal gain/loss covering the gene.

    A small gene inside a focal deletion (e.g. a 4.5 kb gene in a ~212 kb homozygous deletion split over
    adjacent segments) would otherwise be compared with flanks that are deleted too, so the flanks skip the
    whole run of adjacent segments with the same call. Broad (> 3 Mb) and neutral segments leave the flanks
    next to the gene; the sample-wide depth is the reference for those.
    """
    glen = gend - gstart + 1
    size = min(max(3 * glen, CNV_FLANK_MIN), CNV_FLANK_MAX)
    focal = lambda sg: sg["cn_call"] in ("del", "amp") and sg["end"] - sg["start"] + 1 <= CNV_FOCAL_MAX
    altered = [sg for sg in segs if focal(sg)]
    lo = min([gstart] + [sg["start"] for sg in altered])
    hi = max([gend] + [sg["end"] for sg in altered])
    if altered and all_segments:
        calls = {sg["cn_call"] for sg in altered}
        same = sorted((sg for sg in all_segments if focal(sg) and sg["cn_call"] in calls and
                       sg["contig"] == altered[0]["contig"]), key=lambda sg: sg["start"])
        grown = True
        while grown:  # absorb neighbouring segments with the same call, within 20 kb
            grown = False
            for sg in same:
                if sg["end"] < lo and lo - sg["end"] <= 20_000:
                    lo, grown = sg["start"], True
                if sg["start"] > hi and sg["start"] - hi <= 20_000:
                    hi, grown = sg["end"], True
    wins = []
    if lo > 1:
        wins.append((max(1, lo - size), lo - 1))
    if hi < contig_len:
        wins.append((hi + 1, min(contig_len, hi + size)))
    return wins


def baseline_contigs(references: list[str], gene_contigs: list[str]) -> list[str]:
    """Contigs for the sample-wide depth: every autosome (chr1-22 or 1-22), so a gain or loss of the genes' own
    chromosome cannot shift the scale; the genes' contigs only when the reference has no numbered autosomes."""
    auto = [c for c in references if re.fullmatch(r"(chr)?\d+", c)]
    return auto or sorted(set(gene_contigs))


def baseline_depth(bam_path, fasta, contigs: list[str]) -> float:
    """Sample-wide typical depth (MAPQ >= 20): median over random 10 kb windows on these contigs, skipping empty
    windows (gaps, centromeres). GATK's log2 is relative to the sample's overall depth, so this is its scale."""
    import random
    rng = random.Random(0)
    with _open_bam(bam_path, fasta) as bam:
        lens = {c: bam.get_reference_length(c) for c in contigs if bam.get_reference_length(c) > CNV_BASELINE_BIN}
        if not lens:
            return 0.0
        names, weights = list(lens), [lens[c] for c in lens]
        vals = []
        for _ in range(CNV_BASELINE_WINDOWS):
            c = rng.choices(names, weights)[0]
            start = rng.randint(1, lens[c] - CNV_BASELINE_BIN)
            good, _, rlen = _depth_bins(bam, c, start, start + CNV_BASELINE_BIN - 1, CNV_BASELINE_BIN)
            if good[0]:
                vals.append(good[0] * rlen / CNV_BASELINE_BIN)
    vals.sort()
    return vals[len(vals) // 2] if vals else 0.0


AMBIGUOUS_MIN_FRAC = 0.5    # share of reads that are MAPQ < 20 with other hits
AMBIGUOUS_MIN_DEPTH = 0.2   # ... and all-reads depth at least this share of the sample, unless the hits are _alt


def is_ambiguous(amb_frac, alt_frac, total: int, all_depth: float, baseline: float) -> bool:
    """Reads present but unplaceable: most are MAPQ < 20 with other hits, and either those hits are on alternate
    haplotype contigs (_alt), or there are enough of them to be real coverage. A few junk reads left inside a true
    homozygous deletion (mouse or mismapped DNA) are neither."""
    if amb_frac is None or amb_frac < AMBIGUOUS_MIN_FRAC or total < 5:
        return False
    return (alt_frac or 0) >= 0.5 or (bool(baseline) and all_depth >= AMBIGUOUS_MIN_DEPTH * baseline)


def cnv_gene(bam_path, fasta, gene: str, chrom: str, gstart: int, gend: int, segments: list[dict],
             baseline: float | None = None) -> dict:
    """Depth in the gene vs flanks outside its gain/loss segment, the overlapping segment(s), flags, status."""
    glen = gend - gstart + 1
    binsize = max(CNV_BIN_MIN, glen // 50)
    segs = sorted((sg for sg in segments if sg["contig"] == chrom and sg["end"] >= gstart and sg["start"] <= gend),
                  key=lambda sg: -(min(sg["end"], gend) - max(sg["start"], gstart)))
    with _open_bam(bam_path, fasta) as bam:
        clen = bam.get_reference_length(chrom)
        wins = cnv_flank_windows(gstart, gend, segs, clen, segments)
        g_good, g_all, rlen = _depth_bins(bam, chrom, gstart, gend, binsize)
        flank_vals = []
        for a, b in wins:
            fg, _, _ = _depth_bins(bam, chrom, a, b, binsize)
            flank_vals += fg
        # the plotted window: gene, its segment and the flanks, unless that spans more than 3 Mb
        pa, pb = (min([gstart] + [a for a, _ in wins]), max([gend] + [b for _, b in wins]))
        if pb - pa > 3_000_000:
            pad = min(max(3 * glen, CNV_FLANK_MIN), CNV_FLANK_MAX)
            pa, pb = max(1, gstart - pad), min(clen, gend + pad)
        pbin = max(CNV_BIN_MIN, (pb - pa + 1) // 200)
        p_good, p_all, _ = _depth_bins(bam, chrom, pa, pb, pbin)
    to_depth = rlen / binsize
    gene_depth = sum(g_good) / len(g_good) * to_depth if g_good else 0.0
    flank_vals = sorted(v * to_depth for v in flank_vals)
    flank_depth = flank_vals[len(flank_vals) // 2] if flank_vals else 0.0
    ratio = round(gene_depth / flank_depth, 3) if flank_depth else None
    gene_all, gene_good = sum(g_all), sum(g_good)
    frac_low = round(1 - gene_good / gene_all, 3) if gene_all else None
    all_depth = sum(g_all) / len(g_all) * to_depth if g_all else 0.0
    # reads that are present but cannot be placed uniquely: MAPQ < 20 with alternative hits (or MAPQ 0)
    ambiguous = alt_hits = total = 0
    with _open_bam(bam_path, fasta) as bam:
        for r in bam.fetch(chrom, max(0, gstart - 1), gend):
            if r.is_unmapped or r.is_duplicate or r.is_secondary or r.is_supplementary or r.is_qcfail:
                continue
            if not gstart <= r.reference_start + 1 <= gend:
                continue
            total += 1
            if r.mapping_quality < MIN_MAPQ and (r.has_tag("XA") or r.mapping_quality == 0):
                ambiguous += 1
                if r.has_tag("XA") and "_alt," in r.get_tag("XA"):
                    alt_hits += 1
    amb_frac = round(ambiguous / total, 3) if total else None
    starts = [pa + i * pbin for i in range(len(p_good))]
    import math
    if baseline is None:
        with _open_bam(bam_path, fasta) as bam:
            refs = list(bam.references)
        baseline = baseline_depth(bam_path, fasta, baseline_contigs(refs, [chrom]))
    depth_log2 = (round(max(math.log2(gene_depth / baseline), -10.0), 2) if gene_depth > 0 else -10.0) \
        if baseline else None
    depth_log2_all = (round(max(math.log2(all_depth / baseline), -10.0), 2) if all_depth > 0 else -10.0) \
        if baseline else None
    flags = []
    if baseline < 5:
        flags.append("low_depth")
    if len(segs) > 1:
        flags.append("multiple_segments")
    alt_frac = round(alt_hits / total, 3) if total else None
    ambiguous_mapping = is_ambiguous(amb_frac, alt_frac, total, all_depth, baseline)
    if ambiguous_mapping:
        flags.append("ambiguous_mapping")
    elif depth_log2 is not None and depth_log2 < CNV_DEEP_LOG2 and frac_low is not None and frac_low >= 0.5 \
            and gene_all >= 5:
        flags.append("low_mapq_remaining")
    gene_vals = [g * to_depth for g in g_good]
    step = depth_step(gene_vals, gstart, binsize)
    if step and baseline:
        for side in ("left", "right"):
            step[f"{side}_log2"] = round(max(math.log2(step[side] / baseline), -10.0), 2) if step[side] > 0 else -10.0
        flags.append("cn_step")
    parts = segment_parts(gene_vals, gstart, binsize, segs, baseline) if len(segs) > 1 else []
    if parts:     # several segments: each part of the gene is checked against its own segment
        if not all(p["agrees"] for p in parts):
            flags.append("cnv_disagrees")
    elif segs and depth_log2 is not None and not math.isnan(segs[0]["log2"]) \
            and not cnv_log2_agrees(segs[0]["log2"], depth_log2):
        flags.append("cnv_disagrees")
    elif segs and not cnv_visible(segs[0]["cn_call"], depth_log2) and "ambiguous_mapping" not in flags:
        flags.append("cnv_not_visible")
    if not baseline:
        status = "no_coverage"
    elif not segs:
        status = "no_segment"
    else:
        called = segs[0]["cn_call"] in ("del", "amp")
        status = "flagged" if ("cnv_disagrees" in flags or (ambiguous_mapping and called)) else \
            "insufficient" if "cnv_not_visible" in flags else "supported"
    return {"gene": gene, "chrom": chrom, "start": gstart, "end": gend, "segments": segs,
            "gene_depth": round(gene_depth, 1), "baseline_depth": round(baseline, 1), "depth_log2": depth_log2,
            "depth_log2_all": depth_log2_all, "ambiguous_fraction": amb_frac, "alt_fraction": alt_frac, "all_depth": round(all_depth, 1),
            "flank_depth": round(flank_depth, 1), "depth_ratio": ratio,
            "low_mapq_fraction": frac_low, "reads_in_gene": gene_all, "flags": flags, "status": status,
            "flank_windows": wins, "step": step, "parts": parts,
            "plot": {"starts": starts, "good": [g * rlen / pbin for g in p_good], "all": [x * rlen / pbin for x in p_all],
                     "binsize": pbin}, "figures": []}


def cnv_plot(r: dict, out_png: Path, sample: str):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    pl = r["plot"]
    x = [(s + pl["binsize"] / 2) / 1000 for s in pl["starts"]]
    fig, ax = plt.subplots(figsize=(10, 3.6))
    ax.axvspan(r["start"] / 1000, r["end"] / 1000, color="#f2c6b8", alpha=0.5, label=r["gene"])
    ax.plot(x, pl["all"], color="#bbbbbb", lw=1, label="all reads")
    ax.plot(x, pl["good"], color="#1f5f99", lw=1.6, label="MAPQ >= 20")
    for sg in r["segments"]:
        y = r["baseline_depth"] * (2 ** max(min(sg["log2"], 2), -4)) if r["baseline_depth"] else 0
        ax.hlines(y, max(sg["start"], pl["starts"][0]) / 1000, min(sg["end"], pl["starts"][-1] + pl["binsize"]) / 1000,
                  colors="#d85a30", lw=2.2, label=f"segment {sg['cn_call']} (log2 {sg['log2']:.2f})")
    ax.set_xlabel(f"{r['chrom']} position (kb)")
    ax.set_ylabel("read depth (x)")
    ratio = "NA" if r["depth_ratio"] is None else f"{r['depth_ratio']:.2f}"
    ax.axhline(r["baseline_depth"], color="#888780", lw=0.8, ls="--", label="sample-wide depth")
    ax.set_title(f"{sample} · {r['gene']} {r['chrom']}:{r['start']:,}-{r['end']:,} · depth log2 {r['depth_log2']} "
                 f"· local ratio {ratio} · {r['status']}", fontsize=10)
    ax.set_ylim(bottom=0)
    ax.legend(fontsize=7, loc="upper right", ncol=2)
    fig.tight_layout()
    fig.savefig(out_png, dpi=130)
    plt.close(fig)


def take_cnv_screenshots(rows, tumor, genome, snapdir, igv_path, name, log_copy=None,
                         timeout: int = 300, annotation=None) -> tuple[int, str]:
    """IGV coverage/reads view for genes whose window fits (CNV_IGV_MAX); larger genes rely on the depth plot."""
    small = [r for r in rows if (r["end"] - r["start"] + 10_000) <= CNV_IGV_MAX]
    if not small:
        return 0, "all genes are larger than an IGV window; see the depth plots"
    argv, cwd, why = find_igv(igv_path)
    if argv is None:
        return 0, why
    if sys.platform.startswith("linux") and not __import__("os").environ.get("DISPLAY"):
        return 0, "no display available (run with xvfb-run, or use an HPC desktop session)"
    snapdir.mkdir(parents=True, exist_ok=True)
    try:
        igv = IGVSession(argv, cwd, genome, timeout, log_copy)
    except RuntimeError as e:
        return 0, str(e)
    taken, note = 0, ""
    sub = annotation_subset(annotation, [(r["chrom"], r["start"], r["end"]) for r in small],
                            snapdir / "_genes") if annotation else None
    try:
        for r in small:
            a, b = max(1, r["start"] - 5000), r["end"] + 5000
            png = snapdir / f"{re.sub(r'[^A-Za-z0-9._-]+', '_', r['gene'])}_igv.png"
            for c in ["new", f"snapshotDirectory {snapdir}", "preference SAM.DOWNSAMPLE_READS true",
                      f"preference SAM.MAX_VISIBLE_RANGE {CNV_IGV_MAX // 1000 + 20}", "maxPanelHeight 500",
                      f"load {Path(tumor).resolve()} name={name}", *annotation_cmds(sub),
                      f"goto {r['chrom']}:{a}-{b}"]:
                igv.send(c)
            for attempt in range(3):
                time.sleep(1.0 * (attempt + 1))
                igv.send(f"goto {r['chrom']}:{a}-{b}")
                igv.send(f"collapse {name}")   # reads only; the gene track stays expanded
                if sub:
                    igv.send("expand genes")
                igv.send(f"snapshot {png.name}")
                if png.exists() and image_rendered(png):
                    break
            r["_igv_png"] = png
            taken += 1
    except RuntimeError as e:
        note = f"stopped early: {e}"
    finally:
        if note and log_copy:
            log_copy.parent.mkdir(parents=True, exist_ok=True)
            igv.log.flush()
            shutil.copy(igv.home / "igv.log", log_copy)
        igv.close()
    return taken, note


def write_cnv_reports(out: Path, rows: list[dict], meta: dict) -> list[Path]:
    tables = out / "tables"; tables.mkdir(parents=True, exist_ok=True)
    tsv = tables / "cnv_depth.tsv"
    cols = ["gene", "chrom", "start", "end", "segment", "log2", "cn_call", "loh", "gene_depth", "baseline_depth",
            "depth_log2", "flank_depth", "depth_ratio", "low_mapq_fraction", "reads_in_gene", "flags", "status"]
    with open(tsv, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t"); w.writerow(cols)
        for r in rows:
            sg = r["segments"][0] if r["segments"] else None
            w.writerow([r["gene"], r["chrom"], r["start"], r["end"],
                        f"{sg['contig']}:{sg['start']}-{sg['end']}" if sg else "", sg["log2"] if sg else "",
                        sg["cn_call"] if sg else "", int(sg["loh"]) if sg else "", r["gene_depth"], r["baseline_depth"],
                        "" if r["depth_log2"] is None else r["depth_log2"], r["flank_depth"],
                        "" if r["depth_ratio"] is None else r["depth_ratio"], "" if r["low_mapq_fraction"] is None
                        else r["low_mapq_fraction"], r["reads_in_gene"], ";".join(r["flags"]), r["status"]])

    def segtxt(r):
        return "; ".join(f"{sg['cn_call']} log2 {sg['log2']:.2f}{' LOH' if sg['loh'] else ''} "
                         f"({sg['start']:,}-{sg['end']:,})" for sg in r["segments"]) or "none (gap)"

    counts = {k: sum(r["status"] == k for r in rows) for k in ("supported", "flagged", "no_segment", "no_coverage")}
    md = ["# IGV Validator Report: Copy number", "",
          f"**Input**: `{meta['cnv']}` (sample {meta['cnv_sample'] or '-'}) · BAM: {meta['tumor_name']} · "
          f"{len(rows)} gene(s)", f"**Date**: {meta['date']}  ", f"**Skill**: igv-validator {VERSION}", "",
          f"**{counts['supported']} supported, {counts['flagged']} flagged, {counts['no_segment']} without a segment"
          + (f", {counts['no_coverage']} no coverage" if counts["no_coverage"] else "") + "**.", "",
          f"IGV screenshots: {meta['shots']}.", "",
          "A copy-number call is checked by **read depth** (MAPQ >= 20, duplicates removed). **Depth log2** is the "
          "gene's mean depth over the sample-wide typical depth, on the same scale as GATK's log2, so it should "
          "reproduce the segment's log2 (about -1 for a one-copy loss in a diploid genome, very negative for a "
          "homozygous deletion, 0 for neutral). The **local ratio** compares the gene with the regions beside it "
          "(outside a focal gain/loss) and picks up focal changes inside a larger segment. Purity and ploidy shift "
          "these values; **status** summarises agreement with the segment call, it is not a verdict.", "",
          "| Gene | Region | Segment call(s) | Gene depth | Sample depth | Depth log2 | Depth log2, all reads (as IGV "
          "shows) | Local ratio | Ambiguous reads | Flags | Status |", "|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        md.append(f"| {r['gene']} | {r['chrom']}:{r['start']:,}-{r['end']:,} | {segtxt(r)} | {r['gene_depth']}x | "
                  f"{r['baseline_depth']}x | {'NA' if r['depth_log2'] is None else r['depth_log2']} | "
                  f"{'NA' if r['depth_log2_all'] is None else r['depth_log2_all']} | "
                  f"{'NA' if r['depth_ratio'] is None else r['depth_ratio']} | "
                  f"{'NA' if r['ambiguous_fraction'] is None else r['ambiguous_fraction']} | "
                  f"{'; '.join(r['flags']) or 'none'} | **{r['status']}** |")
    md += ["", "## Flags", ""] + [f"- `{k}`: {t}" for k, t in CNV_FLAG_TEXT.items()]
    md += ["", "## Genes", ""]
    for r in rows:
        md += [f"### {r['gene']}: {r['status']}", ""] + [f"![{r['gene']}]({p})" for p in r["figures"]] + [""]
    md += ["---", "", "## Disclaimer", "", f"*{DISCLAIMER}*", ""]
    (out / "report.md").write_text("\n".join(md))
    e = html.escape

    def img(p):
        return "data:image/png;base64," + base64.b64encode((out / p).read_bytes()).decode()

    body = "".join(f"<tr><td><b>{e(r['gene'])}</b></td><td>{e(segtxt(r))}</td><td>{r['gene_depth']}x</td>"
                   f"<td>{r['baseline_depth']}x</td><td>{e(str(r['depth_log2']))}</td><td>{e(str(r['depth_ratio']))}</td>"
                   f"<td>{e('; '.join(r['flags']) or 'none')}</td><td><b>{e(r['status'])}</b></td></tr>" for r in rows)
    figs = "".join(f"<h3>{e(r['gene'])}</h3>" + "".join(f"<img src='{img(p)}' alt='{e(r['gene'])}'>" for p in r["figures"])
                   for r in rows)
    page = (f"<!doctype html><html lang=en><head><meta charset=utf-8><meta name=viewport content='width=device-width,"
            f"initial-scale=1'><title>IGV Validator: Copy Number</title><style>body{{font:15px/1.5 system-ui,sans-serif;"
            f"max-width:1100px;margin:0 auto;padding:24px 16px;background:#fbfaf8;color:#1d1d1b}}table{{border-collapse:"
            f"collapse;width:100%}}td,th{{border-bottom:1px solid #e4e1db;padding:6px 8px;text-align:left}}img{{max-width:"
            f"100%}}@media (prefers-color-scheme:dark){{body{{background:#181817;color:#ecebe7}}}}</style></head><body>"
            f"<h1>IGV Validator: copy number</h1><p>{e(meta['cnv'])} · {e(meta['tumor_name'])} · {e(meta['date'])}</p>"
            f"<table><tr><th>Gene</th><th>Segment call(s)</th><th>Gene depth</th><th>Sample depth</th><th>Depth log2</th>"
            f"<th>Local ratio</th>"
            f"<th>Flags</th><th>Status</th></tr>{body}</table>{figs}<p><i>{e(DISCLAIMER)}</i></p></body></html>")
    (out / "report.html").write_text(page)
    return [out / "report.md", out / "report.html", tsv]


def _inputs(args) -> dict:
    """Absolute input paths, so --summarize --overview can reopen the same BAM and reference."""
    ab = lambda x: str(Path(x).resolve()) if x else None
    return {"tumor": ab(args.tumor), "normal": ab(getattr(args, "normal", None)), "reference": ab(args.reference),
            "igv": igv_location(getattr(args, "igv_path", None))}


def run_cnv(args) -> Path:
    """Copy-number mode: --cnv segments + --regions genes + --tumor BAM."""
    if not args.regions:
        raise InputError("--cnv needs --regions (a BED of the genes to check)")
    if not args.tumor:
        raise InputError("--cnv needs --tumor (the BAM the copy number was called from)")
    for label, p in (("--cnv file", args.cnv), ("--regions file", args.regions)):
        if not Path(p).exists():
            raise InputError(f"{label} not found: {p}")
    _check_inputs(args.cnv, args.tumor, None, args.reference)
    segments, sample = read_cnv_segments(args.cnv, args.cnv_sample)
    genes = [(name, chrom, a, b) for chrom, spans in _read_regions(args.regions).items() for a, b, name in spans]
    with _open_bam(args.tumor, args.reference) as bam:
        missing = sorted({c for _, c, _, _ in genes} - set(bam.references))
    if missing:
        raise InputError(f"contig(s) {', '.join(missing)} are in --regions but not in {args.tumor}")
    out = Path(args.output or f"output/igv_validator_cnv_{datetime.now():%Y%m%d_%H%M%S}")
    if out.exists() and any(out.iterdir()):
        print(f"Warning: {out} is not empty; igv-validator outputs there will be overwritten", file=sys.stderr)
    (out / "figures" / "cnv").mkdir(parents=True, exist_ok=True)
    name = args.tumor_name or Path(args.tumor).stem
    with _open_bam(args.tumor, args.reference) as bam:
        refs = list(bam.references)
    base = baseline_depth(args.tumor, args.reference, baseline_contigs(refs, [c for _, c, _, _ in genes]))
    rows = [cnv_gene(args.tumor, args.reference, g, c, a, b, segments, base) for g, c, a, b in genes]
    for r in rows:
        png = out / "figures" / "cnv" / f"{re.sub(r'[^A-Za-z0-9._-]+', '_', r['gene'])}.png"
        cnv_plot(r, png, sample or name)
        r["figures"].append(str(png.relative_to(out)))
    taken, note = 0, "--no-igv given"
    if not args.no_igv:
        if not args.reference:
            note = "no --reference given (IGV must use the reference the BAM was aligned to)"
        else:
            raw = Path(tempfile.mkdtemp(prefix="igv_cnv_"))
            taken, note = take_cnv_screenshots(rows, args.tumor, Path(args.reference).resolve(), raw, args.igv_path,
                                               name, out / "reproducibility" / "igv.log", args.igv_timeout,
                                               getattr(args, "annotation", None))
            for r in rows:
                if r.pop("_igv_png", None) is not None:
                    src = raw / f"{re.sub(r'[^A-Za-z0-9._-]+', '_', r['gene'])}_igv.png"
                    if src.exists():
                        dst = out / "figures" / "cnv" / src.name
                        sg = r["segments"][0] if r["segments"] else None
                        ratio = "NA" if r["depth_ratio"] is None else f"{r['depth_ratio']:.2f}"
                        caption(src, dst, [f"{r['gene']}  {r['chrom']}:{r['start']:,}-{r['end']:,}  ({name})",
                                           "Segment: " + (f"{sg['cn_call']} log2 {sg['log2']:.2f}" if sg else
                                                          "none covers this gene"),
                                           f"Depth: gene {r['gene_depth']}x vs sample {r['baseline_depth']}x "
                                           f"(log2 {r['depth_log2']}); local ratio {ratio}; status {r['status']}",
                                           "Depth comes from the BAM (MAPQ >= 20), not from the image."])
                        r["figures"].append(str(dst.relative_to(out)))
            shutil.rmtree(raw, ignore_errors=True)
    shots = f"{taken} image(s)" + (f"; {note}" if note else "") if taken else f"skipped: {note}"
    meta = {"cnv": Path(args.cnv).name, "cnv_sample": sample, "tumor_name": name, "shots": shots,
            "date": datetime.now(timezone.utc).strftime("%Y-%m-%d")}
    for r in rows:
        r.pop("plot", None)
    written = write_cnv_reports(out, rows, meta)
    summary = {"mode": "cnv", "genes_checked": len(rows), **{k: sum(r["status"] == k for r in rows)
                                                            for k in ("supported", "flagged", "no_segment",
                                                                      "no_coverage")},
               "screenshots_taken": taken}
    thresholds = {k: globals()[k] for k in ("MIN_MAPQ", "CNV_FLANK_MIN", "CNV_FLANK_MAX", "CNV_BIN_MIN",
                                            "CNV_DEEP_LOG2", "CNV_IGV_MAX", "CNV_LOG2_TOL",
                                            "CNV_FOCAL_MAX", "CNV_BASELINE_WINDOWS", "CNV_BASELINE_BIN")}
    written.append(write_result_json(out, "igv-validator", VERSION, summary,
                                     {"cnv": rows, "thresholds": thresholds, "sample": sample,
                                      "tool": f"copy-number segments ({Path(args.cnv).name})",
                                      "inputs": _inputs(args)},
                                     input_checksum=hashlib.sha256(Path(args.cnv).read_bytes()).hexdigest()))
    write_commands_sh(out, "python skills/igv-validator/igv_validator.py " + shlex.join(sys.argv[1:]))
    write_environment_yml(out, "clawbio-igv-validator", ["pysam>=0.22", "Pillow>=10.1", "matplotlib>=3.8"])
    write_checksums(written + [out / f for r in rows for f in r["figures"]], out, anchor=out)
    return out


# ── Summary across runs ────────────────────────────────────────────────────────

CONCORDANCE = {"supported": "confirmed", "flagged": "questioned", "insufficient": "weak",
               "no_segment": "no call", "no_coverage": "no reads"}
CONCORDANCE_ORDER = ["questioned", "recurrent", "weak", "no reads", "confirmed", "no call"]  # cell colour: first that applies
CONCORDANCE_COLOR = {"questioned": "#f3c7a8", "weak": "#f6e3a1", "no reads": "#e2d9d9", "confirmed": "#c6e8d2", "recurrent": "#dcd6ec",
                     "no call": "#ecebe7"}


def depth_state(log2: float | None) -> str:
    """Words for a depth log2 (gene vs sample-wide depth)."""
    if log2 is None:
        return "unknown"
    if log2 <= CNV_DEEP_LOG2:
        return "deep loss"
    if log2 <= -0.5:
        return "loss"
    if log2 >= 0.4:
        return "gain"
    return "neutral"


def _summary_rows(root: Path, empty_runs: list | None = None) -> list[dict]:
    rows = []
    empty_runs = [] if empty_runs is None else empty_runs
    for rj in sorted(root.rglob("result.json")):
        try:
            d = json.loads(rj.read_text())
        except (OSError, ValueError):
            continue
        if d.get("skill") != "igv-validator":
            continue
        run = rj.parent
        report = str((run / "report.html").relative_to(root))
        data = d.get("data", {})
        tool = data.get("tool", "")
        inputs = data.get("inputs") or {}
        if "cnv" in data:
            sample = data.get("sample") or run.parent.name
            for c in data["cnv"]:
                sg = c["segments"][0] if c["segments"] else None
                ambiguous = "ambiguous_mapping" in c["flags"]
                state = (sg["cn_call"] if sg else "no segment; reads ambiguous" if ambiguous
                         else f"no segment; depth {depth_state(c['depth_log2'])}")
                shows = (f"{c['gene']} covered at {c.get('all_depth', 'NA')}x counting all reads (what IGV "
                         f"displays), {c['gene_depth']}x at MAPQ >= 20; sample average {c['baseline_depth']}x")
                if (c.get("ambiguous_fraction") or 0) >= 0.5:
                    shows += f"; {round(100 * c['ambiguous_fraction'])}% of the reads are ambiguous (MAPQ 0 / alternate hits)"
                rows.append({"sample": sample, "gene": c["gene"], "type": "CNV", "state": state, "tool": tool,
                             "igv_shows": shows, "_bam": inputs.get("tumor"), "_ref": inputs.get("reference"), "_igv": inputs.get("igv"),
                             "_chrom": c["chrom"], "_gstart": c["start"], "_gend": c["end"], "_segs": c["segments"],
                             "_depth": depth_state(c["depth_log2"]), "_depth_log2": c["depth_log2"],
                             "_step": c.get("step"), "_parts": c.get("parts") or [], "_base": c["baseline_depth"],
                             "_ratio": c.get("depth_ratio"),
                             "call": (f"{sg['cn_call']} log2 {sg['log2']:.2f}" + (" LOH" if sg["loh"] else "")) if sg
                             else "no segment",
                             "location": f"{c['chrom']}:{c['start']}-{c['end']}",
                             "reads": f"depth log2 {c['depth_log2']} (gene {c['gene_depth']}x vs sample "
                                      f"{c['baseline_depth']}x); local ratio {c['depth_ratio']}",
                             "caller": "", "flags": ";".join(c["flags"]), "status": c["status"],
                             "concordance": CONCORDANCE.get(c["status"], c["status"]), "report": report,
                             "figures": ";".join(str((run / f).relative_to(root)) for f in c["figures"])})
        if "variants" in data and not data["variants"]:
            empty_runs.append(((data.get("samples") or {}).get("tumor") or run.parent.name, run.name, report))
        for v in data.get("variants", []):
            sample = (data.get("samples") or {}).get("tumor") or run.parent.name
            t = v["tumor"]
            sv = v["kind"] not in ("snv", "insertion", "deletion")
            reads = (f"{t['alt']} reads ({t['split']} split, {t['discordant']} pairs) / {t['depth']}" if sv else
                     f"{t['alt']}/{t['depth']} ({t['vaf_pct']}%)")
            if v.get("normal"):
                n = v["normal"]
                reads += f"; normal {n['alt']}/{n['depth']}"
            c = v.get("caller")
            if sv:
                shows = (f"{t['alt']} reads support the rearrangement ({t['split']} split, {t['discordant']} "
                         f"discordant pairs, counted at both breakends); {t['depth']} reads cover the first breakend" if t["alt"] else
                         f"no read supports the rearrangement ({t['depth']} reads at the breakpoint)")
            elif t["depth"] == 0:
                shows = "no reads at this position"
            else:
                what = v["alt"] if v["kind"] == "snv" else f"the {v['kind']}"
                shows = (f"{t['alt']} of {t['depth']} reads carry {what} ({t['vaf_pct']}%), {t['alt_fwd']} forward / "
                         f"{t['alt_rev']} reverse" if t["alt"] else f"no read of {t['depth']} carries {what}")
            for gene in (v["gene"].split(",") if v["gene"] else ["-"]):
                rows.append({"sample": sample, "gene": gene, "type": "SV" if sv else "SNV/indel", "tool": tool,
                             "igv_shows": shows, "_bam": inputs.get("tumor"), "_ref": inputs.get("reference"), "_igv": inputs.get("igv"),
                             "_chrom": v["chrom"], "_pos": v["pos"], "_chrom2": v.get("chrom2"), "_pos2": v.get("pos2"),
                             "_kind": v["kind"], "_vaf": t.get("vaf_pct"),
                             "state": v["kind"].upper() if sv else "", "call": v["label"],
                             "location": f"{v['chrom']}:{v['pos']}", "reads": reads,
                             "caller": f"{c['alt']}/{c['depth']} ({c['vaf_pct']}%)" if c else "",
                             "flags": ";".join(v["flags"]), "status": v["status"],
                             "concordance": CONCORDANCE.get(v["status"], v["status"]), "report": report,
                             "figures": ";".join(str((run / f).relative_to(root)) for f in v["figures"])})
    return rows


AGREE = {"supported": "yes", "flagged": "no", "insufficient": "unclear", "no_coverage": "unclear"}
WHY = {"recurrent": "same variant in several samples: likely germline or artifact",
       "ambiguous_mapping": "reads are present but map equally well to alternate contigs (MAPQ 0); a MAPQ-filtering "
                            "caller sees a false change",
       "cnv_disagrees": "read depth does not match the called copy number",
       "caller_disagrees": "the caller's read counts differ from the BAM",
       "strand_bias": "all supporting reads on one strand (artifact pattern)",
       "germline_site": "the normal already carries another allele here",
       "normal_support": "supporting reads in the normal",
       "high_depth": "unusually high depth (repeat or mismapped reads)",
       "read_end": "alt bases mostly at read ends (artifact pattern)",
       "low_support": "too few supporting reads", "low_vaf": "under 5% of reads", "low_depth": "low depth"}


SV_SAME_BP = 50   # bp: two SV records whose both ends lie this close are one junction (e.g. Manta + SURVIVOR)


def dedupe_svs(rows: list[dict]) -> list[dict]:
    """Drop SV rows that repeat an earlier row's junction (same sample and gene, both breakends within 50 bp)."""
    kept, seen = [], []
    for r in rows:
        if r["type"] != "SV" or not r.get("_chrom2"):
            kept.append(r)
            continue
        ends = sorted([(r["_chrom"], r["_pos"]), (r["_chrom2"], r["_pos2"])])
        key = (r["sample"], r["gene"].upper())
        dup = any(k == key and all(a[0] == b[0] and abs(a[1] - b[1]) <= SV_SAME_BP for a, b in zip(ends, e))
                  for k, e in seen)
        if not dup:
            seen.append((key, ends))
            kept.append(r)
    return kept


LOW_GERMLINE_VAF = 35.0  # %: an inherited variant sits near 50% (one copy) or 100%, not below this in every sample
RECURRENT_MIN = 3  # samples sharing one SNV/indel before it is treated as germline or artifact


def mark_recurrent(rows: list[dict], min_samples: int = RECURRENT_MIN) -> None:
    """Flag SNVs/indels found at the same position in several samples: independent tumors do not share mutations,
    so the reads carrying them point to a germline variant or an artifact (mapping, PDX mouse reads)."""
    n_samples = len({r["sample"] for r in rows})
    seen, strong = {}, {}
    for r in rows:
        if r["type"] == "SNV/indel":
            seen.setdefault((r["location"], r["call"]), set()).add(r["sample"])
            if r["status"] == "supported":
                strong.setdefault((r["location"], r["call"]), set()).add(r["sample"])
    for r in rows:
        k = (r.get("location"), r["call"])
        # recurrent once enough samples carry it clearly; weak copies of the same variant are then recurrent too
        if r["type"] == "SNV/indel" and len(strong.get(k, ())) >= min_samples and \
                r["status"] in ("supported", "insufficient"):
            r["flags"] = ";".join([f for f in r["flags"].split(";") if f] + ["recurrent"])
            r["concordance"] = "recurrent"
            r["_recurrent"] = f"{len(seen[k])} of {n_samples} samples"
            vafs = [x.get("_vaf") for x in rows if x["type"] == "SNV/indel" and (x.get("location"), x["call"]) == k]
            r["_recurrent_low_vaf"] = bool(vafs) and all(v is not None and v < LOW_GERMLINE_VAF for v in vafs)


REVIEW_TEXT = {
    "ambiguous_mapping": "reads map equally well elsewhere (MAPQ 0)",
    "multiple_segments": "several GATK segments cover the gene",
    "cn_step": "the depth steps inside the gene",
    "cnv_not_visible": "the called change is not visible in the depth",
    "cnv_disagrees": "the caller and the depth disagree",
    "low_depth": "low depth",
    "low_support": "few supporting reads", "low_vaf": "low allele fraction",
    "strand_bias": "all supporting reads on one strand", "read_end": "alt bases mostly at read ends",
    "caller_disagrees": "the caller's counts differ from the BAM", "high_depth": "unusually high depth",
    "normal_support": "supporting reads in the normal",
}
NEAR = [(-0.5, 0.15), (0.4, 0.15), (CNV_DEEP_LOG2, 0.3)]   # depth_state thresholds and how close counts as "near"


def review_reasons(r: dict) -> list[str]:
    """Why a person should look at this row's image before trusting the automatic verdict ([] = no known risk)."""
    flags = [f for f in r["flags"].split(";") if f]
    out = [REVIEW_TEXT[f] for f in flags if f in REVIEW_TEXT]
    if r["type"] == "CNV":
        d = r.get("_depth_log2")
        if d is not None and d > -10 and any(abs(d - t) <= w for t, w in NEAR):
            out.append(f"depth close to a threshold (log2 {d:+.2f})")
        if "_ratio" in r and r["_ratio"] is None and d is not None and r["status"] != "no_coverage":
            out.append("no usable reads around the gene to compare with (duplicated or unmappable flanks)")
    elif r["status"] == "insufficient":
        out.append("weak read support")
    return list(dict.fromkeys(out))


def _check(reasons: list[str], glance: list[str] | None = None) -> str:
    """yes (open the images) / glance (a quick look at the depth plot) / no."""
    if reasons:
        return "yes: " + "; ".join(reasons)
    return "glance: " + "; ".join(glance) if glance else "no"


def _agreement(r: dict) -> dict | None:
    """One line per call the workflow made: does IGV (the reads) agree? Neutral copy number is not a call."""
    if r["type"] == "CNV" and r["state"] not in ("del", "amp") and r["status"] != "flagged":
        return None
    flags = [f for f in r["flags"].split(";") if f]
    agrees = AGREE.get(r["status"], "unclear")
    seg = (r.get("_segs") or [{}])[0]
    small = lambda x: x is not None and abs(x) < CNV_LOG2_TOL
    if r["type"] == "CNV" and agrees == "yes" and small(seg.get("log2")) and small(r.get("_depth_log2", 0)):
        called = f"CNV {r['call']}"
        return {"sample": r["sample"], "gene": r["gene"], "type": r["type"], "tool": r["tool"], "caller_called": called,
                "igv_shows": r["igv_shows"], "igv_agrees": "unclear", "report": r["report"],
                "check_image": _check(review_reasons(r)),
                "why": f"the change is too small to confirm from read depth (a normal gene also reads within "
                       f"{CNV_LOG2_TOL:g} log2 of it)"}
    why = "; ".join(WHY[f] for f in flags if f in WHY) or ("the reads support the call" if agrees == "yes" else "")
    if "recurrent" in flags:
        agrees = "in reads, recurrent"
        why = (f"the reads carry it, but the same variant is in {r['_recurrent']}; independent tumors do not share "
               f"mutations, so it is likely germline or artifact, not a tumor mutation")
        if r.get("_recurrent_low_vaf"):
            why += ("; its allele fractions are too low for an inherited variant (~50% or ~100%), so an artifact "
                    "(e.g. mismapped or mouse reads) is more likely")
    called = f"{r['type']} {r['call']}" if r["type"] != "CNV" else f"CNV {r['call']}"
    return {"sample": r["sample"], "gene": r["gene"], "type": r["type"], "tool": r["tool"], "caller_called": called,
            "igv_shows": r["igv_shows"], "igv_agrees": agrees, "why": why, "report": r["report"],
            "check_image": _check((["IGV does not agree"] if agrees == "no" else []) + review_reasons(r),
                                  ["unclear from the depth; confirm on the depth plot"]
                                  if r["type"] == "CNV" and agrees == "unclear" else None)}


def _no_call_note(r: dict) -> str:
    """Caption for a CNV row with no gain/loss call: say what the reads show instead of assuming agreement."""
    if "ambiguous_mapping" in r["flags"].split(";"):
        return "reads are present but ambiguous (MAPQ 0, alternate contigs); depth cannot judge copy number here"
    if r["state"].startswith("no segment"):
        return "no copy-number call; " + r["state"].split("; ", 1)[1]
    return "no copy-number change; the reads agree"


HEATMAP_FOCAL_BP = 1_000_000   # the heatmap shows a DEL/AMP segment if it is smaller than this ...
HEATMAP_DEEP_LOG2 = 2.0        # ... or deeper than this (log2 < -2 / > 2); broad shallow changes are left out
NAMES = {"DEL": "deletion", "AMP": "gain", "SV": "SV", "SNV": "SNV"}


def _cn_evidence(r: dict, flags: list[str], why: str) -> tuple[str, str, str]:
    """Copy-number part of _evidence (see there), judged from the read depth rather than the caller's number."""
    deep, mild = r["_depth"] == "deep loss", r["_depth"] in ("loss", "gain")
    kind = {"del": "DEL", "amp": "AMP"}.get(r["state"]) or ("AMP" if r["_depth"] == "gain" else "DEL")
    d = r.get("_depth_log2")
    visible = (cnv_visible(r["state"], d) if d is not None else
               r["_depth"] in (("loss", "deep loss") if r["state"] == "del" else ("gain",)))
    if "ambiguous_mapping" in flags:
        return kind, "not supported", f"CNV {r['call']}: {WHY['ambiguous_mapping']}"
    parts = r.get("_parts") or []
    if len(parts) > 1:   # GATK splits the gene: judge each part against its own segment
        words = "; ".join(f"{p['cn_call']} {p['seg_log2']:.2f} from {p['start']:,} (depth {p['depth_log2']:+.2f}"
                          f"{'' if p['agrees'] else ', disagrees'})" for p in parts)
        called = [p for p in parts if p["cn_call"] in ("del", "amp")]
        seen = [p for p in called if p["agrees"] and cnv_visible(p["cn_call"], p["depth_log2"])]
        if any(p["depth_log2"] <= CNV_DEEP_LOG2 for p in parts):
            return "DEL", "depth only", f"GATK splits the gene ({words}); part of it has lost its reads"
        if seen:
            k = "DEL" if seen[0]["cn_call"] == "del" else "AMP"
            return k, "hidden", f"GATK splits the gene ({words}); the change covers part of the gene only"
        if called:
            return kind, "not supported", f"GATK splits the gene ({words}); the depth does not show the change"
        return "", "nothing", ""
    if r["state"] in ("del", "amp"):
        if r["status"] == "flagged" or not visible:
            if deep:   # a broad shallow segment, but the gene itself has lost its reads: a focal loss inside it
                return "DEL", "depth only", (f"GATK {r['call']} over the region, but the read depth in the gene is "
                                             f"log2 {d}: a deeper, focal loss")
            if not visible:
                shown = f" (log2 {d:+.2f})" if d is not None else ""
                return kind, "not supported", (f"GATK {r['call']}, but the read depth in the gene{shown} "
                                               f"is within the normal range")
            return kind, "not supported", (f"GATK {r['call']}, but the read depth in the gene measures log2 "
                                           f"{d:.2f} ({r['_depth']})" if d is not None else
                                           f"CNV {r['call']}: {why or 'depth disagrees'}")
        sg = (r.get("_segs") or [{}])[0]
        size = sg.get("end", 0) - sg.get("start", 0) + 1
        if size < HEATMAP_FOCAL_BP or abs(sg.get("log2", 0)) > HEATMAP_DEEP_LOG2:
            return kind, "real", f"CNV {r['call']}: depth agrees"
        return kind, "hidden", (f"CNV {r['call']} over {size / 1e6:.1f} Mb is real{f' (depth {d:+.2f})' if d is not None else ''} but broad "
                                f"(1 Mb or more) and not deep; the heatmap leaves these out by design (DEL/AMP "
                                f"only under 1 Mb or beyond log2 2)")
    if deep:   # no gain/loss segment, but the reads are gone (a gap between segments, or a missed focal loss)
        return "DEL", "depth only", "no GATK deletion call, but the read depth shows a deep loss"
    if mild:
        return kind, "hidden", (f"the depth shows a mild {r['_depth']} with no GATK call (on chrX, expected "
                                f"for a single X copy if the patient is male); the heatmap leaves these out by design")
    return "", "nothing", ""


def _evidence(r: dict) -> tuple[str, str, str]:
    """(alteration kind, what the reads say, words) for one summary row, in the heatmap's vocabulary.

    Verdicts: real, depth only (no call, but the reads are gone), not supported, weak, recurrent, and hidden
    (true in the reads, but a change the heatmap leaves out by design: broad and shallow copy number)."""
    flags = [f for f in r["flags"].split(";") if f]
    why = "; ".join(WHY[f] for f in flags if f in WHY)
    if r["type"] == "CNV":
        base = _cn_evidence(r, flags, why)
        st = r.get("_step")
        if not st or "ambiguous_mapping" in flags:
            return base
        chrom = r.get("_chrom", "")
        step_text = (f"the depth changes inside the gene at {chrom}:{st['pos']:,} ({st['left']:.0f}x to "
                     f"{st['right']:.0f}x; sample {r.get('_base', '?')}x), a copy-number breakpoint the gene average hides")
        if min(st.get("left_log2", 0), st.get("right_log2", 0)) <= CNV_DEEP_LOG2 and base[:2] != ("DEL", "real"):
            return "DEL", "depth only", "part of the gene has lost its reads: " + step_text
        if base[1] in ("hidden", "nothing") and r["state"].startswith(("no segment", "neutral")):
            return (base[0] or "DEL"), "hidden", step_text + "; the heatmap leaves changes like this out by design"
        return base[0], base[1], f"{base[2]}; {step_text}" if base[2] else step_text
    kind = "SV" if r["type"] == "SV" else "SNV"
    if "recurrent" in flags:
        return kind, "recurrent", f"{r['call']} is in the reads but in {r['_recurrent']}"
    if r["status"] == "supported":
        return kind, "real", f"{r['call']}: {r['igv_shows']}"
    if r["status"] == "flagged":
        return kind, "not supported", f"{r['call']}: {why}"
    return kind, "weak", f"{r['call']}: {why or 'too few reads'}"


def read_heatmap(path: Path) -> list[dict]:
    """Heatmap cells (sample, gene, alteration) from a CSV/TSV such as a heatmap exported from R."""
    p = Path(path)
    if not p.exists():
        raise InputError(f"--heatmap file not found: {p}")
    text = p.read_text()
    rows = list(csv.DictReader(text.splitlines(), delimiter="\t" if "\t" in text.splitlines()[0] else ","))
    if not rows or not {"sample", "gene", "alteration"} <= set(rows[0]):
        raise InputError(f"--heatmap needs columns sample, gene and alteration: {p}")
    return rows


def heatmap_vs_igv(rows: list[dict], heatmap: Path) -> list[dict]:
    """One line per heatmap cell for a gene IGV checked: heatmap label, what IGV found, matches/differs, why.

    Follows the heatmap's rules: one label per cell with SNV/SV before copy number, and copy number only when
    focal or deep. Changes left out by those rules are explained, not counted as differences."""
    checked = {r["gene"].upper() for r in rows}
    out = []
    for cell in read_heatmap(heatmap):
        smp, gene, label = cell["sample"], cell["gene"], (cell["alteration"] or "WT").strip()
        if gene.upper() not in checked:
            continue
        mine = [r for r in rows if r["sample"] == smp and r["gene"].upper() == gene.upper()]
        ev = [_evidence(r) for r in mine]
        ev = [x for x in ev if x[1] != "nothing"]
        wanted = set() if label.upper() == "WT" else set(label.upper().split("+"))
        have = {k for k, v, _ in ev if v in ("real", "depth only")}
        missing = sorted(wanted - have - {"LOH"})
        extra = sorted(have - wanted)
        # a cell shows SNV/SV before copy number, so a copy-number change under an SNV/SV label is not a miss
        hidden_cn = [k for k in extra if k in ("DEL", "AMP") and wanted & {"SNV", "SV"}]
        extra = [k for k in extra if k not in hidden_cn]
        words = lambda k, vs: "; ".join(t for kk, v, t in ev if kk == k and v in vs)
        why = [f"IGV supports the {NAMES.get(k, k)}" for k in sorted(wanted & have)]
        for k in missing:
            bad = words(k, ("not supported", "weak"))
            why.append(f"heatmap {k}, but " + (bad or "no such call was among the calls checked with IGV"))
        for k in extra:
            why.append(f"IGV supports a {NAMES.get(k, k)} the heatmap does not show: {words(k, ('real', 'depth only'))}")
        for k in hidden_cn:
            why.append(f"also a {NAMES[k]} ({words(k, ('real', 'depth only'))}); not shown because a heatmap cell "
                       f"gives SNV/SV priority over copy number")
        rec = [t for k, v, t in ev if v == "recurrent"]
        if rec:
            why.append("; ".join(rec) + ": independent tumors do not share mutations, so it is likely germline or "
                       "artifact" + (", and not showing it is right" if "SNV" not in wanted else ""))
        why += [t for k, v, t in ev if v == "hidden"]
        why += [f"{t}, so leaving it out is right" for k, v, t in ev if v == "not supported" and k not in wanted]
        why += [f"weak: {t}" for k, v, t in ev if v == "weak" and k not in missing]
        if "LOH" in wanted:
            why.append("LOH is not checked (it needs allele counts, not depth)")
        found = "; ".join(f"{k} {v}" for k, v, _ in ev) or "nothing called"
        match = "differs" if missing or extra else ("not checked" if wanted == {"LOH"} else "matches")
        nothing = "nothing was called here and the reads show no change either"
        details = ". ".join(why) or nothing
        # the answer in one line; everything IGV saw in the gene stays in `details`
        if match == "differs":
            short = ". ".join(w for w in why if w.startswith(("heatmap ", "IGV supports a ")))
        elif match == "not checked":
            short = "LOH is not checked (it needs allele counts, not depth)"
        elif wanted:
            short = "; ".join([f"IGV supports the {NAMES.get(k, k)}" for k in sorted(wanted & have)] +
                              [f"the {NAMES[k]} under it is not shown because SNV/SV takes priority in a cell"
                               for k in hidden_cn])
        else:
            reasons = []
            if rec:
                n = re.search(r"in (\d+ of \d+ samples)", rec[0])
                reasons.append(f"the SNV{'s' if len(rec) > 1 else ''} {'are' if len(rec) > 1 else 'is'} "
                               f"recurrent in {n.group(1) if n else 'several samples'} (likely germline or artifact)")
            if any(v == "hidden" for _, v, _ in ev):
                reasons.append("a broad or mild copy-number change, left out by design")
            for k, v, t in ev:
                if v == "not supported":
                    reasons.append(f"the {NAMES.get(k, k)} call is not supported: " +
                                   ("reads are present but ambiguous (MAPQ 0, alternate contigs)"
                                    if "alternate contigs" in t else "the reads disagree with it"))
            if any(v == "weak" for _, v, _ in ev):
                reasons.append("only weak read support for the calls here")
            short = ("WT is right: " + "; ".join(reasons)) if reasons else nothing
        review = (["the verdict differs from the heatmap"] if match == "differs" else []) + \
            [x for r in mine for x in review_reasons(r)] + (["LOH is not checked"] if "LOH" in wanted else [])
        out.append({"sample": smp, "gene": gene, "heatmap": label, "igv_found": found, "match": match,
                    "why": short or details, "details": details,
                    "check_image": _check(list(dict.fromkeys(review)),
                                          ["a copy-number change is left out or not visible; confirm on the depth plot"]
                                          if any(r["type"] == "CNV" and v in ("hidden", "not supported")
                                                 for r, (_, v, _) in zip(mine, [_evidence(x) for x in mine])) else None),
                    "report": next((r.get("report", "") for r in mine), "")})
    return out


OVERVIEW_PAD = 5000  # bp shown on each side of the gene


def overview_plan(rows: list[dict], regions: Path) -> list[dict]:
    """One overview image per (sample, gene in the BED): the gene, its calls marked, and a caption of verdicts."""
    spans = {}
    for chrom, items in _read_regions(regions).items():
        for a, b, name in items:
            c, lo, hi = spans.get(name.upper(), (chrom, a, b))
            spans[name.upper()] = (chrom, min(lo, a), max(hi, b))
    plan = []
    for smp, gene in sorted({(r["sample"], r["gene"].upper()) for r in rows}):
        if gene not in spans:
            continue
        chrom, a, b = spans[gene]
        mine = [r for r in rows if r["sample"] == smp and r["gene"].upper() == gene]
        bam = next((r["_bam"] for r in mine if r.get("_bam")), None)
        if not bam:
            continue
        beds, caption = [], [f"{smp}  ·  {gene}  {chrom}:{a:,}-{b:,}  (all calls in this gene)"]
        for r in mine:
            ag = _agreement(r)
            verdict = f"IGV agrees: {ag['igv_agrees']}" if ag else _no_call_note(r)
            if r["type"] == "CNV":
                for sg in r.get("_segs") or []:
                    beds.append(f"{sg['contig']}\t{sg['start'] - 1}\t{sg['end']}\tGATK {sg['cn_call']} log2 {sg['log2']:.2f}")
                caption.append(f"CNV {r['state']}: {verdict}")
            elif r["type"] == "SV":
                beds.append(f"{r['_chrom']}\t{r['_pos'] - 1}\t{r['_pos']}\tSV {r['_kind']} breakpoint")
                if r.get("_chrom2"):
                    beds.append(f"{r['_chrom2']}\t{r['_pos2'] - 1}\t{r['_pos2']}\tSV {r['_kind']} breakpoint 2")
                    if r["_chrom2"] == r["_chrom"] and r["_kind"] in ("del", "dup", "inv"):
                        lo, hi = sorted((r["_pos"], r["_pos2"]))
                        beds.append(f"{r['_chrom']}\t{lo - 1}\t{hi}\tSV {r['_kind']} span")
                caption.append(f"SV {r['call']}: {verdict}")
            else:
                beds.append(f"{r['_chrom']}\t{r['_pos'] - 1}\t{r['_pos']}\t{r['call']}")
                caption.append(f"SNV/indel {r['call']}: {verdict}")
        plan.append({"sample": smp, "gene": gene, "chrom": chrom, "start": max(1, a - OVERVIEW_PAD),
                     "end": b + OVERVIEW_PAD, "bam": bam, "ref": next((r["_ref"] for r in mine if r.get("_ref")), None),
                     "igv": next((r["_igv"] for r in mine if r.get("_igv")), None),
                     "bed_lines": beds, "caption": caption,
                     "png": re.sub(r"[^A-Za-z0-9._-]+", "_", f"{smp}_{gene}") + ".png"})
    return plan


def take_overviews(plan: list[dict], out: Path, igv_path, timeout: int = 300, annotation=None) -> tuple[dict, str]:
    """Draw the overview images with an isolated IGV; returns {(sample, gene): png path} and a note."""
    done, notes = {}, []
    argv, cwd, why = find_igv(overview_launcher(igv_path, plan))
    if argv is None:
        return done, why
    odir = out / "overview"; odir.mkdir(parents=True, exist_ok=True)
    tmp = Path(tempfile.mkdtemp(prefix="igv_overview_"))
    sub = annotation_subset(annotation, [(p["chrom"], p["start"], p["end"]) for p in plan], tmp / "_genes") \
        if annotation else None
    for ref in sorted({p["ref"] for p in plan if p["ref"]}):
        try:
            igv = IGVSession(argv, cwd, Path(ref), timeout, out / "overview" / "igv.log")
        except RuntimeError as e:
            notes.append(str(e)); continue
        try:
            for p in (x for x in plan if x["ref"] == ref):
                bed = tmp / f"{p['png'][:-4]}_calls.bed"
                bed.write_text("\n".join(p["bed_lines"]) + ("\n" if p["bed_lines"] else ""))
                span_kb = (p["end"] - p["start"]) // 1000 + 20
                locus = f"{p['chrom']}:{p['start']}-{p['end']}"
                for c in ["new", f"snapshotDirectory {tmp}", "preference SAM.DOWNSAMPLE_READS true",
                          f"preference SAM.MAX_VISIBLE_RANGE {span_kb}", "maxPanelHeight 450",
                          f"load {p['bam']} name={p['sample']}"] + ([f"load {bed} name=calls", "expand calls"] if p["bed_lines"] else []) + annotation_cmds(sub) \
                        + [f"goto {locus}"]:
                    igv.send(c)
                png = tmp / p["png"]
                for attempt in range(3):
                    time.sleep(1.5 * (attempt + 1))
                    igv.send(f"goto {locus}"); igv.send(f"collapse {p['sample']}")   # reads only: calls stay expanded
                    if p["bed_lines"]:
                        igv.send("expand calls")
                    if sub:
                        igv.send("expand genes")
                    igv.send(f"snapshot {png.name}")
                    if png.exists() and image_rendered(png):
                        break
                dst = odir / p["png"]
                caption(png, dst, p["caption"][:9] + (["...more calls in summary.tsv"] if len(p["caption"]) > 9 else []))
                done[(p["sample"], p["gene"])] = dst
        except RuntimeError as e:
            notes.append(f"stopped early: {e}")
        finally:
            igv.close()
    shutil.rmtree(tmp, ignore_errors=True)
    return done, "; ".join(notes)


def summarize(root: Path, out: Path | None = None, regions: Path | None = None, igv_path=None,
              timeout: int = 300, heatmap: Path | None = None, annotation: Path | None = None) -> Path:
    """One table and one grid (sample x gene) over every igv-validator run under `root`."""
    root = Path(root)
    if not root.is_dir():
        raise InputError(f"--summarize folder not found: {root}")
    empty_runs = []
    rows = dedupe_svs(_summary_rows(root, empty_runs))
    mark_recurrent(rows)
    if not rows and not empty_runs:
        raise InputError(f"no igv-validator result.json files under {root}; run the checks first")
    out = Path(out) if out else root
    out.mkdir(parents=True, exist_ok=True)
    cols = ["sample", "gene", "type", "tool", "state", "call", "location", "reads", "igv_shows", "caller", "flags",
            "status", "concordance", "report", "figures"]
    with open(out / "summary.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore"); w.writeheader(); w.writerows(rows)
    overviews, ov_note = ({}, "")
    if regions:
        overviews, ov_note = take_overviews(overview_plan(rows, regions), out, igv_path, timeout, annotation)
    agree = [a for a in (_agreement(r) for r in rows) if a]
    cells = heatmap_vs_igv(rows, heatmap) if heatmap else []
    if heatmap:
        with open(out / "heatmap_vs_igv.tsv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=["sample", "gene", "heatmap", "igv_found", "match", "check_image", "why", "details",
                                           "report"],
                               delimiter="\t")
            w.writeheader(); w.writerows(cells)
    acols = ["sample", "gene", "type", "tool", "caller_called", "igv_shows", "igv_agrees", "check_image", "why",
             "report"]
    with open(out / "igv_agreement.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=acols, delimiter="\t"); w.writeheader(); w.writerows(agree)
    samples = sorted({r["sample"] for r in rows})
    genes = sorted({r["gene"] for r in rows})
    e = html.escape
    rel = lambda p: e(str(Path(os.path.relpath(root / p, out))))
    acolor = {"yes": "#c6e8d2", "no": "#f3c7a8", "unclear": "#f6e3a1", "in reads, recurrent": "#dcd6ec"}
    counts = {k: sum(a["igv_agrees"] == k for a in agree) for k in acolor}
    mcolor = {"matches": "#c6e8d2", "differs": "#f3c7a8", "not checked": "#ecebe7"}
    ovl = lambda smp, g: (f" · <a href='{e(str(Path(os.path.relpath(overviews[(smp, g.upper())], out))))}'>overview</a>"
                          if (smp, g.upper()) in overviews else "")
    chk = lambda c: (f"<td style='background:#fbe7a6'><b>yes</b><br><small>{e(c[5:])}</small></td>"
                     if c.startswith("yes") else
                     f"<td style='background:#fdf3d0'>glance<br><small>{e(c[8:])}</small></td>"
                     if c.startswith("glance") else "<td>no</td>")
    warning = ("<div style='border:2px solid #d9a400;background:#fff6d6;padding:10px 14px;margin:12px 0'>"
               "<b>The verdicts below are automatic and can be wrong</b>, especially for copy number (reference "
               "duplications, GC-rich DNA, chromosome ends, noisy samples). Check the image before reporting any "
               "<i>differs</i> row and any row marked <i>Check image: yes</i>; the IGV screenshots and depth plots "
               "are the evidence, the verdict only points you to them.</div>"
               "<div style='border:1px solid #ccc;padding:8px 14px;margin:0 0 12px'><b>How to read this page</b>"
               "<table style='margin-top:6px'><tr><th>You see</th><th>What to do</th></tr>"
               "<tr><td>SNV/SV <i>matches</i>, Check image <i>no</i></td><td>trust it</td></tr>"
               "<tr><td>Deep deletion or clear gain <i>matches</i>, Check image <i>no</i></td><td>trust it</td></tr>"
               "<tr><td>Any <i>differs</i>, or Check image <i>yes</i></td><td>open the IGV image and the depth plot "
               "before reporting; the reason says what to look for</td></tr>"
               "<tr><td>Check image <i>glance</i> (copy number left out, not visible, or unclear)</td>"
               "<td>a quick look at the depth plot: is the blue line where the verdict says?</td></tr>"
               "<tr><td>Anything you will present or publish</td><td>look at the image yourself; never cite the "
               "verdict alone</td></tr></table></div>")
    heat_html = "" if not heatmap else (
        f"<h2>Your heatmap vs IGV</h2><p><b>{sum(c['match'] == 'matches' for c in cells)} match, "
        f"{sum(c['match'] == 'differs' for c in cells)} differ</b> out of {len(cells)} heatmap cell(s) for the genes "
        f"checked. A WT cell can still have a variant in the reads: that matches when the variant is not a tumor "
        f"mutation (for example the same variant in several samples).</p>"
        "<div class=w><table><tr><th>Sample</th><th>Gene</th><th>Heatmap says</th><th>IGV found</th>"
        "<th>Match?</th><th>Check image?</th><th>Why</th><th>Evidence</th></tr>" + "".join(
            f"<tr><td>{e(c['sample'])}</td><td><b>{e(c['gene'])}</b></td><td>{e(c['heatmap'])}</td>"
            f"<td>{e(c['igv_found'])}</td><td style='background:{mcolor.get(c['match'], '#fff')}'><b>"
            f"{e(c['match'])}</b></td>" + chk(c["check_image"]) + f"<td>{e(c['why'])}" +
            (f"<details><summary>details</summary>{e(c['details'])}</details>" if c["details"] != c["why"] else "") +
            "</td><td>" + (f"<a href='{rel(c['report'])}'>report</a>" if c["report"] else "") +
            ovl(c["sample"], c["gene"]) + "</td></tr>"
            for c in sorted(cells, key=lambda c: (c["match"] != "differs", c["gene"], c["sample"]))) + "</table></div>")
    agree_html = warning + heat_html + ((f"<p><i>Overview images: {e(ov_note)}</i></p>" if ov_note else "") +
                  f"<h2>Does IGV agree with the callers?</h2><p><b>{counts['yes']} yes, {counts['no']} no, "
                  f"{counts['unclear']} unclear, {counts['in reads, recurrent']} in reads but recurrent</b> out of {len(agree)} call(s) the workflow made "
                  f"(SNVs/indels, SVs, and copy-number gains/losses; neutral copy number is not a call). "
                  f"<i>Yes</i> means the reads contain what the caller called, not that it is a tumor mutation; "
                  f"<i>in reads, recurrent</i> means the same variant is in {RECURRENT_MIN}+ samples, "
                  f"so it is likely germline or artifact.</p>"
                  "<div class=w><table><tr><th>Sample</th><th>Gene</th><th>Tool</th><th>The caller called</th>"
                  "<th>What IGV shows</th><th>IGV agrees?</th><th>Check image?</th><th>Why</th><th>Evidence</th></tr>" + "".join(
                      f"<tr><td>{e(a['sample'])}</td><td><b>{e(a['gene'])}</b></td><td>{e(a['tool'])}</td>"
                      f"<td>{e(a['caller_called'])}</td><td>{e(a['igv_shows'])}</td>"
                      f"<td style='background:{acolor[a['igv_agrees']]}'><b>{e(a['igv_agrees'])}</b></td>"
                      + chk(a["check_image"]) + f"<td>{e(a['why'])}</td><td><a href='{rel(a['report'])}'>report</a>"
                      + (f" · <a href='{e(str(Path(os.path.relpath(overviews[(a['sample'], a['gene'].upper())], out))))}'>"
                         f"overview</a>" if (a["sample"], a["gene"].upper()) in overviews else "") + "</td></tr>" for a in
                      sorted(agree, key=lambda a: ({"no": 0, "in reads, recurrent": 1, "unclear": 2, "yes": 3}[a["igv_agrees"]], a["gene"],
                                                    a["sample"]))) + "</table></div>")
    grid = "<tr><th>Gene</th>" + "".join(f"<th>{e(s)}</th>" for s in samples) + "</tr>"
    for g in genes:
        grid += f"<tr><th>{e(g)}</th>"
        for smp in samples:
            cell = [r for r in rows if r["gene"] == g and r["sample"] == smp]
            if not cell:
                grid += "<td class=empty>-</td>"
                continue
            # a neutral copy number that the reads agree with is not a finding: shown faintly
            quiet = lambda r: r["type"] == "CNV" and r["state"] in ("neutral", "no segment; depth neutral") \
                and r["concordance"] in ("confirmed", "no call")
            loud = [r for r in cell if not quiet(r)]
            label = lambda r: (f"CNV {r['state']}" if r["type"] == "CNV" else
                               f"SV {r['state']}" if r["type"] == "SV" else r["type"])
            items = "".join(f"<div class='{'quiet' if quiet(r) else ''}'><a href='{rel(r['report'])}'>{e(label(r))}</a>: "
                            f"{e(r['concordance'])}</div>" for r in cell)
            if loud:
                worst = next(k for k in CONCORDANCE_ORDER if any(r["concordance"] == k for r in loud) or k == "no call")
                bg = CONCORDANCE_COLOR.get(worst, "#fff")
            else:
                bg = "#f7f6f3"
            ov = overviews.get((smp, g.upper()))
            if ov:
                items += f"<div><a href='{e(str(Path(os.path.relpath(ov, out))))}'><b>overview image</b></a></div>"
            grid += f"<td style='background:{bg}'>{items}</td>"
        grid += "</tr>"
    detail = "".join(
        f"<tr><td>{e(r['sample'])}</td><td><b>{e(r['gene'])}</b></td><td>{e(r['type'])}</td><td>{e(r['state'])}</td><td>{e(r['call'])}</td>"
        f"<td>{e(r['reads'])}</td><td>{e(r['caller'])}</td><td>{e(r['flags'] or 'none')}</td>"
        f"<td style='background:{CONCORDANCE_COLOR.get(r['concordance'], '#fff')}'>{e(r['concordance'])}</td>"
        f"<td><a href='{rel(r['report'])}'>report</a>"
        + "".join(f" · <a href='{rel(f)}'>image</a>" for f in r["figures"].split(";") if f) + "</td></tr>"
        for r in rows)
    legend = " ".join(f"<span class=chip style='background:{CONCORDANCE_COLOR[k]}'>{k}</span>" for k in CONCORDANCE_ORDER)
    page = (f"<!doctype html><html lang=en><head><meta charset=utf-8><meta name=viewport content='width=device-width,"
            f"initial-scale=1'><title>IGV Validation Summary</title><style>body{{font:14px/1.45 system-ui,sans-serif;"
            f"max-width:1300px;margin:0 auto;padding:24px 16px;background:#fbfaf8;color:#1d1d1b}}table{{border-collapse:"
            f"collapse;margin:8px 0 24px}}th,td{{border:1px solid #e4e1db;padding:6px 8px;text-align:left;"
            f"vertical-align:top}}td.empty{{color:#aaa;text-align:center}}.chip{{padding:2px 8px;border-radius:9px;"
            f"margin-right:4px}}a{{color:inherit}}.w{{overflow-x:auto}}.quiet{{color:#9a988f}}</style></head><body>"
            f"<h1>IGV validation summary</h1><p>{len(rows)} call(s) · {len(samples)} sample(s) · {len(genes)} gene(s) "
            f"· {e(datetime.now(timezone.utc).strftime('%Y-%m-%d'))}</p><p>{legend}</p>"
            + agree_html +
            f"<h2>Details</h2><p>Each call from the workflow is compared with the reads: <b>confirmed</b> = the reads support it "
            f"(status supported); <b>questioned</b> = the reads disagree (e.g. caller_disagrees, cnv_disagrees, "
            f"artifact flags); <b>weak</b> = too few supporting reads; <b>recurrent</b> = the same SNV/indel in several samples (likely germline or artifact); <b>no call</b> = no segment covers the gene "
            f"(depth only; the cell shows what the depth says); <b>no reads</b> = no coverage. Copy number is "
            f"checked for every gene: a neutral state the reads agree with is shown in grey, since it is not a change. "
            f"'Confirmed' means the reads agree with the call, not that a variant is somatic (tumor-only data "
            f"cannot separate germline). These are rule-based summaries, not verdicts.</p>"
            + (f"<p><b>Runs with no calls in the selected genes:</b> "
               + ", ".join(f"<a href='{rel(rep_)}'>{e(smp)} ({e(kind)})</a>" for smp, kind, rep_ in empty_runs)
               + "</p>" if empty_runs else "") +
            f"<h2>Details: samples x genes</h2><div class=w><table>{grid}</table></div>"
            f"<h2>All calls</h2><div class=w><table><tr><th>Sample</th><th>Gene</th><th>Type</th><th>State</th><th>Call</th>"
            f"<th>Reads (BAM)</th><th>Caller</th><th>Flags</th><th>Concordance</th><th>Evidence</th></tr>{detail}"
            f"</table></div><p><i>{e(DISCLAIMER)}</i></p></body></html>")
    (out / "summary.html").write_text(page)
    return out


# ── Main ───────────────────────────────────────────────────────────────────────

def _check_inputs(vcf, tumor, normal, reference):
    files = [("VCF", vcf), ("tumor BAM", tumor)] + ([("normal BAM", normal)] if normal else [])
    for label, p in files:
        if not p or not Path(p).exists():
            raise InputError(f"{label} not found: {p}")
    for label, p in files[1:]:
        p = Path(p)
        idx = [Path(str(p) + s) for s in (".bai", ".csi", ".crai")] + [p.with_suffix(".bai")]
        if not any(i.exists() for i in idx):
            raise InputError(f"no index for the {label} {p} (run: samtools index {p})")
    if reference and not Path(str(reference) + ".fai").exists():
        raise InputError(f"reference {reference} has no .fai index (run: samtools faidx {reference})")


def _check_contigs(variants, bam_path, fasta):
    with _open_bam(bam_path, fasta) as bam:
        names = set(bam.references)
    need = {v.chrom for v in variants} | {v.chrom2 for v in variants if v.chrom2}
    missing = sorted(need - names)
    if missing:
        hint = " (the VCF and BAM disagree on the 'chr' prefix)" if any(
            ("chr" + m in names) or (m.removeprefix("chr") in names) for m in missing) else ""
        raise InputError(f"contig(s) {', '.join(missing)} are in the VCF but not in {bam_path}{hint}")


def _check_reference(variants, bam_path, fasta):
    """The FASTA must be the one the BAMs were aligned to: same contig names and lengths."""
    need = {v.chrom for v in variants} | {v.chrom2 for v in variants if v.chrom2}
    with _open_bam(bam_path, fasta) as bam, pysam.FastaFile(str(fasta)) as fa:
        for c in sorted(need):
            if c not in fa.references:
                raise InputError(f"reference {fasta} has no contig {c}; pass the FASTA the BAMs were aligned to")
            if fa.get_reference_length(c) != bam.get_reference_length(c):
                raise InputError(f"reference {fasta} has {c} of length {fa.get_reference_length(c):,} but the BAM "
                                 f"has {bam.get_reference_length(c):,}; pass the FASTA the BAMs were aligned to")


def _bam_sample(bam_path, fasta) -> str | None:
    with _open_bam(bam_path, fasta) as bam:
        rg = bam.header.to_dict().get("RG", [])
    names = {r.get("SM") for r in rg if r.get("SM")}
    return names.pop() if len(names) == 1 else None


def guess_tumor_sample(samples: list[str]) -> str | None:
    """The one VCF sample whose name says 'tumor', if exactly one does; otherwise None (ask the user)."""
    hits = [x for x in samples if "tumor" in x.lower() or "tumour" in x.lower()]
    return hits[0] if len(hits) == 1 else None


def choose_caller_sample(vcf, tumor_bam, fasta, wanted: str | None) -> tuple[str | None, str]:
    """Which VCF sample column holds the tumor's caller counts: (sample, note)."""
    samples = vcf_samples(vcf)
    if wanted:
        if wanted not in samples:
            raise InputError(f"--caller-sample {wanted} is not in the VCF (samples: {', '.join(samples) or 'none'})")
        return wanted, ""
    sm = _bam_sample(tumor_bam, fasta)
    if sm in samples:
        return sm, ""
    if len(samples) == 1:
        return samples[0], ""
    guess = guess_tumor_sample(samples)
    if guess:
        return guess, f"compared with VCF sample {guess}, chosen by its name (use --caller-sample to change)"
    if not samples:
        return None, "the VCF has no sample columns, so there are no caller counts to compare"
    return None, (f"could not tell which VCF sample is the tumor ({', '.join(samples)}); "
                  "pass --caller-sample to compare caller counts")


def swap_warning(results: dict) -> str:
    """Warn when the 'normal' carries most of the support: --tumor and --normal were probably swapped."""
    if any(r["normal"] is None for r in results.values()):
        return ""  # tumor-only run
    informative = [r for r in results.values() if r["status"] != "no_coverage" and
                   (r["tumor"]["alt"] + r["normal"]["alt"]) >= MIN_SUPPORT]
    flipped = [r for r in informative if r["normal"]["alt"] >= MIN_SUPPORT and
               (r["normal"]["vaf_pct"] or 0) > 2 * (r["tumor"]["vaf_pct"] or 0)]
    if len(informative) >= 2 and len(flipped) * 2 > len(informative):
        return (f"the normal has much more support than the tumor for {len(flipped)} of {len(informative)} variants; "
                "were --tumor and --normal swapped?")
    return ""


def run(args) -> Path:
    if pysam is None:
        raise InputError("pysam is not installed (pip install pysam)")
    if args.demo_cnv:
        args.cnv = args.cnv or DEMO_DIR / "demo_cnv.tsv"
        args.tumor = args.tumor or DEMO_DIR / "demo_tumor.bam"
        args.regions = args.regions or DEMO_DIR / "demo_cnv_genes.bed"
        args.reference = args.reference or DEMO_DIR / "demo_ref.fa"
    if args.cnv:
        return run_cnv(args)
    if args.demo:
        args.vcf = args.vcf or DEMO_DIR / "demo_calls.vcf"
        args.tumor = args.tumor or DEMO_DIR / "demo_tumor.bam"
        if not args.tumor_only:
            args.normal = args.normal or DEMO_DIR / "demo_normal.bam"
        args.reference = args.reference or DEMO_DIR / "demo_ref.fa"
    if args.tumor_only:
        args.normal = None
    if not (args.vcf and args.tumor):
        raise InputError("--vcf and --tumor are required (add --normal for a tumor/normal pair), or use --demo")
    out = Path(args.output or f"output/igv_validator_{datetime.now():%Y%m%d_%H%M%S}")
    _check_inputs(args.vcf, args.tumor, args.normal, args.reference)
    genes = {g.strip() for g in args.genes.split(",") if g.strip()} if args.genes else None
    if args.regions and not Path(args.regions).exists():
        raise InputError(f"--regions file not found: {args.regions}")
    variants, skipped = load_variants(args.vcf, genes, args.variants, args.pass_only, args.regions,
                                      allow_empty=bool(args.regions))
    if len(variants) > args.max_variants:
        for v in variants[args.max_variants:]:
            skipped.append({"id": v.id, "chrom": v.chrom, "pos": v.pos,
                            "reason": f"over --max-variants {args.max_variants}"})
        variants = variants[:args.max_variants]
    _check_contigs(variants, args.tumor, args.reference)
    if args.normal:
        _check_contigs(variants, args.normal, args.reference)
    if args.reference:
        _check_reference(variants, args.tumor, args.reference)
    if out.exists() and any(out.iterdir()):
        print(f"Warning: {out} is not empty; igv-validator outputs there will be overwritten", file=sys.stderr)
    out.mkdir(parents=True, exist_ok=True)

    mode = "tumor_normal" if args.normal else "tumor_only"
    names = (args.tumor_name or Path(args.tumor).stem,
             (args.normal_name or Path(args.normal).stem) if args.normal else None)
    caller_sample, caller_note = choose_caller_sample(args.vcf, args.tumor, args.reference, args.caller_sample)
    results = {}
    for v in variants:
        t = count_variant(v, args.tumor, args.reference)
        n = count_variant(v, args.normal, args.reference) if args.normal else None
        c = v.calls.get(caller_sample) if caller_sample else None
        caller = {"alt": c[0], "depth": c[1], "vaf_pct": _pct(c[0], c[1]), "sample": caller_sample} if c else None
        flags, status = flag_variant(v, t, n, caller)
        results[v.id] = {"variant": v, "tumor": t, "normal": n, "caller": caller, "flags": flags, "status": status,
                         "screenshots": []}
    flag_depth_outliers(results)
    swapped = swap_warning(results)

    taken, note = 0, "--no-igv given"
    if not variants:
        note = "no calls to show"
    elif not args.no_igv:
        if not args.reference:
            note = "no --reference given (IGV must use the reference the BAMs were aligned to)"
        else:
            raw = Path(tempfile.mkdtemp(prefix="igv_shots_"))
            taken, note = take_screenshots(variants, results, args.tumor, args.normal, Path(args.reference).resolve(),
                                           raw, args.igv_path, names, out / "reproducibility" / "igv.log",
                                           args.igv_timeout, getattr(args, "annotation", None))
            figdir = out / "figures" / "igv"
            if taken:
                figdir.mkdir(parents=True, exist_ok=True)
            for vid, r in results.items():
                v, done = r["variant"], []
                for png in r["screenshots"]:
                    lines = [f"{v.id}  {v.gene or ''}  {v.label}",
                             f"{names[0]}: {_support_text(r['tumor'], v.is_sv)}{_strand_text(r['tumor'])}"]
                    if r["normal"] is not None:
                        lines.append(f"{names[1]}: {_support_text(r['normal'], v.is_sv)}{_strand_text(r['normal'])}")
                    if r["caller"]:
                        lines.append(f"Caller reported (VCF {r['caller']['sample']}): {_support_text(r['caller'], False)}")
                    lines.append(f"Status: {r['status']}{' (tumor-only)' if r['normal'] is None else ''}. "
                                 "Counts come from the BAM (pysam), not from the image.")
                    caption(png, figdir / png.name, lines)
                    done.append(figdir / png.name)
                r["screenshots"] = done
            shutil.rmtree(raw, ignore_errors=True)

    rows = []
    for vid, r in results.items():
        v = r["variant"]
        d = {k: val for k, val in asdict(v).items() if k not in ("genes", "calls")}
        rows.append({**d, "label": v.label, "sv": v.is_sv, "tumor": r["tumor"], "normal": r["normal"],
                     "caller": r["caller"], "flags": r["flags"], "status": r["status"],
                     "figures": [str(p.relative_to(out)) for p in r["screenshots"]]})
    meta = {"vcf": Path(args.vcf).name, "tumor_name": names[0], "normal_name": names[1],
            "date": datetime.now(timezone.utc).strftime("%Y-%m-%d"), "screenshots": {"taken": taken, "note": note},
            "swap_warning": swapped, "mode": mode, "caller_sample": caller_sample, "caller_note": caller_note}
    written = write_reports(out, rows, skipped, meta)
    summary = {"variants_checked": len(rows), **{s: sum(r["status"] == s for r in rows)
                                                 for s in ("supported", "flagged", "insufficient", "no_coverage")},
               "not_evaluated": len(skipped), "screenshots_taken": taken, "screenshot_note": note,
               "swap_warning": swapped, "mode": mode, "caller_sample": caller_sample,
               "caller_disagrees": sum("caller_disagrees" in r["flags"] for r in rows)}
    thresholds = {k: globals()[k] for k in ("MIN_MAPQ", "MIN_BQ", "INDEL_WINDOW", "MIN_CLIP", "SPLIT_WINDOW",
                                            "PAIR_WINDOW", "MIN_SV_LEN", "MIN_INV_LEN", "MIN_SUPPORT", "MIN_VAF", "STRAND_MIN",
                                            "READ_END_BP", "LOW_DEPTH", "GERMLINE_OTHER", "CALLER_VAF_DIFF",
                                            "HIGH_DEPTH_X")}
    vcf_sha = hashlib.sha256(Path(args.vcf).read_bytes()).hexdigest()
    written.append(write_result_json(out, "igv-validator", VERSION, summary,
                                     {"variants": [{k: val for k, val in r.items() if k != "sv"} for r in rows],
                                      "not_evaluated": skipped, "thresholds": thresholds,
                                      "samples": {"tumor": names[0], "normal": names[1]},
                                      "tool": vcf_tool(args.vcf), "inputs": _inputs(args)},
                                     input_checksum=vcf_sha))
    write_commands_sh(out, "python skills/igv-validator/igv_validator.py " + shlex.join(sys.argv[1:]))
    write_environment_yml(out, "clawbio-igv-validator", ["pysam>=0.22", "Pillow>=10.1"])
    figs = [f for r in rows for f in r["figures"]]
    write_checksums(written + [out / f for f in figs], out, anchor=out)
    return out


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Validate somatic variant calls against tumor/normal reads, with IGV screenshots.")
    p.add_argument("--vcf", "--input", dest="vcf", help="VCF/VCF.gz of calls (any caller)")
    p.add_argument("--tumor", help="tumor BAM/CRAM (indexed)")
    p.add_argument("--normal", help="matched normal BAM/CRAM (indexed); omit for tumor-only mode")
    p.add_argument("--tumor-only", action="store_true", help="ignore any normal (e.g. with --demo)")
    p.add_argument("--caller-sample", help="VCF sample whose caller counts (AD) to compare (default: the tumor BAM's SM)")
    p.add_argument("--reference", help="indexed FASTA the BAMs were aligned to (screenshots, indel clip rescue, CRAM)")
    p.add_argument("--genes", help="comma-separated gene symbols to check (needs gene names in the VCF)")
    p.add_argument("--variants", help="TSV of chrom pos [ref alt] to check")
    p.add_argument("--regions", help="BED of regions to check (chrom start end [name]); names label the genes, "
                                     "and an SV matches if either breakpoint falls inside")
    p.add_argument("--pass-only", action="store_true", help="only calls with FILTER=PASS")
    p.add_argument("--max-variants", type=int, default=MAX_VARIANTS, help=f"cap (default {MAX_VARIANTS})")
    p.add_argument("--no-igv", action="store_true", help="counts and report only, no screenshots")
    p.add_argument("--igv-timeout", type=int, default=300,
                   help="seconds to wait for IGV to start and load the reference (default 300)")
    p.add_argument("--annotation", help="gene annotation (GTF/GFF3/BED, optionally .gz, e.g. GENCODE for hg38) "
                                         "drawn as a 'genes' track in every IGV screenshot")
    p.add_argument("--igv-path", help="IGV .app bundle (macOS) or igv.sh (Linux); found automatically if omitted")
    p.add_argument("--tumor-name", help="label for the tumor track (default: file name)")
    p.add_argument("--normal-name", help="label for the normal track (default: file name)")
    p.add_argument("--output", help="output directory")
    p.add_argument("--demo", action="store_true", help="run on the bundled synthetic demo data")
    p.add_argument("--cnv", help="copy-number mode: segments TSV (cnv_extracted format or GATK called.seg); "
                                 "needs --tumor and --regions")
    p.add_argument("--cnv-sample", help="sample to take from a multi-sample --cnv file")
    p.add_argument("--demo-cnv", action="store_true", help="run copy-number mode on the synthetic demo")
    p.add_argument("--summarize", help="write summary.html/summary.tsv over every igv-validator run in this folder")
    p.add_argument("--heatmap", help="with --summarize: your heatmap table (CSV/TSV with sample, gene, alteration, e.g. "
                                      "heatmap.csv) to compare cell by cell with IGV")
    p.add_argument("--overview", action="store_true",
                   help="with --summarize and --regions: one IGV image per gene and sample with every call marked")
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    try:
        if getattr(args, "annotation", None) and not Path(args.annotation).exists():
            raise InputError(f"--annotation file not found: {args.annotation}")
        if args.summarize:
            if args.overview and not args.regions:
                raise InputError("--overview needs --regions (the BED of genes to draw)")
            out = summarize(Path(args.summarize), args.output, args.regions if args.overview else None,
                            args.igv_path, args.igv_timeout, args.heatmap, args.annotation)
            print(f"Summary written to {out / 'summary.html'}")
            return 0
        out = run(args)
    except InputError as e:
        print(f"Error: {e}", file=sys.stderr)
        return 2
    print(f"Report written to {out / 'report.md'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
