#!/usr/bin/env python3
"""Generate the synthetic igv-validator demo: reference, tumor/normal BAMs, VCF and truth.

Synthetic demo data for igv-validator. This is NOT real patient data: the reference is
random sequence on two made-up contigs and every read is simulated.

Each planted variant has an exactly known number of supporting reads, written to
demo_truth.json so the tests can check the counters against ground truth:

  V1  demo1:3000  C>T-style SNV, 18 alt reads, both strands            -> supported
  V2  demo1:5000  SNV, 6 alt reads, all reverse strand                  -> strand_bias
  V3  demo1:7000  SNV, 2 alt reads, FILTER=weak_evidence                -> insufficient
  V4  demo1:9000  12 bp deletion, 11 spanning + 2 right- + 2 left-clipped reads
  V5  demo1:12000 <-> demo2:5001 translocation, 6 split reads + 8 discordant pairs
  V6  demo1:15000 SNV, heterozygous in tumor AND normal                 -> normal_support

Copy-number demo (contig demo3, its own RNG so the reads above never change): the tumor has a
homozygous deletion (20,001-25,000: no MAPQ>=20 reads, plus 20 low-MAPQ "mouse-like" reads), a
one-copy loss (35,001-45,000: half coverage) and normal coverage elsewhere. demo_cnv.tsv holds
GATK-style segments in the igv-validator/cnv_extracted format, including one wrong call (a
deletion at 5,000-8,000 where coverage is normal) and a gap (49,001-53,000) with no segment.
demo_cnv_genes.bed names six test genes across those regions. ALTGENE (55,501-57,500) imitates a gene with alternate-haplotype copies:
normal coverage, but every tumor read there has MAPQ 0 and an XA hit on an alternate contig, and the
segment file calls a deep deletion (as a MAPQ-filtering caller would).

The VCF carries caller-style AD counts equal to the truth, except V3, whose AD is deliberately
inflated (12 alt reads claimed, 2 present) so the caller_disagrees check has something to catch.

Run: python make_demo_data.py   (rewrites the files next to this script; seeded, deterministic)
"""
import json
import random
from pathlib import Path

import pysam

HERE = Path(__file__).resolve().parent
RNG = random.Random(20260926)
READ, INSERT, DEPTH_FRAGS, BQ = 100, 350, 380, 30
CONTIGS = {"demo1": 20000, "demo2": 10000}
COMMENT = "Synthetic demo data for igv-validator. This is NOT real patient data."

ref = {c: "".join(RNG.choice("ACGT") for _ in range(n)) for c, n in CONTIGS.items()}
tid = {c: i for i, c in enumerate(CONTIGS)}
# copy-number contig: generated with its own RNG after the SNV/SV contigs, so their reads are unchanged
CNV_CONTIG, CNV_LEN, CNV_DEPTH_FRAGS = "demo3", 60000, 6000
RNG_CNV = random.Random(20260927)
ref[CNV_CONTIG] = "".join(RNG_CNV.choice("ACGT") for _ in range(CNV_LEN))
CONTIGS[CNV_CONTIG] = CNV_LEN
tid[CNV_CONTIG] = len(tid)
HOMDEL, HEMI = (20001, 25000), (35001, 45000)
CNV_SEGMENTS = [  # (start, end, log2, call): igv-validator reads these like GATK called.seg + modelFinal
    (1, 4999, 0.02, "0"), (5000, 8000, -1.5, "-"), (8001, 19999, 0.01, "0"), (20000, 25000, -8.0, "-"),
    (25001, 34999, -0.02, "0"), (35000, 45000, -1.0, "-"), (45001, 49000, 0.03, "0"), (53001, 55000, 0.0, "0"),
    (55001, 58000, -7.0, "-"), (58001, 60000, 0.0, "0")]
ALTREGION = (55001, 58000)
CNV_GENES = [("WRONGDEL", 5500, 7500), ("NEUTRAL", 10000, 15000), ("HOMDEL", 21000, 24000),
             ("HEMI", 36000, 44000), ("GAPGENE", 50000, 52000), ("ALTGENE", 55500, 57500)]  # 1-based inclusive


def other_base(b):
    return RNG.choice([x for x in "ACGT" if x != b])


SNVS = {  # id: (pos 1-based, alt base, tumor alt reads, which strands, gene)
    "V1": (3000, None, 18, "both", "DEMO1"),
    "V2": (5000, None, 6, "reverse", "DEMO2"),
    "V3": (7000, None, 2, "both", "DEMO3"),
}
GERM = ("V6", 15000, "DEMO6")
DEL = ("V4", 9000, 12, "DEMO4")  # anchor pos 1-based, deleted length
BND = ("V5", ("demo1", 12000), ("demo2", 5001), "DEMO5")
SNVS = {k: (p, other_base(ref["demo1"][p - 1]), n, s, g) for k, (p, _, n, s, g) in SNVS.items()}
GERM_ALT = other_base(ref["demo1"][GERM[1] - 1])
WINDOWS = [("demo1", p) for p, *_ in SNVS.values()] + [("demo1", GERM[1]), ("demo1", DEL[1]),
                                                        BND[1], BND[2]]


def seg(qname, chrom, start, seq, cigar, reverse, read1, mate_chrom, mate_start, tlen, extra_flags=0, tags=()):
    a = pysam.AlignedSegment()
    a.query_name, a.query_sequence = qname, seq
    a.flag = 1 | 2 | (16 if reverse else 0) | (32 if not reverse else 0) | (64 if read1 else 128) | extra_flags
    if mate_chrom != chrom:
        a.flag &= ~2  # not a proper pair when mates are on different contigs
    a.reference_id, a.reference_start, a.cigarstring = tid[chrom], start, cigar
    a.mapping_quality = 60
    a.next_reference_id, a.next_reference_start, a.template_length = tid[mate_chrom], mate_start, tlen
    a.query_qualities = pysam.qualitystring_to_array(chr(BQ + 33) * len(seq))
    a.set_tags(list(tags))
    return a


def background(sample, reads):
    """Reference reads around every window: ~60x at the variant position."""
    n = 0
    for chrom, pos in WINDOWS:
        for _ in range(DEPTH_FRAGS):
            s = RNG.randint(pos - 800, pos + 800 - INSERT)
            q = f"{sample}_bg{n}"; n += 1
            r1 = seg(q, chrom, s, ref[chrom][s:s + READ], f"{READ}M", False, True, chrom, s + INSERT - READ, INSERT)
            r2s = s + INSERT - READ
            r2 = seg(q, chrom, r2s, ref[chrom][r2s:r2s + READ], f"{READ}M", True, False, chrom, s, -INSERT)
            reads += [r1, r2]


def covering(reads, chrom, pos1):
    return [r for r in reads if r.reference_id == tid[chrom] and r.reference_start <= pos1 - 1 < r.reference_end
            and "D" not in r.cigarstring and "S" not in r.cigarstring]


def set_base(r, pos1, base):
    i = pos1 - 1 - r.reference_start
    q = r.query_qualities
    s = r.query_sequence
    r.query_sequence = s[:i] + base + s[i + 1:]
    r.query_qualities = q


def plant_snv(reads, pos1, alt, n, strands):
    pool = [r for r in covering(reads, "demo1", pos1) if strands == "both" or r.is_reverse]
    for r in RNG.sample(pool, n):
        set_base(r, pos1, alt)
    return len(covering(reads, "demo1", pos1))


def plant_deletion(sample, reads):
    _, anchor, L, _ = DEL
    hap = ref["demo1"][:anchor] + ref["demo1"][anchor + L:]  # alt haplotype (0-based anchor index = anchor-1)
    out = {"spanning": 0, "clipped": 0}
    plans = [("span", a) for a in RNG.sample(range(20, 80), 11)] + [("rclip", a) for a in (93, 95)] + \
            [("lclip", a) for a in (6, 7)]
    for k, (kind, before) in enumerate(plans):
        # `before` = read bases up to and including the anchor base
        hs = anchor - before  # 0-based start on the alt haplotype
        seq = hap[hs:hs + READ]
        after = READ - before
        if kind == "span":
            start, cigar = hs, f"{before}M{L}D{after}M"; out["spanning"] += 1
        elif kind == "rclip":
            start, cigar = hs, f"{before}M{after}S"; out["clipped"] += 1
        else:
            start, cigar = anchor + L, f"{before}S{after}M"; out["clipped"] += 1
        q = f"{sample}_del{k}"
        m = start + INSERT - READ + 40
        reads.append(seg(q, "demo1", start, seq, cigar, k % 2 == 1, k % 2 == 0, "demo1", m, INSERT))
        reads.append(seg(q, "demo1", m, ref["demo1"][m:m + READ], f"{READ}M", k % 2 == 0, k % 2 == 1, "demo1",
                         start, -INSERT))
    return out


def plant_bnd(sample, reads):
    (c1, p1), (c2, p2) = BND[1], BND[2]
    o1, o2 = p1 - 800, p2 - 1  # hap = c1[p1-800 : p1] + c2[p2-1 : p2-1+800]; junction at hap index 800
    hap = ref[c1][o1:p1] + ref[c2][o2:o2 + 800]
    split, disc = 0, 0
    for k, s in enumerate(RNG.sample(range(470, 501), 6)):  # read2 crosses the junction, >=20 bp each side
        q = f"{sample}_split{k}"
        left = 800 - (s + INSERT - READ); right = READ - left
        r1 = seg(q, c1, o1 + s, hap[s:s + READ], f"{READ}M", False, True, c1, o1 + s + INSERT - READ, INSERT)
        seq2 = hap[s + INSERT - READ:s + INSERT]
        sa_p = f"{c2},{p2},-,{left}S{right}M,60,0;"
        sa_s = f"{c1},{o1 + s + INSERT - READ + 1},-,{left}M{right}S,60,0;"
        r2 = seg(q, c1, o1 + s + INSERT - READ, seq2, f"{left}M{right}S", True, False, c1, o1 + s, -INSERT,
                 tags=[("SA", sa_p)])
        r2s = seg(q, c2, o2, seq2, f"{left}S{right}M", True, False, c1, o1 + s, 0, extra_flags=2048,
                  tags=[("SA", sa_s)])
        reads += [r1, r2, r2s]; split += 1
    for k, s in enumerate(RNG.sample(range(550, 701), 8)):  # read1 on c1, read2 fully on c2
        q = f"{sample}_disc{k}"
        m = o2 + (s + INSERT - READ - 800)
        r1 = seg(q, c1, o1 + s, hap[s:s + READ], f"{READ}M", False, True, c2, m, 0)
        r2 = seg(q, c2, m, hap[s + INSERT - READ:s + INSERT], f"{READ}M", True, False, c1, o1 + s, 0)
        reads += [r1, r2]; disc += 1
    return {"split": split, "discordant": disc, "union": split + disc}


def plant_cnv(sample, reads):
    """Even coverage on demo3, with a homozygous and a one-copy deletion in the tumor."""
    tumor = sample == "demo_tumor"
    n = 0
    for _ in range(CNV_DEPTH_FRAGS):
        s = RNG_CNV.randint(0, CNV_LEN - INSERT)
        a, b = s + 1, s + INSERT  # 1-based fragment span
        if tumor and not (b < HOMDEL[0] or a > HOMDEL[1]):
            continue  # no DNA left in the homozygous deletion
        if tumor and HEMI[0] <= a <= HEMI[1] and RNG_CNV.random() < 0.5:
            continue  # one of two copies lost
        q = f"{sample}_cnv{n}"; n += 1
        r2s = s + INSERT - READ
        reads.append(seg(q, CNV_CONTIG, s, ref[CNV_CONTIG][s:s + READ], f"{READ}M", False, True, CNV_CONTIG, r2s, INSERT))
        reads.append(seg(q, CNV_CONTIG, r2s, ref[CNV_CONTIG][r2s:r2s + READ], f"{READ}M", True, False, CNV_CONTIG, s,
                         -INSERT))
    if tumor:  # the alternate-contig region: reads are there but ambiguous (MAPQ 0, XA to an _alt contig)
        for r in reads:
            if r.reference_id == tid[CNV_CONTIG] and ALTREGION[0] - 1 <= r.reference_start < ALTREGION[1]:
                r.mapping_quality = 0
                r.set_tag("XA", f"{CNV_CONTIG}_alt1,+{r.reference_start - ALTREGION[0] + 2},{READ}M,0;")
    if tumor:  # a few low-MAPQ reads inside the homozygous deletion, as mouse or mismapped reads look
        for k in range(20):
            s = RNG_CNV.randint(HOMDEL[0], HOMDEL[1] - READ)
            r = seg(f"{sample}_lowmq{k}", CNV_CONTIG, s, ref[CNV_CONTIG][s:s + READ], f"{READ}M", k % 2 == 1,
                    True, CNV_CONTIG, s, 0)
            r.mapping_quality = 3
            reads.append(r)


def write_bam(path, reads, sample):
    header = {"HD": {"VN": "1.6", "SO": "unsorted"},
              "SQ": [{"SN": c, "LN": n} for c, n in CONTIGS.items()],
              "RG": [{"ID": sample, "SM": sample}], "CO": [COMMENT]}
    tmp = path.with_suffix(".unsorted.bam")
    with pysam.AlignmentFile(str(tmp), "wb", header=header) as out:
        for r in reads:
            r.set_tag("RG", sample)
            out.write(r)
    pysam.sort("-o", str(path), str(tmp))
    tmp.unlink()
    pysam.index(str(path))


def main():
    (HERE / "demo_ref.fa").write_text("".join(f">{c} {COMMENT}\n" + "\n".join(s[i:i + 60] for i in range(0, len(s), 60))
                                              + "\n" for c, s in ref.items()))
    pysam.faidx(str(HERE / "demo_ref.fa"))
    truth = {}
    for sample in ("demo_tumor", "demo_normal"):
        reads = []
        background(sample, reads)
        tumor = sample == "demo_tumor"
        for vid, (pos, alt, n, strands, _) in SNVS.items():
            depth = plant_snv(reads, pos, alt, n if tumor else 0, strands)
            truth.setdefault(vid, {})[sample] = {"alt": n if tumor else 0, "depth": depth}
        cov = covering(reads, "demo1", GERM[1])
        k = len(cov) // 2
        for r in RNG.sample(cov, k):
            set_base(r, GERM[1], GERM_ALT)
        truth.setdefault("V6", {})[sample] = {"alt": k, "depth": len(cov)}
        if tumor:
            d = plant_deletion(sample, reads)
            truth["V4"] = {sample: {"alt": d["spanning"] + d["clipped"], "alt_without_clip_rescue": d["spanning"]},
                           "demo_normal": {"alt": 0}}
            truth["V5"] = {sample: plant_bnd(sample, reads), "demo_normal": {"union": 0}}
        plant_cnv(sample, reads)
        write_bam(HERE / f"{sample}.bam", reads, sample)
    cols = ["sample", "contig", "start", "end", "num_cr", "log2", "call", "maf", "cn_call", "loh", "log2_capped"]
    cn = {"-": "del", "+": "amp", "0": "neutral"}
    (HERE / "demo_cnv.tsv").write_text("\t".join(cols) + "\n" + "".join(
        "\t".join(map(str, ["demo_tumor", CNV_CONTIG, a, b, (b - a + 1) // 100, l2, c, "NaN", cn[c], 0, max(l2, -10)]))
        + "\n" for a, b, l2, c in CNV_SEGMENTS))
    (HERE / "demo_cnv_genes.bed").write_text(f"# {COMMENT}\n" + "".join(
        f"{CNV_CONTIG}\t{a - 1}\t{b}\t{g}\n" for g, a, b in CNV_GENES))
    truth["CNV"] = {"HOMDEL": {"ratio_max": 0.1, "low_mapq": True}, "HEMI": {"ratio_min": 0.35, "ratio_max": 0.65},
                    "NEUTRAL": {"ratio_min": 0.8, "ratio_max": 1.2}, "WRONGDEL": {"disagrees": True},
                    "GAPGENE": {"no_segment": True}, "ALTGENE": {"ambiguous": True}}

    anchor, L = DEL[1], DEL[2]
    (c1, p1), (c2, p2) = BND[1], BND[2]
    b1, b2 = ref[c1][p1 - 1], ref[c2][p2 - 1]
    rows = [(c, p, v, ref["demo1"][p - 1], a, "weak_evidence" if v == "V3" else "PASS", f"GENE={g}")
            for v, (p, a, _, _, g) in SNVS.items() for c in ["demo1"]]
    rows.append(("demo1", anchor, "V4", ref["demo1"][anchor - 1:anchor + L], ref["demo1"][anchor - 1], "PASS",
                 f"GENE={DEL[3]}"))
    rows.append((c1, p1, "V5_1", b1, f"{b1}[{c2}:{p2}[", "PASS", f"SVTYPE=BND;MATEID=V5_2;GENE={BND[3]}"))
    rows.append((c2, p2, "V5_2", b2, f"]{c1}:{p1}]{b2}", "PASS", f"SVTYPE=BND;MATEID=V5_1;GENE={BND[3]}"))
    rows.append(("demo1", GERM[1], "V6", ref["demo1"][GERM[1] - 1], GERM_ALT, "PASS", f"GENE={GERM[2]}"))
    rows.sort(key=lambda r: (list(CONTIGS).index(r[0]), r[1]))
    lines = ["##fileformat=VCFv4.2", f"##comment={COMMENT}"]
    lines += [f"##contig=<ID={c},length={n}>" for c, n in CONTIGS.items()]
    lines += ['##FILTER=<ID=PASS,Description="All filters passed">',
              '##FILTER=<ID=weak_evidence,Description="Synthetic low-evidence call">',
              '##INFO=<ID=GENE,Number=1,Type=String,Description="Synthetic gene symbol">',
              '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
              '##INFO=<ID=MATEID,Number=1,Type=String,Description="ID of mate breakend">',
              '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
              '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Caller-reported allelic depths (synthetic)">',
              "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tdemo_normal\tdemo_tumor"]
    def ad(vid, sample):
        """Caller-style AD (ref,alt) from the truth; V3's tumor count is deliberately inflated."""
        t = truth.get(vid, {}).get(sample)
        if not t or "depth" not in t:
            return "."
        alt = 12 if (vid, sample) == ("V3", "demo_tumor") else t["alt"]
        return f"{t['depth'] - t['alt']},{alt}"
    truth["V4"]["demo_tumor"]["depth"] = 60  # caller-style depth for the deletion (15 alt of 60)
    truth["V4"]["demo_normal"]["depth"] = 63
    lines += [f"{c}\t{p}\t{i}\t{r}\t{a}\t.\t{f}\t{inf}\tGT:AD\t0/0:{ad(i, 'demo_normal')}\t0/1:{ad(i, 'demo_tumor')}"
              for c, p, i, r, a, f, inf in rows]
    (HERE / "demo_calls.vcf").write_text("\n".join(lines) + "\n")
    (HERE / "demo_variants.tsv").write_text(
        f"# {COMMENT}\nchrom\tpos\tref\talt\ndemo1\t{SNVS['V1'][0]}\t{ref['demo1'][SNVS['V1'][0] - 1]}\t{SNVS['V1'][1]}\n"
        f"demo1\t{anchor}\t{ref['demo1'][anchor - 1:anchor + L]}\t{ref['demo1'][anchor - 1]}\n")
    (HERE / "demo_truth.json").write_text(json.dumps({"comment": COMMENT, "variants": truth}, indent=2) + "\n")
    print("wrote", ", ".join(sorted(p.name for p in HERE.iterdir() if p.name != Path(__file__).name)))


if __name__ == "__main__":
    main()
