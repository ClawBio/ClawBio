#!/usr/bin/env python3
"""Build the human (GRCh38) reference files BJ-Expression 2.0.0 expects, from
public sources, so a run does not depend on BioSkryb's S3 bucket.

Inputs (all public): the NCBI GRCh38 no-alt analysis-set genome, the Ensembl 104
GTF, and the Ensembl 104 cDNA and ncRNA FASTAs.

Formats were read from the pipeline's own scripts, not guessed:
  gtf_file      contigs must match the BAMs, so they are renamed chr-style
  tx2gene       5 columns, no header, whitespace separated:
                transcript_id, transcript_biotype, gene_id, gene_biotype, gene_symbol
                (transcript_id must equal the Salmon FASTA ID exactly)
  genebody_ref  12-column BED (geneBodyCoverageIsai.py skips anything else)
  bed_files_dir <biotype>_transcripts.bed files (insert-size step, GRCh38 only)

The result mirrors the layout of `genomes_base`, so pass the output root as
`--genomes-base`. STAR and Salmon indexes are built separately (they need a lot
of RAM) with the commands this module produces.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
import shutil
import sys
from collections import Counter
from pathlib import Path

PRIMARY = (*map(str, range(1, 23)), "X", "Y", "MT")
ATTR = re.compile(r'(\w+) "([^"]*)"')
SAFE_BIOTYPE = re.compile(r"[A-Za-z0-9_]+")

# Layout under genomes_base, as in BioSkryb's genomes.config for GRCh38.
SEQ = Path("genomes/Homo_sapiens/NCBI/GRCh38/Sequence")
STAR_DIR = SEQ / "star_chr_eqnames_wgs"
SALMON_DIR = SEQ / "salmon_chr_eqnames_wgs"
GTF_NAME = "Homo_sapiens.GRCh38.104.chr.gtf"
TX2GENE_PATH = SEQ / "salmon" / "Ensemble_GRCh38_cDNA-ncRNA.complete.tx2gene"
GENOME_MD5_NAME = "GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"


class BuilderError(Exception):
    """A user-facing failure."""


# ──────────────────────────────────────────────
# GTF
# ──────────────────────────────────────────────

def ucsc_name(contig: str) -> str | None:
    """Ensembl primary contig -> UCSC/NCBI-style name; None if not primary."""
    if contig == "MT":
        return "chrM"
    return f"chr{contig}" if contig in PRIMARY else None


def parse_attributes(field: str) -> dict[str, str]:
    return dict(ATTR.findall(field))


def _open_text(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def rename_gtf(lines, fasta_contigs: set[str] | None = None):
    """Rename contigs chr-style; drop non-primary contigs and (optionally)
    contigs absent from the genome FASTA. Returns (lines, stats)."""
    out: list[str] = []
    dropped: Counter = Counter()
    not_in_fasta: Counter = Counter()
    kept = 0
    for raw in lines:
        ln = raw.rstrip("\n")
        if ln.startswith("#"):
            out.append(ln)
            continue
        fields = ln.split("\t")
        new = ucsc_name(fields[0])
        if new is None:
            dropped[fields[0]] += 1
            continue
        if fasta_contigs is not None and new not in fasta_contigs:
            not_in_fasta[new] += 1
            continue
        fields[0] = new
        out.append("\t".join(fields))
        kept += 1
    return out, {"kept": kept, "dropped_by_contig": dict(dropped),
                 "not_in_fasta": dict(not_in_fasta)}


def transcripts(gtf_lines, with_versions: bool = True) -> list[dict]:
    """One record per transcript that has exons, in a canonical order."""
    tx: dict[str, dict] = {}
    exons: dict[str, list] = {}
    for raw in gtf_lines:
        if raw.startswith("#"):
            continue
        f = raw.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        feature = f[2]
        if feature not in ("transcript", "exon"):
            continue
        a = parse_attributes(f[8])
        key = a.get("transcript_id")
        if not key:
            continue
        if feature == "exon":
            exons.setdefault(key, []).append((int(f[3]), int(f[4])))
            continue
        tver, gver = a.get("transcript_version"), a.get("gene_version")
        gene = a.get("gene_id", "")
        gene_full = f"{gene}.{gver}" if with_versions and gver else gene
        tx[key] = {
            "chrom": f[0], "strand": f[6],
            "tx_id": f"{key}.{tver}" if with_versions and tver else key,
            "gene_id": gene_full,
            "symbol": a.get("gene_name") or gene_full,
            "gene_biotype": a.get("gene_biotype") or "unknown",
            "tx_biotype": a.get("transcript_biotype") or "unknown",
        }
    records = []
    for key, rec in tx.items():
        ex = sorted(exons.get(key, []))
        if not ex:
            continue
        records.append({**rec, "exons": ex, "start": ex[0][0], "end": ex[-1][1]})
    return sorted(records, key=lambda r: (r["chrom"], r["start"], r["end"], r["tx_id"]))


def _clean(value: str) -> str:
    """read.table splits on whitespace and treats # and quotes specially."""
    return re.sub(r"[\s#\"']+", "_", value) or "NA"


def tx2gene_lines(records: list[dict]) -> list[str]:
    rows = sorted({(r["tx_id"], r["tx_biotype"], r["gene_id"], r["gene_biotype"], r["symbol"])
                   for r in records})
    return ["\t".join(_clean(c) for c in row) for row in rows]


def bed12(rec: dict) -> str:
    start0 = rec["start"] - 1
    sizes = [e - s + 1 for s, e in rec["exons"]]
    offsets = [s - 1 - start0 for s, _ in rec["exons"]]
    return "\t".join(map(str, [
        rec["chrom"], start0, rec["end"], rec["tx_id"], 0, rec["strand"],
        start0, rec["end"], 0, len(sizes),
        ",".join(map(str, sizes)), ",".join(map(str, offsets))]))


def write_biotype_beds(records: list[dict], outdir: Path) -> dict[str, int]:
    by_biotype: dict[str, list[str]] = {}
    for r in records:
        if not SAFE_BIOTYPE.fullmatch(r["tx_biotype"]):
            raise BuilderError(f"unsafe transcript biotype name {r['tx_biotype']!r}; "
                               "it becomes a file name")
        by_biotype.setdefault(r["tx_biotype"], []).append(bed12(r))
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    for biotype, lines in by_biotype.items():
        (outdir / f"{biotype}_transcripts.bed").write_text("\n".join(lines) + "\n")
    return {b: len(v) for b, v in by_biotype.items()}


# ──────────────────────────────────────────────
# FASTA
# ──────────────────────────────────────────────

def fasta_ids(path: Path):
    with _open_text(Path(path)) as fh:
        for ln in fh:
            if ln.startswith(">"):
                yield ln[1:].split()[0] if ln[1:].strip() else ""


def decoy_names(path: Path) -> list[str]:
    return list(fasta_ids(path))


def filter_fasta(paths: list[Path], keep_ids: set[str], out: Path) -> dict:
    """Write records whose ID is in keep_ids, each ID once, in input order."""
    seen: set[str] = set()
    stats = {"kept": 0, "dropped": 0, "duplicates": 0}
    with open(out, "w") as w:
        for p in paths:
            write = False
            with _open_text(Path(p)) as fh:
                for ln in fh:
                    if ln.startswith(">"):
                        rid = ln[1:].split()[0] if ln[1:].strip() else ""
                        if rid in seen:
                            stats["duplicates"] += 1
                            write = False
                        elif rid in keep_ids:
                            seen.add(rid)
                            stats["kept"] += 1
                            write = True
                        else:
                            stats["dropped"] += 1
                            write = False
                    if write:
                        w.write(ln if ln.endswith("\n") else ln + "\n")
    return stats


def missing_from_tx2gene(ids, known: set[str]) -> list[str]:
    return sorted(i for i in ids if i not in known)


# ──────────────────────────────────────────────
# Index commands and checks
# ──────────────────────────────────────────────

def star_index_command(genome_fa: Path, gtf: Path, outdir: Path, threads: int,
                       sjdb_overhang: int = 100, log_prefix: Path | None = None) -> list[str]:
    cmd = ["STAR", "--runMode", "genomeGenerate", "--runThreadN", str(threads),
           "--genomeDir", str(outdir), "--genomeFastaFiles", str(genome_fa),
           "--sjdbGTFfile", str(gtf), "--sjdbOverhang", str(sjdb_overhang)]
    if log_prefix is not None:
        cmd += ["--outFileNamePrefix", str(log_prefix)]
    return cmd


def salmon_index_command(gentrome: Path, decoys: Path, outdir: Path, threads: int) -> list[str]:
    return ["salmon", "index", "-t", str(gentrome), "-d", str(decoys),
            "-i", str(outdir), "-p", str(threads)]


def md5_matches(path: Path, expected: str) -> bool:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest() == expected.lower()


def _expected_md5(md5_file: Path, name: str) -> str:
    for ln in Path(md5_file).read_text().splitlines():
        parts = ln.split()
        if len(parts) == 2 and Path(parts[1]).name == name:
            return parts[0]
    raise BuilderError(f"{name} not listed in {md5_file}")


# ──────────────────────────────────────────────
# Housekeeping-gene BED (genebody_ref)
# ──────────────────────────────────────────────

HOUSEKEEPING_BED = "hg38.HouseKeepingGenes.bed"  # name the pipeline's config uses


def parse_gene_list(text: str) -> list[tuple[str, str]]:
    """Two columns, Ensembl gene ID then symbol. Versions are stripped: the list
    comes from a different annotation release than the GTF, so versions differ."""
    genes: list[tuple[str, str]] = []
    seen: set[str] = set()
    for n, raw in enumerate(text.splitlines(), start=1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) != 2 or not re.fullmatch(r"ENSG\d+(\.\d+)?", parts[0]):
            raise BuilderError(f"gene list line {n}: expected 'ENSG... SYMBOL', got {raw!r}")
        gid = parts[0].split(".")[0]
        if gid not in seen:
            seen.add(gid)
            genes.append((gid, parts[1]))
    return genes


def housekeeping_bed(records: list[dict], genes: list[tuple[str, str]]):
    """12-column BED of every transcript of the listed genes, matched on the
    unversioned gene ID. Symbol disagreements and misses are reported."""
    by_gene: dict[str, list[dict]] = {}
    by_symbol: dict[str, set[str]] = {}
    for r in records:
        gid = r["gene_id"].split(".")[0]
        by_gene.setdefault(gid, []).append(r)
        by_symbol.setdefault(r["symbol"], set()).add(gid)
    wanted = {gid for gid, _ in genes if gid in by_gene}
    unmatched, mismatch, symbol_only = [], [], []
    for gid, sym in genes:
        if gid in by_gene:
            gtf_sym = by_gene[gid][0]["symbol"]
            if gtf_sym != sym:
                mismatch.append([gid, sym, gtf_sym])
        else:
            unmatched.append(gid)
            hit = sorted(by_symbol.get(sym, ()))
            if hit:
                symbol_only.append([gid, sym, hit[0]])
    lines = [bed12(r) for r in records if r["gene_id"].split(".")[0] in wanted]
    if not lines:
        raise BuilderError("no housekeeping genes matched the annotation by gene ID; "
                           "an empty BED would silently produce no gene-body coverage")
    return lines, {
        "requested": len(genes), "matched_genes": sorted(wanted), "transcripts": len(lines),
        "per_gene_transcripts": {g: len(by_gene[g]) for g in sorted(wanted)},
        "unmatched": unmatched, "symbol_mismatch": mismatch, "symbol_only": symbol_only}


def write_housekeeping(out: Path, genes_file: Path) -> dict:
    out = Path(out)
    prep = out / "prepare_report.json"
    gtf = out / STAR_DIR / GTF_NAME
    if not prep.is_file() or not gtf.is_file():
        raise BuilderError(f"{out} has no prepared references; run 'prepare' first")
    versions = json.loads(prep.read_text())["transcript_ids_have_versions"]
    genes = parse_gene_list(Path(genes_file).read_text())
    with open(gtf) as fh:
        recs = transcripts(fh, with_versions=versions)
    lines, report = housekeeping_bed(recs, genes)
    bed = out / STAR_DIR / HOUSEKEEPING_BED
    bed.write_text("\n".join(lines) + "\n")
    report["genes_file_sha256"] = hashlib.sha256(Path(genes_file).read_bytes()).hexdigest()
    report["bed"] = str(bed)
    (out / "housekeeping_report.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


# ──────────────────────────────────────────────
# Cell-typing reference (celltype_ref)
# ──────────────────────────────────────────────

CELLTYPE_DIR = Path("dev-resources/rnaseq_celltype_qbio")
CELLTYPE_NAME = "db_celltype_qbio.RDS"  # the script runs readRDS() on exactly this name
# scrnaseq_celltype_refactored_pipeline.R indexes the list by these fixed names;
# a list without them makes that step fail (the pipeline then ignores the failure).
CELLTYPE_EXPECTED = ("HPCA_singleR", "GTEx_singleR", "GTExmeds_singleR",
                     "TCGA_STN_singleR", "TCGA_PT_singleR")
CELLTYPE_IMAGE = "public.ecr.aws/v2u5x5f1/basejumper:rnaseq_celltype_0.4"


def missing_celltype_entries(names) -> list[str]:
    have = set(names)
    return [n for n in CELLTYPE_EXPECTED if n not in have]


def inspect_celltype_names(apptainer: str, sandbox: Path, rds: Path,
                           bind: Path | None = None) -> list[str] | None:
    """Entry names of the list in the RDS, read inside the cell-typing image
    (it has R). None if it could not be read."""
    import subprocess
    bind = Path(bind) if bind else Path.home()
    r = ('x <- readRDS(commandArgs(TRUE)[1]); '
         'cat(if (is.list(x)) names(x) else character(0), sep = "\\n")')
    try:
        proc = subprocess.run(
            [apptainer, "exec", "--writable", "--bind", f"{bind}:{bind}", str(sandbox),
             "Rscript", "-e", r, str(rds)], capture_output=True, text=True, timeout=600)
    except (OSError, subprocess.TimeoutExpired):
        return None
    if proc.returncode != 0:
        return None
    return [ln for ln in proc.stdout.splitlines() if ln.strip()]


def place_celltype_ref(out: Path, rds: Path, force: bool = False,
                       names: list[str] | None = None) -> dict:
    """Copy an RDS to where the pipeline will look for it, under the name its
    script reads. A copy (not a link) so containers that cannot see the source
    path, such as /orgbucket, still find it."""
    rds, out = Path(rds), Path(out)
    if not rds.is_file():
        raise BuilderError(f"cell-typing RDS not found: {rds}")
    dest = out / CELLTYPE_DIR / CELLTYPE_NAME
    digest = hashlib.sha256(rds.read_bytes()).hexdigest()
    if dest.exists():
        if hashlib.sha256(dest.read_bytes()).hexdigest() != digest and not force:
            raise BuilderError(f"{dest} already holds a different file; pass --force to replace it")
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(rds, dest)
    missing = missing_celltype_entries(names) if names is not None else None
    return {
        "source": str(rds), "placed": str(dest), "sha256": digest,
        "inspected": names is not None,
        "entries_found": list(names) if names is not None else None,
        "expected_entries_missing": missing,
        "compatible_with_pipeline_script": (not missing) if names is not None else None,
    }


# ──────────────────────────────────────────────
# build-indexes: STAR and Salmon inside the sandboxes, with a RAM guard
# ──────────────────────────────────────────────

# The pipeline aligns with STAR 2.7.6a and quantifies with Salmon 1.6.0, so the
# indexes are built with exactly those versions (the sandboxes from
# apptainer_sandboxes.py).
STAR_IMAGE = "quay.io/biocontainers/star:2.7.6a--0"
SALMON_IMAGE = "quay.io/biocontainers/salmon:1.6.0--h84f40af_0"
# Free-memory floors, in GiB. CONSERVATIVE ESTIMATES, not measurements: a human
# STAR index is commonly quoted at about 30 GB, Salmon with decoys at well under
# that. An OOM kill restarts the container and loses background work, so we
# refuse to start without headroom. Override with --min-free-gb-star/-salmon.
DEFAULT_MIN_FREE_GB = {"star": 40.0, "salmon": 16.0}
STAR_FILES = ("Genome", "SA", "SAindex", "genomeParameters.txt", "chrName.txt")


def available_memory_gb(meminfo_text: str) -> float:
    m = re.search(r"^MemAvailable:\s+(\d+)\s*kB", meminfo_text, re.M)
    if not m:
        raise BuilderError("MemAvailable not found in the memory info")
    return int(m.group(1)) / 1024 / 1024


def _cache_name(image: str) -> str:
    """The sandbox directory name; one rule, owned by apptainer_sandboxes.py."""
    import importlib.util
    path = Path(__file__).resolve().parent / "apptainer_sandboxes.py"
    spec = importlib.util.spec_from_file_location("_bj_apptainer_sandboxes", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.cache_name(image)


def _star_complete(d: Path) -> list[str]:
    return [n for n in STAR_FILES if not (d / n).is_file()]


def _salmon_complete(d: Path) -> list[str]:
    ok = (d / "info.json").is_file() or (d / "versionInfo.json").is_file()
    return [] if ok else ["info.json or versionInfo.json"]


def build_indexes(out: Path, sandboxes: Path, threads: int, only: str | None = None,
                  sjdb_overhang: int = 100, min_free_gb: dict | None = None,
                  bind: Path | None = None, meminfo_text: str | None = None) -> dict:
    import subprocess
    out, sandboxes = Path(out), Path(sandboxes)
    work = out / "_work"
    need = [out / "prepare_report.json", out / STAR_DIR / GTF_NAME,
            work / "genome.fa", work / "gentrome.fa", work / "decoys.txt"]
    if not all(p.is_file() for p in need):
        raise BuilderError(f"{out} has no prepared references; run 'prepare' first")
    which = ["star", "salmon"] if only is None else [only]
    apptainer = shutil.which("apptainer") or shutil.which("singularity")
    if not apptainer:
        raise BuilderError("apptainer not found on PATH")

    targets = {"star": (out / STAR_DIR / "index_ensemble", _star_complete, STAR_IMAGE),
               "salmon": (out / SALMON_DIR / "index", _salmon_complete, SALMON_IMAGE)}
    report: dict[str, str] = {}
    todo = []
    for k in which:
        outdir, complete, image = targets[k]
        box = sandboxes / _cache_name(image)
        if not box.is_dir():
            raise BuilderError(f"sandbox for {image} not found at {box}; run "
                               "apptainer_sandboxes.py first")
        if outdir.is_dir() and any(outdir.iterdir()):
            missing = complete(outdir)
            if missing:
                raise BuilderError(f"{outdir} exists but is incomplete (missing "
                                   f"{', '.join(missing)}); delete it and re-run")
            report[k] = "skipped"
        else:
            todo.append(k)

    floors = {**DEFAULT_MIN_FREE_GB, **(min_free_gb or {})}
    text = meminfo_text if meminfo_text is not None else Path("/proc/meminfo").read_text()
    avail = available_memory_gb(text)
    for k in todo:
        if avail < floors[k]:
            raise BuilderError(
                f"only {avail:.0f} GiB of memory is available; the {k} index needs about "
                f"{floors[k]:.0f} GiB (an estimate). Ask for more RAM (resize the container), "
                "then re-run; an out-of-memory kill would restart the container.")

    bind = Path(bind) if bind else Path.home()
    for k in todo:
        outdir, complete, image = targets[k]
        if k == "star":
            tool = star_index_command(work / "genome.fa", out / STAR_DIR / GTF_NAME, outdir,
                                      threads, sjdb_overhang, log_prefix=work / "star_")
        else:
            tool = salmon_index_command(work / "gentrome.fa", work / "decoys.txt", outdir, threads)
        outdir.parent.mkdir(parents=True, exist_ok=True)
        cmd = [apptainer, "exec", "--writable", "--bind", f"{bind}:{bind}",
               str(sandboxes / _cache_name(image)), *tool]
        proc = subprocess.run(cmd, cwd=work, capture_output=True, text=True)
        if proc.returncode != 0:
            raise BuilderError(f"{k} index build failed (exit {proc.returncode}): "
                               f"{proc.stderr.strip()[-400:]}")
        missing = complete(outdir)
        if missing:
            raise BuilderError(f"{k} index is incomplete after the build (missing "
                               f"{', '.join(missing)})")
        report[k] = "built"
    return report


# ──────────────────────────────────────────────
# prepare
# ──────────────────────────────────────────────

def prepare(downloads: Path, out: Path, versions: str = "auto") -> dict:
    d, out = Path(downloads), Path(out)
    need = ["genome.fna.gz", "ncbi_md5checksums.txt", "ensembl104.gtf.gz",
            "ensembl104.cdna.fa.gz", "ensembl104.ncrna.fa.gz"]
    missing = [n for n in need if not (d / n).is_file()]
    if missing:
        raise BuilderError(f"missing download(s) in {d}: {', '.join(missing)}")
    if not md5_matches(d / "genome.fna.gz", _expected_md5(d / "ncbi_md5checksums.txt", GENOME_MD5_NAME)):
        raise BuilderError("genome.fna.gz does not match NCBI's md5; re-download it")

    work = out / "_work"
    work.mkdir(parents=True, exist_ok=True)
    star_dir, salmon_dir = out / STAR_DIR, out / SALMON_DIR
    star_dir.mkdir(parents=True, exist_ok=True)
    (out / TX2GENE_PATH).parent.mkdir(parents=True, exist_ok=True)

    # genome: decompress once, collect contig names
    contigs: set[str] = set()
    with gzip.open(d / "genome.fna.gz", "rt") as src, open(work / "genome.fa", "w") as dst:
        for ln in src:
            if ln.startswith(">"):
                contigs.add(ln[1:].split()[0])
            dst.write(ln)

    # transcript ID style follows the FASTA the Salmon index will be built from
    first_id = next(fasta_ids(d / "ensembl104.cdna.fa.gz"), "")
    use_versions = (("." in first_id) if versions == "auto" else versions == "yes")

    with gzip.open(d / "ensembl104.gtf.gz", "rt") as fh:
        renamed, gtf_stats = rename_gtf(fh, contigs)
    (star_dir / GTF_NAME).write_text("\n".join(renamed) + "\n")

    recs = transcripts(renamed, with_versions=use_versions)
    tx_lines = tx2gene_lines(recs)
    (out / TX2GENE_PATH).write_text("\n".join(tx_lines) + "\n")
    biotypes = write_biotype_beds(recs, star_dir / "biotypes")

    tx_ids = {ln.split("\t", 1)[0] for ln in tx_lines}
    fa_stats = filter_fasta([d / "ensembl104.cdna.fa.gz", d / "ensembl104.ncrna.fa.gz"],
                            tx_ids, work / "transcripts.fa")
    (work / "decoys.txt").write_text("\n".join(decoy_names(work / "genome.fa")) + "\n")
    with open(work / "gentrome.fa", "wb") as g:  # transcripts first, then genome
        for part in ("transcripts.fa", "genome.fa"):
            with open(work / part, "rb") as p:
                shutil.copyfileobj(p, g, 1 << 24)
    kept_ids = set(fasta_ids(work / "transcripts.fa"))
    report = {
        "genome_md5_verified": True,
        "transcript_ids_have_versions": use_versions,
        "genome_contigs": len(contigs),
        "gtf": gtf_stats,
        "transcripts_with_exons": len(recs),
        "tx2gene_rows": len(tx_lines),
        "biotype_files": biotypes,
        "transcript_fasta": fa_stats,
        "gtf_transcripts_without_fasta_sequence": len(tx_ids - kept_ids),
        "layout": {"gtf": str(out / STAR_DIR / GTF_NAME), "tx2gene": str(out / TX2GENE_PATH),
                   "biotypes_dir": str(star_dir / "biotypes"),
                   "star_index_dir": str(star_dir / "index_ensemble"),
                   "salmon_index_dir": str(salmon_dir / "index"),
                   "work": str(work)},
    }
    (out / "prepare_report.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="bj-reference-builder", description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("prepare", help="Rename GTF, write tx2gene and biotype BEDs, filter transcripts, make decoys")
    p.add_argument("--downloads", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True, help="Root to pass as --genomes-base")
    p.add_argument("--versions", choices=("auto", "yes", "no"), default="auto",
                   help="Put .version on transcript/gene IDs (auto = follow the Ensembl FASTA)")
    c = sub.add_parser("celltype", help="Wire an RDS in as the pipeline's cell-typing reference")
    c.add_argument("--rds", type=Path, required=True)
    c.add_argument("--out", type=Path, required=True, help="The root 'prepare' wrote")
    c.add_argument("--sandboxes", type=Path, default=Path.home() / "bj-expression-sandboxes")
    c.add_argument("--force", action="store_true", help="Replace a different file already there")
    c.add_argument("--no-inspect", action="store_true",
                   help="Do not read the list's entry names in the container")
    b = sub.add_parser("build-indexes", help="Build the STAR and Salmon indexes in the sandboxes (needs lots of RAM)")
    b.add_argument("--out", type=Path, required=True, help="The root 'prepare' wrote")
    b.add_argument("--sandboxes", type=Path, default=Path.home() / "bj-expression-sandboxes")
    b.add_argument("--threads", type=int, default=8)
    b.add_argument("--only", choices=("star", "salmon"))
    b.add_argument("--sjdb-overhang", type=int, default=100)
    b.add_argument("--min-free-gb-star", type=float, default=DEFAULT_MIN_FREE_GB["star"])
    b.add_argument("--min-free-gb-salmon", type=float, default=DEFAULT_MIN_FREE_GB["salmon"])
    h = sub.add_parser("housekeeping", help="Write the 12-column housekeeping-gene BED (genebody_ref)")
    h.add_argument("--genes", type=Path, required=True,
                   help="Two columns: Ensembl gene ID (version optional), symbol")
    h.add_argument("--out", type=Path, required=True, help="The root 'prepare' wrote")
    args = ap.parse_args(argv)
    if args.cmd == "celltype":
        try:
            rep = place_celltype_ref(args.out, args.rds, args.force)
            if not args.no_inspect:
                exe = shutil.which("apptainer") or shutil.which("singularity")
                box = args.sandboxes / _cache_name(CELLTYPE_IMAGE)
                if exe and box.is_dir():
                    # inspect the placed copy: the source may be on a mount the container cannot see
                    names = inspect_celltype_names(exe, box, Path(rep["placed"]))
                    if names is not None:
                        rep = place_celltype_ref(args.out, args.rds, args.force, names)
        except BuilderError as exc:
            print(f"bj-reference-builder: error: {exc}", file=sys.stderr)
            return 1
        (args.out / "celltype_report.json").write_text(json.dumps(rep, indent=2) + "\n")
        print(f"Placed: {rep['placed']}")
        if not rep["inspected"]:
            print("Entry names were not inspected, so compatibility with the pipeline's "
                  "cell-typing script is unknown.")
        elif rep["compatible_with_pipeline_script"]:
            print("All five entries the pipeline's script needs are present.")
        else:
            print("WARNING: the pipeline's cell-typing script needs these entries, which are "
                  f"missing: {', '.join(rep['expected_entries_missing'])}. That step will fail "
                  "and the pipeline will ignore the failure after two retries.")
        return 0
    if args.cmd == "build-indexes":
        try:
            rep = build_indexes(args.out, args.sandboxes, args.threads, args.only,
                                args.sjdb_overhang,
                                {"star": args.min_free_gb_star, "salmon": args.min_free_gb_salmon})
        except BuilderError as exc:
            print(f"bj-reference-builder: error: {exc}", file=sys.stderr)
            return 1
        print(json.dumps(rep, indent=2))
        return 0
    if args.cmd == "housekeeping":
        try:
            rep = write_housekeeping(args.out, args.genes)
        except BuilderError as exc:
            print(f"bj-reference-builder: error: {exc}", file=sys.stderr)
            return 1
        print(json.dumps({k: rep[k] for k in ("requested", "transcripts", "unmatched",
                                              "symbol_mismatch", "symbol_only")}, indent=2))
        return 0
    try:
        report = prepare(args.downloads, args.out, args.versions)
    except BuilderError as exc:
        print(f"bj-reference-builder: error: {exc}", file=sys.stderr)
        return 1
    print(json.dumps({k: report[k] for k in ("transcripts_with_exons", "tx2gene_rows",
                                             "transcript_fasta", "transcript_ids_have_versions")}, indent=2))
    print(f"Report: {args.out / 'prepare_report.json'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
