#!/usr/bin/env python3
# Vendored from UKDRI/informatics_data_skills @ 7cc3e6e (geo/scripts/geo.py).
# Copyright (c) 2026 UK Dementia Research Institute. Licensed MIT.
# Kept close to upstream so changes there remain easy to merge; the ClawBio
# CLI contract lives in geo_fetch.py, not here.
"""Query metadata and download data from NCBI Gene Expression Omnibus (GEO).

Metadata is retrieved through NCBI E-utilities (db=gds); data files are pulled
from the GEO FTP tree at https://ftp.ncbi.nlm.nih.gov/geo/ .

Standard library only (urllib) -- no pip install required.

Examples:
    python geo_fetch.py metadata GSE2553
    python geo_fetch.py samples GSE2553
    python geo_fetch.py files GSE2553
    python geo_fetch.py download GSE2553 --matrix --out ./geo_out
    python geo_fetch.py download GSE2553 --suppl --out ./geo_out
    python geo_fetch.py search "breast cancer RNA-seq" --organism "Homo sapiens" --limit 20
    python geo_fetch.py metadata-table GSE2553 --out metadata.tsv
    python geo_fetch.py runtable GSE110009 --out .          # SraRunTable.csv + SRR_Acc_List.txt
    python geo_fetch.py samplesheet GSE110009 --assay bulk --from-runtable SraRunTable.csv --fastq-dir ./fastq
"""
import argparse
import json
import os
import sys
import time
import urllib.parse
import urllib.request

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
FTP_BASE = "https://ftp.ncbi.nlm.nih.gov/geo"
USER_AGENT = "geo-skill/1.0 (https://www.ncbi.nlm.nih.gov/geo/)"

# NCBI accepts an optional email + api_key to lift E-utilities rate limits.
# These are read from the environment but, per the "No unsolicited credentials"
# rule, are NEVER sent unless the user explicitly opts in with
# --use-ncbi-credentials — the mere presence of the env vars is not consent.
NCBI_EMAIL = os.environ.get("NCBI_EMAIL", "")
NCBI_API_KEY = os.environ.get("NCBI_API_KEY", "")
_USE_CREDENTIALS = False  # flipped on only by --use-ncbi-credentials


def _eutil_params(extra):
    p = dict(extra)
    if not _USE_CREDENTIALS:
        return p  # opt-in required; ignore any env-var credentials
    if NCBI_EMAIL:
        p["email"] = NCBI_EMAIL
        p["tool"] = "geo-skill"
    if NCBI_API_KEY:
        p["api_key"] = NCBI_API_KEY
    return p


def http_get(url, retries=3):
    last = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=120) as r:
                return r.read()
        except Exception as e:  # noqa: BLE001
            last = e
            time.sleep(1.5 * (attempt + 1))
    raise SystemExit(f"GET failed for {url}: {last}")


def get_json(base, params):
    url = base + "?" + urllib.parse.urlencode(_eutil_params(params))
    return json.loads(http_get(url).decode("utf-8", "replace"))


def geo_ftp_dir(accession):
    """Return the FTP subdirectory + type dir for a GEO accession.

    GSE2553 -> ('series', 'GSE2nnn'), GSM12345 -> ('samples', 'GSM12nnn').
    """
    accession = accession.strip().upper()
    prefix = accession[:3]
    num = accession[3:]
    type_dir = {"GSE": "series", "GSM": "samples", "GPL": "platforms"}.get(prefix)
    if not type_dir:
        raise SystemExit(f"Unsupported accession for FTP download: {accession}")
    stub = f"{prefix}{num[:-3]}nnn" if len(num) > 3 else f"{prefix}nnn"
    return type_dir, stub


def resolve_uid(accession):
    """Resolve a GEO accession (GSE/GSM/GPL/GDS) to a db=gds UID."""
    res = get_json(
        f"{EUTILS}/esearch.fcgi",
        {"db": "gds", "term": f"{accession}[ACCN]", "retmode": "json", "retmax": "20"},
    )
    ids = res.get("esearchresult", {}).get("idlist", [])
    if not ids:
        raise SystemExit(f"No GEO record found for accession {accession}")
    # Prefer the UID whose summary accession matches exactly.
    summ = get_json(
        f"{EUTILS}/esummary.fcgi",
        {"db": "gds", "id": ",".join(ids), "retmode": "json"},
    ).get("result", {})
    for uid in summ.get("uids", []):
        if summ[uid].get("accession", "").upper() == accession.upper():
            return uid, summ[uid]
    uid = ids[0]
    return uid, summ.get(uid, {})


def cmd_metadata(args):
    _, rec = resolve_uid(args.accession)
    if args.json:
        print(json.dumps(rec, indent=2))
        return
    fields = [
        "accession", "title", "taxon", "gdstype", "gpl", "n_samples",
        "pdat", "summary",
    ]
    for f in fields:
        val = rec.get(f, "")
        print(f"{f:12s}: {val}")


def cmd_samples(args):
    _, rec = resolve_uid(args.accession)
    samples = rec.get("samples", [])
    if args.json:
        print(json.dumps(samples, indent=2))
        return
    if not samples:
        print("No sample list in summary (query the SOFT file for full sample metadata).")
        return
    for s in samples:
        print(f"{s.get('accession', '')}\t{s.get('title', '')}")


def list_ftp_files(url):
    """Parse an Apache-style FTP HTML directory listing into filenames."""
    import re
    html = http_get(url).decode("utf-8", "replace")
    names = re.findall(r'href="([^"?/][^"]*)"', html)
    # Keep only real entries in this directory: no absolute URLs (footer links),
    # no protocol/query/anchor junk, no parent-dir links.
    return [n for n in names if "://" not in n and ":" not in n and n not in ("..", ".")]


def cmd_files(args):
    type_dir, stub = geo_ftp_dir(args.accession)
    base = f"{FTP_BASE}/{type_dir}/{stub}/{args.accession.upper()}"
    out = {}
    for sub in ("matrix", "soft", "miniml", "suppl"):
        try:
            files = list_ftp_files(f"{base}/{sub}/")
            if files:
                out[sub] = [f"{base}/{sub}/{f}" for f in files]
        except SystemExit:
            continue
    if args.json:
        print(json.dumps(out, indent=2))
        return
    for sub, urls in out.items():
        print(f"[{sub}]")
        for u in urls:
            print(f"  {u}")
    if not out:
        print("No files found on GEO FTP for this accession.")


def download(url, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    name = url.rstrip("/").split("/")[-1]
    dest = os.path.join(out_dir, name)
    print(f"  downloading {name} ...", file=sys.stderr)
    data = http_get(url)
    with open(dest, "wb") as fh:
        fh.write(data)
    print(f"  saved {dest} ({len(data)} bytes)", file=sys.stderr)
    return dest


def cmd_download(args):
    type_dir, stub = geo_ftp_dir(args.accession)
    base = f"{FTP_BASE}/{type_dir}/{stub}/{args.accession.upper()}"
    subs = []
    if args.matrix:
        subs.append("matrix")
    if args.soft:
        subs.append("soft")
    if args.miniml:
        subs.append("miniml")
    if args.suppl:
        subs.append("suppl")
    if not subs:
        subs = ["matrix"]  # sensible default
    for sub in subs:
        try:
            files = list_ftp_files(f"{base}/{sub}/")
        except SystemExit:
            print(f"[{sub}] not available", file=sys.stderr)
            continue
        for f in files:
            download(f"{base}/{sub}/{f}", os.path.join(args.out, args.accession.upper(), sub))


# ---- nf-core/scrnaseq samplesheet ----
# GEO does not host raw reads; they live in SRA/ENA. We resolve the linked
# SRA study / BioProject from the GEO record, then pull FASTQ links from ENA.
import csv
import re

ENA_PORTAL = "https://www.ebi.ac.uk/ena/portal/api"
SAMPLESHEET_COLS = ["sample", "fastq_1", "fastq_2"]


# --- output field hygiene: control chars (incl. CR/LF/tab) -> _, so no
# value can break the line or cell it is written into. ---
_CTRL_RE = re.compile(r"[\x00-\x1f\x7f-\x9f]")  # CR, LF, tab, other control chars


def _safe_field(v):
    """Value safe for one CSV/TSV cell: control chars -> _, runs of spaces collapsed."""
    return re.sub(r" +", " ", _CTRL_RE.sub("_", v or "")).strip()


def clean_sample(name):
    """nf-core sample id: keep the conservative safe set [A-Za-z0-9._-];
    whitespace and every other character -> _."""
    s = re.sub(r"[^A-Za-z0-9._-]+", "_", (name or "").strip())
    return s.strip("._-") or "sample"


def finalize_sample_ids(rows):
    """Clean each row's 'sample' in place; warn on rewrite, error on collision.

    A collision (two *distinct* originals normalizing to one id) is fatal: the
    pipeline concatenates rows sharing a `sample`, so it would silently pool
    different biological samples.
    """
    seen, warned = {}, set()
    for r in rows:
        orig = r.get("sample", "") or ""
        clean = clean_sample(orig)
        if clean != orig and orig not in warned:
            print(f"warning: sample {orig!r} -> {clean!r}", file=sys.stderr)
            warned.add(orig)
        if clean in seen and seen[clean] != orig:
            raise SystemExit(
                f"sample-name collision: {orig!r} and {seen[clean]!r} both normalize "
                f"to {clean!r}. Rename one or choose a different --group-by so distinct "
                "samples stay distinct (the pipeline concatenates rows sharing a sample).")
        seen.setdefault(clean, orig)
        r["sample"] = clean
    return rows


def to_https(url):
    if url.startswith("ftp://"):
        return "https://" + url[len("ftp://"):]
    if not url.startswith("http"):
        return "https://" + url
    return url


def fastq_pair(urls):
    urls = [u for u in urls if u]

    def b(u):
        return u.rsplit("/", 1)[-1].lower()

    r1 = [u for u in urls if any(t in b(u) for t in ("_r1", "_1.fastq", "_1.fq"))]
    r2 = [u for u in urls if any(t in b(u) for t in ("_r2", "_2.fastq", "_2.fq"))]
    if r1 and r2:
        return r1[0], r2[0]
    non_idx = [u for u in urls if not any(t in b(u) for t in ("_i1", "_i2", "_r3", "_3.fastq"))]
    if len(non_idx) >= 2:
        return non_idx[0], non_idx[1]
    if len(non_idx) == 1:
        return non_idx[0], ""
    if urls:
        return urls[0], (urls[1] if len(urls) > 1 else "")
    return "", ""


# ---- local FASTQ names: fasterq-dump output / cellranger symlinks ----
# (these are the names fasterq-dump itself writes; anything consuming them
#  -- including a future sra-tools skill -- must agree with this shape)
#
# The Illumina/bcl2fastq shape that cellranger's own FASTQ parser requires — a bare
# _R1.fastq.gz is NOT matched by it, while nf-core/scrnaseq accepts this form too.
CELLRANGER_R1 = "_S1_L001_R1_001.fastq.gz"
CELLRANGER_R2 = "_S1_L001_R2_001.fastq.gz"

_READ_MAP_RE = re.compile(r"([1-9]),([1-9])")


def parse_read_map(value):
    """Validate '<r1>,<r2>' 1-based read positions; None when not given."""
    if not value:
        return None
    m = _READ_MAP_RE.fullmatch(value.strip())
    if not m:
        raise SystemExit(
            f"--read-map must be '<r1>,<r2>' with 1-based positions (e.g. 3,4), got "
            f"{value!r}. For a 10x run whose technical reads are stored as separate "
            "files: 3 files (single index) -> 2,3; 4 files (dual index) -> 3,4.")
    r1, r2 = int(m.group(1)), int(m.group(2))
    if r1 == r2:
        raise SystemExit(f"--read-map R1 and R2 must differ, got {value!r}.")
    return r1, r2


def fastq_names_cellranger(run, fastq_dir):
    """Paths of the symlinks written by `sra job-scripts --cellranger-links`."""
    return (os.path.join(fastq_dir, run + CELLRANGER_R1),
            os.path.join(fastq_dir, run + CELLRANGER_R2))


def fastq_names_sra(run, layout, fastq_dir, read_map=None):
    """fasterq-dump output names: <run>_1.fastq.gz/_2.fastq.gz, or <run>.fastq.gz
    for a SINGLE layout. `read_map` selects other suffixes, e.g. (3, 4).

    These names are *constructed*, so no filename heuristic runs here. With
    --include-technical a 10x run yields 3-4 files whose _1/_2 are index/barcode
    reads rather than the cDNA pair — pass a read map, or use naming 'cellranger'.
    """
    if read_map:
        r1, r2 = read_map
        return (os.path.join(fastq_dir, f"{run}_{r1}.fastq.gz"),
                os.path.join(fastq_dir, f"{run}_{r2}.fastq.gz"))
    if (layout or "").strip().upper().startswith("PAIRED"):
        return (os.path.join(fastq_dir, f"{run}_1.fastq.gz"),
                os.path.join(fastq_dir, f"{run}_2.fastq.gz"))
    return os.path.join(fastq_dir, f"{run}.fastq.gz"), ""


def local_fastq_names(run, layout, fastq_dir, naming, read_map=None):
    """--fastq-dir dispatch: 'cellranger' symlink names, else fasterq-dump names."""
    if naming == "cellranger":
        if not (layout or "").strip().upper().startswith("PAIRED"):
            print(f"warning: {run} layout is not PAIRED; the cellranger link job "
                  "skips single-end runs, so this path may not exist", file=sys.stderr)
        return fastq_names_cellranger(run, fastq_dir)
    return fastq_names_sra(run, layout, fastq_dir, read_map)


def pick_by_read_map(urls, read_map, run=""):
    """Select (fastq_1, fastq_2) from a run's archive URLs by 1-based position.

    Bypasses fastq_pair: for a 10x run whose technical reads are separate files the
    heuristic drops the cDNA read.
    """
    ordered = sorted((u for u in urls if u), key=lambda u: u.rsplit("/", 1)[-1])
    r1, r2 = read_map
    if max(r1, r2) > len(ordered):
        raise SystemExit(
            f"--read-map {r1},{r2} needs at least {max(r1, r2)} files but run "
            f"{run or '?'} has {len(ordered)}. A study mixing chemistries needs one "
            "run of sheets per read layout.")
    return ordered[r1 - 1], ordered[r2 - 1]


def check_fastq_opts(args):
    """Guardrails for the path/naming flags; returns the parsed read map.

    The path modes are a closed, mutually
    exclusive choice.
    """
    fastq_dir = getattr(args, "fastq_dir", None)
    naming = getattr(args, "fastq_naming", None)
    if naming and not fastq_dir:
        raise SystemExit("--fastq-naming only applies with --fastq-dir DIR.")
    if fastq_dir and getattr(args, "local_dir", None):
        raise SystemExit(
            "--fastq-dir and --local-dir are mutually exclusive: --fastq-dir names the "
            "flat fasterq-dump output, --local-dir the per-run `download` layout.")
    if naming == "cellranger":
        if args.assay != "scrna":
            raise SystemExit(
                "--fastq-naming cellranger describes a 10x barcode/cDNA read pair; use "
                "--assay scrna, or --fastq-naming sra for bulk data.")
        if getattr(args, "read_map", None):
            print("warning: --read-map is ignored with --fastq-naming cellranger — the "
                  "link job already resolved which reads are R1/R2", file=sys.stderr)
            return None
    return parse_read_map(getattr(args, "read_map", None))


def rows_from_runtable(path, fastq_dir, group_by, naming="sra", read_map=None):
    """Build samplesheet rows from a local SraRunTable.csv (E-utilities runinfo CSV).

    fastq columns are fasterq-dump filenames keyed by the `Run` accession; single vs
    paired is read from the `LibraryLayout` column (SINGLE -> SRR.fastq.gz; PAIRED ->
    SRR_1.fastq.gz/SRR_2.fastq.gz — the `sra` skill's dump job gzips its output).
    No network access.
    """
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        cols = reader.fieldnames or []
        run_col = "Run" if "Run" in cols else (cols[0] if cols else "Run")
        layout_col = "LibraryLayout" if "LibraryLayout" in cols else None
        # --group-by if it names a runtable column, else fall back to a sample-ish column.
        sample_col = group_by if group_by in cols else next(
            (c for c in ("SampleName", "Sample", run_col) if c in cols), run_col)
        rows = []
        for r in reader:
            run = (r.get(run_col) or "").strip()
            if not run:
                continue
            layout = r.get(layout_col, "") if layout_col else ""
            f1, f2 = local_fastq_names(run, layout, fastq_dir, naming, read_map)
            # raw sample here; finalize_sample_ids cleans + collision-checks later
            rows.append({"sample": (r.get(sample_col) or run),
                         "fastq_1": f1, "fastq_2": f2})
    return rows


def resolve_sra_project(accession):
    """Find the SRA study / BioProject linked to a GEO accession."""
    url = ("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc="
           + accession + "&targ=self&form=text&view=brief")
    text = http_get(url).decode("utf-8", "replace")
    for pat in (r"PRJ[A-Z]+\d+", r"SRP\d+|ERP\d+|DRP\d+", r"SRX\d+|ERX\d+"):
        m = re.search(pat, text)
        if m:
            return m.group(0)
    return None


def ena_runs(project):
    fields = ("run_accession,experiment_accession,sample_accession,sample_alias,"
              "sample_title,library_layout,fastq_ftp")
    url = (f"{ENA_PORTAL}/filereport?accession={project}"
           f"&result=read_run&fields={fields}&format=tsv&limit=0")
    text = http_get(url).decode("utf-8", "replace")
    lines = [ln for ln in text.splitlines() if ln.strip()]
    if len(lines) < 2:
        return []
    header = lines[0].split("\t")
    return [dict(zip(header, ln.split("\t"))) for ln in lines[1:]]


def write_samplesheet(rows, out, assay, strandedness="auto"):
    """Write an nf-core sample sheet.

    assay 'scrna' → nf-core/scrnaseq columns (sample,fastq_1,fastq_2);
    assay 'bulk'  → nf-core/rnaseq columns (+ strandedness).
    """
    cols = list(SAMPLESHEET_COLS) + (["strandedness"] if assay == "bulk" else [])
    rows = sorted(rows, key=lambda x: (x["sample"], x["fastq_1"]))
    with open(out, "w", newline="") as fh:
        # LF, not csv.writer's default CRLF: a trailing \r ends up inside the last
        # field's value when the pipeline reads the row.
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(cols)
        for r in rows:
            row = dict(r)
            if assay == "bulk":
                row["strandedness"] = strandedness
            w.writerow([_safe_field(row.get(c, "")) for c in cols])
    return len(rows)


def cmd_samplesheet(args):
    read_map = check_fastq_opts(args)
    if args.from_runtable:
        if not args.fastq_dir:
            raise SystemExit("--from-runtable requires --fastq-dir DIR "
                             "(samplesheet fastq paths must include a directory).")
        rows = rows_from_runtable(args.from_runtable, args.fastq_dir, args.group_by,
                                  args.fastq_naming, read_map)
        if not rows:
            raise SystemExit(f"No run rows found in {args.from_runtable}.")
        finalize_sample_ids(rows)
        n = write_samplesheet(rows, args.out, args.assay, args.strandedness)
        print(f"Wrote {n} row(s) for {len({r['sample'] for r in rows})} sample(s) to {args.out} "
              f"(nf-core/{'rnaseq' if args.assay == 'bulk' else 'scrnaseq'}, "
              f"SRA fastq names from {args.from_runtable})", file=sys.stderr)
        return
    project = resolve_sra_project(args.accession)
    if not project:
        raise SystemExit(
            f"Could not find a linked SRA/BioProject for {args.accession}. "
            "Check the series relations, or query ENA directly with the SRP/PRJNA accession.")
    print(f"{args.accession} -> SRA project {project}", file=sys.stderr)
    runs = [r for r in ena_runs(project) if r.get("fastq_ftp")]
    if not runs:
        raise SystemExit(
            f"No public FASTQ found in ENA for {project} (linked to {args.accession}). "
            "The raw data may be under controlled access (e.g. dbGaP) or not mirrored to ENA.")
    rows = []
    for r in runs:
        run = r.get("run_accession", "run")
        if args.fastq_dir:
            # names are constructed from the run accession; no heuristic, no URLs
            f1, f2 = local_fastq_names(run, r.get("library_layout"), args.fastq_dir,
                                       args.fastq_naming, read_map)
        else:
            urls = [to_https(l) for l in r.get("fastq_ftp", "").split(";") if l]
            if read_map:
                f1, f2 = pick_by_read_map(urls, read_map, run)
            else:
                f1, f2 = fastq_pair(urls)
            if args.local_dir:
                f1 = os.path.join(args.local_dir, run, f1.rsplit("/", 1)[-1]) if f1 else ""
                f2 = os.path.join(args.local_dir, run, f2.rsplit("/", 1)[-1]) if f2 else ""
        sample = r.get(args.group_by) or r.get("sample_accession") or r.get("run_accession")
        rows.append({"sample": sample or "", "fastq_1": f1, "fastq_2": f2})
    finalize_sample_ids(rows)
    n = write_samplesheet(rows, args.out, args.assay, args.strandedness)
    print(f"Wrote {n} row(s) for {len({r['sample'] for r in rows})} sample(s) to {args.out} "
          f"(nf-core/{'rnaseq' if args.assay == 'bulk' else 'scrnaseq'})", file=sys.stderr)


# ---- SRA Run Selector files (SraRunTable.csv + SRR_Acc_List.txt) ----
import io


def sra_runinfo_csv(project):
    """Fetch the SRA `runinfo` CSV for a study/project via E-utilities.

    This is the same table the SRA Run Selector serves as SraRunTable.csv. We
    esearch db=sra (using history so large studies work) then efetch rettype=runinfo.
    """
    res = get_json(
        f"{EUTILS}/esearch.fcgi",
        {"db": "sra", "term": project, "retmode": "json",
         "usehistory": "y", "retmax": "0"},
    )
    er = res.get("esearchresult", {})
    webenv, qk = er.get("webenv"), er.get("querykey")
    if not webenv or not qk or er.get("count", "0") == "0":
        return ""
    url = f"{EUTILS}/efetch.fcgi?" + urllib.parse.urlencode(_eutil_params(
        {"db": "sra", "query_key": qk, "WebEnv": webenv,
         "rettype": "runinfo", "retmode": "text"}))
    return http_get(url).decode("utf-8", "replace")


def cmd_runtable(args):
    """Download the SRA Run Selector files for a GEO accession.

    Resolves the linked SRA study / BioProject, then writes SraRunTable.csv (the
    full runinfo table) and SRR_Acc_List.txt (one run accession per line).
    """
    project = resolve_sra_project(args.accession)
    if not project:
        raise SystemExit(
            f"Could not find a linked SRA/BioProject for {args.accession}. "
            "The series may have no sequencing data in SRA.")
    print(f"{args.accession} -> SRA project {project}", file=sys.stderr)
    text = sra_runinfo_csv(project)
    if not text.strip():
        raise SystemExit(
            f"No SRA runinfo returned for {project} (linked to {args.accession}). "
            "The raw data may be under controlled access or not in SRA.")

    os.makedirs(args.out, exist_ok=True)
    run_table = os.path.join(args.out, "SraRunTable.csv")
    with open(run_table, "w", newline="") as fh:
        fh.write(text if text.endswith("\n") else text + "\n")

    rows = [r for r in csv.reader(io.StringIO(text)) if r]
    header = rows[0] if rows else []
    run_idx = header.index("Run") if "Run" in header else 0
    accs = [r[run_idx].strip() for r in rows[1:] if len(r) > run_idx and r[run_idx].strip()]
    acc_list = os.path.join(args.out, "SRR_Acc_List.txt")
    with open(acc_list, "w") as fh:
        fh.write("\n".join(accs) + ("\n" if accs else ""))

    print(f"Wrote {run_table} ({len(accs)} run(s)) and {acc_list}", file=sys.stderr)


# ---- harmonized metadata.tsv (shared shape across skills) ----
import html

METADATA_COLS = ["sample", "replicate", "species", "sex", "age",
                 "condition", "genotype", "treatment", "tissue"]
NA = "NA"

FIELD_KEYS = {
    "species": ["organism", "scientific_name", "species", "organism scientific name"],
    "sex": ["sex", "gender"],
    "age": ["age", "developmental stage", "dev stage", "age at collection", "age at sampling"],
    "condition": ["disease", "disease state", "condition", "phenotype",
                  "health state", "clinical information", "diagnosis"],
    "genotype": ["genotype", "genotype/variation", "variation", "strain",
                 "strain/background", "background"],
    "treatment": ["treatment", "agent", "compound", "stimulus",
                  "perturbation", "dose", "treatment protocol"],
    "tissue": ["tissue", "organism part", "tissue type", "tissue region", "source tissue"],
}
_REPLICATE_KEYS = ("replicate", "biological replicate", "technical replicate",
                   "replicate number")
_IDENTITY_KEYS = ("source name", "sample", "sample name", "sample_accession",
                  "run_accession", "assay name", "name", "title")
_NULLS = {"", "na", "n/a", "none", "null", "unknown", "not applicable",
          "not available", "not collected", "not provided", "missing", "--"}


def _norm_key(k):
    k = (k or "").strip().lower()
    m = re.match(r"(?:characteristics|comment|factor\s*value|factorvalue)\s*\[(.+)\]$", k)
    return m.group(1).strip() if m else k


def _clean_val(v):
    v = re.sub(r" +", " ", _CTRL_RE.sub("_", html.unescape(v or ""))).strip()
    return "" if v.lower() in _NULLS else v


def _col_name(k):
    """Column-safe name for a promoted extra characteristic."""
    return re.sub(r"[^0-9a-z]+", "_", k.strip().lower()).strip("_") or "field"


def harmonize_row(sample, replicate, attrs):
    """Map source-native (name -> value) annotations onto the metadata columns.

    Core fields (METADATA_COLS) are always present; missing ones become 'NA'. Every
    other characteristic is promoted to its own column so nothing is dropped.
    """
    norm = {}
    for k, v in (attrs or {}).items():
        nk, cv = _norm_key(k), _clean_val(v)
        if nk and cv and nk not in norm:
            norm[nk] = cv
    row = {"sample": _clean_val(sample) or NA,
           "replicate": _clean_val(replicate) or NA}
    used = set(_REPLICATE_KEYS) | set(_IDENTITY_KEYS)
    for field, keys in FIELD_KEYS.items():
        used.update(keys)  # every synonym is consumed, not just the matched one
        val = ""
        for cand in keys:
            if norm.get(cand):
                val = norm[cand]
                break
        row[field] = val or NA
    for k, v in norm.items():
        if k not in used:
            row.setdefault(_col_name(k), v)
    return row


def write_metadata_tsv(rows, out):
    """Write rows as TSV; header = core columns + union of any extra columns."""
    header = list(METADATA_COLS)
    for r in rows:
        for k in r:
            if k not in header:
                header.append(k)
    with open(out, "w", newline="") as fh:
        # LF, not csv.writer's default CRLF, so the TSV is byte-identical
        # across platforms.
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        for r in rows:
            w.writerow([r.get(c, NA) for c in header])
    return len(rows), len(header)


def geo_sample_soft(gsm):
    url = ("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc="
           + gsm + "&targ=self&form=text&view=quick")
    return http_get(url).decode("utf-8", "replace")


def parse_geo_sample(text):
    """Pull organism + `!Sample_characteristics_ch*` (tag: value) from SOFT text."""
    attrs = {}
    organism = ""
    for ln in text.splitlines():
        ln = ln.strip()
        if "=" not in ln:
            continue
        key, _, val = ln.partition("=")
        key, val = key.strip(), val.strip()
        if key.startswith("!Sample_organism") and val:
            organism = val
        elif key.startswith("!Sample_characteristics") and val:
            if ":" in val:
                tag, tv = val.split(":", 1)
                attrs.setdefault(tag.strip(), tv.strip())
            else:
                attrs.setdefault(val, "")
    if organism:
        attrs.setdefault("organism", organism)
    return attrs


def cmd_metadata_table(args):
    """Write a harmonized metadata.tsv (one row per sample/GSM).

    Per-sample characteristics are read from each GSM's SOFT record; each GSM is a
    row, with `replicate` taken from a replicate characteristic when present else 1.
    """
    _, rec = resolve_uid(args.accession)
    samples = rec.get("samples", [])
    if not samples:
        raise SystemExit(
            f"No sample list in the GEO summary for {args.accession}; "
            "cannot build a per-sample metadata table.")
    rows = []
    for s in samples:
        gsm = s.get("accession", "")
        attrs = {}
        if gsm:
            try:
                attrs = parse_geo_sample(geo_sample_soft(gsm))
            except SystemExit:
                attrs = {}
        attrs.setdefault("organism", rec.get("taxon", ""))
        replicate = ""
        for k, v in attrs.items():
            if _norm_key(k) in _REPLICATE_KEYS and v:
                replicate = v
                break
        rows.append(harmonize_row(gsm or s.get("title", ""), replicate or "1", attrs))
    n, ncols = write_metadata_tsv(rows, args.out)
    print(f"Wrote {n} sample row(s) x {ncols} column(s) for {args.accession} to {args.out}",
          file=sys.stderr)


def cmd_search(args):
    term = args.query
    if args.organism:
        term += f' AND "{args.organism}"[Organism]'
    if args.type:
        term += f" AND {args.type}[ETYP]"
    res = get_json(
        f"{EUTILS}/esearch.fcgi",
        {"db": "gds", "term": term, "retmode": "json", "retmax": str(args.limit)},
    )
    ids = res.get("esearchresult", {}).get("idlist", [])
    count = res.get("esearchresult", {}).get("count", "0")
    if not ids:
        print(f"0 results for: {term}")
        return
    summ = get_json(
        f"{EUTILS}/esummary.fcgi",
        {"db": "gds", "id": ",".join(ids), "retmode": "json"},
    ).get("result", {})
    if args.json:
        print(json.dumps([summ[u] for u in summ.get("uids", [])], indent=2))
        return
    print(f"{count} total hits; showing {len(ids)}:\n")
    for uid in summ.get("uids", []):
        r = summ[uid]
        print(f"{r.get('accession', ''):12s} {r.get('gdstype', ''):20s} n={r.get('n_samples', '?')}  {r.get('title', '')}")


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    # Shared opt-in: send NCBI credentials only when the user explicitly asks.
    cred = argparse.ArgumentParser(add_help=False)
    cred.add_argument("--use-ncbi-credentials", action="store_true",
                      help="opt in to send NCBI_EMAIL / NCBI_API_KEY (read from the "
                           "environment) with E-utilities requests to lift NCBI rate "
                           "limits; OFF by default even when those env vars are set")

    m = sub.add_parser("metadata", parents=[cred], help="Series/dataset summary via E-utilities")
    m.add_argument("accession")
    m.add_argument("--json", action="store_true")
    m.set_defaults(func=cmd_metadata)

    s = sub.add_parser("samples", parents=[cred], help="List samples (GSMs) in a series")
    s.add_argument("accession")
    s.add_argument("--json", action="store_true")
    s.set_defaults(func=cmd_samples)

    f = sub.add_parser("files", parents=[cred], help="List downloadable files on GEO FTP")
    f.add_argument("accession")
    f.add_argument("--json", action="store_true")
    f.set_defaults(func=cmd_files)

    d = sub.add_parser("download", parents=[cred], help="Download files from GEO FTP")
    d.add_argument("accession")
    d.add_argument("--out", default="./geo_out")
    d.add_argument("--matrix", action="store_true", help="series matrix (default)")
    d.add_argument("--soft", action="store_true", help="SOFT family file")
    d.add_argument("--miniml", action="store_true", help="MINiML file")
    d.add_argument("--suppl", action="store_true", help="supplementary/raw files")
    d.set_defaults(func=cmd_download)

    q = sub.add_parser("search", parents=[cred], help="Free-text search of GEO DataSets")
    q.add_argument("query")
    q.add_argument("--organism", help='e.g. "Homo sapiens"')
    q.add_argument("--type", help="entry type filter, e.g. gse, gds")
    q.add_argument("--limit", type=int, default=20)
    q.add_argument("--json", action="store_true")
    q.set_defaults(func=cmd_search)

    mt = sub.add_parser("metadata-table", parents=[cred],
                        help="Write a harmonized metadata.tsv (one row per sample/GSM)")
    mt.add_argument("accession", help="GSE accession")
    mt.add_argument("--out", default="metadata.tsv")
    mt.set_defaults(func=cmd_metadata_table)

    rt = sub.add_parser("runtable", parents=[cred],
                        help="Download SRA Run Selector files (SraRunTable.csv + SRR_Acc_List.txt)")
    rt.add_argument("accession", help="GSE (or GSM) accession")
    rt.add_argument("--out", default=".", help="output directory (default current dir)")
    rt.set_defaults(func=cmd_runtable)

    ss = sub.add_parser("samplesheet", parents=[cred],
                        help="Build an nf-core samplesheet.csv via the linked SRA/ENA data (--assay scrna|bulk)")
    ss.add_argument("accession", help="GSE (or GSM) accession")
    ss.add_argument("--assay", choices=["scrna", "bulk"], required=True,
                    help="scrna -> nf-core/scrnaseq (sample,fastq_1,fastq_2); "
                         "bulk -> nf-core/rnaseq (adds a strandedness column)")
    ss.add_argument("--strandedness", choices=["auto", "forward", "reverse", "unstranded"],
                    default="auto",
                    help="value for the rnaseq strandedness column (only --assay bulk; default auto)")
    ss.add_argument("--out", default="samplesheet.csv")
    ss.add_argument("--group-by", default="sample_accession",
                    help="ENA run field for the `sample` column "
                         "(sample_accession, sample_alias, sample_title, experiment_accession)")
    ss.add_argument("--local-dir",
                    help="write local paths <dir>/<run>/<file> instead of ENA URLs")
    ss.add_argument("--from-runtable", metavar="PATH",
                    help="build fastq names offline from a local SraRunTable.csv: "
                         "SRR.fastq.gz (SINGLE) or SRR_1.fastq.gz/SRR_2.fastq.gz "
                         "(PAIRED), per the LibraryLayout column. Skips ENA. "
                         "Requires --fastq-dir.")
    ss.add_argument("--fastq-dir", metavar="DIR",
                    help="the flat fasterq-dump output dir (the `sra` skill's "
                         "--fastq-dir); fastq names are built from the run accession "
                         "instead of ENA URLs. Required with --from-runtable; mutually "
                         "exclusive with --local-dir")
    ss.add_argument("--fastq-naming", choices=["sra", "cellranger"],
                    help="with --fastq-dir: 'sra' -> <run>_1.fastq.gz/_2.fastq.gz "
                         "(default), 'cellranger' -> the "
                         "<run>_S1_L001_R{1,2}_001.fastq.gz symlinks written by "
                         "`sra job-scripts --cellranger-links` (--assay scrna only)")
    ss.add_argument("--read-map", metavar="R1,R2",
                    help="declare which reads are the cDNA pair, 1-based (e.g. 3,4). "
                         "Needed for a 10x run whose technical reads are separate files: "
                         "3 files -> 2,3; 4 files (dual index) -> 3,4. Without it the "
                         "R1/R2 filename heuristic runs and would drop the cDNA read")
    ss.set_defaults(func=cmd_samplesheet)

    args = p.parse_args(argv)
    if getattr(args, "use_ncbi_credentials", False):
        global _USE_CREDENTIALS
        _USE_CREDENTIALS = True
    args.func(args)


if __name__ == "__main__":
    main()
