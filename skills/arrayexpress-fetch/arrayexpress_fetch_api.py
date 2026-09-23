#!/usr/bin/env python3
# Vendored from UKDRI/informatics_data_skills @ 7cc3e6e (arrayexpress/scripts/arrayexpress.py).
# Copyright (c) 2026 UK Dementia Research Institute. Licensed MIT.
# Kept close to upstream so changes there remain easy to merge; the ClawBio
# CLI contract lives in arrayexpress_fetch.py, not here.
"""Query metadata and download data from ArrayExpress (functional genomics).

ArrayExpress is now the functional-genomics collection inside BioStudies, so this
uses the BioStudies REST API (https://www.ebi.ac.uk/biostudies/api/v1) with
ArrayExpress-aware helpers for MAGE-TAB (IDF/SDRF) and raw/processed data. File
bytes come from the host /studies/{accession}/info advertises as `httpLink`,
which is not the same tree for every collection -- see file_url().

Standard library only (urllib) -- no pip install required.

Examples:
    python arrayexpress_fetch.py --command metadata --accession E-MTAB-11448
    python arrayexpress_fetch.py --command files --accession E-MTAB-11448
    python arrayexpress_fetch.py --command sdrf --accession E-MTAB-11448
    python arrayexpress_fetch.py --command download --accession E-MTAB-11448 --magetab
    python arrayexpress_fetch.py --command download --accession E-MTAB-11448 --processed
    python arrayexpress_fetch.py --command search --query "single cell heart" --limit 20
    python arrayexpress_fetch.py --command metadata-table --accession E-MTAB-11448
"""
import argparse
import json
import os
import sys
import time
import urllib.parse
import urllib.request

API = "https://www.ebi.ac.uk/biostudies/api/v1"
FILES = "https://www.ebi.ac.uk/biostudies/files"
COLLECTION = "ArrayExpress"
USER_AGENT = "arrayexpress-skill/1.0"


def http_get(url, retries=3):
    last = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=180) as r:
                return r.read()
        except Exception as e:  # noqa: BLE001
            last = e
            time.sleep(1.5 * (attempt + 1))
    raise SystemExit(f"GET failed for {url}: {last}")


def get_json(url):
    return json.loads(http_get(url).decode("utf-8", "replace"))


def study_json(accession):
    return get_json(f"{API}/studies/{accession}")


def attrs_to_dict(attributes):
    return {a.get("name"): a.get("value") for a in (attributes or []) if isinstance(a, dict)}


def walk_files(node):
    if isinstance(node, dict):
        if node.get("type") == "file" and "path" in node:
            yield node
        for v in node.values():
            yield from walk_files(v)
    elif isinstance(node, list):
        for v in node:
            yield from walk_files(v)


def list_files(accession):
    d = study_json(accession)
    return list(walk_files(d.get("section", {})))


def classify(file_rec):
    """Classify a file entry as idf/sdrf/processed/raw/other.

    Prefer the file's own MAGE-TAB attributes (Type/Description), then fall back
    to the filename.
    """
    path = file_rec["path"] if isinstance(file_rec, dict) else str(file_rec)
    p = path.lower()
    if p.endswith(".idf.txt"):
        return "idf"
    if p.endswith(".sdrf.txt"):
        return "sdrf"
    hint = ""
    if isinstance(file_rec, dict):
        hint = " ".join(attrs_to_dict(file_rec.get("attributes")).values()).lower()
    if "processed" in hint or "processed" in p or \
            p.endswith((".counts.txt", ".tsv", ".mtx", ".h5", ".h5ad", ".rds")):
        return "processed"
    if "raw" in hint or p.endswith((".fastq.gz", ".fq.gz", ".cel", ".cel.gz", ".bam", ".cram")):
        return "raw"
    return "other"


def cmd_metadata(args):
    d = study_json(args.accession)
    if args.json:
        print(json.dumps(d, indent=2))
        return
    sec = d.get("section", {})
    sa = attrs_to_dict(sec.get("attributes"))
    ta = attrs_to_dict(d.get("attributes"))
    print(f"accession   : {d.get('accno','')}")
    print(f"title       : {sa.get('Title', ta.get('Title',''))}")
    print(f"release     : {ta.get('ReleaseDate','')}")
    if "Organism" in sa:
        print(f"organism    : {sa['Organism']}")
    if "Description" in sa:
        print(f"\ndescription : {sa['Description']}")
    files = list(walk_files(sec))
    kinds = {}
    for f in files:
        k = classify(f)
        kinds[k] = kinds.get(k, 0) + 1
    print(f"\nfiles       : {len(files)} total  " +
          "  ".join(f"{k}={v}" for k, v in sorted(kinds.items())))


def cmd_files(args):
    files = list_files(args.accession)
    if args.json:
        print(json.dumps(files, indent=2))
        return
    print(f"{len(files)} file(s) for {args.accession}:\n")
    for f in files:
        print(f"{classify(f):10s} {str(f.get('size','?')):>14s}  {f['path']}")


_INFO_CACHE = {}


def study_info(accession):
    """Cached /studies/{accession}/info, which carries the real download base."""
    if accession not in _INFO_CACHE:
        try:
            _INFO_CACHE[accession] = get_json(f"{API}/studies/{accession}/info")
        except SystemExit:
            _INFO_CACHE[accession] = {}
    return _INFO_CACHE[accession]


def file_url(accession, path):
    """Build the download URL for one attached file.

    DIVERGES FROM UPSTREAM, deliberately, and identically to biostudies-fetch.
    Upstream hardcodes `https://www.ebi.ac.uk/biostudies/files/{acc}/{path}`.
    That route is NOT broken -- re-checked live 2026-09-22, it answers 302 and
    redirects to the very URL this function builds. Two reasons to resolve it
    directly anyway:

    1. There is no single base tree to hardcode. /studies/{accession}/info
       advertises a different `httpLink` per collection:
           E-MTAB-10030 -> ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/030/E-MTAB-10030
           S-BSST2074   -> ftp.ebi.ac.uk/pub/databases/biostudies/S-BSST/074/S-BSST2074
       Upstream's constant works only because the redirect papers over this.
    2. The redirect is slow and unreliable for these payloads. Measured
       2026-09-22: 0.3 s direct vs 12.4 s via www, and a 24 GB file timed out
       entirely on the www route.

    Files live under `Files/` beneath the advertised base. The legacy path stays
    as a fallback for when /info is unreachable or omits httpLink.
    """
    quoted = urllib.parse.quote(path)
    base = (study_info(accession) or {}).get("httpLink")
    if base:
        return f"{base.rstrip('/')}/Files/{quoted}"
    return f"{FILES}/{accession}/{quoted}"


def cmd_sdrf(args):
    files = list_files(args.accession)
    sdrf = [f for f in files if classify(f) == "sdrf"]
    if not sdrf:
        raise SystemExit("No SDRF file found for this accession.")
    data = http_get(file_url(args.accession, sdrf[0]["path"]))
    sys.stdout.write(data.decode("utf-8", "replace"))


def download_file(accession, path, out_dir):
    dest = os.path.join(out_dir, path)
    os.makedirs(os.path.dirname(dest) or ".", exist_ok=True)
    if os.path.exists(dest):
        print(f"  exists, skipping {path}", file=sys.stderr)
        return
    print(f"  downloading {path} ...", file=sys.stderr)
    req = urllib.request.Request(file_url(accession, path), headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=600) as r, open(dest + ".part", "wb") as fh:
        while True:
            chunk = r.read(1 << 20)
            if not chunk:
                break
            fh.write(chunk)
    os.replace(dest + ".part", dest)
    print(f"  saved {dest}", file=sys.stderr)


def cmd_download(args):
    files = list_files(args.accession)
    wanted = set()
    if args.magetab:
        wanted |= {"idf", "sdrf"}
    if args.processed:
        wanted.add("processed")
    if args.raw:
        wanted.add("raw")
    if not wanted:
        wanted = {"idf", "sdrf", "processed", "raw", "other"}  # everything
    sel = [f for f in files if classify(f) in wanted]
    if not sel:
        raise SystemExit("No matching files for the selected categories.")
    out = os.path.join(args.out, args.accession)
    for f in sel:
        download_file(args.accession, f["path"], out)
    print(f"Downloaded {len(sel)} file(s) to {out}", file=sys.stderr)


# ---- nf-core/scrnaseq samplesheet (built from the SDRF) ----
import csv
import io
import re

ENA_PORTAL = "https://www.ebi.ac.uk/ena/portal/api"
SAMPLESHEET_COLS = ["sample", "fastq_1", "fastq_2"]


# --- output field hygiene ---
# Rule: control characters, CR/LF/tab included, become "_", and files are
# written with LF endings. One stray character in a source annotation would
# otherwise split a row and corrupt the whole TSV/CSV.
_CTRL_RE = re.compile(r"[\x00-\x1f\x7f-\x9f]")  # CR, LF, tab, other control chars


def _safe_field(v):
    """Value safe for one CSV/TSV cell: control chars -> _, runs of spaces collapsed."""
    return re.sub(r" +", " ", _CTRL_RE.sub("_", v or "")).strip()


def clean_sample(name):
    """nf-core sample id: keep the conservative safe set [A-Za-z0-9._-];
    whitespace and every other character -> _.

    Deliberately conservative: these ids become filenames and channel keys in
    a pipeline, so anything outside the safe set is replaced rather than
    escaped."""
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
# Rule: these names must match what fasterq-dump actually wrote, which is why
# the naming mode is explicit rather than guessed. The same sample-sheet
# contract is shared with the other archive skills, so do not diverge here.
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
            print(f"warning: {run} looks single-end; the cellranger link job skips "
                  "single-end runs, so this path may not exist", file=sys.stderr)
        return fastq_names_cellranger(run, fastq_dir)
    return fastq_names_sra(run, layout, fastq_dir, read_map)


def pick_by_read_map(urls, read_map, run=""):
    """Select (fastq_1, fastq_2) from a run's archive URLs by 1-based position.

    Bypasses fastq_pair: for a 10x run whose technical reads are separate files
    the filename heuristic drops the cDNA read, so --read-map must be given
    explicitly. This is the highest-risk behaviour in the skill -- a silently
    wrong read assignment produces a plausible but meaningless quantification.
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

    --local-dir, --fastq-dir and the default remote-URL mode are a closed,
    mutually exclusive set: exactly one decides where the sample sheet's paths
    point. Accepting two would silently produce a sheet mixing local and remote
    paths.
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


def get_sdrf_text(accession):
    files = list_files(accession)
    sdrf = [f for f in files if classify(f) == "sdrf"]
    if not sdrf:
        raise SystemExit(f"No SDRF file found for {accession}.")
    return http_get(file_url(accession, sdrf[0]["path"])).decode("utf-8", "replace")


def ena_run_fastq(run):
    url = (f"{ENA_PORTAL}/filereport?accession={run}"
           f"&result=read_run&fields=fastq_ftp&format=tsv&limit=0")
    text = http_get(url).decode("utf-8", "replace")
    links = []
    for ln in text.splitlines()[1:]:
        links.extend(x for x in ln.split("\t")[-1].split(";") if x)
    return [to_https(x) for x in links]


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


def sdrf_col(header, name):
    """Index of `name` in an SDRF header, case-insensitively."""
    for i, h in enumerate(header):
        if h.lower() == name.lower():
            return i
    return None


def group_sdrf_rows(rows, header):
    """Collapse SDRF lines into one group per sample x replicate.

    A conventional bulk paired-end SDRF puts **each FASTQ on its own line**, so
    a 12-sample study has 24 lines. Both the samplesheet and the metadata table
    are per sample x replicate, so both group here rather than each rolling
    their own -- they disagreed before, and the metadata table shipped one row
    per FASTQ.

    The key is `Comment[ENA_RUN]` falling back to `Source Name`: keying on the
    run keeps a sample sequenced across several runs distinguishable, which is
    what "x replicate" means, while still merging the mates of one run.

    Returns [(key, {"sample", "run", "rows"[]})] in first-seen order.
    """
    si = sdrf_col(header, "Source Name")
    ri = sdrf_col(header, "Comment[ENA_RUN]")
    groups, order = {}, []
    for idx, row in enumerate(rows):
        if not any(c.strip() for c in row):
            continue
        sample = row[si].strip() if si is not None and si < len(row) else ""
        run = row[ri].strip() if ri is not None and ri < len(row) else ""
        key = run or sample or f"row{idx}"
        g = groups.get(key)
        if g is None:
            g = {"sample": sample, "run": run, "rows": []}
            groups[key] = g
            order.append(key)
        g["rows"].append(row)
        if sample and not g["sample"]:
            g["sample"] = sample
    return [(k, groups[k]) for k in order]


def cmd_samplesheet(args):
    read_map = check_fastq_opts(args)
    rows = list(csv.reader(io.StringIO(get_sdrf_text(args.accession)), delimiter="\t"))
    header = [h.strip() for h in rows[0]]

    def col(name):
        return sdrf_col(header, name)

    si = col("Source Name")
    fi = [i for i, h in enumerate(header) if h.lower() == "comment[fastq_uri]"]
    if si is None:
        raise SystemExit("SDRF has no 'Source Name' column; cannot build a samplesheet.")

    groups = {}
    order = []
    for key, g in group_sdrf_rows(rows[1:], header):
        g = {"sample": g["sample"], "run": g["run"],
             "uris": [to_https(row[i].strip()) for row in g["rows"] for i in fi
                      if i < len(row) and row[i].strip()]}
        groups[key] = g
        order.append(key)

    out_rows = []
    for key in order:
        g = groups[key]
        urls = g["uris"]
        if not urls and g["run"]:
            urls = ena_run_fastq(g["run"])  # fall back to ENA for this run
        if not urls:
            continue
        if args.fastq_dir:
            # Names are constructed from the run accession, so there must be one.
            if not g["run"]:
                raise SystemExit(
                    f"--fastq-dir needs a run accession to build filenames, but SDRF "
                    f"group {key!r} has no Comment[ENA_RUN]. Drop --fastq-dir, or add "
                    "the run accessions.")
            # The SDRF has no layout column; the URI count is the available signal.
            layout = "PAIRED" if len(urls) >= 2 else "SINGLE"
            f1, f2 = local_fastq_names(g["run"], layout, args.fastq_dir,
                                       args.fastq_naming, read_map)
            out_rows.append({"sample": (g["sample"] or g["run"] or ""),
                             "fastq_1": f1, "fastq_2": f2})
            continue
        if read_map:
            f1, f2 = pick_by_read_map(urls, read_map, g["run"])
        else:
            f1, f2 = fastq_pair(urls)
        if args.local_dir and g["run"]:
            f1 = os.path.join(args.local_dir, g["run"], f1.rsplit("/", 1)[-1]) if f1 else ""
            f2 = os.path.join(args.local_dir, g["run"], f2.rsplit("/", 1)[-1]) if f2 else ""
        # raw sample here; finalize_sample_ids cleans + collision-checks later
        out_rows.append({"sample": (g["sample"] or g["run"] or ""),
                         "fastq_1": f1, "fastq_2": f2})

    if not out_rows:
        raise SystemExit(
            "No FASTQ files found in the SDRF (this may be an array study with no "
            "sequencing data, or reads are only in ENA under controlled access).")
    finalize_sample_ids(out_rows)
    n = write_samplesheet(out_rows, args.out, args.assay, args.strandedness)
    print(f"Wrote {n} row(s) for {len({r['sample'] for r in out_rows})} "
          f"sample(s) to {args.out} "
          f"(nf-core/{'rnaseq' if args.assay == 'bulk' else 'scrnaseq'})", file=sys.stderr)


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


def _is_attr_col(header):
    return bool(re.match(r"(?:characteristics|factor\s*value|factorvalue)\s*\[",
                         header.strip().lower()))


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
        # LF, not csv.writer's default CRLF: a CR in a TSV cell breaks the file
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        for r in rows:
            w.writerow([r.get(c, NA) for c in header])
    return len(rows), len(header)


def ebi_biosample_attrs(accession):
    """Fetch characteristics from the EBI BioSamples API for a BioSample id.

    Returns {} for non-BioSample identifiers (only SAME*/SAMEA*/SAMN*/SAMD*).
    """
    if not re.match(r"^SAM[NED]", (accession or "").strip(), re.I):
        return {}
    try:
        data = json.loads(
            http_get("https://www.ebi.ac.uk/biosamples/samples/" + accession)
            .decode("utf-8", "replace"))
    except SystemExit:
        return {}
    attrs = {}
    for name, vals in (data.get("characteristics") or {}).items():
        if name.strip().lower().startswith(
                ("ena-", "ena ", "arrayexpress-", "insdc", "external id",
                 "ncbi submission", "sra accession", "submitter", "submission",
                 "broker")):
            continue  # skip archive bookkeeping, keep biological characteristics
        if isinstance(vals, list) and vals and isinstance(vals[0], dict):
            text = vals[0].get("text", "")
            if text:
                attrs[name] = text
    return attrs


def merge_biosample(sample_id, attrs):
    """Fill missing/extra fields from EBI BioSamples when sample_id is a BioSample."""
    bs = ebi_biosample_attrs(sample_id)
    if not bs:
        return attrs
    merged = dict(attrs)
    for k, v in bs.items():
        merged.setdefault(k, v)
    return merged


def cmd_metadata_table(args):
    """Write a harmonized metadata.tsv, one row per sample x replicate.

    Characteristics[...] and FactorValue[...] columns become the harmonized fields;
    `Source Name` is the sample and a technical/biological replicate column (or the
    ENA run, or 1) the replicate.

    Grouping is not optional. A bulk paired-end SDRF has one line per FASTQ, so
    iterating lines directly emitted two identical rows per sample -- 24 rows
    for the 12 samples of E-MTAB-5688 -- and made two BioSamples lookups where
    one was needed. `group_sdrf_rows` is the same collapse the samplesheet uses.
    """
    table = list(csv.reader(io.StringIO(get_sdrf_text(args.accession)), delimiter="\t"))
    if not table:
        raise SystemExit("SDRF is empty.")
    header = [h.strip() for h in table[0]]
    rows = []
    for _key, group in group_sdrf_rows(table[1:], header):
        sample = group["sample"]
        replicate = ""
        attrs = {}
        bios = ""
        # Merge the group's lines: they describe one sample x replicate and
        # differ only in their file columns, but take the first non-empty value
        # of each attribute rather than assuming line 1 carries them all.
        for raw in group["rows"]:
            for i, h in enumerate(header):
                v = raw[i].strip() if i < len(raw) else ""
                if not v:
                    continue
                if not replicate and _norm_key(h) in _REPLICATE_KEYS:
                    replicate = v
                if _is_attr_col(h):
                    attrs.setdefault(h, v)
            # ArrayExpress stores the BioSample id in the sample or Comment[BioSD_SAMPLE].
            if not bios and re.match(r"^SAM[NED]", sample, re.I):
                bios = sample
            if not bios:
                for i in range(len(header)):
                    cell = raw[i].strip() if i < len(raw) else ""
                    if re.match(r"^SAM[NED][A-Z]?\d+$", cell, re.I):
                        bios = cell
                        break
        if bios:
            attrs = merge_biosample(bios, attrs)
        rows.append(harmonize_row(sample, replicate or group["run"] or "1", attrs))
    if not rows:
        raise SystemExit("SDRF has no data rows.")
    n, ncols = write_metadata_tsv(rows, args.out)
    print(f"Wrote {n} sample x replicate row(s) x {ncols} column(s) for "
          f"{args.accession} to {args.out}", file=sys.stderr)


def cmd_search(args):
    params = {"query": args.query, "collection": COLLECTION,
              "pageSize": args.limit, "page": 1}
    res = get_json(f"{API}/search?" + urllib.parse.urlencode(params))
    if args.json:
        print(json.dumps(res, indent=2))
        return
    hits = res.get("hits", [])
    print(f"{res.get('totalHits','?')} total hits; showing {len(hits)}:\n")
    for h in hits:
        print(f"{h.get('accession',''):16s} {h.get('title','')}")


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    m = sub.add_parser("metadata", help="Study metadata + file breakdown")
    m.add_argument("accession")
    m.add_argument("--json", action="store_true")
    m.set_defaults(func=cmd_metadata)

    f = sub.add_parser("files", help="List files, classified (idf/sdrf/raw/processed)")
    f.add_argument("accession")
    f.add_argument("--json", action="store_true")
    f.set_defaults(func=cmd_files)

    s = sub.add_parser("sdrf", help="Print the SDRF (sample/data relationship) table")
    s.add_argument("accession")
    s.set_defaults(func=cmd_sdrf)

    d = sub.add_parser("download", help="Download files by category")
    d.add_argument("accession")
    d.add_argument("--magetab", action="store_true", help="IDF + SDRF")
    d.add_argument("--processed", action="store_true", help="processed data matrices")
    d.add_argument("--raw", action="store_true", help="raw data (FASTQ/CEL/BAM)")
    d.add_argument("--out", default="./ae_out")
    d.set_defaults(func=cmd_download)

    q = sub.add_parser("search", help="Search the ArrayExpress collection")
    q.add_argument("query")
    q.add_argument("--limit", type=int, default=20)
    q.add_argument("--json", action="store_true")
    q.set_defaults(func=cmd_search)

    mt = sub.add_parser("metadata-table",
                        help="Write a harmonized metadata.tsv (one row per sample x replicate)")
    mt.add_argument("accession")
    mt.add_argument("--out", default="metadata.tsv")
    mt.set_defaults(func=cmd_metadata_table)

    ss = sub.add_parser("samplesheet",
                        help="Build an nf-core samplesheet.csv from the SDRF (--assay scrna|bulk)")
    ss.add_argument("accession")
    ss.add_argument("--assay", choices=["scrna", "bulk"], required=True,
                    help="scrna -> nf-core/scrnaseq (sample,fastq_1,fastq_2); "
                         "bulk -> nf-core/rnaseq (adds a strandedness column)")
    ss.add_argument("--strandedness", choices=["auto", "forward", "reverse", "unstranded"],
                    default="auto",
                    help="value for the rnaseq strandedness column (only --assay bulk; default auto)")
    ss.add_argument("--out", default="samplesheet.csv")
    ss.add_argument("--local-dir",
                    help="write local paths <dir>/<run>/<file> instead of FASTQ URLs")
    ss.add_argument("--fastq-dir", metavar="DIR",
                    help="write paths into the flat fasterq-dump output dir (the `sra` "
                         "skill's --fastq-dir) instead of FASTQ URLs; names are built "
                         "from Comment[ENA_RUN], which must be present. Mutually "
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
    args.func(args)


if __name__ == "__main__":
    main()
