#!/usr/bin/env python3
# Vendored from UKDRI/informatics_data_skills @ 7cc3e6e (pride/scripts/pride.py).
# Copyright (c) 2026 UK Dementia Research Institute. Licensed MIT.
# Kept close to upstream so changes there remain easy to merge; the ClawBio
# CLI contract lives in pride_fetch.py, not here.
"""Query metadata and download data from the PRIDE Archive (proteomics).

Uses the PRIDE Archive REST API v3
(https://www.ebi.ac.uk/pride/ws/archive/v3). File downloads use the public FTP
locations advertised in each file record.

Standard library only (urllib) -- no pip install required.

Examples:
    python pride_fetch.py metadata PXD000001
    python pride_fetch.py files PXD000001
    python pride_fetch.py files PXD000001 --ext raw
    python pride_fetch.py download PXD000001 --ext mzid --out ./pride_out
    python pride_fetch.py search "phosphoproteome" --limit 20
    python pride_fetch.py samplesheet PXD000561                 # minimal-valid <acc>.sdrf.tsv
    python pride_fetch.py samplesheet PXD000001 --from generate --acquisition dda
    python pride_fetch.py download-script PXD000561 --ext raw --out dl.sh
    python pride_fetch.py metadata-table PXD000561 --out metadata.tsv
"""
import argparse
import csv
import io
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request

API = "https://www.ebi.ac.uk/pride/ws/archive/v3"
USER_AGENT = "pride-skill/1.0"

# Minimal valid SDRF for the quantms / quantmsdiann (DIA-NN) pipeline, in order.
# https://github.com/bigbio/quantmsdiann/blob/main/docs/usage.md#minimal-valid-metadata-example
MINIMAL_SDRF_COLUMNS = [
    "source name",
    "characteristics[organism]",
    "characteristics[organism part]",
    "characteristics[disease]",
    "characteristics[biological replicate]",
    "assay name",
    "technology type",
    "comment[technical replicate]",
    "comment[data file]",
    "comment[file uri]",
    "comment[fraction identifier]",
    "comment[label]",
    "comment[instrument]",
    "comment[proteomics data acquisition method]",
    "comment[cleavage agent details]",
    "comment[modification parameters]",
    "comment[precursor mass tolerance]",
    "comment[fragment mass tolerance]",
    "factor value[condition]",
]

# Defaults for minimal columns that are not per-file (placeholders the user should
# review). Values follow the documented minimal example.
SDRF_DEFAULTS = {
    "characteristics[organism]": "not available",
    "characteristics[organism part]": "not available",
    "characteristics[disease]": "not available",
    "characteristics[biological replicate]": "1",
    "technology type": "proteomic profiling by mass spectrometry",
    "comment[technical replicate]": "1",
    "comment[fraction identifier]": "1",
    "comment[label]": "label free sample",
    "comment[instrument]": "not available",
    "comment[proteomics data acquisition method]": "data-independent acquisition",
    "comment[cleavage agent details]": "NT=Trypsin;AC=MS:1001251",
    "comment[modification parameters]": "NT=Carbamidomethyl;AC=UNIMOD:4;MT=Fixed;PP=Anywhere;TA=C",
    "comment[precursor mass tolerance]": "10 ppm",
    "comment[fragment mass tolerance]": "20 ppm",
    "factor value[condition]": "not available",
}

# Mass-spec data files accepted by quantms/quantmsdiann (+ compressed variants).
MS_FILE_EXTS = (
    ".raw", ".raw.gz",
    ".mzml", ".mzml.gz",
    ".d", ".d.zip", ".d.tar", ".d.tar.gz",
    ".wiff", ".dia",
)


def http_get(url, retries=3):
    last = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT, "Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=180) as r:
                return r.read()
        except Exception as e:  # noqa: BLE001
            last = e
            time.sleep(1.5 * (attempt + 1))
    raise SystemExit(f"GET failed for {url}: {last}")


def get_json(path, params=None):
    url = f"{API}{path}"
    if params:
        url += "?" + urllib.parse.urlencode(params)
    return json.loads(http_get(url).decode("utf-8", "replace"))


def cmd_metadata(args):
    rec = get_json(f"/projects/{args.accession}")
    if args.json:
        print(json.dumps(rec, indent=2))
        return
    print(f"accession   : {rec.get('accession','')}")
    print(f"title       : {rec.get('title','')}")
    print(f"submitted   : {rec.get('submissionDate','')}   published: {rec.get('publicationDate','')}")
    print(f"type        : {rec.get('submissionType','')}")
    print(f"organisms   : {', '.join(_names(rec.get('organisms', [])))}")
    print(f"instruments : {', '.join(_names(rec.get('instruments', [])))}")
    print(f"diseases    : {', '.join(_names(rec.get('diseases', [])))}")
    print(f"keywords    : {', '.join(rec.get('keywords', []))}")
    if rec.get("doi"):
        print(f"doi         : {rec['doi']}")
    print(f"\ndescription : {rec.get('projectDescription','')}")


def _names(items):
    out = []
    for it in items:
        if isinstance(it, dict):
            out.append(it.get("name") or it.get("value") or str(it))
        else:
            out.append(str(it))
    return out


def iter_files(accession):
    """Yield all file records for a project, following pagination."""
    page = 0
    while True:
        batch = get_json(f"/projects/{accession}/files",
                         {"pageSize": 100, "page": page})
        if not batch:
            break
        for f in batch:
            yield f
        if len(batch) < 100:
            break
        page += 1


def ftp_url(file_rec):
    for loc in file_rec.get("publicFileLocations", []):
        if isinstance(loc, dict) and "FTP" in (loc.get("name") or ""):
            return loc.get("value")
    # fall back to first available location
    locs = file_rec.get("publicFileLocations", [])
    return locs[0].get("value") if locs else None


def cmd_files(args):
    files = list(iter_files(args.accession))
    if args.ext:
        files = [f for f in files if _fname(f).lower().endswith(args.ext.lower())]
    if args.json:
        print(json.dumps(files, indent=2))
        return
    print(f"{len(files)} file(s) for {args.accession}:\n")
    for f in files:
        cat = (f.get("fileCategory") or {}).get("value", "")
        print(f"{_fname(f)}\t[{cat}]\t{ftp_url(f) or ''}")


def _fname(file_rec):
    url = ftp_url(file_rec) or ""
    return url.rstrip("/").split("/")[-1] or file_rec.get("accession", "")


def download_file(url, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    dl = url.replace("ftp://", "https://") if url.startswith("ftp://") else url
    name = dl.rstrip("/").split("/")[-1]
    dest = os.path.join(out_dir, name)
    if os.path.exists(dest):
        print(f"  exists, skipping {name}", file=sys.stderr)
        return dest
    print(f"  downloading {name} ...", file=sys.stderr)
    req = urllib.request.Request(dl, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=600) as r, open(dest + ".part", "wb") as fh:
        while True:
            chunk = r.read(1 << 20)
            if not chunk:
                break
            fh.write(chunk)
    os.replace(dest + ".part", dest)
    print(f"  saved {dest}", file=sys.stderr)
    return dest


def cmd_download(args):
    files = list(iter_files(args.accession))
    if args.ext:
        files = [f for f in files if _fname(f).lower().endswith(args.ext.lower())]
    if not files:
        raise SystemExit("No matching files.")
    out = os.path.join(args.out, args.accession)
    n = 0
    for f in files:
        url = ftp_url(f)
        if url:
            download_file(url, out)
            n += 1
    print(f"Downloaded {n} file(s) to {out}", file=sys.stderr)


def _https(url):
    return url.replace("ftp://", "https://") if url.startswith("ftp://") else url


# --- output field hygiene: control chars (incl. CR/LF/tab) -> _, so no
# value can break the line or cell it is written into. ---
_CTRL_RE = re.compile(r"[\x00-\x1f\x7f-\x9f]")  # CR, LF, tab, other control chars


def _safe_field(v):
    """Value safe for one line/cell: control chars -> _, runs of spaces collapsed."""
    return re.sub(r" +", " ", _CTRL_RE.sub("_", v or "")).strip()


# Transfer-robustness defaults are shared with the other archive skills, so a
# PRIDE script behaves like an ENA or GEO one. PRIDE keeps its own emitter
# because it is file-list driven and supports --unzip.
from clawbio.common.download_script import (  # noqa: E402
    CONNECT_TIMEOUT, CURL_PROBE, RETRIES, RETRY_DELAY, STALL_BYTES, STALL_SECONDS,
)


def _dl_cmd(tool, url, outdir):
    name = _safe_field(url.rstrip("/").split("/")[-1])
    dest = f'"{outdir}/{name}"'
    if tool == "curl":
        # -f fail on HTTP errors, -sS quiet but still report errors, -L follow
        # redirects, --retry with a delay, --connect-timeout plus the
        # --speed-limit/--speed-time pair so a stalled transfer aborts and
        # retries instead of hanging until walltime, and -C - to resume.
        # Proteomics RAW files run to tens of GB, so resume matters most here.
        return (f'curl -fsSL --retry {RETRIES} --retry-delay {RETRY_DELAY} '
                f'$RETRY_ALL --connect-timeout {CONNECT_TIMEOUT} '
                f'--speed-limit {STALL_BYTES} --speed-time {STALL_SECONDS} '
                f'-C - --create-dirs -o {dest} "{url}"')
    # wget: -q fully quiet, -c to resume. Verified against GNU Wget 1.25.0 on
    # 2026-09-22, re-verified 2026-09-23, that -c with -O resumes (206 Partial
    # Content) rather than restarting. NOT -nc, which is --no-clobber and would
    # skip an existing partial file instead of finishing it. Proteomics RAW
    # files run to tens of GB, so resume matters most here.
    return (f'wget -q --tries={RETRIES} --waitretry={RETRY_DELAY} '
            f'--timeout={CONNECT_TIMEOUT} -c -O {dest} "{url}"')


def _slurm_header(args):
    sb = [
        "#!/bin/bash",
        "#",
        f"#SBATCH --job-name={args.job_name}",
    ]
    sb.append(f"#SBATCH --partition={args.partition}" if args.partition
              else "# #SBATCH --partition=<your_partition>   # set for your cluster")
    sb.append(f"#SBATCH --account={args.account}" if args.account
              else "# #SBATCH --account=<your_account>        # set if required")
    sb += [
        f"#SBATCH --cpus-per-task={args.cpus}",
        f"#SBATCH --mem={args.mem}",
        f"#SBATCH --time={args.time}",
        "#SBATCH --output=slurm-%j.out",
        "#SBATCH --error=slurm-%j.err",
    ]
    if args.email:
        sb += ["#SBATCH --mail-type=END,FAIL", f"#SBATCH --mail-user={args.email}"]
    else:
        sb += ["# #SBATCH --mail-type=END,FAIL", "# #SBATCH --mail-user=<you@example.org>"]
    return "\n".join(sb)


def cmd_download_script(args):
    """Generate a bash + SLURM script that downloads a project's data files.

    PRIDE data files are typically vendor .raw files or .zip archives (which for
    Bruker timsTOF data contain .d directories). One quiet curl/wget command is
    emitted per file; use --unzip to also extract .zip archives after download.
    """
    files = list(iter_files(args.accession))
    if args.ext:
        files = [f for f in files if _fname(f).lower().endswith(args.ext.lower())]
    entries = []
    for f in files:
        url = ftp_url(f)
        if not url:
            continue
        url = _https(url)
        if _CTRL_RE.search(url) or '"' in url:  # never embed an unsafe URL in the shell
            print(f"warning: skipping file with an unsafe URL: {url!r}", file=sys.stderr)
            continue
        cat = (f.get("fileCategory") or {}).get("value", "")
        entries.append((_safe_field(_fname(f)), url, _safe_field(cat)))
    if not entries:
        raise SystemExit("No matching files for this project (check --ext).")

    parts = []
    if not args.no_slurm:
        parts.append(_slurm_header(args))
        parts.append("")
    else:
        parts.append("#!/bin/bash")
    parts += ["set -euo pipefail", "", f'OUTDIR="{args.outdir}"', 'mkdir -p "$OUTDIR"', ""]
    if args.tool == "curl":
        # Required: _dl_cmd references $RETRY_ALL, and `set -u` makes an
        # unbound variable a hard error.
        parts += [CURL_PROBE]
    n_zip = 0
    for name, url, cat in entries:
        parts.append(f"# {cat or 'file'}: {name}")
        parts.append(_dl_cmd(args.tool, url, "$OUTDIR"))
        if args.unzip and name.lower().endswith(".zip"):
            parts.append(f'unzip -q -o "$OUTDIR/{name}" -d "$OUTDIR"')
            n_zip += 1
        parts.append("")
    parts.append(f'echo "Downloaded {len(entries)} file(s) to $OUTDIR"')
    parts.append("")
    with open(args.out, "w") as fh:
        fh.write("\n".join(parts))
    mode = os.stat(args.out).st_mode
    os.chmod(args.out, mode | ((mode & 0o444) >> 2))  # mirror read bits into execute
    msg = (f"Wrote {args.out}: {len(entries)} download command(s) using {args.tool}")
    if args.unzip and n_zip:
        msg += f" (+{n_zip} unzip step(s))"
    print(msg + ".", file=sys.stderr)
    if not args.no_slurm:
        print(f"Submit with: sbatch {args.out}   (or run directly: ./{args.out})", file=sys.stderr)


def _norm_col(col):
    return re.sub(r"\s+", " ", col.strip().lower())


def sdrf_out_path(out):
    """Enforce the .sdrf.tsv extension (quantms rejects .sdrf/.tsv/.csv)."""
    if out.endswith(".sdrf.tsv"):
        return out
    return re.sub(r"\.(sdrf\.tsv|sdrf|tsv|csv)$", "", out, flags=re.IGNORECASE) + ".sdrf.tsv"


def _cv_term(item):
    name = (item or {}).get("name", "")
    acc = (item or {}).get("accession", "")
    if name and acc:
        return f"NT={name};AC={acc}"
    return name or "not available"


def is_ms_file(name):
    return name.lower().endswith(MS_FILE_EXTS)


def minimal_defaults(accession, acquisition):
    """Defaults for the shared minimal columns, enriched from project metadata."""
    d = dict(SDRF_DEFAULTS)
    d["comment[proteomics data acquisition method]"] = (
        "data-independent acquisition" if acquisition == "dia" else "data-dependent acquisition")
    try:
        meta = get_json(f"/projects/{accession}")
    except SystemExit:
        return d
    if meta.get("organisms"):
        name = re.sub(r"\s*\(.*\)$", "", meta["organisms"][0].get("name", "")).strip()
        if name:
            d["characteristics[organism]"] = name
    if meta.get("diseases"):
        d["characteristics[disease]"] = meta["diseases"][0].get("name") or d["characteristics[disease]"]
    if meta.get("instruments"):
        d["comment[instrument]"] = _cv_term(meta["instruments"][0])
    return d


def generate_minimal_rows(accession, files, defaults, local_dir=None):
    rows = []
    for i, f in enumerate(files, start=1):
        url = ftp_url(f)
        fname = _fname(f)
        row = dict(defaults)
        row["source name"] = f"Sample {i}"
        row["assay name"] = f"run {i}"
        row["comment[data file]"] = os.path.join(local_dir, fname) if local_dir else fname
        row["comment[file uri]"] = _https(url) if url else ""
        rows.append(row)
    return rows


def write_minimal_tsv(rows, out):
    with open(out, "w", newline="") as fh:
        # LF, not csv.writer's default CRLF: a trailing \r ends up inside the last
        # field's value when the pipeline reads the row.
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(MINIMAL_SDRF_COLUMNS)
        for r in rows:
            w.writerow([_safe_field(r.get(c, "")) for c in MINIMAL_SDRF_COLUMNS])


def complete_existing_sdrf(text, out, defaults, local_dir=None):
    """Write the submitter SDRF, appending any missing minimal columns with defaults."""
    rows = list(csv.reader(io.StringIO(text), delimiter="\t"))
    if not rows:
        raise SystemExit("Downloaded SDRF is empty.")
    header = rows[0]
    present = {_norm_col(c) for c in header}
    missing = [c for c in MINIMAL_SDRF_COLUMNS if _norm_col(c) not in present]
    df_idx = next((i for i, c in enumerate(header)
                   if _norm_col(c) == _norm_col("comment[data file]")), None)
    with open(out, "w", newline="") as fh:
        # LF, not csv.writer's default CRLF (see write_minimal_tsv)
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow([_safe_field(c) for c in header + missing])
        extra = [defaults.get(c, "not available") for c in missing]
        for row in rows[1:]:
            if not any(c.strip() for c in row):
                continue
            if local_dir and df_idx is not None and df_idx < len(row) and row[df_idx].strip():
                row = list(row)
                row[df_idx] = os.path.join(local_dir, os.path.basename(row[df_idx].strip()))
            w.writerow([_safe_field(c) for c in row + extra])
    return missing


def validate_minimal(path):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
    present = {_norm_col(c) for c in header}
    return [c for c in MINIMAL_SDRF_COLUMNS if _norm_col(c) not in present]


def cmd_samplesheet(args):
    """Write a minimal-valid SDRF sample sheet for the project.

    The output conforms to the quantms/quantmsdiann minimal SDRF definition
    (19 required columns, tab-delimited, .sdrf.tsv extension). If the project has a
    submitter SDRF it is used and any missing minimal columns are appended with
    defaults; otherwise a minimal SDRF is generated from the data files + metadata.
    """
    out = sdrf_out_path(args.out or f"{args.accession}.sdrf.tsv")
    if args.out and out != args.out:
        print(f"note: output renamed to {out} (quantms requires the .sdrf.tsv extension)",
              file=sys.stderr)

    urls = get_json(f"/files/sdrf/{args.accession}")
    if args.source == "pride" and not urls:
        raise SystemExit(
            f"No submitter SDRF in PRIDE for {args.accession}. "
            "Use --from generate to build a minimal SDRF from the data files.")
    defaults = minimal_defaults(args.accession, args.acquisition)

    if urls and args.source in ("auto", "pride"):
        text = http_get(_https(urls[0])).decode("utf-8", "replace")
        missing = complete_existing_sdrf(text, out, defaults, args.local_dir)
        print(f"Wrote PRIDE submitter SDRF -> {out}", file=sys.stderr)
        if missing:
            print(f"Completed {len(missing)} missing minimal column(s) with defaults: "
                  f"{', '.join(missing)}", file=sys.stderr)
        if len(urls) > 1:
            print(f"note: {len(urls)} SDRF files exist for this project; used the first.",
                  file=sys.stderr)
    else:
        if args.source == "auto":
            print("No submitter SDRF in PRIDE; generating a minimal SDRF from "
                  "data files + metadata.", file=sys.stderr)
        files = [f for f in iter_files(args.accession) if is_ms_file(_fname(f))]
        if not files:
            raise SystemExit(
                "No MS data files (.raw/.mzML/.d/.wiff) found; cannot generate a minimal SDRF.")
        rows = generate_minimal_rows(args.accession, files, defaults, args.local_dir)
        write_minimal_tsv(rows, out)
        print(f"Generated minimal SDRF for {len(rows)} run(s) -> {out}", file=sys.stderr)

    still_missing = validate_minimal(out)
    n_ok = len(MINIMAL_SDRF_COLUMNS) - len(still_missing)
    print(f"Minimal SDRF columns present: {n_ok}/{len(MINIMAL_SDRF_COLUMNS)}.", file=sys.stderr)
    if still_missing:
        print(f"  WARNING still missing: {', '.join(still_missing)}", file=sys.stderr)
    print("Review placeholder values (acquisition method, instrument, tolerances, "
          "enzyme, modifications, organism part, factor value) before running quantms.",
          file=sys.stderr)


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
        # LF, not csv.writer's default CRLF, so the TSV is byte-identical
        # across platforms.
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
    """Write a harmonized metadata.tsv (one row per sample x replicate).

    Uses the submitter SDRF when present (characteristics[...] / factor value[...]
    columns, `source name` as sample, technical/biological replicate as replicate);
    otherwise emits a single project-level row from PRIDE metadata.
    """
    try:
        sdrf_urls = get_json(f"/files/sdrf/{args.accession}")
    except SystemExit:
        sdrf_urls = []
    rows = []
    if sdrf_urls:
        text = http_get(_https(sdrf_urls[0])).decode("utf-8", "replace")
        table = list(csv.reader(io.StringIO(text), delimiter="\t"))
        if not table:
            raise SystemExit("Downloaded SDRF is empty.")
        header = [h.strip() for h in table[0]]
        for raw in table[1:]:
            if not any(c.strip() for c in raw):
                continue
            cells = {header[i]: (raw[i] if i < len(raw) else "") for i in range(len(header))}
            sample = next((v for k, v in cells.items() if k.strip().lower() == "source name"), "")
            replicate = ""
            for k, v in cells.items():
                if _norm_key(k) in _REPLICATE_KEYS and v.strip():
                    replicate = v.strip()
                    break
            attrs = {k: v for k, v in cells.items() if _is_attr_col(k)}
            attrs = merge_biosample(sample, attrs)
            rows.append(harmonize_row(sample, replicate or "1", attrs))
        src = "submitter SDRF"
    else:
        meta = get_json(f"/projects/{args.accession}")
        attrs = {}
        if meta.get("organisms"):
            attrs["organism"] = re.sub(r"\s*\(.*\)$", "",
                                       meta["organisms"][0].get("name", "")).strip()
        if meta.get("diseases"):
            attrs["disease"] = meta["diseases"][0].get("name", "")
        if meta.get("instruments"):
            attrs["instrument"] = ", ".join(_names(meta["instruments"]))
        rows.append(harmonize_row(args.accession, "1", attrs))
        src = "project metadata (no SDRF)"
    n, ncols = write_metadata_tsv(rows, args.out)
    print(f"Wrote {n} row(s) x {ncols} column(s) for {args.accession} from {src} to {args.out}",
          file=sys.stderr)


def cmd_search(args):
    params = {"keyword": args.query, "pageSize": args.limit, "page": 0}
    res = get_json("/search/projects", params)
    # response may be a list or a paged object; normalise
    items = res if isinstance(res, list) else res.get("_embedded", {}).get("projects", res.get("content", []))
    if args.json:
        print(json.dumps(items, indent=2))
        return
    print(f"Showing up to {args.limit} results for '{args.query}':\n")
    for r in items:
        print(f"{r.get('accession',''):12s} {r.get('title','')}")


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    m = sub.add_parser("metadata", help="Project metadata for a PXD accession")
    m.add_argument("accession")
    m.add_argument("--json", action="store_true")
    m.set_defaults(func=cmd_metadata)

    f = sub.add_parser("files", help="List files in a project")
    f.add_argument("accession")
    f.add_argument("--ext", help="filter by extension, e.g. raw, mzid, mzML, mgf")
    f.add_argument("--json", action="store_true")
    f.set_defaults(func=cmd_files)

    d = sub.add_parser("download", help="Download project files")
    d.add_argument("accession")
    d.add_argument("--ext", help="only download files with this extension")
    d.add_argument("--out", default="./pride_out")
    d.set_defaults(func=cmd_download)

    q = sub.add_parser("search", help="Keyword search of PRIDE projects")
    q.add_argument("query")
    q.add_argument("--limit", type=int, default=20)
    q.add_argument("--json", action="store_true")
    q.set_defaults(func=cmd_search)

    ss = sub.add_parser("samplesheet",
                        help="Write a minimal-valid SDRF sample sheet (.sdrf.tsv) for quantms/quantmsdiann")
    ss.add_argument("accession")
    ss.add_argument("--out", help="output path (default <accession>.sdrf.tsv; extension is enforced)")
    ss.add_argument("--from", dest="source", choices=["auto", "pride", "generate"], default="auto",
                    help="auto: use submitter SDRF if present else generate; "
                         "pride: submitter SDRF only; generate: build from data files + metadata")
    ss.add_argument("--acquisition", choices=["dia", "dda"], default="dia",
                    help="value for comment[proteomics data acquisition method] (default dia)")
    ss.add_argument("--local-dir",
                    help="write comment[data file] as local paths <dir>/<file> "
                         "(matching download-script output) instead of bare filenames")
    ss.set_defaults(func=cmd_samplesheet)

    mt = sub.add_parser("metadata-table",
                        help="Write a harmonized metadata.tsv (one row per sample x replicate)")
    mt.add_argument("accession")
    mt.add_argument("--out", default="metadata.tsv")
    mt.set_defaults(func=cmd_metadata_table)

    ds = sub.add_parser("download-script",
                        help="Generate a bash + SLURM script to download project files (.raw/.zip)")
    ds.add_argument("accession")
    ds.add_argument("--ext", help="only include files with this extension, e.g. raw, zip, d.zip")
    ds.add_argument("--tool", choices=["curl", "wget"], default="curl")
    ds.add_argument("--outdir", default="pride_data", help="download destination directory")
    ds.add_argument("--out", default="download_pride.sh", help="generated script path")
    ds.add_argument("--unzip", action="store_true",
                    help="append `unzip` steps for .zip archives (e.g. Bruker .d directories)")
    ds.add_argument("--no-slurm", action="store_true", help="omit the SLURM header (plain bash)")
    ds.add_argument("--job-name", default="pride_download")
    ds.add_argument("--partition")
    ds.add_argument("--account")
    ds.add_argument("--cpus", default="1")
    ds.add_argument("--mem", default="4G")
    ds.add_argument("--time", default="24:00:00")
    ds.add_argument("--email")
    ds.set_defaults(func=cmd_download_script)

    args = p.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
