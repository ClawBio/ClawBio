#!/usr/bin/env python3
# Vendored from UKDRI/informatics_data_skills @ 7cc3e6e (biostudies/scripts/biostudies.py).
# Copyright (c) 2026 UK Dementia Research Institute. Licensed MIT.
# Kept close to upstream so changes there remain easy to merge; the ClawBio
# CLI contract lives in biostudies_fetch.py, not here.
"""Query metadata and download data from EMBL-EBI BioStudies.

Uses the BioStudies REST API (https://www.ebi.ac.uk/biostudies/api/v1). Files are
downloaded from https://www.ebi.ac.uk/biostudies/files/{accession}/{path}.

Standard library only (urllib) -- no pip install required.

Examples:
    python biostudies_fetch.py --command metadata S-BSST123
    python biostudies_fetch.py files    S-BSST123
    python biostudies_fetch.py download S-BSST123 --out ./bs_out
    python biostudies_fetch.py search "spatial transcriptomics" --limit 20
    python biostudies_fetch.py search "cancer" --collection BioImages --limit 10
    python biostudies_fetch.py metadata-table S-BSST123 --out metadata.tsv
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
USER_AGENT = "biostudies-skill/1.0"


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


def cmd_metadata(args):
    d = study_json(args.accession)
    if args.json:
        print(json.dumps(d, indent=2))
        return
    print(f"accession : {d.get('accno','')}")
    top = attrs_to_dict(d.get("attributes"))
    sec = d.get("section", {})
    secattr = attrs_to_dict(sec.get("attributes"))
    for key in ("Title", "ReleaseDate", "AttachTo"):
        if key in top:
            print(f"{key:10s}: {top[key]}")
    for key in ("Title", "Description", "Organism"):
        if key in secattr:
            print(f"{key:10s}: {secattr[key]}")
    print(f"\nfiles     : {count_files(sec)} file entries "
          f"(run `files {args.accession}` to list)")


def walk_files(node):
    """Recursively yield file dicts (type == 'file') from a study section tree."""
    if isinstance(node, dict):
        if node.get("type") == "file" and "path" in node:
            yield node
        for v in node.values():
            yield from walk_files(v)
    elif isinstance(node, list):
        for v in node:
            yield from walk_files(v)


def count_files(section):
    return sum(1 for _ in walk_files(section))


def list_files(accession):
    d = study_json(accession)
    return list(walk_files(d.get("section", {})))


def cmd_files(args):
    files = list_files(args.accession)
    if args.json:
        print(json.dumps(files, indent=2))
        return
    print(f"{len(files)} file(s) for {args.accession}:\n")
    for f in files:
        size = f.get("size", "?")
        desc = attrs_to_dict(f.get("attributes")).get("Description", "")
        print(f"{f['path']}\t{size}\t{desc}")


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

    DIVERGES FROM UPSTREAM. Upstream hardcoded
    `https://www.ebi.ac.uk/biostudies/files/{accession}/{path}`, which now
    returns 404 for every study (checked 2026-09-21 against S-BSST2074,
    E-MTAB-2037, E-MTAB-13429 and E-MTAB-11467). The live API advertises the
    real base as `httpLink` on /studies/{accession}/info, with files under
    `Files/`, and the BioStudies web UI uses that same link. Resolve it from
    there and keep the legacy path only as a fallback.
    """
    quoted = urllib.parse.quote(path)
    base = (study_info(accession) or {}).get("httpLink")
    if base:
        return f"{base.rstrip('/')}/Files/{quoted}"
    return f"{FILES}/{accession}/{quoted}"


def download_file(accession, path, out_dir):
    url = file_url(accession, path)
    dest = os.path.join(out_dir, path)
    os.makedirs(os.path.dirname(dest) or ".", exist_ok=True)
    if os.path.exists(dest):
        print(f"  exists, skipping {path}", file=sys.stderr)
        return dest
    print(f"  downloading {path} ...", file=sys.stderr)
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
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
    files = list_files(args.accession)
    if args.match:
        files = [f for f in files if args.match.lower() in f["path"].lower()]
    if not files:
        raise SystemExit("No matching files.")
    out = os.path.join(args.out, args.accession)
    for f in files:
        download_file(args.accession, f["path"], out)
    print(f"Downloaded {len(files)} file(s) to {out}", file=sys.stderr)


# ---- harmonized metadata.tsv (shared shape across skills) ----
import csv
import html
import re

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


# Clean output fields: control chars (incl. CR/LF/tab) -> _, so that no cell
# can break the TSV it is written into.
_CTRL_RE = re.compile(r"[\x00-\x1f\x7f-\x9f]")


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
        # Clean output fields: LF, not csv.writer's default CRLF, so the TSV is
        # byte-identical across platforms.
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        for r in rows:
            w.writerow([r.get(c, NA) for c in header])
    return len(rows), len(header)


def _first_attr(attrs, candidates):
    for k, v in attrs.items():
        if _norm_key(k) in candidates and _clean_val(v):
            return v
    return ""


def collect_sample_nodes(section):
    """Best-effort: find sample-like subsection nodes in the PageTab tree.

    A node counts as a sample if its `type` mentions "sample", or its attributes
    include an organism field. The root `section` is excluded (its organism is the
    study-level organism, not a sample); results are deduped by object identity.
    """
    found, seen = [], set()

    def walk(node, is_root):
        if isinstance(node, dict):
            typ = str(node.get("type", "")).lower()
            attrs = attrs_to_dict(node.get("attributes"))
            norm = {_norm_key(k) for k in attrs}
            if not is_root and id(node) not in seen and ("sample" in typ or "organism" in norm):
                seen.add(id(node))
                found.append(node)
            for v in node.values():
                walk(v, False)
        elif isinstance(node, list):
            for v in node:
                walk(v, False)

    walk(section, True)
    return found


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
    """Write a harmonized metadata.tsv (one row per sample subsection).

    BioStudies studies are heterogeneous: this locates sample-like subsections and
    maps their attributes to the common schema. When no per-sample structure is
    present, a single study-level row is emitted from the section attributes.
    """
    d = study_json(args.accession)
    sec = d.get("section", {})
    rows = []
    for node in collect_sample_nodes(sec):
        attrs = attrs_to_dict(node.get("attributes"))
        sample = node.get("accno") or _first_attr(attrs, ("name", "sample", "source name", "title"))
        replicate = _first_attr(attrs, set(_REPLICATE_KEYS))
        attrs = merge_biosample(sample, attrs)
        rows.append(harmonize_row(sample, replicate or "1", attrs))
    if not rows:
        rows.append(harmonize_row(d.get("accno", ""), "1", attrs_to_dict(sec.get("attributes"))))
        note = " (no per-sample structure; used study-level attributes)"
    else:
        note = ""
    n, ncols = write_metadata_tsv(rows, args.out)
    print(f"Wrote {n} row(s) x {ncols} column(s) for {args.accession} to {args.out}{note}",
          file=sys.stderr)


def cmd_search(args):
    params = {"query": args.query, "pageSize": args.limit, "page": 1}
    if args.collection:
        params["collection"] = args.collection
    url = f"{API}/search?" + urllib.parse.urlencode(params)
    res = get_json(url)
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

    m = sub.add_parser("metadata", help="Study metadata for an accession")
    m.add_argument("accession")
    m.add_argument("--json", action="store_true")
    m.set_defaults(func=cmd_metadata)

    f = sub.add_parser("files", help="List files attached to a study")
    f.add_argument("accession")
    f.add_argument("--json", action="store_true")
    f.set_defaults(func=cmd_files)

    d = sub.add_parser("download", help="Download study files")
    d.add_argument("accession")
    d.add_argument("--match", help="only files whose path contains this substring")
    d.add_argument("--out", default="./bs_out")
    d.set_defaults(func=cmd_download)

    mt = sub.add_parser("metadata-table",
                        help="Write a harmonized metadata.tsv (one row per sample subsection)")
    mt.add_argument("accession")
    mt.add_argument("--out", default="metadata.tsv")
    mt.set_defaults(func=cmd_metadata_table)

    q = sub.add_parser("search", help="Search BioStudies")
    q.add_argument("query")
    q.add_argument("--collection", help="restrict to a collection, e.g. BioImages, ArrayExpress")
    q.add_argument("--limit", type=int, default=20)
    q.add_argument("--json", action="store_true")
    q.set_defaults(func=cmd_search)

    args = p.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
