#!/usr/bin/env python3
"""
article-data-fetcher — ClawBio skill
Downloads genomics data files (VCF, FASTA, H5AD, CSV, BAM, etc.) deposited
by authors in public repositories (GEO, ENA, Zenodo, Figshare, Dryad, OSF)
for a given article DOI or PubMed ID.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional
from urllib.parse import urlparse

import requests
from tqdm import tqdm

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

SUPPORTED_EXTENSIONS = {
    "vcf", "vcf.gz", "bcf",
    "fasta", "fa", "fna", "fastq", "fastq.gz",
    "bam", "bai", "cram",
    "h5ad", "h5", "loom",
    "csv", "tsv", "txt", "xlsx",
    "json", "yaml",
    "bed", "gff", "gtf",
    "gz", "zip", "tar.gz",
    "mtx", "mtx.gz",
}

DEMO_ACCESSION = "GSE100866"  # CITE-seq CBMC dataset — has csv.gz supplementary files

ENTREZ_BASE = "https://edata.ncbi.nlm.nih.gov/geo/query/acc.cgi"
CROSSREF_BASE = "https://api.crossref.org/works"
ENA_PORTAL = "https://www.ebi.ac.uk/ena/portal/api/filereport"
ZENODO_BASE = "https://zenodo.org/api/records"
FIGSHARE_BASE = "https://api.figshare.com/v2/articles"

ACCESSION_PATTERNS = {
    "geo": re.compile(r"\bGSE\d+\b", re.IGNORECASE),
    "sra": re.compile(r"\b(PRJNA|SRP|ERP|DRP)\d+\b", re.IGNORECASE),
    "arrayexpress": re.compile(r"\bE-[A-Z]+-\d+\b", re.IGNORECASE),
    "zenodo": re.compile(r"10\.5281/zenodo\.(\d+)|zenodo\.org/records?/(\d+)", re.IGNORECASE),
    "zenodo_community": re.compile(r"zenodo\.org/communities/([\w-]+)", re.IGNORECASE),
    "figshare": re.compile(r"10\.6084/m9\.figshare\.(\d+)|figshare\.com/articles/[^/]+/(\d+)", re.IGNORECASE),
    "dryad": re.compile(r"10\.5061/dryad\.[\w]+", re.IGNORECASE),
    "osf": re.compile(r"osf\.io/([\w]{5})", re.IGNORECASE),
}

SIZE_WARNING_BYTES = 10 * 1024 ** 3  # 10 GB


# ---------------------------------------------------------------------------
# Resolution helpers
# ---------------------------------------------------------------------------

def resolve_doi(doi: str) -> dict:
    """Fetch article metadata from Crossref."""
    url = f"{CROSSREF_BASE}/{doi}"
    r = requests.get(url, timeout=15, headers={"User-Agent": "ClawBio/1.0"})
    r.raise_for_status()
    data = r.json()["message"]
    return {
        "title": (data.get("title") or ["Unknown"])[0],
        "doi": doi,
        "abstract": data.get("abstract", ""),
        "links": [l.get("URL", "") for l in data.get("link", [])],
    }


def resolve_pmid(pmid: str) -> dict:
    """Fetch article metadata from NCBI Entrez, including DOI if available."""
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
        f"?db=pubmed&id={pmid}&retmode=json"
    )
    r = requests.get(url, timeout=15)
    r.raise_for_status()
    result = r.json()["result"][pmid]
    # Extract DOI from articleids list
    doi = None
    for aid in result.get("articleids", []):
        if aid.get("idtype") == "doi":
            doi = aid.get("value")
            break
    return {
        "title": result.get("title", "Unknown"),
        "pmid": pmid,
        "doi": doi,
        "abstract": "",
        "links": [],
    }


def resolve_pmcid(pmcid: str) -> dict:
    """
    Fetch article metadata + full-text from PubMed Central.
    Fetches the HTML page directly — more reliable than efetch XML for newer articles
    and captures all href links including Zenodo communities, GitHub, OSF, etc.
    """
    numeric_id = re.sub(r"[^\d]", "", pmcid)

    # 1. Get article title from Entrez summary
    title = "Unknown"
    try:
        summary_url = (
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
            f"?db=pmc&id={numeric_id}&retmode=json"
        )
        r = requests.get(summary_url, timeout=15)
        if r.status_code == 200:
            result = r.json().get("result", {}).get(numeric_id, {})
            title = result.get("title", "Unknown")
    except requests.RequestException:
        pass

    # 2. Fetch the HTML page — captures all external links including community URLs
    page_text = ""
    links: list[str] = []
    doi = None
    try:
        page_url = f"https://pmc.ncbi.nlm.nih.gov/articles/{pmcid}/"
        r = requests.get(page_url, timeout=20, headers={"User-Agent": "ClawBio/1.0"})
        if r.status_code == 200:
            page_text = r.text
            # Extract all href values so patterns can match full URLs
            links = re.findall(r'href="([^"]+)"', page_text)
            # Extract DOI from <meta name="citation_doi"> or og:url
            doi_meta = re.search(r'<meta[^>]+name="citation_doi"[^>]+content="([^"]+)"', page_text)
            if doi_meta:
                doi = doi_meta.group(1).strip()
            if not doi:
                doi_meta = re.search(r'<meta[^>]+content="(10\.[^"]+)"[^>]+name="citation_doi"', page_text)
                if doi_meta:
                    doi = doi_meta.group(1).strip()
            # Also try to find a bare DOI in the page text as last resort
            if not doi:
                m = re.search(r'\b(10\.\d{4,}/\S+)', page_text)
                if m:
                    doi = m.group(1).rstrip('.,;)')
    except requests.RequestException:
        pass

    return {
        "title": title,
        "pmcid": pmcid,
        "doi": doi,
        "abstract": page_text + " " + " ".join(links),
        "links": links,
    }


def parse_identifier(identifier: str) -> tuple[str, str]:
    """Return (kind, value) — kind is 'doi', 'pmid', 'pmcid', or 'url'."""
    identifier = identifier.strip()
    # Unwrap doi.org resolver URLs → bare DOI
    doi_url_match = re.match(r"https?://(?:dx\.)?doi\.org/(10\..+)", identifier)
    if doi_url_match:
        return ("doi", doi_url_match.group(1))
    if identifier.startswith("10."):
        return ("doi", identifier)
    if re.match(r"^PMC\d+$", identifier, re.IGNORECASE):
        return ("pmcid", identifier.upper())
    if re.match(r"^(PMID:?)?\d{7,9}$", identifier, re.IGNORECASE):
        pmid = re.sub(r"[^\d]", "", identifier)
        return ("pmid", pmid)
    if identifier.startswith("http"):
        # Detect PMCIDs embedded in URLs e.g. /pmc/articles/PMC1234567/
        m = re.search(r"PMC\d+", identifier, re.IGNORECASE)
        if m:
            return ("pmcid", m.group(0).upper())
        return ("url", identifier)
    raise ValueError(
        f"Cannot recognise identifier '{identifier}'. "
        "Provide a DOI (10.xxxx/… or https://doi.org/…), PMID, PMCID (e.g. PMC1234567), or repository URL."
    )


def _wrap_filename(fname: str, width: int) -> list[str]:
    """Wrap a filename into lines of at most `width` chars, breaking at - or . where possible."""
    if len(fname) <= width:
        return [fname]
    lines = []
    while len(fname) > width:
        # Find the rightmost natural break (- or .) before the width limit
        cut = None
        for sep in ("-", "."):
            pos = fname.rfind(sep, 0, width)
            if pos > width // 2:  # only accept if not too far left
                if cut is None or pos > cut:
                    cut = pos + 1  # break after the separator
        if cut is None:
            cut = width  # no good break found — hard cut
        lines.append(fname[:cut])
        fname = fname[cut:]
    if fname:
        lines.append(fname)
    return lines


def print_file_listing(files: list[dict], title: str = "") -> None:
    """Print a single-line-per-file table; long filenames wrap onto continuation lines."""
    if title:
        print(title)
    col_fname = 42
    col_size  = 10
    col_type  = 28
    divider = f"  {'─'*4}  {'─'*col_fname}  {'─'*col_size}  {'─'*col_type}  {'─'*6}"
    header  = f"  {'#':>4}  {'FILENAME':<{col_fname}}  {'SIZE':>{col_size}}  {'DATA TYPE':<{col_type}}  SOURCE"
    print(f"\n{header}")
    print(divider)
    for i, f in enumerate(files, 1):
        size_str = _fmt_size(f.get("size_bytes"))
        dtype    = f.get("data_type", "data file")
        repo     = f['repository']
        chunks   = _wrap_filename(f['filename'], col_fname)
        # First line: number + first chunk + size + type + source
        print(f"  [{i:>2}]  {chunks[0]:<{col_fname}}  {size_str:>{col_size}}  {dtype:<{col_type}}  {repo}")
        # Continuation lines: remaining filename chunks, aligned under filename column
        for chunk in chunks[1:]:
            print(f"        {chunk}")
    print()


def prompt_and_select_files(files: list[dict]) -> list[dict]:
    """
    Interactively ask the user which files to download.
    Accepts:
      - Numbers or ranges:  1,3,5-7
      - File extensions:    vcf,fasta,h5ad
      - 'all'
    Returns the selected subset of files.
    """
    print(
        "Which files would you like to download?\n"
        "  • By number / range  →  1,3,5-7\n"
        "  • By extension       →  vcf,fasta,h5ad\n"
        "  • Everything         →  all\n"
        "  • Cancel             →  none\n"
    )
    answer = input("→ ").strip().lower()
    if not answer or answer == "none":
        print("Download cancelled.")
        return []

    if answer == "all":
        return files

    selected: list[dict] = []

    # Split on commas, then decide whether each token is a number/range or an extension
    tokens = [t.strip() for t in answer.split(",") if t.strip()]
    number_tokens: list[str] = []
    ext_tokens: list[str] = []

    for token in tokens:
        # A number range looks like "3" or "5-7"
        if re.match(r"^\d+(-\d+)?$", token):
            number_tokens.append(token)
        else:
            ext_tokens.append(token)

    # Resolve number/range tokens → file indices (1-based)
    indices: set[int] = set()
    for token in number_tokens:
        if "-" in token:
            start, end = token.split("-", 1)
            indices.update(range(int(start), int(end) + 1))
        else:
            indices.add(int(token))

    for idx in sorted(indices):
        if 1 <= idx <= len(files):
            selected.append(files[idx - 1])
        else:
            print(f"  ⚠️  No file at position {idx} — skipped")

    # Resolve extension tokens
    if ext_tokens:
        ext_set = set(ext_tokens)
        for f in files:
            if file_matches_types(f["filename"], ext_set) and f not in selected:
                selected.append(f)

    return selected


# ---------------------------------------------------------------------------
# Accession discovery
# ---------------------------------------------------------------------------
# Data type inference
# ---------------------------------------------------------------------------

# (keyword fragments → data type label)
_DATA_TYPE_RULES: list[tuple[list[str], str]] = [
    (["count", "umi", "raw", "matrix", "expression", "mtx"], "count matrix"),
    (["norm", "normalised", "normalized", "tpm", "fpkm", "rpkm", "cpm"], "normalised expression"),
    (["clr", "clr-transformed", "adt"], "protein / ADT expression"),
    (["metadata", "meta", "sample_info", "sampleinfo", "barcode", "cell_info"], "sample / cell metadata"),
    (["cluster", "annotation", "label", "cell_type", "celltype"], "cluster / cell-type annotation"),
    (["embedding", "umap", "tsne", "pca", "latent"], "embedding / dimensionality reduction"),
    (["variant", "vcf", "snp", "indel", "mutation", "somatic", "germline"], "variant calls"),
    (["fastq", "fasta", "fa", "fna", "seq"], "raw sequences"),
    (["bam", "cram", "sam", "alignment", "mapped"], "read alignments"),
    (["bed", "peak", "atac", "chip", "cut&run"], "genomic intervals / peaks"),
    (["gff", "gtf", "annotation", "genome"], "genome annotation"),
    (["h5ad", "h5", "anndata", "seurat", "rds", "loom"], "single-cell object"),
    (["deg", "differential", "de_", "lfc", "fold_change", "deseq", "edger"], "differential expression"),
    (["supplement", "readme", "readme", "description"], "supplementary / readme"),
    (["network", "ppi", "interaction", "edge", "node"], "interaction network"),
    (["pathway", "gsea", "enrichment", "go_", "kegg"], "pathway / enrichment"),
    (["phenotype", "clinical", "cohort", "patient", "subject"], "clinical / phenotype data"),
]

def infer_data_type(filename: str, description: str = "") -> str:
    """Guess what kind of data a file contains from its name and optional description."""
    text = (filename + " " + description).lower()
    for keywords, label in _DATA_TYPE_RULES:
        if any(k in text for k in keywords):
            return label
    return "data file"


# ---------------------------------------------------------------------------

def discover_accessions(text: str) -> dict[str, list[str]]:
    """Scan free text for repository accession numbers."""
    found: dict[str, list[str]] = {}
    for repo, pattern in ACCESSION_PATTERNS.items():
        matches = list({m.group(0) for m in pattern.finditer(text)})
        if matches:
            found[repo] = matches
    return found


def fetch_geo_files(accession: str) -> list[dict]:
    """Return file listing for a GEO series, including sizes and data type tags."""
    ftp_base = (
        f"https://ftp.ncbi.nlm.nih.gov/geo/series/"
        f"{accession[:-3]}nnn/{accession}/suppl/"
    )
    # Also fetch GEO series page for file descriptions
    geo_descriptions: dict[str, str] = {}
    try:
        series_url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}&targ=self&form=text&view=quick"
        sr = requests.get(series_url, timeout=10)
        if sr.status_code == 200:
            # GEO SOFT text has lines like: !Series_supplementary_file = ftp://...filename.csv.gz
            for line in sr.text.splitlines():
                if "supplementary_file" in line.lower():
                    fname = line.split("/")[-1].strip()
                    # Look for a description line nearby (not always present, best-effort)
                    geo_descriptions[fname] = line
    except requests.RequestException:
        pass

    try:
        r = requests.get(ftp_base, timeout=10)
        if r.status_code == 200:
            # GEO FTP Apache listing format:
            # <a href="filename.csv.gz">filename.csv.gz</a>  2017-01-01 00:00  945K
            rows = re.findall(
                r'href="([^"/][^"]*)">[^<]+</a>\s+[\d\-]+ [\d:]+\s+([\d.]+[KMGT]?)\s',
                r.text,
            )
            files = []
            for fname, size_str in rows:
                if "?" in fname or fname.startswith("http"):
                    continue
                size_bytes = _parse_size(size_str)
                desc = geo_descriptions.get(fname, "")
                files.append({
                    "filename": fname,
                    "url": ftp_base + fname,
                    "repository": "GEO",
                    "accession": accession,
                    "size_bytes": size_bytes,
                    "data_type": infer_data_type(fname, desc),
                })
            if files:
                return files
            # Fallback: filenames only, no sizes
            filenames = re.findall(r'href="([^"/][^"]*\.\w+)"', r.text)
            return [
                {
                    "filename": f,
                    "url": ftp_base + f,
                    "repository": "GEO",
                    "accession": accession,
                    "size_bytes": None,
                    "data_type": infer_data_type(f),
                }
                for f in filenames
                if not f.startswith("?") and not f.startswith("http")
            ]
    except requests.RequestException:
        pass
    return []


def _parse_size(size_str: str) -> Optional[int]:
    """Convert Apache size strings like '1.2M', '340K', '3.1G' to bytes."""
    size_str = size_str.strip()
    if not size_str or size_str == "-":
        return None
    multipliers = {"K": 1024, "M": 1024**2, "G": 1024**3, "T": 1024**4}
    suffix = size_str[-1].upper()
    if suffix in multipliers:
        try:
            return int(float(size_str[:-1]) * multipliers[suffix])
        except ValueError:
            return None
    try:
        return int(size_str)
    except ValueError:
        return None


def _fmt_size(size_bytes: Optional[int]) -> str:
    """Format a byte count as KB / MB / GB automatically."""
    if size_bytes is None:
        return "unknown"
    if size_bytes >= 1024**3:
        return f"{size_bytes / 1024**3:.1f} GB"
    if size_bytes >= 1024**2:
        return f"{size_bytes / 1024**2:.1f} MB"
    if size_bytes >= 1024:
        return f"{size_bytes / 1024:.1f} KB"
    return f"{size_bytes} B"


def fetch_ena_files(accession: str) -> list[dict]:
    """Return file listing for an ENA/SRA project."""
    params = {
        "accession": accession,
        "result": "read_run",
        "fields": "run_accession,fastq_ftp,submitted_ftp,sra_ftp",
        "format": "json",
        "limit": 100,
    }
    try:
        r = requests.get(ENA_PORTAL, params=params, timeout=15)
        if r.status_code == 200:
            files = []
            for run in r.json():
                for field in ("fastq_ftp", "submitted_ftp"):
                    for url in (run.get(field) or "").split(";"):
                        url = url.strip()
                        if url:
                            fname = url.split("/")[-1]
                            files.append({
                                "filename": fname,
                                "url": "https://" + url,
                                "repository": "ENA",
                                "accession": accession,
                                "size_bytes": None,
                                "data_type": infer_data_type(fname),
                            })
            return files
    except requests.RequestException:
        pass
    return []


def fetch_zenodo_files(record_id: str) -> list[dict]:
    """Return file listing for a Zenodo record."""
    try:
        r = requests.get(f"{ZENODO_BASE}/{record_id}", timeout=10)
        if r.status_code == 200:
            data = r.json()
            description = data.get("metadata", {}).get("description", "")
            return [
                {
                    "filename": f["key"],
                    "url": f["links"]["self"],
                    "repository": "Zenodo",
                    "accession": f"zenodo.{record_id}",
                    "size_bytes": f.get("size"),
                    "checksum": f.get("checksum", "").replace("md5:", ""),
                    "data_type": infer_data_type(f["key"], description),
                }
                for f in data.get("files", [])
            ]
    except requests.RequestException:
        pass
    return []


def fetch_zenodo_community_files(community_acc: str) -> list[dict]:
    """Return file listings for all records in a Zenodo community."""
    # community_acc may be the full match e.g. "zenodo.org/communities/lungcancer"
    slug = community_acc.split("/")[-1]
    files: list[dict] = []
    try:
        url = f"https://zenodo.org/api/communities/{slug}/records"
        r = requests.get(url, params={"size": 25}, timeout=15)
        if r.status_code == 200:
            hits = r.json().get("hits", {}).get("hits", [])
            print(f"  📦 Zenodo community '{slug}' — {len(hits)} deposit(s) containing the following files:")
            for hit in hits:
                record_id = str(hit.get("id", ""))
                title = hit.get("metadata", {}).get("title", record_id)
                description = hit.get("metadata", {}).get("description", "")
                n_files = len(hit.get("files", []))
                print(f"     • {title[:70]}  ({n_files} file{'s' if n_files != 1 else ''})")
                for f in hit.get("files", []):
                    files.append({
                        "filename": f.get("key", f.get("filename", "unknown")),
                        "url": f.get("links", {}).get("self", ""),
                        "repository": f"Zenodo:{slug}",
                        "accession": f"zenodo.{record_id}",
                        "size_bytes": f.get("size"),
                        "checksum": f.get("checksum", "").replace("md5:", ""),
                        "data_type": infer_data_type(f.get("key", ""), description),
                    })
    except requests.RequestException:
        pass
    return files


def search_zenodo_by_article(doi: Optional[str], title: str) -> list[dict]:
    """Search Zenodo for deposits referencing the article DOI or title."""
    files: list[dict] = []
    queries: list[str] = []
    if doi:
        queries.append(f'related_identifier.identifier:"{doi}"')
    if title and title != "Unknown":
        short_title = re.sub(r'[<>"\']', '', title[:60])
        queries.append(f'"{short_title}"')
    for query in queries:
        try:
            print(f"  🔎 Zenodo: {query}")
            r = requests.get(
                "https://zenodo.org/api/records",
                params={"q": query, "size": 10, "sort": "mostrecent"},
                timeout=15, headers={"User-Agent": "ClawBio/1.0"},
            )
            if r.status_code != 200:
                continue
            hits = r.json().get("hits", {}).get("hits", [])
            if not hits:
                continue
            print(f"     Found {len(hits)} record(s) on Zenodo")
            for hit in hits:
                record_id = str(hit.get("id", ""))
                description = hit.get("metadata", {}).get("description", "")
                rec_title = hit.get("metadata", {}).get("title", record_id)
                n = len(hit.get("files", []))
                print(f"     • {rec_title[:70]}  ({n} file{'s' if n!=1 else ''})")
                for f in hit.get("files", []):
                    entry = {
                        "filename": f.get("key", "unknown"),
                        "url": f.get("links", {}).get("self", ""),
                        "repository": "Zenodo",
                        "accession": f"zenodo.{record_id}",
                        "size_bytes": f.get("size"),
                        "checksum": f.get("checksum", "").replace("md5:", ""),
                        "data_type": infer_data_type(f.get("key", ""), description),
                    }
                    if not any(x["url"] == entry["url"] for x in files):
                        files.append(entry)
            if files:
                break
        except requests.RequestException:
            continue
    return files


def search_geo_by_article(doi: Optional[str], pmid: Optional[str], title: str) -> list[dict]:
    """Search NCBI GEO DataSets for studies linked to the article."""
    files: list[dict] = []
    # Build query terms
    terms: list[str] = []
    if pmid:
        terms.append(f"{pmid}[PMID]")
    if doi:
        terms.append(f"{doi}[DOI]")
    if title and title != "Unknown":
        short = re.sub(r'[<>"\'\[\]]', '', title[:80])
        terms.append(f'"{short}"[Title]')
    for term in terms:
        try:
            print(f"  🔎 GEO: {term}")
            r = requests.get(
                "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
                params={"db": "gds", "term": term, "retmax": 10, "retmode": "json"},
                timeout=15,
            )
            if r.status_code != 200:
                continue
            ids = r.json().get("esearchresult", {}).get("idlist", [])
            if not ids:
                continue
            # Convert GEO DataSets UIDs → GSE accessions via esummary
            sr = requests.get(
                "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
                params={"db": "gds", "id": ",".join(ids), "retmode": "json"},
                timeout=15,
            )
            if sr.status_code != 200:
                continue
            for uid, doc in sr.json().get("result", {}).items():
                if uid == "uids":
                    continue
                acc = doc.get("accession", "")
                if acc.startswith("GSE"):
                    print(f"     • {acc}: {doc.get('title', '')[:70]}")
                    files.extend(fetch_geo_files(acc))
            if files:
                break
        except requests.RequestException:
            continue
    return files


def search_ena_by_article(doi: Optional[str], pmid: Optional[str]) -> list[dict]:
    """Search ENA for studies linked to the article DOI or PMID."""
    files: list[dict] = []
    refs = []
    if doi:
        refs.append(("study_doi", doi))
    if pmid:
        refs.append(("pubmed_id", pmid))
    for field, value in refs:
        try:
            print(f"  🔎 ENA: {field}={value}")
            r = requests.get(
                "https://www.ebi.ac.uk/ena/portal/api/search",
                params={
                    "query": f'{field}="{value}"',
                    "result": "study",
                    "fields": "study_accession,study_title",
                    "format": "json",
                    "limit": 10,
                },
                timeout=15,
            )
            if r.status_code != 200:
                continue
            studies = r.json()
            for s in studies:
                acc = s.get("study_accession", "")
                print(f"     • {acc}: {s.get('study_title', '')[:70]}")
                files.extend(fetch_ena_files(acc))
            if files:
                break
        except requests.RequestException:
            continue
    return files


def search_dryad_by_doi(doi: Optional[str]) -> list[dict]:
    """Search Dryad for datasets linked to the article DOI."""
    if not doi:
        return []
    files: list[dict] = []
    try:
        print(f"  🔎 Dryad: {doi}")
        r = requests.get(
            "https://datadryad.org/api/v2/search",
            params={"q": doi, "per_page": 10},
            timeout=15,
            headers={"Accept": "application/json"},
        )
        if r.status_code != 200:
            return []
        for ds in r.json().get("_embedded", {}).get("stash:datasets", []):
            identifier = ds.get("identifier", "")
            title = ds.get("title", identifier)
            download_uri = ds.get("_links", {}).get("stash:download", {}).get("href", "")
            print(f"     • {title[:70]}")
            if download_uri:
                url = f"https://datadryad.org{download_uri}" if download_uri.startswith("/") else download_uri
                files.append({
                    "filename": f"{identifier.replace('/', '_')}.zip",
                    "url": url,
                    "repository": "Dryad",
                    "accession": identifier,
                    "size_bytes": None,
                    "data_type": infer_data_type(title),
                })
    except requests.RequestException:
        pass
    return files


def search_datacite_by_doi(doi: str) -> list[dict]:
    """
    Query DataCite for datasets that cite or are related to the article DOI.
    DataCite covers Zenodo, Dryad, Figshare, and many institutional repositories.
    Returns accessions/URLs for further resolution where possible.
    """
    files: list[dict] = []
    try:
        print(f"  🔎 DataCite (meta-index): {doi}")
        r = requests.get(
            "https://api.datacite.org/dois",
            params={
                "query": f'relatedIdentifiers.relatedIdentifier:"{doi}"',
                "resource-type-id": "dataset",
                "page[size]": 10,
            },
            timeout=15,
            headers={"Accept": "application/json"},
        )
        if r.status_code != 200:
            return []
        for item in r.json().get("data", []):
            attrs = item.get("attributes", {})
            ds_doi = attrs.get("doi", "")
            ds_title = (attrs.get("titles") or [{}])[0].get("title", ds_doi)
            publisher = attrs.get("publisher", "")
            url = attrs.get("url", "")
            print(f"     • [{publisher}] {ds_title[:60]}")
            # Try to resolve known repositories from the dataset DOI/URL
            if "zenodo" in url.lower():
                m = re.search(r"zenodo\.org/records?/(\d+)", url, re.IGNORECASE)
                if m:
                    files.extend(fetch_zenodo_files(m.group(1)))
                    continue
            if "dryad" in url.lower() or "dryad" in publisher.lower():
                if url:
                    files.append({
                        "filename": f"{ds_doi.replace('/', '_')}.zip",
                        "url": url,
                        "repository": "Dryad (DataCite)",
                        "accession": ds_doi,
                        "size_bytes": None,
                        "data_type": infer_data_type(ds_title),
                    })
                continue
            # Generic: just surface the landing page URL so user can visit it
            if url:
                files.append({
                    "filename": ds_title[:80],
                    "url": url,
                    "repository": publisher or "DataCite",
                    "accession": ds_doi,
                    "size_bytes": None,
                    "data_type": infer_data_type(ds_title),
                })
    except requests.RequestException:
        pass
    return files


# ---------------------------------------------------------------------------
# Download helpers
# ---------------------------------------------------------------------------

def file_matches_types(filename: str, types: set[str]) -> bool:
    if "all" in types:
        return True
    name_lower = filename.lower()
    for t in types:
        if name_lower.endswith(f".{t}"):
            return True
    return False


def md5_file(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def unzip_file(path: Path) -> Path:
    """Decompress a .gz file in-place, returning the path to the unzipped file."""
    if path.suffix != ".gz":
        return path
    out_path = path.with_suffix("")  # strip .gz
    print(f"  📦 Unzipping → {out_path.name}")
    with gzip.open(path, "rb") as f_in, open(out_path, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
    path.unlink()  # remove the .gz original
    return out_path


def download_file(url: str, dest: Path, expected_md5: Optional[str] = None, unzip: bool = True) -> bool:
    """Stream-download a file with a progress bar. Returns True on success."""
    try:
        r = requests.get(url, stream=True, timeout=30)
        if r.status_code in (401, 403):
            print(f"  ⛔ Access denied (HTTP {r.status_code}): {url}")
            print("     This file may be behind a paywall.")
            return False
        r.raise_for_status()
        total = int(r.headers.get("content-length", 0))
        dest.parent.mkdir(parents=True, exist_ok=True)
        with open(dest, "wb") as f, tqdm(
            total=total, unit="B", unit_scale=True, desc=dest.name, leave=False
        ) as bar:
            for chunk in r.iter_content(chunk_size=65536):
                f.write(chunk)
                bar.update(len(chunk))
        if expected_md5:
            actual = md5_file(dest)
            if actual != expected_md5:
                print(f"  ❌ Checksum mismatch for {dest.name} — file deleted")
                dest.unlink()
                return False
        if unzip and dest.suffix == ".gz":
            unzip_file(dest)
        return True
    except requests.RequestException as exc:
        print(f"  ❌ Download failed: {exc}")
        return False


# ---------------------------------------------------------------------------
# Main orchestration
# ---------------------------------------------------------------------------

def run(
    identifier: str,
    file_types: set[str],
    output_dir: Path,
    non_interactive: bool = False,
    unzip: bool = True,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    # 1. Resolve article
    kind, value = parse_identifier(identifier)
    print(f"\n🔍 Resolving article: {identifier}")
    try:
        if kind == "doi":
            meta = resolve_doi(value)
        elif kind == "pmid":
            meta = resolve_pmid(value)
        elif kind == "pmcid":
            meta = resolve_pmcid(value)
        else:
            meta = {"title": identifier, "doi": identifier, "abstract": "", "links": []}
    except Exception as exc:
        print(f"⚠️  Could not resolve article metadata: {exc}")
        meta = {"title": identifier, "abstract": "", "links": []}

    print(f"📄 Article: {meta['title']}")

    # 2. Discover accessions from abstract + links
    search_text = meta.get("abstract", "") + " ".join(meta.get("links", []))
    # Also try fetching the Crossref page for a data availability note
    accessions = discover_accessions(search_text)

    # If identifier itself is an accession, treat it directly
    if kind == "url":
        for repo, pattern in ACCESSION_PATTERNS.items():
            m = pattern.search(value)
            if m:
                accessions.setdefault(repo, []).append(m.group(0))

    if not accessions:
        print("\n⚠️  No repository accessions found in the article text.")
        print("   Searching external repositories for linked datasets…\n")
        # Use the DOI from the resolved metadata, regardless of what the user provided
        doi  = meta.get("doi") or (value if kind == "doi" else None)
        pmid = meta.get("pmid") or (value if kind == "pmid" else None)
        if doi:
            print(f"   (using DOI: {doi})")

        all_files = []
        if doi:
            all_files.extend(search_datacite_by_doi(doi))
        all_files.extend(search_zenodo_by_article(doi, meta["title"]))
        all_files.extend(search_geo_by_article(doi, pmid, meta["title"]))
        all_files.extend(search_ena_by_article(doi, pmid))
        all_files.extend(search_dryad_by_doi(doi))

        # Deduplicate by URL
        seen: set[str] = set()
        deduped = []
        for f in all_files:
            if f["url"] not in seen:
                seen.add(f["url"])
                deduped.append(f)
        all_files = deduped

        if not all_files:
            print(
                "\n⚠️  Nothing found across Zenodo, GEO, ENA, Dryad, or DataCite.\n"
                "   Please check the paper's Data Availability Statement and provide\n"
                "   the accession number directly (e.g. GSE123456, PRJNA654321)."
            )
            return
    else:
        # 3. Fetch file listings from discovered accessions
        all_files: list[dict] = []
        for repo, acc_list in accessions.items():
            for acc in acc_list:
                print(f"📂 Fetching file list from {repo.upper()} ({acc})…")
                if repo == "geo":
                    all_files.extend(fetch_geo_files(acc))
                elif repo in ("sra", "arrayexpress"):
                    all_files.extend(fetch_ena_files(acc))
                elif repo == "zenodo":
                    record_id = re.search(r"(\d+)$", acc)
                    if record_id:
                        all_files.extend(fetch_zenodo_files(record_id.group(1)))
                elif repo == "zenodo_community":
                    all_files.extend(fetch_zenodo_community_files(acc))

    if not all_files:
        print("⚠️  No downloadable files discovered. The repository may be private or the data not yet public.")
        return

    # 4. Present available files
    print_file_listing(all_files, title=f"\nFound {len(all_files)} file(s):")

    # 5. Confirm file types (unless non-interactive)
    if not non_interactive:
        selected = prompt_and_select_files(all_files)
        if not selected:
            print("No files selected. Aborting.")
            return
    else:
        selected = [f for f in all_files if file_matches_types(f["filename"], file_types)] if file_types else all_files

    if not selected:
        print(f"\n⚠️  No files match the selection.")
        return

    # Warn if total size is large
    known_sizes = [f["size_bytes"] for f in selected if f.get("size_bytes")]
    if known_sizes:
        total_bytes = sum(known_sizes)
        if total_bytes > SIZE_WARNING_BYTES:
            print(
                f"\n⚠️  WARNING: Total known download size is "
                f"{_fmt_size(total_bytes)}. Proceed? [y/N] "
            )
            if not non_interactive:
                answer = input("→ ").strip().lower()
                if answer != "y":
                    print("Aborted.")
                    return

    # 6. Download
    print(f"\n⬇️  Downloading {len(selected)} file(s) to {output_dir}/\n")
    manifest_files = []
    for f in selected:
        acc_dir = output_dir / f["accession"]
        dest = acc_dir / f["filename"]
        print(f"  ⬇  {f['filename']}")
        ok = download_file(f["url"], dest, expected_md5=f.get("checksum"), unzip=unzip)
        manifest_files.append({
            **f,
            "local_path": str(dest),
            "downloaded": ok,
            "downloaded_at": datetime.now(timezone.utc).isoformat(),
        })
        if ok:
            print(f"  ✅ {f['filename']}")

    # 7. Write manifest
    manifest = {
        "article": identifier,
        "title": meta["title"],
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "output_dir": str(output_dir),
        "file_types_requested": list(file_types),
        "files": manifest_files,
    }
    manifest_path = output_dir / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2))

    # 8. Write report
    downloaded = [f for f in manifest_files if f["downloaded"]]
    skipped = [f for f in all_files if not file_matches_types(f["filename"], file_types)]
    report_lines = [
        "# Article Data Fetcher — Download Report",
        "",
        f"**Article**: {meta['title']}",
        f"**Identifier**: {identifier}",
        f"**Date**: {datetime.now(timezone.utc).strftime('%Y-%m-%d')}",
        "",
        f"## Files Downloaded ({len(downloaded)})",
        "",
    ]
    for f in downloaded:
        size_str = _fmt_size(f.get("size_bytes"))
        report_lines.append(f"- ✅ `{f['filename']}` ({size_str}) — {f['repository']}")
    if skipped:
        report_lines += ["", f"## Files Skipped — type not requested ({len(skipped)})", ""]
        for f in skipped:
            report_lines.append(f"- ⏭  `{f['filename']}` — {f['repository']}")
    report_lines += [
        "",
        "---",
        "*ClawBio is a research tool. Verify data integrity before use in analysis.*",
    ]
    report_path = output_dir / "report.md"
    report_path.write_text("\n".join(report_lines))

    print(f"\n✅ Done. {len(downloaded)}/{len(selected)} files downloaded.")
    print(f"   Manifest : {manifest_path}")
    print(f"   Report   : {report_path}")


# ---------------------------------------------------------------------------
# Demo mode
# ---------------------------------------------------------------------------

def run_demo(output_dir: Path, unzip: bool = True) -> None:
    print(f"🧬 Running demo with GEO accession {DEMO_ACCESSION}…")
    output_dir.mkdir(parents=True, exist_ok=True)
    files = fetch_geo_files(DEMO_ACCESSION)
    if not files:
        print("⚠️  Could not reach GEO FTP in demo mode. Check your internet connection.")
        return

    # Show the same file listing as the real flow
    print_file_listing(files, title=f"\nFound {len(files)} file(s) in {DEMO_ACCESSION}:")

    # Ask the user which files they want — same as the real flow
    selected = prompt_and_select_files(files)
    if not selected:
        print("No files selected. Aborting.")
        return

    print(f"\n⬇️  Downloading {len(selected)} file(s) to {output_dir}/\n")
    manifest_files = []
    for f in selected:
        dest = output_dir / f["accession"] / f["filename"]
        print(f"  ⬇  {f['filename']}")
        ok = download_file(f["url"], dest, unzip=unzip)
        manifest_files.append({**f, "local_path": str(dest), "downloaded": ok})
        if ok:
            print(f"  ✅ {f['filename']}")

    manifest_path = output_dir / "manifest.json"
    manifest_path.write_text(json.dumps({"demo": True, "accession": DEMO_ACCESSION, "files": manifest_files}, indent=2))
    report_path = output_dir / "report.md"
    report_path.write_text(
        f"# Demo Report\n\nAccession: {DEMO_ACCESSION}\n"
        f"Files downloaded: {sum(1 for f in manifest_files if f['downloaded'])}\n"
    )
    print(f"\n✅ Demo complete.\n   Manifest : {manifest_path}\n   Report   : {report_path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Download genomics data files from repositories linked to a paper."
    )
    parser.add_argument("--id", dest="identifier", help="DOI, PMID, or repository URL")
    parser.add_argument(
        "--types",
        default="",
        help="Comma-separated file extensions to download (e.g. vcf,fasta,h5ad) or 'all'",
    )
    parser.add_argument(
        "--output", default="./downloads", help="Output directory (default: ./downloads)"
    )
    parser.add_argument("--demo", action="store_true", help="Run in demo mode")
    parser.add_argument(
        "--no-unzip",
        action="store_true",
        help="Keep .gz files compressed instead of unzipping after download",
    )
    args = parser.parse_args()

    output_dir = Path(args.output)
    unzip = not args.no_unzip

    if args.demo:
        if args.output == "./downloads":
            output_dir = Path("skills/article-data-fetcher/demo")
        run_demo(output_dir, unzip=unzip)
        return

    if not args.identifier:
        parser.error("--id is required unless --demo is used")

    types: set[str] = set()
    if args.types:
        raw = args.types.strip().lower()
        types = {"all"} if raw == "all" else {t.strip() for t in raw.split(",") if t.strip()}

    run(identifier=args.identifier, file_types=types, output_dir=output_dir, unzip=unzip)


if __name__ == "__main__":
    main()
