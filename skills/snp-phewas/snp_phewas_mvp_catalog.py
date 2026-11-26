"""The Million Veteran Program outcome catalog, built from the GWAS Catalog's own records.

MVP's genome-wide PheWAS (Verma A et al., Science 2024;385:eadj1182,
doi:10.1126/science.adj1182, PMID 39024449) is served by the NHGRI-EBI GWAS Catalog as
one study per (trait, population): 7,537 studies on the 2026-10-19 snapshot. The Catalog
record carries no MVP trait id, so everything here is read from the record itself:

  * `disease_trait` -> the trait kind and code. A phecode study ends in `(PheCode 447)`;
    a quantitative study ends in `(..., inv-norm transformed)` and names its summary
    statistic (mean / maximum / minimum of the EHR values, or a single survey value); any
    other binary study is a survey item.
  * `discovery_ancestry` -> the population (EUR / AFR / AMR / EAS), or META for the
    Catalog's multi-population meta-analyses.
  * `initial_sample_size` -> cases and controls (binary) or individuals (quantitative).

Which studies can be read by position is not in the record: only the studies the Catalog
has HARMONIZED have a sorted, bgzipped, tabix-indexed file (`harmonised_list.txt`). Every
other study is still a catalog row, with `access_status` saying why it cannot be read
(`raw_only_unindexed`: the Catalog hosts an unsorted plain-gzip file only;
`no_summary_statistics`: nothing deposited). A lookup returns those as
`source_file_not_readable` rows, never as absences.

Two halves, so the committed catalog is reproducible offline:

  `refresh_snapshot(out_dir)`  (network, run by hand) fetches the three source snapshots
                               into `catalogs/mvp_source/`: the study records (REST v2), the MVP
                               lines of `harmonised_list.txt`, and each harmonized file's
                               metadata (`<acc>.h.tsv.gz-meta.yaml`: md5, date, assembly).
  `build_rows(...)`            (offline) turns the pinned snapshots into catalog rows, with
                               the write-time assertions.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
import re
from collections import Counter
from typing import Any, Iterable

from snp_phewas_crosswalk import normalize_phecode, quant_key  # noqa: F401

PUBMED_ID = 39024449
STUDIES_URL = ("https://www.ebi.ac.uk/gwas/rest/api/v2/studies?pubmed_id=39024449"
               "&page={page}&size=200")
HARMONISED_LIST_URL = ("https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
                       "harmonised_list.txt")
DATA_BASE_URL = "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"

# The fields of a REST v2 study record the snapshot keeps (the rest is presentation).
STUDY_FIELDS = ("accession_id", "disease_trait", "initial_sample_size", "discovery_ancestry",
                "efo_traits", "full_summary_stats_available", "full_summary_stats", "cohort",
                "terms_of_license", "pubmed_id")

META_COLUMNS = ("accession", "data_file_name", "data_file_md5sum", "date_metadata_last_modified",
                "genome_assembly", "is_sorted", "is_harmonised", "file_type")

ACCESS_HARMONISED = "harmonised"
ACCESS_RAW_ONLY = "raw_only_unindexed"
ACCESS_NO_SUMSTATS = "no_summary_statistics"
ACCESS_STATUSES = (ACCESS_HARMONISED, ACCESS_RAW_ONLY, ACCESS_NO_SUMSTATS)
ACCESS_STATUS_LABELS = {
    ACCESS_HARMONISED: "harmonized by the GWAS Catalog (sorted, position-indexed file)",
    ACCESS_RAW_ONLY: ("the GWAS Catalog holds this study without a position index (it has "
                      "not harmonized it), so it cannot be read by position"),
    ACCESS_NO_SUMSTATS: "deposited without summary statistics",
}

KIND_PHECODE = "phecode"
KIND_QUANT = "quantitative"
KIND_SURVEY = "survey_binary"
TRAIT_KIND_LABELS = {
    KIND_PHECODE: "disease (phecode, from the health record)",
    KIND_QUANT: "measurement (lab, vital sign or survey number)",
    KIND_SURVEY: "survey answer",
}

POPULATIONS = ("EUR", "AFR", "AMR", "EAS")
META = "META"
# The Catalog's ancestry wording for MVP's genetically inferred groups.
_ANCESTRY = {
    "European (U.S.)": "EUR",
    "African American or Afro-Caribbean (U.S.)": "AFR",
    "Hispanic or Latin American (U.S.)": "AMR",
    "East Asian (U.S.)": "EAS",
}

_PHECODE_RE = re.compile(r"^(?P<label>.*?)\s*\(PheCode\s+(?P<code>[0-9.]+)\)\s*$")
_QUANT_RE = re.compile(r"^(?P<base>.*?)\s*\((?P<inner>[^()]*transformed)\)\s*$")
_STATS = {"mean": "mean", "maximum": "maximum", "minimum": "minimum",
          "minimm": "minimum"}  # `minimm`: a typo in the Catalog's HDL-C minimum record
QUANT_STAT_WORDS = {"mean": "mean of EHR values", "maximum": "maximum of EHR values",
                    "minimum": "minimum of EHR values", "single": ""}
_N_RE = re.compile(r"^(?P<n>\d[\d,]*)\s+.*?\b(?P<what>cases|controls|individuals)\s*$")

CATALOG_COLUMNS: tuple[str, ...] = (
    "outcome_id", "population", "meta_populations", "trait_text", "trait_kind", "phecode",
    "phecode_display", "quant_base", "quant_qualifier", "quant_stat", "outcome_label",
    "category_native", "efo_ids", "n", "n_cases", "n_controls", "access_status", "filename",
    "filename_tabix", "md5_hex", "harmonised_date",
)


class MvpCatalogError(ValueError):
    """A study record that the catalog rules cannot place; never silently skipped."""


# --- parsing one record -------------------------------------------------------------------

def parse_trait(text: str) -> dict[str, str]:
    """{trait_kind?, phecode, phecode_display, quant_base, quant_qualifier, quant_stat,
    label} from the Catalog's `disease_trait`. The kind of a non-phecode, non-quantitative
    trait is decided by the sample size (binary = survey item), not here."""
    t = " ".join(text.split())
    m = _PHECODE_RE.match(t)
    if m:
        code = m.group("code")
        return {"trait_kind": KIND_PHECODE, "phecode": normalize_phecode(code),
                "phecode_display": code, "quant_base": "", "quant_qualifier": "",
                "quant_stat": "", "label": m.group("label")}
    m = _QUANT_RE.match(t)
    if m:
        toks = [s.strip() for s in m.group("inner").split(",")]
        stat, quals = "single", []
        for tk in toks[:-1]:  # the last token is the transformation
            if tk in _STATS:
                stat = _STATS[tk]
            else:
                quals.append(tk)
        base = m.group("base")
        words = [*quals, QUANT_STAT_WORDS[stat]]
        label = f"{base} ({', '.join(w for w in words if w)}{', ' if any(words) else ''}inverse-normal transformed)"
        return {"trait_kind": KIND_QUANT, "phecode": "", "phecode_display": "",
                "quant_base": base, "quant_qualifier": ", ".join(quals), "quant_stat": stat,
                "label": label}
    return {"trait_kind": "", "phecode": "", "phecode_display": "", "quant_base": "",
            "quant_qualifier": "", "quant_stat": "", "label": t}


def parse_population(discovery_ancestry: list[str]) -> tuple[str, str]:
    """(population, meta_populations) from the Catalog's discovery ancestry entries."""
    pops = []
    for a in discovery_ancestry:
        name = a.split(" ", 1)[1] if a[:1].isdigit() else a
        if name not in _ANCESTRY:
            raise MvpCatalogError(f"unknown MVP discovery ancestry {a!r}")
        pops.append(_ANCESTRY[name])
    if not pops:
        raise MvpCatalogError("a study with no discovery ancestry")
    if len(pops) == 1:
        return pops[0], ""
    return META, "+".join(p for p in POPULATIONS if p in pops)


def parse_sample_size(text: str) -> tuple[int | None, int | None, int]:
    """(n_cases, n_controls, n) summed over the populations of `initial_sample_size`."""
    cases = controls = indiv = 0
    seen = Counter()
    for seg in text.split(", "):  # thousands separators carry no space
        m = _N_RE.match(seg.strip())
        if not m:
            raise MvpCatalogError(f"unparsable sample size segment {seg!r} in {text!r}")
        n = int(m.group("n").replace(",", ""))
        seen[m.group("what")] += 1
        if m.group("what") == "cases":
            cases += n
        elif m.group("what") == "controls":
            controls += n
        else:
            indiv += n
    if seen["individuals"] and (seen["cases"] or seen["controls"]):
        raise MvpCatalogError(f"mixed binary and quantitative sample size {text!r}")
    if seen["individuals"]:
        return None, None, indiv
    if not (seen["cases"] and seen["controls"]):
        raise MvpCatalogError(f"a binary sample size without cases and controls {text!r}")
    return cases, controls, cases + controls


# --- the offline build --------------------------------------------------------------------

def _md5(b: bytes) -> str:
    return hashlib.md5(b).hexdigest()


def read_studies(jsonl_gz: bytes) -> list[dict[str, Any]]:
    return [json.loads(line) for line in gzip.decompress(jsonl_gz).decode("utf-8").splitlines()
            if line.strip()]


def read_harmonised_paths(subset_gz: bytes) -> dict[str, str]:
    """{accession: relative path of the harmonized file} from the MVP lines of the list."""
    out: dict[str, str] = {}
    for line in gzip.decompress(subset_gz).decode("utf-8").splitlines():
        p = line.strip()
        if not p:
            continue
        p = p[2:] if p.startswith("./") else p
        acc = p.split("/")[1]
        if acc in out:
            raise MvpCatalogError(f"{acc} listed twice in harmonised_list.txt")
        out[acc] = p
    return out


def read_meta(meta_gz: bytes) -> dict[str, dict[str, str]]:
    rows = csv.DictReader(io.StringIO(gzip.decompress(meta_gz).decode("utf-8")), delimiter="\t")
    return {r["accession"]: r for r in rows}


def build_rows(studies: list[dict[str, Any]], harmonised: dict[str, str],
               meta: dict[str, dict[str, str]]) -> list[dict[str, str]]:
    """One catalog row per study, with the write-time assertions:

    * `outcome_id` (the GCST accession) unique and non-null on every row, and the rows equal
      the study list as SETS;
    * `filename`, `filename_tabix` and `md5_hex` set exactly when `access_status` is
      harmonised, and the harmonized file is GRCh38 (else the build refuses rather than
      read it at the wrong position);
    * every harmonized path names a study in the list (a stray line fails loudly).
    """
    accs = [s["accession_id"] for s in studies]
    if len(set(accs)) != len(accs) or not all(accs):
        raise MvpCatalogError("the study list repeats or omits an accession")
    stray = set(harmonised) - set(accs)
    if stray:
        raise MvpCatalogError(f"harmonized paths for studies not in the list: {sorted(stray)[:5]}")
    rows = []
    for s in sorted(studies, key=lambda x: x["accession_id"]):
        acc = s["accession_id"]
        pop, meta_pops = parse_population(s["discovery_ancestry"])
        n_cases, n_controls, n = parse_sample_size(s["initial_sample_size"])
        tr = parse_trait(s["disease_trait"])
        binary = n_controls is not None
        if tr["trait_kind"] == KIND_QUANT and binary:
            raise MvpCatalogError(f"{acc}: a transformed trait with cases and controls")
        if tr["trait_kind"] == KIND_PHECODE and not binary:
            raise MvpCatalogError(f"{acc}: a phecode trait without cases and controls")
        if not tr["trait_kind"]:
            if not binary:
                raise MvpCatalogError(f"{acc}: an untyped quantitative trait {s['disease_trait']!r}")
            tr["trait_kind"] = KIND_SURVEY
            tr["label"] = f"{tr['label']} (survey answer)"
        if acc in harmonised:
            access = ACCESS_HARMONISED
            m = meta.get(acc)
            if m is None or not m.get("data_file_md5sum"):
                raise MvpCatalogError(f"{acc}: harmonised but no metadata md5 in the snapshot")
            # `is_sorted` is NOT required: 12 of the 4,442 harmonized MVP files' metadata say
            # false, all 12 carry a .tbi (HEAD-checked 2026-10-19), and tabix refuses to
            # index a file whose positions are out of order, so the flag is stale there.
            if m["genome_assembly"] != "GRCh38":
                raise MvpCatalogError(f"{acc}: harmonized file is {m['genome_assembly']}, "
                                      f"not GRCh38; not readable at a GRCh38 position")
            filename = harmonised[acc]
            if not filename.endswith(m["data_file_name"]):
                raise MvpCatalogError(f"{acc}: metadata names {m['data_file_name']}, list names {filename}")
            md5, date = m["data_file_md5sum"], m["date_metadata_last_modified"]
        else:
            access = (ACCESS_RAW_ONLY if s.get("full_summary_stats_available")
                      else ACCESS_NO_SUMSTATS)
            filename, md5, date = "", "", ""
        efo = s.get("efo_traits") or []
        rows.append({
            "outcome_id": acc, "population": pop, "meta_populations": meta_pops,
            "trait_text": " ".join(s["disease_trait"].split()),
            "trait_kind": tr["trait_kind"], "phecode": tr["phecode"],
            "phecode_display": tr["phecode_display"], "quant_base": tr["quant_base"],
            "quant_qualifier": tr["quant_qualifier"], "quant_stat": tr["quant_stat"],
            "outcome_label": tr["label"],
            "category_native": (f"PheCode {tr['phecode']}" if tr["trait_kind"] == KIND_PHECODE
                                else TRAIT_KIND_LABELS[tr["trait_kind"]]),
            "efo_ids": ";".join(e["efo_id"] for e in efo),
            "n": str(n), "n_cases": "" if n_cases is None else str(n_cases),
            "n_controls": "" if n_controls is None else str(n_controls),
            "access_status": access, "filename": filename,
            "filename_tabix": f"{filename}.tbi" if filename else "",
            "md5_hex": md5, "harmonised_date": date,
        })
    assert_catalog(rows, set(accs))
    return rows


def assert_catalog(rows: list[dict[str, str]], study_ids: set[str]) -> None:
    """Raised (not `assert`ed) write-time invariants of the MVP catalog."""
    ids = [r["outcome_id"] for r in rows]
    if not all(ids) or len(set(ids)) != len(ids) or set(ids) != study_ids:
        raise MvpCatalogError(
            f"catalog rows do not partition the {len(study_ids)} studies "
            f"({len(set(ids))} distinct of {len(ids)}; {len(study_ids - set(ids))} missing)")
    for r in rows:
        if r["access_status"] not in ACCESS_STATUSES:
            raise MvpCatalogError(f"{r['outcome_id']}: unknown access status {r['access_status']!r}")
        readable = r["access_status"] == ACCESS_HARMONISED
        for col in ("filename", "filename_tabix", "md5_hex"):
            if bool(r[col]) != readable:
                raise MvpCatalogError(f"{r['outcome_id']}: {col} must be set exactly when "
                                      f"harmonised (status {r['access_status']})")
        if r["trait_kind"] not in TRAIT_KIND_LABELS:
            raise MvpCatalogError(f"{r['outcome_id']}: unknown trait kind {r['trait_kind']!r}")


def write_catalog(rows: Iterable[dict[str, str]], path: str) -> int:
    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=list(CATALOG_COLUMNS), delimiter="\t",
                       lineterminator="\n")
    w.writeheader()
    n = 0
    for r in rows:
        w.writerow(r)
        n += 1
    with open(path, "wb") as fh:
        # filename="" and mtime=0: the gzip header carries neither the output path nor the
        # time, so a rebuild from the same snapshots is byte-identical wherever it is written
        with gzip.GzipFile(filename="", fileobj=fh, mode="wb", mtime=0) as gz:
            gz.write(buf.getvalue().encode("utf-8"))
    return n


# --- the network snapshot (run by hand) ---------------------------------------------------

def _gz(text: str) -> bytes:
    bio = io.BytesIO()
    with gzip.GzipFile(fileobj=bio, mode="wb", mtime=0) as gz:
        gz.write(text.encode("utf-8"))
    return bio.getvalue()


def refresh_snapshot(out_dir: str, *, workers: int = 8, log=print) -> dict[str, Any]:  # pragma: no cover - network
    """Fetch the three source snapshots into `out_dir`; returns their md5s for the pin."""
    import os
    from concurrent.futures import ThreadPoolExecutor

    import requests
    import yaml

    from snp_phewas_tabix import get_with_retry

    sess = requests.Session()

    def page_of(page: int) -> dict[str, Any]:
        r = get_with_retry(sess, STUDIES_URL.format(page=page), attempts=5, backoff_s=5.0,
                           timeout=300)
        r.raise_for_status()
        return r.json()

    first = page_of(0)
    total_pages, total = first["page"]["totalPages"], first["page"]["totalElements"]
    pages = [first]
    with ThreadPoolExecutor(max_workers=4) as ex:  # one page of 200 takes ~20 s to serve
        for i, d in enumerate(ex.map(page_of, range(1, total_pages)), 2):
            pages.append(d)
            log(f"study list page {i}/{total_pages}")
    studies = [{k: s.get(k) for k in STUDY_FIELDS}
               for d in pages for s in d["_embedded"]["studies"]]
    if len(studies) != total or len({s["accession_id"] for s in studies}) != total:
        raise MvpCatalogError(f"read {len(studies)} study records of {total} "
                              f"({len({s['accession_id'] for s in studies})} distinct)")
    studies.sort(key=lambda s: s["accession_id"])
    accs = {s["accession_id"] for s in studies}
    hl = get_with_retry(sess, HARMONISED_LIST_URL, attempts=5, backoff_s=5.0, timeout=600)
    hl.raise_for_status()
    full = hl.content
    lines = sorted(ln for ln in full.decode("utf-8").splitlines()
                   if len(ln.split("/")) > 2 and ln.split("/")[2] in accs)
    paths = read_harmonised_paths(_gz("\n".join(lines) + "\n"))

    def meta_of(acc: str) -> dict[str, str]:
        url = DATA_BASE_URL + paths[acc] + "-meta.yaml"
        for attempt in range(4):
            try:
                resp = sess.get(url, timeout=60)
                if resp.status_code == 200:
                    m = yaml.safe_load(resp.text)
                    return {"accession": acc, "data_file_name": m.get("data_file_name", ""),
                            "data_file_md5sum": m.get("data_file_md5sum", ""),
                            "date_metadata_last_modified": str(m.get("date_metadata_last_modified", "")),
                            "genome_assembly": str(m.get("genome_assembly", "")),
                            "is_sorted": str(m.get("is_sorted", "")).lower(),
                            "is_harmonised": str(m.get("is_harmonised", "")).lower(),
                            "file_type": str(m.get("file_type", ""))}
            except requests.RequestException:
                pass  # a transient error: retry, then raise below
        raise MvpCatalogError(f"{acc}: harmonized metadata could not be read from {url}")

    metas: list[dict[str, str]] = []
    with ThreadPoolExecutor(max_workers=workers) as ex:
        for i, m in enumerate(ex.map(meta_of, sorted(paths)), 1):
            metas.append(m)
            if i % 250 == 0 or i == len(paths):
                log(f"harmonized metadata {i}/{len(paths)}")
    os.makedirs(out_dir, exist_ok=True)
    blobs = {
        "studies_jsonl_gz": _gz("".join(json.dumps(s, sort_keys=True) + "\n" for s in studies)),
        "harmonised_subset_gz": _gz("\n".join(lines) + "\n"),
        "harmonised_meta_gz": _gz("".join(["\t".join(META_COLUMNS) + "\n"]
                                          + ["\t".join(m[c] for c in META_COLUMNS) + "\n"
                                             for m in metas])),
    }
    names = {"studies_jsonl_gz": "gwas_catalog_studies_pmid39024449.jsonl.gz",
             "harmonised_subset_gz": "harmonised_list_mvp.txt.gz",
             "harmonised_meta_gz": "harmonised_meta_mvp.tsv.gz"}
    out: dict[str, Any] = {"harmonised_list_md5": _md5(full),
                           "harmonised_list_last_modified": hl.headers.get("Last-Modified", ""),
                           "n_studies": len(studies), "n_harmonised": len(paths)}
    for k, blob in blobs.items():
        with open(os.path.join(out_dir, names[k]), "wb") as fh:
            fh.write(blob)
        out[f"{names[k]}_md5"] = _md5(blob)
    return out
