#!/usr/bin/env python3
"""Choose the QTL and the GWAS for a locuscompare run from an Open Targets credible-set page.

An Open Targets Platform credible-set page (for example
https://platform.opentargets.org/credible-set/2f6755df0a7737586c9ff36fa2643c61)
lists the molecular QTLs that colocalise with one GWAS signal. locuscompare draws
one pair at a time, so somebody has to pick the QTL, and the page's GWAS study has
to be readable by region (the GWAS Catalog harmonised file with its tabix index).
This helper makes both choices explicit and repeatable:

    python ot_credible_set.py options <studyLocusId or URL> [--gene SYMBOL] [--qtl ROW_OR_ID]
    python ot_credible_set.py config  <studyLocusId or URL> --qtl ROW_OR_ID [--gwas GCST...] --out config.yaml

`options` prints the colocalising QTLs (table A) and a recommendation with its
reasons. Once a QTL is chosen (or the page lists only one) it checks the page's
GWAS file; if that file is definitely not readable by region, it lists other GWAS
studies that colocalise with the same QTL (table B) and recommends one. `config`
writes the locuscompare config for the chosen pair, recording every substitution
(a different GWAS study, an older release of the QTL data) as provenance and as
plain-language caveats that the render prints.

Rules, all in code below:

* Table A folds rows that measure the same thing (same gene, tissue, study and QTL
  kind) and keeps the one with the highest colocalisation probability (H4). Rows of
  different kinds (gene-level eQTL vs exon-level eQTL, say) are never folded.
* Order: the gene named with --gene first, then H4 (high to low), then sample size.
* At most 15 rows are printed per table; the number not printed is stated.
* Every row Open Targets returned ends in exactly one state: shown, hidden by the
  display cap, or excluded with a named reason; the counts must add up to what
  Open Targets reported, or the helper stops with an error.
* A GWAS file check that fails (timeout, network error, server error) is
  "unknown", never "no": it does not trigger a substitution and is never
  recommended as a substitute.

Standard library only. Data: Open Targets Platform GraphQL API (data CC0 1.0) and
the GWAS Catalog FTP directory listings.
"""
from __future__ import annotations

import argparse
import csv
import datetime as _dt
import json
import os
import re
import sys
import urllib.error
import urllib.request
from dataclasses import asdict, dataclass, field
from pathlib import Path

OT_GRAPHQL_URL = "https://api.platform.opentargets.org/api/v4/graphql"
GWAS_FTP_BASE = "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics"
PROBE_TIMEOUT_S = 20
OT_TIMEOUT_S = 60
OT_PAGE_SIZE = 500
DISPLAY_CAP = 15
WINDOW_BP = 1_000_000
EQTL_RELEASE = "r7"

MOLQTL_STUDY_TYPES = ["eqtl", "pqtl", "sqtl", "tuqtl", "sceqtl", "scpqtl", "scsqtl", "sctuqtl"]

_STUDY_LOCUS_RE = re.compile(r"\b([0-9a-f]{32})\b")
_QTD_RE = re.compile(r"QTD\d{6}")
_GCST_RE = re.compile(r"GCST\d+")

# ----- GWAS file check (one GET of the study's harmonised/ listing) -----------------
PROBE_INDEXED = "indexed"
PROBE_NO_INDEX = "no_index"
PROBE_NO_HARMONISED_DIR = "no_harmonised_dir"
PROBE_NOT_GCST = "not_gcst"
PROBE_UNKNOWN = "unknown"
PROBE_NOT_PROBED = "not_probed"
PROBE_LABELS = {
    PROBE_INDEXED: "yes: harmonised file with a tabix index",
    PROBE_NO_INDEX: "no: no indexed harmonised file under the name the fetcher reads",
    PROBE_NO_HARMONISED_DIR: "no: the study has no harmonised folder",
    PROBE_NOT_GCST: "no: not a GWAS Catalog study",
    PROBE_UNKNOWN: "unknown: the check failed",
    PROBE_NOT_PROBED: "not checked",
}
# Only these say the file is definitely not readable by region. `unknown` and
# `not_probed` are not evidence of absence.
PROBE_DEFINITE_NO = frozenset({PROBE_NO_INDEX, PROBE_NO_HARMONISED_DIR, PROBE_NOT_GCST})

# ----- Can the bundled fetchers read the QTL? -----------------------------------------
EXPOSURE_EQTL_CATALOGUE = "eqtl_catalogue"
EXPOSURE_UKB_PPP = "ukb_ppp"
EXPOSURE_NOT_FETCHABLE = "not_fetchable"

# ----- Row states ----------------------------------------------------------------------
STATUS_SHOWN = "shown"
STATUS_HIDDEN = "hidden_by_cap"
EXCLUDE_REASONS = {
    "collapsed": "folded into a shown or hidden row measuring the same thing with a higher H4",
    "is_page_study": "the page's own GWAS study, so not a substitute for it",
}

# ----- QTL kind -------------------------------------------------------------------------
# eQTL Catalogue quantification code -> QTL kind. Must equal `_QUANT_METHOD_MODALITY`
# in locuscompare_region_render.py (a test compares the two).
_QUANT_METHOD_MODALITY: dict[str, str] = {
    "ge": "eqtl",
    "microarray": "eqtl",
    "afc": "eqtl",  # eQTL Catalogue writes "aFC"; codes are matched lower-cased
    "exon": "exon_qtl",
    "tx": "tx_qtl",
    "txrev": "tuqtl",
    "leafcutter": "sqtl",
    "majiq": "sqtl",
    "aptamer": "pqtl",
}
KIND_LABELS = {
    "eqtl": "eQTL (gene)",
    "exon_qtl": "exon eQTL",
    "tx_qtl": "transcript eQTL",
    "tuqtl": "transcript-usage QTL",
    "sqtl": "sQTL",
    "pqtl": "pQTL",
    "sceqtl": "single-cell eQTL",
    "qtl": "QTL (kind not stated)",
}
QUANT_PHRASES = {
    "ge": "expression",
    "microarray": "expression (microarray)",
    "afc": "expression",
    "exon": "exon expression",
    "tx": "transcript expression",
    "txrev": "transcript usage",
    "leafcutter": "splicing",
    "majiq": "splicing",
    "aptamer": "plasma protein (aptamer)",
}


class SelectionError(Exception):
    """A request the helper refuses; the message says why. Exit code 2."""


class PartitionError(RuntimeError):
    """Rows were lost or double-counted between Open Targets and the tables."""


# =========================================================================== network
def _post(query: str, variables: dict) -> dict:
    """POST one GraphQL query to Open Targets and return the parsed response."""
    body = json.dumps({"query": query, "variables": variables}).encode()
    req = urllib.request.Request(
        OT_GRAPHQL_URL, data=body, headers={"Content-Type": "application/json"}
    )
    with urllib.request.urlopen(req, timeout=OT_TIMEOUT_S) as resp:
        return json.loads(resp.read().decode())


def _get(url: str, timeout: float) -> tuple[int, str]:
    """GET a URL; return (HTTP status, body). Network failures raise."""
    try:
        with urllib.request.urlopen(url, timeout=timeout) as resp:
            return resp.status, resp.read().decode("utf-8", "replace")
    except urllib.error.HTTPError as exc:
        return exc.code, ""


_STUDY_FIELDS = """
  id studyType projectId traitFromSource publicationFirstAuthor publicationDate
  nSamples hasSumstats summarystatsLocation condition
  target { id approvedSymbol }
  biosample { biosampleId biosampleName }
  diseases { id name }
"""

CREDIBLE_SET_QUERY = """
query ($id: String!) { credibleSet(studyLocusId: $id) {
  studyLocusId studyType
  variant { id rsIds }
  locus(page: {index: 0, size: 500}) { count rows { variant { id } posteriorProbability } }
  study { %s }
} }""" % _STUDY_FIELDS

COLOC_QUERY = """
query ($id: String!, $types: [StudyTypeEnum!], $index: Int!, $size: Int!) {
 credibleSet(studyLocusId: $id) {
  colocalisation(studyTypes: $types, page: {index: $index, size: $size}) {
    count
    rows { h4 h3 clpp colocalisationMethod numberColocalisingVariants betaRatioSignAverage
      otherStudyLocus { studyLocusId variant { id } study { %s } } }
  } } }""" % _STUDY_FIELDS


class SourceError(RuntimeError):
    """Open Targets could not be read (network, API error, empty answer). Exit code 3."""


def _data(resp: dict) -> dict:
    if resp.get("errors"):
        raise SourceError(f"Open Targets API error: {resp['errors'][0].get('message')}")
    if resp.get("data") is None:
        # An empty answer is not evidence that the credible set does not exist.
        raise SourceError("Open Targets returned no data; try again")
    return resp["data"]


def fetch_credible_set(study_locus_id: str) -> dict:
    cs = _data(_post(CREDIBLE_SET_QUERY, {"id": study_locus_id})).get("credibleSet")
    if not cs:
        raise SelectionError(f"Open Targets has no credible set with id {study_locus_id}")
    return cs


def fetch_colocalisations(study_locus_id: str, study_types: list[str]) -> tuple[int, list[dict]]:
    """All colocalisation rows of one credible set, every page. Returns (count, rows)."""
    rows: list[dict] = []
    count = None
    index = 0
    while True:
        cs = _data(_post(COLOC_QUERY, {"id": study_locus_id, "types": study_types,
                                       "index": index, "size": OT_PAGE_SIZE})).get("credibleSet")
        if not cs:
            raise SelectionError(f"Open Targets has no credible set with id {study_locus_id}")
        page = cs["colocalisation"]
        count = page["count"] if count is None else count
        rows.extend(page["rows"])
        if not page["rows"] or len(rows) >= count:
            break
        index += 1
    return count, rows


# =========================================================================== probes
@dataclass
class Probe:
    status: str
    detail: str = ""


def harmonised_dir_url(accession: str) -> str:
    n = int(accession[4:])
    lo = ((n - 1) // 1000) * 1000 + 1
    hi = lo + 999
    width = len(accession) - 4
    bucket = f"GCST{lo:0{width}d}-GCST{hi:0{width}d}"
    return f"{GWAS_FTP_BASE}/{bucket}/{accession}/harmonised/"


def probe_gwas(accession: str) -> Probe:
    """Is `<accession>.h.tsv.gz` with its `.tbi` in the study's harmonised/ folder?"""
    if not _GCST_RE.fullmatch(accession or ""):
        return Probe(PROBE_NOT_GCST)
    try:
        status, body = _get(harmonised_dir_url(accession), PROBE_TIMEOUT_S)
    except Exception as exc:  # noqa: BLE001 - any failure is "unknown", never "no"
        return Probe(PROBE_UNKNOWN, type(exc).__name__)
    if status == 404:
        return Probe(PROBE_NO_HARMONISED_DIR)
    if status != 200:
        return Probe(PROBE_UNKNOWN, f"HTTP {status}")
    names = set(re.findall(r'href="([^"]+)"', body))
    if f"{accession}.h.tsv.gz" in names and f"{accession}.h.tsv.gz.tbi" in names:
        return Probe(PROBE_INDEXED)
    # A definite "no index" needs a real directory listing. A 200 carrying anything else
    # (an empty body, a throttling or error page) says nothing about the file.
    if not any(n.endswith((".gz", ".tbi", ".yaml", ".txt", ".log")) for n in names):
        return Probe(PROBE_UNKNOWN, "HTTP 200 without a directory listing")
    return Probe(PROBE_NO_INDEX)


# =========================================================================== eQTL index
def _index_candidates() -> list[Path]:
    here = Path(__file__).resolve().parent
    env = os.environ.get("EQTL_CATALOGUE_INDEX")
    return ([Path(env)] if env else []) + [
        here.parent / "eqtl-catalogue-region-fetch" / "data" / "dataset_index_r7.tsv",
        here.parent / "eqtl_catalogue_region_fetch" / "data" / "dataset_index_r7.tsv",
    ]


def load_eqtl_index() -> dict[str, dict]:
    """eQTL Catalogue r7 dataset table, keyed by dataset id (QTD...)."""
    for path in _index_candidates():
        if path.is_file():
            with path.open(newline="") as fh:
                return {r["dataset_id"]: r for r in csv.DictReader(fh, delimiter="\t")}
    raise SelectionError(
        "cannot find the eQTL Catalogue dataset table (dataset_index_r7.tsv) next to this "
        "skill; set EQTL_CATALOGUE_INDEX to its path"
    )


# =========================================================================== rows
@dataclass
class QtlRow:
    study_locus_id: str
    study_id: str
    study_type: str
    project_id: str
    target_id: str
    symbol: str
    biosample_id: str
    tissue: str
    n_samples: int | None
    h4: float
    sign: float | None
    method: str
    location: str
    trait_from_source: str
    quant: str
    kind: str
    qtd: str | None
    index_row: dict | None
    exposure: str
    older_release: bool
    status: str = ""
    row_number: int | None = None
    folded: list[str] = field(default_factory=list)

    @property
    def sort_n(self) -> int:
        if self.index_row and str(self.index_row.get("sample_size", "")).isdigit():
            return int(self.index_row["sample_size"])
        return self.n_samples or 0


@dataclass
class GwasRow:
    study_locus_id: str
    study_id: str
    trait: str
    author: str
    year: str
    n_samples: int | None
    disease_ids: list[str]
    disease_names: list[str]
    h4: float
    sign: float | None
    method: str
    same_disease: bool
    lead: str = ""
    probe: Probe = field(default_factory=lambda: Probe(PROBE_NOT_PROBED))
    status: str = ""
    row_number: int | None = None


def _quant_from_study(study_id: str, index_row: dict | None) -> str:
    """The quantification code, read from the study id's own tokens
    (gtex-v10_ge_lung_..., hipsci_exon_ipsc_...), else from the dataset table."""
    for tok in (study_id or "").lower().split("_"):
        if tok in _QUANT_METHOD_MODALITY:
            return tok
    if index_row:
        return (index_row.get("quant_method") or "").strip().lower()
    return ""


def _kind(quant: str, study_type: str, project_id: str) -> str:
    if project_id.upper().startswith("UKB_PPP"):
        return "pqtl"
    code = _QUANT_METHOD_MODALITY.get(quant, "qtl")
    if code == "eqtl" and (study_type or "").lower().startswith("sc"):
        code = "sceqtl"
    return code


def qtl_row(raw: dict, index: dict[str, dict]) -> QtlRow:
    other = raw["otherStudyLocus"]
    s = other["study"]
    loc = s.get("summarystatsLocation") or ""
    m = _QTD_RE.search(loc)
    qtd = m.group(0) if m else None
    idx = index.get(qtd) if qtd else None
    project = s.get("projectId") or ""
    if project.upper().startswith("UKB_PPP"):
        exposure = EXPOSURE_UKB_PPP
    elif idx is not None:
        exposure = EXPOSURE_EQTL_CATALOGUE
    else:
        exposure = EXPOSURE_NOT_FETCHABLE
    older = bool(idx) and project.lower().startswith("gtex-v10") and idx.get("study_label") == "GTEx"
    quant = _quant_from_study(s["id"], idx)
    target = s.get("target") or {}
    bios = s.get("biosample") or {}
    return QtlRow(
        study_locus_id=other["studyLocusId"], study_id=s["id"], study_type=s.get("studyType") or "",
        project_id=project, target_id=target.get("id") or "", symbol=target.get("approvedSymbol") or "",
        biosample_id=bios.get("biosampleId") or "", tissue=bios.get("biosampleName") or "",
        n_samples=s.get("nSamples"), h4=float(raw.get("h4") or 0.0), sign=raw.get("betaRatioSignAverage"),
        method=raw.get("colocalisationMethod") or "", location=loc,
        trait_from_source=s.get("traitFromSource") or "", quant=quant,
        kind=_kind(quant, s.get("studyType") or "", project), qtd=qtd, index_row=idx,
        exposure=exposure, older_release=older,
    )


def gwas_row(raw: dict, page_disease_ids: set[str]) -> GwasRow:
    other = raw["otherStudyLocus"]
    s = other["study"]
    dis = s.get("diseases") or []
    ids = [d["id"] for d in dis]
    return GwasRow(
        study_locus_id=other["studyLocusId"], study_id=s["id"], trait=s.get("traitFromSource") or "",
        author=s.get("publicationFirstAuthor") or "", year=(s.get("publicationDate") or "")[:4],
        n_samples=s.get("nSamples"), disease_ids=ids, disease_names=[d["name"] for d in dis],
        h4=float(raw.get("h4") or 0.0), sign=raw.get("betaRatioSignAverage"),
        method=raw.get("colocalisationMethod") or "", same_disease=bool(page_disease_ids & set(ids)),
        lead=(other.get("variant") or {}).get("id") or "",
    )


def check_partition(rows: list, count: int, what: str) -> dict:
    """Every row in exactly one state, and all of Open Targets' rows accounted for."""
    tally = {"shown": 0, "hidden_by_cap": 0, "excluded": 0}
    for r in rows:
        if r.status == STATUS_SHOWN:
            tally["shown"] += 1
        elif r.status == STATUS_HIDDEN:
            tally["hidden_by_cap"] += 1
        elif r.status.startswith("excluded:") and r.status[9:] in EXCLUDE_REASONS:
            tally["excluded"] += 1
        else:
            raise PartitionError(f"{what}: row {r.study_id} has no valid state ({r.status!r})")
    if sum(tally.values()) != len(rows) or len(rows) != count:
        raise PartitionError(
            f"{what}: Open Targets reported {count} rows, {len(rows)} were read, and "
            f"{tally} were accounted for"
        )
    return tally


def _cap(ordered: list) -> None:
    for i, r in enumerate(ordered):
        if i < DISPLAY_CAP:
            r.status, r.row_number = STATUS_SHOWN, i + 1
        else:
            r.status = STATUS_HIDDEN


def build_qtl_table(rows: list[QtlRow], gene: str | None = None) -> list[QtlRow]:
    """Fold same-measurement rows, order, cap. Returns the kept rows in display order."""
    groups: dict[tuple, list[QtlRow]] = {}
    for r in rows:
        # An unstated kind is never folded with anything: its study id joins the key.
        key = (r.target_id, r.biosample_id, r.project_id, r.kind,
               r.study_id if r.kind == "qtl" else "")
        groups.setdefault(key, []).append(r)
    kept: list[QtlRow] = []
    for members in groups.values():
        members.sort(key=lambda r: (-r.h4, -r.sort_n, r.study_id))
        best = members[0]
        for other in members[1:]:
            other.status = "excluded:collapsed"
            best.folded.append(other.study_id)
        kept.append(best)
    g = (gene or "").upper()
    kept.sort(key=lambda r: (0 if g and r.symbol.upper() == g else 1, -r.h4, -r.sort_n, r.study_id))
    _cap(kept)
    return kept


def build_gwas_table(rows: list[GwasRow], page_study_id: str) -> list[GwasRow]:
    best: dict[str, GwasRow] = {}
    for r in rows:
        if r.study_id == page_study_id:
            r.status = "excluded:is_page_study"
            continue
        cur = best.get(r.study_id)
        if cur is None or r.h4 > cur.h4:
            if cur is not None:
                cur.status = "excluded:collapsed"
            best[r.study_id] = r
        else:
            r.status = "excluded:collapsed"
    kept = sorted(best.values(),
                  key=lambda r: (0 if r.same_disease else 1, -r.h4, -(r.n_samples or 0), r.study_id))
    _cap(kept)
    return kept


# =========================================================================== selection
@dataclass
class Recommendation:
    row: object | None
    reasons: list[str]


def recommend_qtl(kept: list[QtlRow], gene: str | None, among: list[QtlRow] | None = None) -> Recommendation:
    """Recommend from the shown rows, or from `among` (rows matching what the user asked for)."""
    shown = among if among is not None else [r for r in kept if r.status == STATUS_SHOWN]
    readable = [r for r in shown if r.exposure != EXPOSURE_NOT_FETCHABLE]
    if not readable:
        return Recommendation(None, [
            "none of the QTLs shown can be read by the bundled fetchers; supply the QTL's "
            "summary statistics as a file (sumstats_path) instead"])
    g = (gene or "").upper()
    reasons: list[str] = []
    pool = readable
    if g:
        named = [r for r in readable if r.symbol.upper() == g]
        if named:
            pool = named
            reasons.append(f"it is a {named[0].symbol} QTL, the gene you named")
        else:
            reasons.append(f"no readable QTL for {gene} is listed, so this recommends another gene")
    pick = min(pool, key=lambda r: (-r.h4, -r.sort_n, r.study_id))
    reasons.append(f"highest colocalisation probability (H4 {pick.h4:.4f}) among the "
                   f"{len(pool)} QTLs {'listed' if pool is not readable else 'shown'} that the "
                   f"bundled fetchers can read" + (f" for {pick.symbol}" if pool is not readable else ""))
    if pick.folded:
        reasons.append(f"{len(pick.folded)} other row(s) for the same gene, tissue, study and QTL "
                       f"kind were folded into this one")
    if pick.older_release:
        reasons.append(_older_release_sentence(pick))
    if pick.kind not in ("eqtl", "pqtl", "sceqtl", "qtl"):
        reasons.append(_pooled_feature_sentence(pick))
    return Recommendation(pick, reasons)


def _older_release_sentence(r: QtlRow) -> str:
    idx_n = (r.index_row or {}).get("sample_size", "?")
    return (f"Open Targets computed this colocalisation on GTEx v10 (n={r.n_samples}); the bundled "
            f"fetcher reads GTEx v8 from eQTL Catalogue {EQTL_RELEASE} (n={idx_n}), so the plot "
            f"shows the same gene and tissue from an older release")


def _pooled_feature_sentence(r: QtlRow) -> str:
    return (f"this is a {KIND_LABELS[r.kind]}; the config selects the dataset's rows by gene "
            f"({r.target_id}), which pools every {KIND_LABELS[r.kind]} feature of {r.symbol} in "
            f"the window; the feature Open Targets colocalised is {r.trait_from_source}")


def recommend_gwas(kept: list[GwasRow]) -> Recommendation:
    shown = [r for r in kept if r.status == STATUS_SHOWN]
    same = [r for r in shown if r.same_disease and r.probe.status == PROBE_INDEXED]
    if same:
        pick = same[0]
        return Recommendation(pick, [
            "same Open Targets disease as the page study",
            f"highest H4 ({pick.h4:.4f}) among the same-disease studies whose file is indexed",
        ])
    other = [r for r in shown if r.probe.status == PROBE_INDEXED]
    if other:
        pick = other[0]
        return Recommendation(pick, [
            "none of the same-disease studies shown has an indexed harmonised file",
            f"DIFFERENT TRAIT: {pick.trait} is not the page's disease; highest H4 "
            f"({pick.h4:.4f}) among the indexed studies shown",
        ])
    return Recommendation(None, [
        "no study shown has an indexed harmonised file (an unchecked or failed check is not "
        "counted as readable); supply the GWAS summary statistics as a file (sumstats_path)"])


def select_qtl(token: str, kept: list[QtlRow], all_rows: list[QtlRow]) -> QtlRow:
    token = token.strip()
    if token.isdigit():
        for r in kept:
            if r.status == STATUS_SHOWN and r.row_number == int(token):
                return r
        raise SelectionError(f"no shown QTL row numbered {token}")
    for r in all_rows:
        if token in (r.study_locus_id, r.study_id):
            return r
    raise SelectionError(f"{token} is not one of the QTLs that colocalise with this credible set")


# Words a user may use for a QTL kind, beyond its label and quantification phrase.
_KIND_WORDS = {
    "eqtl": "eqtl expression gene ge",
    "exon_qtl": "eqtl exon",
    "tx_qtl": "eqtl transcript tx",
    "tuqtl": "tuqtl transcript usage txrev",
    "sqtl": "sqtl splicing",
    "pqtl": "pqtl protein",
    "sceqtl": "sceqtl eqtl single-cell single cell",
    "qtl": "qtl",
}


def _is_exact_qtl_ref(token: str) -> bool:
    """A row number or a 32-hex credible-set id: refers to one row or none, never a search."""
    t = token.strip().lower()
    return t.isdigit() or bool(re.fullmatch(r"[0-9a-f]{32}", t))


def match_qtls(text: str, kept: list[QtlRow]) -> list[QtlRow]:
    """Rows whose description contains every word of `text` (case-insensitive), e.g.
    "PCSK9 lung" or "lung exon". Searched over all kept rows, shown or hidden by the cap."""
    words = [w for w in re.split(r"[\s,;]+", text.lower()) if w]
    out = []
    for r in kept:
        hay = " ".join([
            r.symbol, r.target_id, r.tissue, r.biosample_id, r.study_id, r.project_id, r.quant,
            KIND_LABELS[r.kind], QUANT_PHRASES.get(r.quant, ""), _KIND_WORDS.get(r.kind, ""),
        ]).lower()
        if words and all(w in hay for w in words):
            out.append(r)
    return out


def parse_study_locus_id(text: str) -> str:
    m = _STUDY_LOCUS_RE.search((text or "").lower())
    if not m:
        raise SelectionError(f"not a credible-set id or URL: {text!r} (expected 32 hex characters)")
    return m.group(1)


# =========================================================================== analysis
@dataclass
class Analysis:
    study_locus_id: str
    credible_set: dict
    qtl_count: int
    qtl_rows: list[QtlRow]
    qtl_kept: list[QtlRow]
    qtl_tally: dict
    qtl_rec: Recommendation
    selected: QtlRow | None = None
    page_probe: Probe | None = None
    gwas_count: int | None = None
    gwas_rows: list[GwasRow] = field(default_factory=list)
    gwas_kept: list[GwasRow] = field(default_factory=list)
    gwas_tally: dict | None = None
    gwas_rec: Recommendation | None = None
    match_text: str | None = None
    matches: list[QtlRow] = field(default_factory=list)
    match_rec: Recommendation | None = None


def analyse(target: str, gene: str | None = None, qtl: str | None = None, probe: bool = True,
            show_substitutes: bool = False) -> Analysis:
    slid = parse_study_locus_id(target)
    cs = fetch_credible_set(slid)
    if (cs.get("studyType") or "").lower() != "gwas":
        raise SelectionError(
            f"credible set {slid} is a {cs.get('studyType')} signal, not a GWAS one; open the "
            "GWAS credible set it colocalises with and start from there")
    index = load_eqtl_index()
    count, raw = fetch_colocalisations(slid, MOLQTL_STUDY_TYPES)
    rows = [qtl_row(r, index) for r in raw]
    kept = build_qtl_table(rows, gene)
    tally = check_partition(rows, count, "molecular QTL colocalisations")
    a = Analysis(slid, cs, count, rows, kept, tally, recommend_qtl(kept, gene))

    if qtl:
        try:
            a.selected = select_qtl(qtl, kept, rows)
        except SelectionError:
            if _is_exact_qtl_ref(qtl):
                raise
            # Not a row number or id: treat it as a description and narrow to what matches.
            a.match_text = qtl
            a.matches = match_qtls(qtl, kept)
            if not a.matches:
                raise SelectionError(
                    f"no QTL here matches {qtl!r}; run options without --qtl to see every row")
            if len(a.matches) == 1:
                a.selected = a.matches[0]
            else:
                a.match_rec = recommend_qtl(kept, gene, among=a.matches)
    elif len(kept) == 1:
        a.selected = kept[0]
    if a.selected is None:
        return a
    if a.selected.exposure == EXPOSURE_NOT_FETCHABLE:
        # Nothing on the GWAS side is worth checking for a QTL the fetchers cannot read.
        return a

    page = cs["study"]
    a.page_probe = probe_gwas(page["id"]) if probe else Probe(PROBE_NOT_PROBED)
    if a.page_probe.status in PROBE_DEFINITE_NO or show_substitutes:
        page_dis = {d["id"] for d in page.get("diseases") or []}
        gcount, graw = fetch_colocalisations(a.selected.study_locus_id, ["gwas"])
        a.gwas_count = gcount
        a.gwas_rows = [gwas_row(r, page_dis) for r in graw]
        a.gwas_kept = build_gwas_table(a.gwas_rows, page["id"])
        a.gwas_tally = check_partition(a.gwas_rows, gcount, "GWAS colocalisations of the chosen QTL")
        for r in a.gwas_kept:
            if r.status == STATUS_SHOWN:
                r.probe = probe_gwas(r.study_id) if probe else Probe(PROBE_NOT_PROBED)
        a.gwas_rec = recommend_gwas(a.gwas_kept)
    return a


# =========================================================================== printing
def _direction(sign) -> str:
    if sign is None:
        return "?"
    return "same" if sign > 0 else ("opposite" if sign < 0 else "mixed")


def _exposure_cell(r: QtlRow) -> str:
    if r.exposure == EXPOSURE_UKB_PPP:
        return "yes: UKB-PPP"
    if r.exposure == EXPOSURE_EQTL_CATALOGUE:
        cell = f"yes: eQTL Catalogue {r.qtd}"
        if r.older_release:
            cell += f" (older release: GTEx v8, n={r.index_row.get('sample_size')})"
        return cell
    return "no: not in the bundled fetchers"


def _table(header: list[str], body: list[list[str]]) -> str:
    widths = [max(len(str(c)) for c in col) for col in zip(header, *body)] if body else [len(h) for h in header]
    fmt = "  ".join("{:<%d}" % w for w in widths)
    return "\n".join(fmt.format(*row).rstrip() for row in [header] + body)


def _study_desc(page: dict) -> str:
    who = page.get("publicationFirstAuthor") or ""
    yr = (page.get("publicationDate") or "")[:4]
    cite = f" ({who} {yr})".rstrip() if who else ""
    return f"{page['id']}, {page.get('traitFromSource')}{cite}"


def render_text(a: Analysis) -> str:
    cs, page = a.credible_set, a.credible_set["study"]
    var = cs["variant"]
    rs = (var.get("rsIds") or [""])[0]
    out = [
        f"Open Targets credible set {a.study_locus_id}",
        f"  GWAS: {_study_desc(page)}, n={page.get('nSamples')}",
        f"  lead variant: {var['id']}" + (f" ({rs})" if rs else ""),
        "",
        f"Molecular QTLs that colocalise with this signal: {a.qtl_count} listed by Open Targets, "
        f"{a.qtl_tally['shown']} shown.",
    ]
    body = [[str(r.row_number), r.symbol, KIND_LABELS[r.kind], r.tissue, r.project_id,
             f"{r.h4:.4f}", _direction(r.sign), str(r.n_samples), _exposure_cell(r)]
            for r in a.qtl_kept if r.status == STATUS_SHOWN]
    out.append(_table(["#", "gene", "kind", "tissue", "study", "H4", "direction", "n",
                       "readable by the bundled fetchers"], body))
    if a.qtl_tally["excluded"]:
        out.append(f"{a.qtl_tally['excluded']} more row(s) were folded into a row above (same gene, "
                   "tissue, study and QTL kind; the highest H4 is kept).")
    if a.qtl_tally["hidden_by_cap"]:
        out.append(f"{a.qtl_tally['hidden_by_cap']} more not shown (lower in the order); "
                   "--json writes every row.")
    out.append("")
    rec = a.qtl_rec
    if rec.row is not None:
        r = rec.row
        out.append(f"Recommended QTL: row {r.row_number}, {r.symbol} {KIND_LABELS[r.kind]} in "
                   f"{r.tissue} ({r.study_id})")
    else:
        out.append("No QTL to recommend.")
    out += [f"  - {x}" for x in rec.reasons]

    if a.matches and a.selected is None:
        out += ["", f"{len(a.matches)} QTLs match {a.match_text!r}; ask the user which one:"]
        mbody = [[str(r.row_number or "-"), r.symbol, KIND_LABELS[r.kind], r.tissue, r.project_id,
                  f"{r.h4:.4f}", _exposure_cell(r)] for r in a.matches]
        out.append(_table(["#", "gene", "kind", "tissue", "study", "H4",
                           "readable by the bundled fetchers"], mbody))
        rec = a.match_rec
        if rec.row is not None:
            out.append(f"Recommended among these: {rec.row.symbol} {KIND_LABELS[rec.row.kind]} in "
                       f"{rec.row.tissue} ({rec.row.study_id})")
            out += [f"  - {x}" for x in rec.reasons]
    if a.selected is None:
        out += ["", "Next: choose a QTL, then run",
                f"  python ot_credible_set.py options {a.study_locus_id} --qtl <row number or study id>"]
        return "\n".join(out)

    s = a.selected
    out += ["", f"Chosen QTL: {s.symbol} {KIND_LABELS[s.kind]} in {s.tissue} ({s.study_id})"]
    if a.page_probe is None:
        # Only an unreadable QTL stops before the GWAS side is looked at.
        out += [f"{s.study_id} cannot be read by the bundled fetchers, so no plot can be made "
                "from it here; supply its summary statistics as a file (sumstats_path) or "
                "choose another QTL."]
        return "\n".join(out)
    p = a.page_probe
    out.append(f"Page GWAS file check ({page['id']}): {PROBE_LABELS[p.status]}"
               + (f" ({p.detail})" if p.detail else ""))
    if p.status in (PROBE_UNKNOWN, PROBE_NOT_PROBED):
        out.append(f"The page study {page['id']} is kept: its file could not be checked, which is "
                   "not evidence that it is unreadable. The render may still work; if the GWAS "
                   "fetch fails, rerun options with --show-substitutes.")
    if a.gwas_rec is not None:
        if p.status in PROBE_DEFINITE_NO:
            out.append(f"The page study {page['id']} cannot be read by region, so these GWAS studies "
                       "that colocalise with the chosen QTL could stand in for it:")
        else:
            out.append("GWAS studies that colocalise with the chosen QTL (shown on request):")
        out.append(f"{a.gwas_count} listed by Open Targets, {a.gwas_tally['shown']} shown.")
        gbody = [[str(r.row_number), r.study_id, r.trait[:45], f"{r.author} {r.year}".strip(),
                  "same disease" if r.same_disease else "different trait", f"{r.h4:.4f}",
                  _direction(r.sign), str(r.n_samples),
                  PROBE_LABELS[r.probe.status] + (f" ({r.probe.detail})" if r.probe.detail else "")]
                 for r in a.gwas_kept if r.status == STATUS_SHOWN]
        out.append(_table(["#", "study", "trait", "first author", "vs page", "H4", "direction", "n",
                           "readable by region"], gbody))
        n_page = sum(1 for r in a.gwas_rows if r.status == "excluded:is_page_study")
        n_fold = sum(1 for r in a.gwas_rows if r.status == "excluded:collapsed")
        if n_page:
            out.append(f"{n_page} row(s) are the page's own study and are not offered.")
        if n_fold:
            out.append(f"{n_fold} row(s) repeat a study above with a lower H4 and were folded.")
        if a.gwas_tally["hidden_by_cap"]:
            out.append(f"{a.gwas_tally['hidden_by_cap']} more not shown (lower in the order, not "
                       "checked); --json writes every row.")
        g = a.gwas_rec
        if g.row is not None:
            out.append(f"Recommended GWAS: {g.row.study_id}, {g.row.trait}"
                       + (f" ({g.row.author} {g.row.year})" if g.row.author else ""))
            out += [f"  - {x}" for x in g.reasons]
            out.append(f"  - substitution: the plot would show {g.row.study_id}, not the page "
                       f"study {page['id']}")
        else:
            out.append("No GWAS to recommend.")
            out += [f"  - {x}" for x in g.reasons]
    gw = ""
    if a.gwas_rec and a.gwas_rec.row is not None and p.status in PROBE_DEFINITE_NO:
        gw = f" --gwas {a.gwas_rec.row.study_id}"
    out += ["", "Next: confirm the pair with the user, then run",
            f"  python ot_credible_set.py config {a.study_locus_id} --qtl {s.study_locus_id}{gw} "
            "--out config.yaml"]
    return "\n".join(out)


def to_json(a: Analysis) -> dict:
    def q(r: QtlRow) -> dict:
        d = asdict(r)
        d.pop("index_row")
        d["kind_label"] = KIND_LABELS[r.kind]
        return d

    def g(r: GwasRow) -> dict:
        d = asdict(r)
        d["probe_label"] = PROBE_LABELS[r.probe.status]
        return d

    return {
        "study_locus_id": a.study_locus_id,
        "credible_set": a.credible_set,
        "qtl": {"count": a.qtl_count, "tally": a.qtl_tally, "rows": [q(r) for r in a.qtl_rows],
                "recommended": a.qtl_rec.row.study_locus_id if a.qtl_rec.row else None,
                "reasons": a.qtl_rec.reasons},
        "selected_qtl": a.selected.study_locus_id if a.selected else None,
        "page_probe": asdict(a.page_probe) if a.page_probe else None,
        "gwas": None if a.gwas_rec is None else {
            "count": a.gwas_count, "tally": a.gwas_tally, "rows": [g(r) for r in a.gwas_rows],
            "recommended": a.gwas_rec.row.study_id if a.gwas_rec.row else None,
            "reasons": a.gwas_rec.reasons},
    }


# =========================================================================== config
def _y(v) -> str:
    """One YAML scalar. JSON strings are valid YAML double-quoted scalars."""
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, int):
        return str(v)
    if isinstance(v, float):
        # YAML 1.1 reads "1e-05" as a string; a float needs a dot before the exponent.
        t = repr(v)
        return t.replace("e", ".0e") if "e" in t and "." not in t else t
    return json.dumps("" if v is None else str(v))


def _exposure_label(r: QtlRow) -> str:
    phrase = QUANT_PHRASES.get(r.quant, "QTL")
    if r.exposure == EXPOSURE_UKB_PPP:
        return f"{r.symbol} plasma protein (UKB-PPP {r.project_id.split('_')[-1]})"
    idx = r.index_row or {}
    tissue = idx.get("tissue_label") or r.tissue
    study = (idx.get("study_label") or r.project_id).replace("_", " ")
    if study == "GTEx":
        study = "GTEx v8"
    return f"{r.symbol} {phrase} - {tissue} ({study} via eQTL Catalogue {EQTL_RELEASE})"


def _page_members(cs: dict) -> list[str]:
    return [r["variant"]["id"] for r in ((cs.get("locus") or {}).get("rows") or [])]


def build_config(a: Analysis, gwas: str | None, today: str | None = None,
                 lead: str | None = None) -> str:
    s = a.selected
    if s is None and a.matches:
        raise SelectionError(
            f"{len(a.matches)} QTLs match {a.match_text!r} ("
            + ", ".join(r.study_id for r in a.matches)
            + "); ask the user which one and pass its study id with --qtl")
    if s is None:
        raise SelectionError("choose a QTL with --qtl (a row number from options, or a study id)")
    if s.exposure == EXPOSURE_NOT_FETCHABLE:
        raise SelectionError(f"{s.study_id} cannot be read by the bundled fetchers; supply its "
                             "summary statistics as a file (sumstats_path) instead")
    cs, page = a.credible_set, a.credible_set["study"]
    substitutions: list[dict] = []
    caveats: list[str] = []
    p = a.page_probe or Probe(PROBE_NOT_PROBED)

    if gwas and gwas != page["id"]:
        # A study can colocalise through more than one credible set: keep its best H4,
        # the row the substitutes table shows.
        allowed = {}
        for r in sorted(a.gwas_rows, key=lambda r: r.h4):
            allowed[r.study_id] = r
        if not allowed:
            count, raw = fetch_colocalisations(s.study_locus_id, ["gwas"])
            page_dis = {d["id"] for d in page.get("diseases") or []}
            allowed = {r.study_id: r for r in (gwas_row(x, page_dis) for x in raw)}
        if gwas not in allowed:
            raise SelectionError(f"{gwas} does not colocalise with {s.study_id} in Open Targets; "
                                 "pick one from the substitutes table")
        g = allowed[gwas]
        gp = probe_gwas(gwas) if p.status != PROBE_NOT_PROBED else Probe(PROBE_NOT_PROBED)
        if gp.status in PROBE_DEFINITE_NO:
            raise SelectionError(f"{gwas} cannot be read by region ({PROBE_LABELS[gp.status]})")
        outcome_check = gp.status
        if gp.status in (PROBE_UNKNOWN, PROBE_NOT_PROBED):
            caveats.append(
                f"Whether {gwas} can be read by region is unknown (file check: "
                f"{PROBE_LABELS[gp.status]}{'; ' + gp.detail if gp.detail else ''}).")
        reason = (f"the page study cannot be read by region ({PROBE_LABELS[p.status][4:]})"
                  if p.status in PROBE_DEFINITE_NO else "chosen by the user")
        substitutions.append({"side": "outcome", "requested": page["id"], "used": gwas,
                              "reason": reason})
        caveats.append(
            f"The outcome is {gwas} ({g.trait}), not the page study {page['id']} "
            f"({page.get('traitFromSource')}): {reason}."
            + ("" if g.same_disease else " It is a different trait from the page's disease."))
        outcome = {"id": gwas, "trait": g.trait, "h4": g.h4, "sign": g.sign,
                   "study_locus_id": g.study_locus_id}
    else:
        if p.status in PROBE_DEFINITE_NO:
            raise SelectionError(
                f"the page study {page['id']} cannot be read by region "
                f"({PROBE_LABELS[p.status]}); pass --gwas with a substitute from "
                f"`options {a.study_locus_id} --qtl {s.study_locus_id}`")
        if p.status in (PROBE_UNKNOWN, PROBE_NOT_PROBED):
            caveats.append(
                f"Whether the page study {page['id']} can be read by region is unknown (file check: "
                f"{PROBE_LABELS[p.status]}{'; ' + p.detail if p.detail else ''}); it was kept.")
        outcome_check = p.status
        outcome = {"id": page["id"], "trait": page.get("traitFromSource"), "h4": s.h4,
                   "sign": s.sign, "study_locus_id": a.study_locus_id}
    if not _GCST_RE.fullmatch(outcome["id"]):
        raise SelectionError(f"{outcome['id']} is not a GWAS Catalog study; the bundled fetcher reads "
                             "GWAS Catalog accessions only")
    if s.older_release:
        substitutions.append({"side": "exposure", "requested": s.study_id,
                              "used": f"{s.qtd} (eQTL Catalogue {EQTL_RELEASE}, GTEx v8)",
                              "reason": "the bundled fetcher reads eQTL Catalogue r7, which holds GTEx v8"})
        caveats.append(_older_release_sentence(s)[0].upper() + _older_release_sentence(s)[1:] + ".")
    if s.kind not in ("eqtl", "pqtl", "sceqtl", "qtl"):
        caveats.append(_pooled_feature_sentence(s)[0].upper() + _pooled_feature_sentence(s)[1:] + ".")

    # The lead. Default: the page's lead; `--lead` may name any member of the page's
    # credible set (a recorded substitution). Whether a lead is usable (in both datasets
    # and the LD panel) is only known at render time, so the members go into the config
    # as `lead.candidates` (PIP high to low) and the renderer falls back to the first
    # usable one, stating the swap. Delete `candidates` from the config to pin the lead.
    var = cs["variant"]
    lead_id, lead_reason = var["id"], ""
    members = _page_members(cs)
    if lead:
        if lead not in members:
            raise SelectionError(f"{lead} is not in the page's credible set ({', '.join(members)})")
        lead_id, lead_reason = lead, "chosen by the user from the page's credible set"
    locus = cs.get("locus") or {}
    cand_rows = sorted(
        (r for r in (locus.get("rows") or []) if r.get("posteriorProbability") is not None),
        key=lambda r: -float(r["posteriorProbability"]))
    if locus.get("count") is not None and int(locus["count"]) > len(locus.get("rows") or []):
        caveats.append(
            f"Open Targets lists {locus['count']} members of this credible set; the "
            f"{len(locus.get('rows') or [])} it returned are the lead candidates.")
    if lead_id != var["id"]:
        substitutions.append({"side": "lead", "requested": var["id"], "used": lead_id,
                              "reason": lead_reason})
        caveats.append(f"The lead requested is {lead_id}, not the page's lead {var['id']}: "
                       f"{lead_reason}. The render may still use another credible-set member "
                       f"if {lead_id} turns out not to be usable, and says so if it does.")
    chrom, pos = lead_id.split("_")[:2]
    rs = (var.get("rsIds") or [None])[0] if lead_id == var["id"] else None
    today = today or _dt.date.today().isoformat()
    L = ["schema_version: \"1.0\"", "",
         f"# Written by ot_credible_set.py from https://platform.opentargets.org/credible-set/{a.study_locus_id}",
         "", "lead:", f"  variant_id: {_y(lead_id)}"]
    if rs:
        L.append(f"  rs_id: {_y(rs)}")
    L += [f"  chromosome: {_y(chrom)}", f"  position_bp: {int(pos)}", f"  window_bp: {WINDOW_BP}"]
    if cand_rows:
        L.append("  candidates:  # credible-set members, PIP high to low; delete to pin the lead")
        L += [f"    - {{variant_id: {_y(r['variant']['id'])}, pip: {float(r['posteriorProbability'])}}}"
              for r in cand_rows]
    L += ["",
          "exposure:", f"  trait_label: {_y(_exposure_label(s))}", f"  gene_symbol: {_y(s.symbol)}",
          "  fetch:"]
    if s.exposure == EXPOSURE_UKB_PPP:
        L += ["    source: ukb_ppp", f"    protein_label: {_y(s.symbol)}",
              f"    ancestry: {_y(s.project_id.split('_')[-1])}"]
    else:
        L += ["    source: eqtl_catalogue", f"    dataset_id: {s.qtd}",
              f"    molecular_trait_id: {s.target_id}"]
    L += ["", "outcome:", f"  trait_label: {_y(outcome['trait'])}", "  fetch:",
          "    source: gwas_catalog", f"    accession: {outcome['id']}", "",
          "ld:", "  source: 1000g_phase3_grch38", "  super_pop: EUR", "",
          "gene_track:", "  source: gencode_v39", "",
          "harmonisation:", "  exclude_palindromic: true", "",
          "provenance:", "  source: Open Targets Platform", f"  retrieved: {_y(today)}",
          f"  credible_set: {_y(a.study_locus_id)}", f"  page_study: {_y(page['id'])}",
          f"  page_study_file_check: {_y(p.status)}",
          f"  outcome_study_file_check: {_y(outcome_check)}",
          f"  qtl_study: {_y(s.study_id)}", f"  qtl_credible_set: {_y(s.study_locus_id)}",
          f"  qtl_kind: {_y(s.kind)}",
          f"  coloc_qtl_vs_page_h4: {s.h4}", f"  coloc_qtl_vs_page_beta_ratio_sign: {_y(s.sign)}",
          f"  outcome_study: {_y(outcome['id'])}", f"  coloc_qtl_vs_outcome_h4: {outcome['h4']}",
          f"  coloc_qtl_vs_outcome_beta_ratio_sign: {_y(outcome['sign'])}"]
    if substitutions:
        L.append("  substitutions:")
        for sub in substitutions:
            L.append(f"    - side: {_y(sub['side'])}")
            L += [f"      {k}: {_y(sub[k])}" for k in ("requested", "used", "reason")]
    else:
        L.append("  substitutions: []")
    if caveats:
        L += ["", "caveats:"] + [f"  - {_y(c)}" for c in caveats]
    return "\n".join(L) + "\n"


# =========================================================================== CLI
def _parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog="ot_credible_set.py",
        description="Choose the QTL and GWAS for a locuscompare run from an Open Targets credible set.")
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("options", "config"):
        p = sub.add_parser(name)
        p.add_argument("credible_set", help="studyLocusId (32 hex characters) or the page URL")
        p.add_argument("--gene", help="list this gene's QTLs first and prefer it")
        p.add_argument("--qtl", help="row number from options, a QTL studyId / studyLocusId, or a "
                                     "description such as 'PCSK9 lung' (several matches are listed)")
        p.add_argument("--no-probe", action="store_true", help="skip the GWAS file checks")
        p.add_argument("--show-substitutes", action="store_true",
                       help="list substitute GWAS studies even when the page study is not known to be unreadable")
        if name == "options":
            p.add_argument("--json", help="also write every row and its state to this file")
        else:
            p.add_argument("--gwas", help="GWAS Catalog accession to use as the outcome")
            p.add_argument("--lead", help="centre the plot on this member of the page's credible set")
            p.add_argument("--out", required=True, help="where to write config.yaml")
    return ap


class OutputError(RuntimeError):
    """The output file could not be written. Exit code 4."""


def _write(path: Path, text: str) -> None:
    try:
        path.write_text(text)
    except OSError as exc:  # not an Open Targets failure; say what it is
        raise OutputError(f"could not write {path}: {exc}") from exc


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        if args.cmd == "config" and not args.qtl:
            raise SelectionError("config needs --qtl: run options first and choose a QTL with the user")
        a = analyse(args.credible_set, gene=args.gene, qtl=args.qtl, probe=not args.no_probe,
                    show_substitutes=args.show_substitutes)
        if args.cmd == "options":
            print(render_text(a))
            if args.json:
                _write(Path(args.json), json.dumps(to_json(a), indent=1, default=str))
            return 0
        text = build_config(a, args.gwas, lead=args.lead)
        _write(Path(args.out), text)
        print(f"wrote {args.out}")
        for line in text.split("caveats:\n", 1)[1:]:
            print("caveats:\n" + line.rstrip())
        return 0
    except SelectionError as exc:
        print(f"refused: {exc}", file=sys.stderr)
        return 2
    except OutputError as exc:
        print(str(exc), file=sys.stderr)
        return 4
    except (SourceError, OSError) as exc:  # URLError and timeouts are OSErrors
        print(f"could not read Open Targets: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 3
    except PartitionError as exc:
        print(f"internal check failed: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
