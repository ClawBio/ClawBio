#!/usr/bin/env python3
"""Fetch a single named API source for a resolved variant.

Cache-first: if {skill_dir}/data/demo_{rsid}.json exists and contains
api_results.{api}, write that. Otherwise hit the live API.

Network paths are best-effort; on failure they emit a stub
{"source": api, "status": "error", "message": ...} so downstream merge
still has a record. This keeps the DAG honest without making any single
source a hard failure.
"""
from __future__ import annotations

import argparse
import json
import sys
import urllib.error
import urllib.request
import urllib.parse
from pathlib import Path
from typing import Any


VALID_APIS = [
    "gwas_catalog",
    "open_targets",
    "open_targets_credsets",
    "pheweb_ukb",
    "finngen",
    "pheweb_bbj",
    "gtex",
    "eqtl_catalogue",
]


def _http_get_json(url: str, timeout: float = 15.0) -> Any:
    req = urllib.request.Request(url, headers={"Accept": "application/json", "User-Agent": "clawbio2-gwas-lookup"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.load(r)


def _http_post_json(url: str, body: dict, timeout: float = 15.0) -> Any:
    data = json.dumps(body).encode()
    req = urllib.request.Request(
        url,
        data=data,
        headers={"Accept": "application/json", "Content-Type": "application/json", "User-Agent": "clawbio2-gwas-lookup"},
    )
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.load(r)


def from_cache(skill_dir: Path, rsid: str, api: str) -> dict | None:
    cache = skill_dir / "data" / f"demo_{rsid}.json"
    if not cache.is_file():
        return None
    bundle = json.loads(cache.read_text())
    return bundle.get("api_results", {}).get(api)


# ---- Live fetchers (best-effort; minimal — return raw payload under `result`) ----


def fetch_gwas_catalog(variant: dict) -> dict:
    rsid = variant["rsid"]
    url = f"https://www.ebi.ac.uk/gwas/rest/api/singleNucleotidePolymorphisms/{rsid}/associations"
    data = _http_get_json(url)
    assocs = data.get("_embedded", {}).get("associations", [])

    def _pval(a: dict) -> float:
        v = a.get("pvalue") or a.get("pValue")
        try:
            return float(v) if v is not None else 1.0
        except (TypeError, ValueError):
            return 1.0

    top = sorted(assocs, key=_pval)[:10]

    # The renderer reads `traits` as a string list but the GWAS REST API
    # returns trait names only via _links.efoTraits.href — resolve here so
    # live mode emits the same shape as the cached fixture.
    trait_cache: dict[str, list[str]] = {}
    for a in top:
        href = a.get("_links", {}).get("efoTraits", {}).get("href")
        if not href:
            a["traits"] = []
            continue
        if href not in trait_cache:
            try:
                tdata = _http_get_json(href)
                trait_cache[href] = [
                    t.get("trait") for t in tdata.get("_embedded", {}).get("efoTraits", [])
                    if t.get("trait")
                ]
            except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, json.JSONDecodeError):
                trait_cache[href] = []
        a["traits"] = trait_cache[href]

    return {
        "source": "gwas_catalog",
        "status": "ok",
        "rsid": rsid,
        "total_associations": len(assocs),
        "raw": top + assocs[10:50],
    }


_OPEN_TARGETS_DEPRECATION = (
    "Open Targets Genetics was sunset in 2024; migration to the merged "
    "Open Targets platform is pending"
)


def fetch_open_targets(variant: dict) -> dict:
    return {
        "source": "open_targets",
        "status": "skipped",
        "rsid": variant.get("rsid", ""),
        "message": _OPEN_TARGETS_DEPRECATION,
    }


def fetch_open_targets_credsets(variant: dict) -> dict:
    return {
        "source": "open_targets_credsets",
        "status": "skipped",
        "rsid": variant.get("rsid", ""),
        "message": _OPEN_TARGETS_DEPRECATION,
    }


def fetch_pheweb_ukb(variant: dict) -> dict:
    vid = variant.get("variant_ids", {}).get("pheweb", "")
    if not vid:
        return {"source": "pheweb_ukb", "status": "skipped", "message": "no pheweb id"}
    url = f"https://pheweb.org/UKB-TOPMed/api/variant/{urllib.parse.quote(vid)}"
    data = _http_get_json(url)
    return {"source": "pheweb_ukb", "status": "ok", "rsid": variant["rsid"], "data": data}


def fetch_finngen(variant: dict) -> dict:
    # FinnGen r13+ gates /api/variant behind Google account auth (/get_authorized
    # redirect); the public r12 freeze is offline. Without credentials there is
    # no usable live path, so emit the same skipped-stub shape as Open Targets.
    return {
        "source": "finngen",
        "status": "skipped",
        "rsid": variant.get("rsid", ""),
        "message": "FinnGen API now requires authenticated access (r13+ /get_authorized); public anonymous variant lookup is unavailable",
    }


def fetch_pheweb_bbj(variant: dict) -> dict:
    vid = variant.get("variant_ids", {}).get("bbj", "")
    if not vid:
        return {"source": "pheweb_bbj", "status": "skipped", "message": "no bbj id"}
    url = f"https://pheweb.jp/api/variant/{urllib.parse.quote(vid)}"
    data = _http_get_json(url)
    return {"source": "pheweb_bbj", "status": "ok", "rsid": variant["rsid"], "data": data}


def fetch_gtex(variant: dict) -> dict:
    vid = variant.get("variant_ids", {}).get("gtex", "")
    if not vid:
        return {"source": "gtex", "status": "skipped", "message": "no gtex variant id"}
    url = f"https://gtexportal.org/api/v2/association/singleTissueEqtl?variantId={urllib.parse.quote(vid)}"
    data = _http_get_json(url)
    return {"source": "gtex", "status": "ok", "rsid": variant["rsid"], "data": data}


def fetch_eqtl_catalogue(variant: dict) -> dict:
    rsid = variant.get("rsid", "")
    if not rsid:
        return {"source": "eqtl_catalogue", "status": "skipped", "message": "no rsid"}
    url = f"https://www.ebi.ac.uk/eqtl/api/v3/associations?rsid={urllib.parse.quote(rsid)}&size=50"
    data = _http_get_json(url)
    return {"source": "eqtl_catalogue", "status": "ok", "rsid": rsid, "data": data}


FETCHERS = {
    "gwas_catalog": fetch_gwas_catalog,
    "open_targets": fetch_open_targets,
    "open_targets_credsets": fetch_open_targets_credsets,
    "pheweb_ukb": fetch_pheweb_ukb,
    "finngen": fetch_finngen,
    "pheweb_bbj": fetch_pheweb_bbj,
    "gtex": fetch_gtex,
    "eqtl_catalogue": fetch_eqtl_catalogue,
}


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--api", required=True, choices=VALID_APIS)
    p.add_argument("--variant", required=True, type=Path, help="path to variant.json (output of resolve step)")
    p.add_argument("--skill-dir", required=True, type=Path)
    p.add_argument("--out", required=True, type=Path)
    args = p.parse_args()

    variant = json.loads(args.variant.read_text())
    rsid = variant.get("rsid", "")

    cached = from_cache(args.skill_dir, rsid, args.api)
    if cached is not None:
        result = cached
        source = "cache"
    else:
        try:
            result = FETCHERS[args.api](variant)
            source = "live"
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, json.JSONDecodeError) as e:
            result = {"source": args.api, "status": "error", "rsid": rsid, "message": f"{type(e).__name__}: {e}"}
            source = "error_stub"

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2))
    status = result.get("status", "?")
    print(f"[fetch] api={args.api} rsid={rsid} source={source} status={status}")


if __name__ == "__main__":
    main()
