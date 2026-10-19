"""rsID -> Open Targets variant ids (GRCh38 `chr_pos_ref_alt`), and the weight that ranks
the alleles one rsID names.

Two GraphQL queries against the Open Targets Platform API. The effect allele of this skill
is DEFINED as the Open Targets `alt` of the chosen variant id, so the resolver is part of the
method, not a convenience: another resolver (Ensembl `/variation`, dbSNP) could pick a
different allele at a multi-allelic site and silently change every sign in the output.

* `search_variant_ids`: every variant id the search returns for a query, in its own order,
  so the caller can say which allele it answered for when an rsID names several variants.
* `gwas_weight`: (GWAS credible sets, highest gnomAD allele frequency) for one variant id,
  which ranks the alleles at one position. The search lists rs7903146's rare C>G (no GWAS
  credible set) before its common C>T (888 credible sets, frequency 0.29 in non-Finnish
  Europeans; both read live 2026-10-19), and reading the first id made every outcome in
  every biobank an allele mismatch.
"""

from __future__ import annotations

import time
from typing import Any

import requests

OT_GRAPHQL = "https://api.platform.opentargets.org/api/v4/graphql"
DEFAULT_TIMEOUT_S = 60.0
RETRY_ATTEMPTS = 3
RETRY_BACKOFF_S = 0.5
GWAS_STUDY_TYPE = "gwas"

VARIANT_SEARCH_QUERY = """
query VariantSearch($queryString: String!, $size: Int!) {
  search(queryString: $queryString, entityNames: ["variant"], page: {index: 0, size: $size}) {
    hits { id entity }
  }
}
"""

VARIANT_GWAS_WEIGHT_QUERY = """
query VariantGwasWeight($variantId: String!, $studyTypes: [StudyTypeEnum!]) {
  variant(variantId: $variantId) {
    id
    alleleFrequencies { populationName alleleFrequency }
    credibleSets(studyTypes: $studyTypes, page: { index: 0, size: 1 }) { count }
  }
}
"""


class OpenTargetsError(RuntimeError):
    """The Open Targets API could not answer (network, HTTP error, GraphQL error, bad body)."""


def _post(query: str, variables: dict[str, Any], *, endpoint: str = OT_GRAPHQL,
          session: requests.Session | None = None,
          timeout_s: float = DEFAULT_TIMEOUT_S) -> dict[str, Any]:
    """The reply's `data`, or `OpenTargetsError`. A reply with no data (a GraphQL failure
    can arrive as HTTP 200 with `"data": null`) is an outage, never an empty answer."""
    sess = session or requests.Session()
    body = {"query": query, "variables": variables}
    last: Exception | None = None
    for i in range(RETRY_ATTEMPTS):
        if i:
            time.sleep(RETRY_BACKOFF_S * (2 ** (i - 1)))
        try:
            r = sess.post(endpoint, json=body, timeout=timeout_s,
                          headers={"Content-Type": "application/json"})
        except requests.RequestException as e:
            last = e
            continue
        if r.status_code in (429, 500, 502, 503, 504):
            last = OpenTargetsError(f"HTTP {r.status_code}")
            continue
        if r.status_code != 200:
            raise OpenTargetsError(f"HTTP {r.status_code}")
        try:
            payload = r.json()
        except ValueError as e:
            raise OpenTargetsError("non-JSON body") from e
        if payload.get("errors"):
            raise OpenTargetsError("GraphQL error: " + str(payload["errors"])[:200])
        if payload.get("data") is None:
            raise OpenTargetsError("the reply carried no data")
        return payload["data"]
    raise OpenTargetsError(type(last).__name__ if isinstance(last, requests.RequestException)
                           else str(last))


def search_variant_ids(query: str, *, size: int = 10, endpoint: str = OT_GRAPHQL,
                       session: requests.Session | None = None,
                       timeout_s: float = DEFAULT_TIMEOUT_S) -> list[str]:
    """Every variant id the Open Targets search returns for `query`, in its order.

    An empty list means Open Targets answered and holds no such variant; a failure to
    answer raises `OpenTargetsError` (an outage is never reported as "no variant").
    """
    data = _post(VARIANT_SEARCH_QUERY,
                 {"queryString": (query or "").strip(), "size": size},
                 endpoint=endpoint, session=session, timeout_s=timeout_s)
    hits = ((data.get("search") or {}).get("hits")) or []
    return [h["id"] for h in hits if h.get("entity") == "variant" and h.get("id")]


def gwas_weight(variant_id: str, *, endpoint: str = OT_GRAPHQL,
                session: requests.Session | None = None,
                timeout_s: float = DEFAULT_TIMEOUT_S) -> tuple[int, float]:
    """(GWAS credible sets, highest gnomAD allele frequency) for one variant id; (0, 0.0)
    for an id Open Targets does not know. A failure to answer raises `OpenTargetsError`."""
    data = _post(VARIANT_GWAS_WEIGHT_QUERY,
                 {"variantId": variant_id, "studyTypes": [GWAS_STUDY_TYPE]},
                 endpoint=endpoint, session=session, timeout_s=timeout_s)
    v = data.get("variant")
    if not v:
        return 0, 0.0
    n = int(((v.get("credibleSets") or {}).get("count")) or 0)
    af = max((float(a.get("alleleFrequency") or 0.0)
              for a in (v.get("alleleFrequencies") or [])), default=0.0)
    return n, af
