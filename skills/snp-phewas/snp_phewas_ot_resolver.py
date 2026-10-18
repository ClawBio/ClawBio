"""rsID -> Open Targets variant ids (GRCh38 `chr_pos_ref_alt`), in Open Targets' rank order.

One GraphQL `search` query against the Open Targets Platform API. The effect allele of
this skill is DEFINED as the Open Targets `alt` of the top hit, so the resolver is part of
the method, not a convenience: another resolver (Ensembl `/variation`, dbSNP) could pick a
different allele at a multi-allelic site and silently change every sign in the output.

All hits are returned so the caller can say which allele it answered for when an rsID
names several variants at one position.
"""

from __future__ import annotations

import time

import requests

OT_GRAPHQL = "https://api.platform.opentargets.org/api/v4/graphql"
DEFAULT_TIMEOUT_S = 60.0
RETRY_ATTEMPTS = 3
RETRY_BACKOFF_S = 0.5

VARIANT_SEARCH_QUERY = """
query VariantSearch($queryString: String!, $size: Int!) {
  search(queryString: $queryString, entityNames: ["variant"], page: {index: 0, size: $size}) {
    hits { id entity }
  }
}
"""


class OpenTargetsError(RuntimeError):
    """The Open Targets API could not answer (network, HTTP error, GraphQL error, bad body)."""


def search_variant_ids(query: str, *, size: int = 10, endpoint: str = OT_GRAPHQL,
                       session: requests.Session | None = None,
                       timeout_s: float = DEFAULT_TIMEOUT_S) -> list[str]:
    """Every variant id the Open Targets search returns for `query`, best first.

    An empty list means Open Targets answered and holds no such variant; a failure to
    answer raises `OpenTargetsError` (an outage is never reported as "no variant").
    """
    sess = session or requests.Session()
    body = {"query": VARIANT_SEARCH_QUERY,
            "variables": {"queryString": (query or "").strip(), "size": size}}
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
        hits = (((payload.get("data") or {}).get("search") or {}).get("hits")) or []
        return [h["id"] for h in hits if h.get("entity") == "variant" and h.get("id")]
    raise OpenTargetsError(type(last).__name__ if isinstance(last, requests.RequestException)
                           else str(last))
