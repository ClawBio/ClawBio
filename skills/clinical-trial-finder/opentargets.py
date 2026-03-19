#!/usr/bin/env python3
"""OpenTargets Platform API client — gene-to-disease associations.

GraphQL endpoint: https://api.platform.opentargets.org/api/v4/graphql
No authentication required. Public API, no rate limiting for reasonable use.

Domain decisions (documented in SKILL.md):
- Score threshold ≥ 0.6 to exclude speculative associations
- Max 5 diseases per gene to keep trial queries focused
- MONDO/Orphanet IDs used as-is; no mapping to MeSH (OpenTargets handles this)
"""

import json
import urllib.error
import urllib.request
from typing import NamedTuple

OT_GRAPHQL = "https://api.platform.opentargets.org/api/v4/graphql"

DEFAULT_MIN_SCORE = 0.6
DEFAULT_MAX_DISEASES = 5


class Disease(NamedTuple):
    id: str
    name: str
    score: float


def resolve_gene(symbol: str) -> tuple[str, str]:
    """Resolve a gene symbol to (ensembl_id, approved_name).

    Raises ValueError if the symbol is not found in OpenTargets.
    """
    query = """
    query ResolveGene($q: String!) {
      search(queryString: $q, entityNames: ["target"], page: {index: 0, size: 1}) {
        hits {
          id
          object {
            ... on Target {
              approvedSymbol
              approvedName
            }
          }
        }
      }
    }
    """
    data = _graphql(query, {"q": symbol})
    hits = data.get("search", {}).get("hits", [])
    if not hits:
        raise ValueError(f"Gene not found in OpenTargets: {symbol!r}")
    hit = hits[0]
    return hit["id"], hit["object"]["approvedName"]


def get_diseases(
    ensembl_id: str,
    min_score: float = DEFAULT_MIN_SCORE,
    max_results: int = DEFAULT_MAX_DISEASES,
) -> list[Disease]:
    """Return top diseases associated with a gene, filtered by score.

    Fetches 3x max_results and filters by min_score so the final list
    respects both the score threshold and the count limit.
    """
    query = """
    query GeneDiseases($id: String!, $size: Int!) {
      target(ensemblId: $id) {
        associatedDiseases(page: {index: 0, size: $size}) {
          rows {
            disease { id name }
            score
          }
        }
      }
    }
    """
    data = _graphql(query, {"id": ensembl_id, "size": max_results * 3})
    rows = data.get("target", {}).get("associatedDiseases", {}).get("rows", [])
    diseases = [
        Disease(id=r["disease"]["id"], name=r["disease"]["name"], score=r["score"])
        for r in rows
        if r["score"] >= min_score
    ]
    return diseases[:max_results]


def _graphql(query: str, variables: dict) -> dict:
    """POST a GraphQL query to OpenTargets. Returns the data payload."""
    payload = json.dumps({"query": query, "variables": variables}).encode()
    req = urllib.request.Request(
        OT_GRAPHQL,
        data=payload,
        headers={"Content-Type": "application/json"},
        method="POST",
    )
    try:
        with urllib.request.urlopen(req, timeout=15) as resp:
            result = json.loads(resp.read().decode())
    except urllib.error.URLError as exc:
        raise RuntimeError(f"OpenTargets Platform unreachable: {exc}") from exc
    except json.JSONDecodeError as exc:
        raise RuntimeError(f"OpenTargets returned malformed JSON: {exc}") from exc

    if "errors" in result:
        msgs = [e.get("message", "unknown") for e in result["errors"]]
        raise RuntimeError(f"OpenTargets GraphQL error: {'; '.join(msgs)}")

    return result.get("data", {})
