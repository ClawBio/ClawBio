"""scoring.py — Turn an OLS4 search response into candidates ranked by string similarity.

OLS4's `/api/search` does not return a numeric relevance score in its default
JSON (no `score` field on the docs), so this module computes a
deterministic string-similarity ratio between the input value and each
candidate's label/synonyms. This is lexical similarity, not biological or probabilistic
confidence. It keeps every value traceable to (a) the raw
OLS4 response and (b) a fixed, auditable formula — never an opaque or
model-guessed number.
"""

from __future__ import annotations

import difflib

SYNONYM_FIELDS = (
    "exact_synonyms",
    "related_synonyms",
    "narrow_synonyms",
    "broad_synonyms",
)


def normalise(value: str) -> str:
    """Case/whitespace-fold a raw metadata value.

    Used both as the cache/fixture key and as the basis for scoring, so the
    same normalisation must be used everywhere a value is looked up.
    """
    return " ".join(str(value).strip().lower().split())


def candidate_names(doc: dict) -> list[str]:
    """All label + synonym strings OLS4 returned for one candidate document."""
    names = [doc.get("label", "")]
    for field in SYNONYM_FIELDS:
        names.extend(doc.get(field) or [])
    return [n for n in names if n]


def string_similarity(query: str, doc: dict) -> float:
    """Best string-similarity ratio (0-1) between `query` and any of the
    candidate's label/synonyms. difflib.SequenceMatcher.ratio() is used
    because it is deterministic, needs no extra dependency, and is easy for a
    human to reproduce by hand when auditing a flagged row."""
    norm_query = normalise(query)
    best = 0.0
    for name in candidate_names(doc):
        ratio = difflib.SequenceMatcher(None, norm_query, normalise(name)).ratio()
        if ratio > best:
            best = ratio
    return round(best, 4)


def rank_candidates(query: str, docs: list[dict], top_n: int = 3) -> list[dict]:
    """Compute string similarity for every OLS4 doc against `query` and return the top `top_n`,
    sorted by similarity descending. Ties keep OLS4's original relevance order
    (Python's sort is stable)."""
    scored = []
    for doc in docs:
        scored.append(
            {
                "ontology_id": doc.get("obo_id") or doc.get("short_form") or "",
                "label": doc.get("label", ""),
                "string_similarity": string_similarity(query, doc),
                "iri": doc.get("iri", ""),
            }
        )
    scored.sort(key=lambda c: c["string_similarity"], reverse=True)
    return scored[:top_n]


def best_candidate(candidates: list[dict], min_similarity: float) -> tuple[dict | None, bool]:
    """Return (top candidate or None, flagged).

    `flagged` is True whenever there is no top candidate, or its string similarity is
    below `min_similarity`. This function never silently returns a weak match
    pick as if it were trustworthy — the caller must check `flagged`.
    """
    if not candidates:
        return None, True
    top = candidates[0]
    return top, top["string_similarity"] < min_similarity
