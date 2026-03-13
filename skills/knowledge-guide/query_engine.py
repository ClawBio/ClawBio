"""Query engine for the knowledge-guide skill.

Resolves free-text queries, topic lookups, tool lookups, and concept
lookups against the GTN cache. Scoring follows the same weighted-keyword
pattern used in galaxy_bridge/tool_recommender.py.
"""
from __future__ import annotations

import re

# Stopwords to strip from free-text queries
_STOPWORDS = frozenset({
    "a", "an", "the", "is", "are", "was", "were", "be", "been", "being",
    "have", "has", "had", "do", "does", "did", "will", "would", "could",
    "should", "may", "might", "shall", "can", "need", "must",
    "in", "on", "at", "to", "for", "of", "with", "by", "from", "as",
    "into", "through", "during", "before", "after", "above", "below",
    "between", "out", "off", "over", "under", "again", "further", "then",
    "and", "but", "or", "nor", "not", "so", "yet", "both", "either",
    "neither", "each", "every", "all", "any", "few", "more", "most",
    "other", "some", "such", "no", "only", "own", "same",
    "what", "which", "who", "whom", "this", "that", "these", "those",
    "how", "why", "when", "where",
    "i", "me", "my", "we", "our", "you", "your", "he", "she", "it",
    "they", "them", "its", "his", "her", "their",
    "about", "does", "matter",
})


def _tokenize(text: str) -> list[str]:
    """Lowercase, strip punctuation, remove stopwords."""
    words = re.findall(r"[a-z0-9]+(?:-[a-z0-9]+)*", text.lower())
    return [w for w in words if w not in _STOPWORDS]


def _flatten_tutorials(cache: dict) -> list[dict]:
    """Extract all tutorials from the cache into a flat list."""
    tutorials = []
    for topic in cache.get("topics", []):
        for tut in topic.get("tutorials", []):
            tutorials.append(tut)
    return tutorials


def search_tutorials(
    query: str,
    cache: dict,
    max_results: int = 5,
) -> list[dict]:
    """Free-text search across all cached tutorials.

    Scores each tutorial by keyword matches in title (weight 4),
    objectives (weight 3), topic title (weight 2), and tool names (weight 1).

    Returns top-N results sorted by score descending.
    """
    tokens = _tokenize(query)
    if not tokens:
        return []

    tutorials = _flatten_tutorials(cache)
    scored: list[tuple[float, dict]] = []

    for tut in tutorials:
        score = 0.0
        title_lower = tut.get("title", "").lower()
        objectives = tut.get("objectives", [])
        if not isinstance(objectives, list):
            objectives = []
        objectives_text = " ".join(str(o) for o in objectives).lower()
        topic_lower = tut.get("topic", "").lower()
        tools_text = " ".join(str(t) for t in tut.get("tools", [])).lower()

        # Full phrase match (highest signal)
        query_lower = query.lower().strip()
        if query_lower in title_lower:
            score += 10.0
        if query_lower in objectives_text:
            score += 6.0

        # Per-token scoring
        for token in tokens:
            if token in title_lower:
                score += 4.0
            if token in objectives_text:
                score += 3.0
            if token in topic_lower:
                score += 2.0
            if token in tools_text:
                score += 1.0

        if score > 0:
            scored.append((score, tut))

    scored.sort(key=lambda x: x[0], reverse=True)

    results = []
    for score, tut in scored[:max_results]:
        results.append({
            "title": tut["title"],
            "name": tut.get("name", ""),
            "topic": tut.get("topic", ""),
            "time_estimation": tut.get("time_estimation", ""),
            "level": tut.get("level", ""),
            "objectives": tut.get("objectives", []),
            "url": tut.get("url", ""),
            "relevance_score": round(score, 1),
        })
    return results


def lookup_topic(topic_id: str, cache: dict) -> dict | None:
    """Direct topic lookup by ID. Returns topic dict or None."""
    for topic in cache.get("topics", []):
        if topic["name"] == topic_id:
            return topic
    return None


def lookup_tool(tool_name: str, cache: dict) -> list[dict]:
    """Reverse lookup: tool name → tutorials that use it.

    Searches the tool_index keys (partial match) and also scans
    tutorial tool lists for direct name matches.
    """
    tool_lower = tool_name.lower().strip()
    results = []
    seen = set()

    # 1. Check tool_index (keyed by repo path like devteam/fastqc/fastqc)
    for tool_key, tutorials in cache.get("tool_index", {}).items():
        if tool_lower in tool_key.lower():
            for tut in tutorials:
                key = tut.get("path", tut.get("title", ""))
                if key not in seen:
                    seen.add(key)
                    results.append(tut)

    # 2. Scan tutorial tool lists for direct matches
    for topic in cache.get("topics", []):
        for tut in topic.get("tutorials", []):
            for tool in tut.get("tools", []):
                if tool_lower in str(tool).lower():
                    key = f"{tut.get('topic', '')}/{tut.get('name', '')}"
                    if key not in seen:
                        seen.add(key)
                        results.append({
                            "path": key,
                            "title": tut["title"],
                            "topic": tut.get("topic", ""),
                            "url": tut.get("url", ""),
                            "time_estimation": tut.get("time_estimation", ""),
                            "level": tut.get("level", ""),
                            "objectives": tut.get("objectives", []),
                        })

    return results


def lookup_concept(
    concept: str,
    cache: dict,
    recommendations: dict | None = None,
) -> list[dict]:
    """Fuzzy concept lookup — checks skill_recommendations.json first,
    then falls back to free-text search.

    Returns list of tutorial matches.
    """
    concept_lower = concept.lower().strip()
    results = []
    seen = set()

    # 1. Check recommendations for matching concepts
    if recommendations:
        for skill_name, entry in recommendations.items():
            for c in entry.get("concepts", []):
                if concept_lower in c.lower() or c.lower() in concept_lower:
                    for tut in entry.get("tutorials", []):
                        if tut.get("id", "") not in seen:
                            seen.add(tut["id"])
                            results.append(tut)

    # 2. Fall back to free-text search
    if not results:
        results = search_tutorials(concept, cache)

    return results
