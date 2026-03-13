"""GTN API client — fetches topics, tutorials, and tool mappings from the
Galaxy Training Network (https://training.galaxyproject.org).

All functions return parsed JSON dicts. No caching here — that's handled
by gtn_cache_builder.py.
"""
from __future__ import annotations

import requests

GTN_BASE = "https://training.galaxyproject.org/training-material/api"


def fetch_topics() -> list[dict]:
    """Fetch all GTN topics (/api/topics.json)."""
    resp = requests.get(f"{GTN_BASE}/topics.json", timeout=30)
    resp.raise_for_status()
    return resp.json()


def fetch_topic_detail(topic_id: str) -> dict:
    """Fetch full detail for a single topic (/api/topics/{id}.json).

    Returns dict with 'name', 'title', 'summary', 'materials' (list of
    tutorial metadata including title, time_estimation, objectives, tools).
    """
    resp = requests.get(f"{GTN_BASE}/topics/{topic_id}.json", timeout=30)
    resp.raise_for_status()
    return resp.json()


def fetch_tool_tutorial_map() -> dict:
    """Fetch the tool→tutorial reverse index (/api/top-tools.json).

    Returns dict keyed by tool repo path (e.g. 'devteam/fastqc/fastqc'),
    each value has 'tool_id' (list of [full_id, version]) and 'tutorials'
    (list of [path, title, topic, url]).
    """
    resp = requests.get(f"{GTN_BASE}/top-tools.json", timeout=60)
    resp.raise_for_status()
    return resp.json()


def fetch_tutorial_content(topic_id: str, tutorial_id: str) -> dict | None:
    """Fetch full tutorial content for deep pull.

    Returns dict with 'name', 'title', 'content' (HTML/markdown body),
    or None if the tutorial is not found.
    """
    url = f"{GTN_BASE}/topics/{topic_id}/tutorials/{tutorial_id}/tutorial.json"
    try:
        resp = requests.get(url, timeout=30)
        resp.raise_for_status()
        return resp.json()
    except requests.HTTPError:
        return None
