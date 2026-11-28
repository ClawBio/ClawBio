"""Live smoke test for EUCTR search in clinical-trial-finder.

Gated by @pytest.mark.live and RUN_LIVE_TESTS=1 so the offline suite stays
network-free. Run locally with:

    RUN_LIVE_TESTS=1 uv run pytest skills/clinical-trial-finder/tests/test_live_euctr.py
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import euctr  # noqa: E402

pytestmark = pytest.mark.live


@pytest.mark.skipif(
    os.getenv("RUN_LIVE_TESTS") != "1",
    reason="live tests gated on RUN_LIVE_TESTS=1",
)
def test_live_euctr_search_returns_parsed_trials() -> None:
    trials = euctr.fetch_euctr("breast cancer", max_results=5)
    assert 1 <= len(trials) <= 5
    seen_ids = set()
    for t in trials:
        assert t["nct_id"].startswith("EUCTR")
        assert t["nct_id"] not in seen_ids
        seen_ids.add(t["nct_id"])
        assert t["title"]
        assert not t["title"].startswith("EU trial ")
        assert t["source"] == "euctr"
