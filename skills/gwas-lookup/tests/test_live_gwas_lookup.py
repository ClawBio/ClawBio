"""Live smoke test for gwas-lookup: every source answers for rs3798220.

The offline suite runs on pre-fetched fixtures, so it cannot see an upstream
API changing its schema (Open Targets, HTTP 400) or retiring (eQTL Catalogue,
HTTP 410). This test can. Gated by @pytest.mark.live and RUN_LIVE_TESTS=1 so
the offline suite stays network-free.

Run locally with:
    RUN_LIVE_TESTS=1 pytest skills/gwas-lookup/tests/test_live_gwas_lookup.py
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import gwas_lookup  # noqa: E402

pytestmark = pytest.mark.live


@pytest.mark.skipif(
    os.getenv("RUN_LIVE_TESTS") != "1",
    reason="live tests gated on RUN_LIVE_TESTS=1",
)
def test_live_every_source_answers(tmp_path) -> None:
    gwas_lookup.run_lookup("rs3798220", tmp_path / "out", use_cache=False, make_figures=False)
    summary = json.loads((tmp_path / "out" / "result.json").read_text())["summary"]
    assert summary["apis_failed"] == {}, summary["apis_failed"]
