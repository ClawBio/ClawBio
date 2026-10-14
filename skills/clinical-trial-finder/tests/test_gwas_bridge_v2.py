"""gwas_bridge against the GWAS Catalog REST API v2, no network.

Uses the raw v2 responses recorded for gwas-lookup on 2026-10-08
(skills/gwas-lookup/tests/fixtures/gwas_catalog_v2_*.json):
rs7903146 has 353 associations over two pages of 200.
"""

from __future__ import annotations

import importlib.util
import io
import json
import sys
import urllib.error
import urllib.parse
from pathlib import Path

import pytest

SKILL_DIR = Path(__file__).resolve().parent.parent
GWAS_FIXTURES = SKILL_DIR.parent / "gwas-lookup" / "tests" / "fixtures"

spec = importlib.util.spec_from_file_location("gwas_bridge", SKILL_DIR / "gwas_bridge.py")
gwas_bridge = importlib.util.module_from_spec(spec)
sys.modules["gwas_bridge"] = gwas_bridge
spec.loader.exec_module(gwas_bridge)


def _fixture(name: str) -> dict:
    return json.loads((GWAS_FIXTURES / name).read_text())


PAGES = {
    0: _fixture("gwas_catalog_v2_rs7903146_page0.json"),
    1: _fixture("gwas_catalog_v2_rs7903146_page1.json"),
}


class _Resp(io.BytesIO):
    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


def _serve(monkeypatch, pages=PAGES, fail_page=None, fail_code=None):
    seen: list[str] = []

    def fake_urlopen(req, timeout=None):
        url = req.full_url if hasattr(req, "full_url") else req
        seen.append(url)
        query = urllib.parse.parse_qs(urllib.parse.urlparse(url).query)
        page = int(query.get("page", ["0"])[0])
        if fail_code is not None and (fail_page is None or page == fail_page):
            raise urllib.error.HTTPError(url, fail_code, "err", {}, io.BytesIO(b""))
        body = pages[page]
        return _Resp(body.encode() if isinstance(body, str) else json.dumps(body).encode())

    monkeypatch.setattr(gwas_bridge.urllib.request, "urlopen", fake_urlopen)
    monkeypatch.setattr(gwas_bridge, "_try_cached_result", lambda *a, **k: None)
    return seen


def test_api_base_is_v2():
    assert gwas_bridge._GWAS_API == "https://www.ebi.ac.uk/gwas/rest/api/v2"


def test_queries_v2_associations_with_rs_id_filter(monkeypatch):
    seen = _serve(monkeypatch)
    gwas_bridge.resolve_rsid("rs7903146")
    parsed = urllib.parse.urlparse(seen[0])
    query = urllib.parse.parse_qs(parsed.query)
    assert parsed.path == "/gwas/rest/api/v2/associations"
    assert query["rs_id"] == ["rs7903146"]
    assert query["page"] == ["0"]
    assert "singleNucleotidePolymorphisms" not in seen[0]


def test_reads_every_page_not_just_the_first(monkeypatch):
    seen = _serve(monkeypatch)
    gwas_bridge.resolve_rsid("rs7903146")
    pages = [urllib.parse.parse_qs(urllib.parse.urlparse(u).query)["page"][0] for u in seen]
    assert pages == ["0", "1"]


def test_extracts_traits_and_mapped_genes(monkeypatch):
    _serve(monkeypatch)
    result = gwas_bridge.resolve_rsid("rs7903146", max_traits=5)
    assert result["rsid"] == "rs7903146"
    assert "Type 2 Diabetes Mellitus" in result["traits"]
    assert len(result["traits"]) == 5
    # Disease traits rank before measurement traits.
    assert not any("Measurement" in t for t in result["traits"][:1])
    assert result["genes"][0] == "TCF7L2"
    assert set(result["genes"]) == {"TCF7L2", "MTNR1B"}


def test_http_410_raises_never_returns_empty(monkeypatch):
    _serve(monkeypatch, fail_code=410)
    with pytest.raises(ValueError, match="HTTP 410"):
        gwas_bridge.resolve_rsid("rs7903146")


def test_http_error_on_later_page_raises_not_partial(monkeypatch):
    _serve(monkeypatch, fail_page=1, fail_code=503)
    with pytest.raises(ValueError, match="HTTP 503"):
        gwas_bridge.resolve_rsid("rs7903146")


def test_zero_hits_raises_no_associations(monkeypatch):
    empty = {"_links": {}, "page": {"size": 200, "totalElements": 0, "totalPages": 0, "number": 0}}
    _serve(monkeypatch, pages={0: empty})
    with pytest.raises(ValueError, match="No genome-wide significant"):
        gwas_bridge.resolve_rsid("rs999999999999")


@pytest.mark.parametrize(
    "payload",
    [
        "[]",
        "not json",
        {"_embedded": {"associations": []}},
        {"page": {"totalElements": 1, "totalPages": 1}, "_embedded": {"associations": "x"}},
    ],
)
def test_malformed_payload_raises(monkeypatch, payload):
    _serve(monkeypatch, pages={0: payload})
    with pytest.raises(ValueError, match="(?i)malformed"):
        gwas_bridge.resolve_rsid("rs1")


def test_missing_p_value_is_skipped_not_crash(monkeypatch):
    page = json.loads(json.dumps(PAGES[0]))
    page["page"]["totalPages"] = 1
    page["_embedded"]["associations"][0]["p_value"] = None
    _serve(monkeypatch, pages={0: page})
    result = gwas_bridge.resolve_rsid("rs7903146")
    assert result["traits"]
