"""GWAS Catalog REST API v2 client: parsing, paging and error reporting.

Fixtures are raw live v2 responses recorded on 2026-10-08:
  gwas_catalog_v2_rs7903146_page{0,1}.json
      /v2/associations?rs_id=rs7903146&size=200&page={0,1}&sort=p_value&direction=asc
      totalElements 353, totalPages 2 (200 + 153 records)
  gwas_catalog_v2_rs2187668.json
      same query for rs2187668, totalElements 16, one page
  gwas_catalog_legacy_410.txt
      body of the retired v1 endpoint, which answers HTTP 410
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest
import requests

SKILL_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(SKILL_DIR))

from gwas_lookup_api import gwas_catalog  # noqa: E402

FIXTURES = Path(__file__).parent / "fixtures"


def _fixture(name: str) -> dict:
    return json.loads((FIXTURES / name).read_text())


RS7903146_PAGES = {
    0: _fixture("gwas_catalog_v2_rs7903146_page0.json"),
    1: _fixture("gwas_catalog_v2_rs7903146_page1.json"),
}
RS2187668 = _fixture("gwas_catalog_v2_rs2187668.json")


class FakeClient:
    """Serves recorded pages by the requested page number; records calls."""

    def __init__(self, pages: dict[int, dict]):
        self.pages = pages
        self.calls: list[tuple[str, dict]] = []

    def get(self, endpoint, params=None):
        params = dict(params or {})
        self.calls.append((endpoint, params))
        return self.pages[params.get("page", 0)]


def _install(monkeypatch, client):
    monkeypatch.setattr(gwas_catalog, "_make_client", lambda *a, **k: client)
    return client


def test_base_url_is_v2():
    assert gwas_catalog.BASE_URL == "https://www.ebi.ac.uk/gwas/rest/api/v2"


def test_queries_v2_associations_filtered_by_rs_id(monkeypatch):
    client = _install(monkeypatch, FakeClient({0: RS2187668}))
    gwas_catalog.get_associations("rs2187668", max_hits=100)
    endpoint, params = client.calls[0]
    assert endpoint == "associations"
    assert params["rs_id"] == "rs2187668"
    assert params["page"] == 0
    # Strongest associations first, so max_hits keeps the most significant.
    assert params["sort"] == "p_value"
    assert params["direction"] == "asc"
    assert params["size"] <= gwas_catalog.PAGE_SIZE


def test_parses_v2_fields_into_unchanged_schema(monkeypatch):
    _install(monkeypatch, FakeClient({0: RS2187668}))
    result = gwas_catalog.get_associations("rs2187668", max_hits=100)
    assert result["status"] == "ok"
    assert result["source"] == "gwas_catalog"
    assert result["rsid"] == "rs2187668"
    assert result["total_associations"] == 16
    assert result["total_available"] == 16
    assert result["truncated"] is False

    first = result["associations"][0]
    legacy_keys = {
        "pvalue", "pvalue_mlog", "pvalue_exponent", "risk_allele",
        "risk_frequency", "or_beta", "beta_num", "beta_direction",
        "beta_unit", "ci", "traits", "study_accession",
    }
    assert legacy_keys <= set(first)
    assert first["pvalue"] == 8e-93
    assert first["pvalue_mlog"] == 8
    assert first["pvalue_exponent"] == -93
    assert first["risk_allele"] == "rs2187668-A"
    assert first["risk_frequency"] == "0.1297"
    assert first["or_beta"] == 4.32
    assert first["beta_num"] is None
    assert first["ci"] == "[NR]"
    assert first["traits"] == ["membranous glomerulonephritis"]
    assert first["efo_ids"] == ["MONDO_0005376"]
    assert first["study_accession"] == "GCST000984"
    assert first["mapped_genes"] == ["HLA-DQA1"]


def test_parses_beta_string_into_number_unit_direction(monkeypatch):
    _install(monkeypatch, FakeClient({0: RS2187668}))
    result = gwas_catalog.get_associations("rs2187668", max_hits=100)
    beta_row = next(a for a in result["associations"] if a["study_accession"] == "GCST90470245")
    assert beta_row["beta_num"] == pytest.approx(0.12451027)
    assert beta_row["beta_unit"] == "unit"
    assert beta_row["beta_direction"] == "increase"
    assert beta_row["or_beta"] is None


@pytest.mark.parametrize(
    "raw, expected",
    [
        ("0.03939 unit increase", (0.03939, "unit", "increase")),
        ("1.2 kg/m2 decrease", (1.2, "kg/m2", "decrease")),
        ("0.5 increase", (0.5, None, "increase")),
        ("-", (None, None, None)),
        ("", (None, None, None)),
        (None, (None, None, None)),
    ],
)
def test_parse_beta(raw, expected):
    assert gwas_catalog._parse_beta(raw) == expected


def test_pages_until_all_records_when_max_hits_exceeds_one_page(monkeypatch):
    client = _install(monkeypatch, FakeClient(RS7903146_PAGES))
    result = gwas_catalog.get_associations("rs7903146", max_hits=1000)
    assert [c[1]["page"] for c in client.calls] == [0, 1]
    assert result["total_associations"] == 353
    assert result["total_available"] == 353
    assert result["truncated"] is False
    assert len(result["associations"]) == 353


def test_max_hits_spanning_pages_is_honoured_exactly(monkeypatch):
    _install(monkeypatch, FakeClient(RS7903146_PAGES))
    result = gwas_catalog.get_associations("rs7903146", max_hits=250)
    assert len(result["associations"]) == 250
    assert result["total_associations"] == 250
    assert result["total_available"] == 353
    assert result["truncated"] is True


def test_max_hits_within_first_page_stops_after_one_request(monkeypatch):
    client = _install(monkeypatch, FakeClient(RS7903146_PAGES))
    result = gwas_catalog.get_associations("rs7903146", max_hits=20)
    assert len(client.calls) == 1
    assert len(result["associations"]) == 20
    assert result["truncated"] is True


def test_short_page_before_reported_total_is_flagged_not_hidden(monkeypatch):
    """If the API stops early the shortfall is reported, never silently dropped."""
    page0 = RS7903146_PAGES[0]
    empty = {"_links": {}, "page": {"size": 200, "totalElements": 353, "totalPages": 2, "number": 1}}
    _install(monkeypatch, FakeClient({0: page0, 1: empty}))
    result = gwas_catalog.get_associations("rs7903146", max_hits=1000)
    assert len(result["associations"]) == 200
    assert result["total_available"] == 353
    assert result["truncated"] is True


def test_zero_hits_is_ok_with_zero(monkeypatch):
    empty = {"_links": {}, "page": {"size": 200, "totalElements": 0, "totalPages": 0, "number": 0}}
    _install(monkeypatch, FakeClient({0: empty}))
    result = gwas_catalog.get_associations("rs999999999999")
    assert result["status"] == "ok"
    assert result["associations"] == []
    assert result["total_associations"] == 0


def _response(status: int, body: str) -> requests.Response:
    resp = requests.Response()
    resp.status_code = status
    resp._content = body.encode()
    resp.url = "https://www.ebi.ac.uk/gwas/rest/api/v2/associations"
    resp.reason = "Gone" if status == 410 else "Error"
    return resp


def test_http_410_is_status_error_never_zero_associations(monkeypatch):
    body = (FIXTURES / "gwas_catalog_legacy_410.txt").read_text()
    monkeypatch.setattr(requests.Session, "get", lambda self, *a, **k: _response(410, body))
    result = gwas_catalog.get_associations("rs7903146", use_cache=False)
    assert result["status"] == "error"
    assert "410" in result["message"]
    assert "associations" not in result


def test_http_error_on_second_page_is_error_not_partial(monkeypatch):
    class FailingSecondPage(FakeClient):
        def get(self, endpoint, params=None):
            if (params or {}).get("page") == 1:
                raise requests.HTTPError("503 Service Unavailable")
            return super().get(endpoint, params)

    _install(monkeypatch, FailingSecondPage(RS7903146_PAGES))
    result = gwas_catalog.get_associations("rs7903146", max_hits=1000)
    assert result["status"] == "error"
    assert "503" in result["message"]


@pytest.mark.parametrize(
    "payload",
    [
        [],
        "not json object",
        {"_embedded": {"associations": []}},  # no page block
        {"page": {"totalElements": 1}, "_embedded": {"associations": "x"}},
        {"page": {"totalElements": 1}, "_embedded": {"associations": ["x"]}},
    ],
)
def test_malformed_payload_is_status_error(monkeypatch, payload):
    _install(monkeypatch, FakeClient({0: payload}))
    result = gwas_catalog.get_associations("rs1")
    assert result["status"] == "error"
    assert "malformed" in result["message"].lower()


def test_provenance_records_v2_base_url(tmp_path):
    from gwas_lookup_core.report import write_reproducibility

    write_reproducibility(tmp_path, "python gwas_lookup.py --rsid rs7903146", [])
    versions = json.loads((tmp_path / "reproducibility" / "api_versions.json").read_text())
    assert versions["apis"]["gwas_catalog"] == "https://www.ebi.ac.uk/gwas/rest/api/v2"
