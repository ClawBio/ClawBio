"""FinnGen R12 PheWAS client: parsing results, counts, and error reporting.

Fixture is a real live FinnGen R12 response recorded on 2026-11-26:
  finngen_r12_rs3798220.json
      https://r12.finngen.fi/api/variant/6:160540105-T-C
      Subset of real R12 `results` rows preserving the exact response schema
      (`regions`, `results`, `variant`, with `n_case` and `n_control`).
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

SKILL_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(SKILL_DIR))

from gwas_lookup_api import finngen  # noqa: E402

FIXTURES = Path(__file__).parent / "fixtures"
RS3798220 = json.loads((FIXTURES / "finngen_r12_rs3798220.json").read_text())


class FakeClient:
    def __init__(self, payload):
        self.payload = payload
        self.calls: list[str] = []

    def get(self, endpoint, params=None):
        self.calls.append(endpoint)
        return self.payload


def _install(monkeypatch, payload):
    client = FakeClient(payload)
    monkeypatch.setattr(finngen, "_make_client", lambda *a, **k: client)
    return client


def test_reads_r12_results_and_case_control_counts(monkeypatch):
    client = _install(monkeypatch, RS3798220)
    result = finngen.get_phewas("6", 160540105, "T", "C")

    assert client.calls == ["api/variant/6:160540105-T-C"]
    assert result["source"] == "finngen"
    assert result["status"] == "ok"
    assert result["variant"] == "6:160540105-T-C"
    assert result["total_associations"] == len(RS3798220["results"])
    assert result["total_associations"] > 0

    first = result["associations"][0]
    expected = RS3798220["results"][0]
    assert first["phenocode"] == expected["phenocode"]
    assert first["phenostring"] == expected["phenostring"]
    assert first["category"] == expected["category"]
    assert first["pval"] == pytest.approx(expected["pval"])
    assert first["beta"] == pytest.approx(expected["beta"])
    assert first["sebeta"] == pytest.approx(expected["sebeta"])
    assert first["maf"] == pytest.approx(expected["maf"])
    assert first["num_cases"] == expected["n_case"]
    assert first["num_controls"] == expected["n_control"]
    assert first["num_cases"] == 4677
    assert first["num_controls"] == 480316


@pytest.mark.parametrize(
    "payload",
    [
        [],
        "not json object",
        {},
        {"phenos": [{"phenocode": "LEGACY"}]},
        {"results": None},
        {"results": "not a list"},
        {"results": []},
        {"results": ["not a dict"]},
    ],
)
def test_missing_or_empty_results_returns_error(monkeypatch, payload):
    _install(monkeypatch, payload)
    result = finngen.get_phewas("6", 160540105, "T", "C")
    assert result["source"] == "finngen"
    assert result["status"] == "error"
    assert "associations" not in result
