"""
test_protocols_io.py — Test suite for protocols.io bridge skill

Run with: pytest skills/protocols-io/tests/test_protocols_io.py -v

Uses pre-cached demo data and mocked API responses — no network required.
"""

import json
import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from protocols_io import (
    DISCLAIMER,
    _load_demo_json,
    _parse_protocol_id,
    _RateLimiter,
    _sha256,
    _strip_html,
    _write_reproducibility,
    download_protocol_pdf,
    format_search_results,
    format_protocol_detail,
    format_steps,
    search_protocols,
    get_protocol,
    get_protocol_steps,
    load_tokens,
    main,
    save_tokens,
)


# ---------------------------------------------------------------------------
# Demo data loading
# ---------------------------------------------------------------------------


def test_load_demo_search_results():
    """Demo search JSON loads and has expected structure."""
    data = _load_demo_json("demo_search_results.json")
    assert "items" in data
    assert "pagination" in data
    assert len(data["items"]) == 5


def test_load_demo_protocol():
    """Demo protocol JSON loads and has expected structure."""
    data = _load_demo_json("demo_protocol.json")
    assert "payload" in data
    payload = data["payload"]
    assert payload["id"] == 12001
    assert payload["title"] == "Total RNA Extraction from Tissue (TRIzol)"
    assert len(payload["steps"]) == 8
    assert len(payload["materials"]) == 5


def test_load_demo_missing_file():
    """Loading a nonexistent demo file exits with error."""
    with pytest.raises(SystemExit):
        _load_demo_json("nonexistent_file.json")


# ---------------------------------------------------------------------------
# Formatters — search results
# ---------------------------------------------------------------------------


def test_format_search_results_header():
    """Search report includes query and result count."""
    data = _load_demo_json("demo_search_results.json")
    md = format_search_results(data, "RNA extraction")
    assert '# Protocols.io Search: "RNA extraction"' in md
    assert "5 results found" in md


def test_format_search_results_items():
    """Each search result renders with title, creator, DOI."""
    data = _load_demo_json("demo_search_results.json")
    md = format_search_results(data, "RNA extraction")
    assert "Tree Mapping for Leaf Collection" in md
    assert "Gene calling with Prodigal" in md
    assert "Lysis Buffer" in md
    assert "Vellend Cabo" in md
    assert "dx.doi.org" in md


def test_format_search_results_urls():
    """Search results include protocols.io URLs."""
    data = _load_demo_json("demo_search_results.json")
    md = format_search_results(data, "test")
    assert "https://www.protocols.io/view/" in md


def test_format_search_results_disclaimer():
    """Search report includes the safety disclaimer."""
    data = _load_demo_json("demo_search_results.json")
    md = format_search_results(data, "test")
    assert DISCLAIMER in md


def test_format_search_results_empty():
    """Empty search results render gracefully."""
    data = {"items": [], "pagination": {"total_results": 0}}
    md = format_search_results(data, "nothing")
    assert "0 results found" in md
    assert DISCLAIMER in md


# ---------------------------------------------------------------------------
# Formatters — protocol detail
# ---------------------------------------------------------------------------


def test_format_protocol_detail_title():
    """Protocol detail includes title and metadata."""
    data = _load_demo_json("demo_protocol.json")
    md = format_protocol_detail(data)
    assert "# Total RNA Extraction from Tissue (TRIzol)" in md
    assert "Demo Researcher" in md
    assert "dx.doi.org" in md


def test_format_protocol_detail_sections():
    """Protocol detail includes all expected sections."""
    data = _load_demo_json("demo_protocol.json")
    md = format_protocol_detail(data)
    assert "## Description" in md
    assert "## Guidelines" in md
    assert "## Before You Start" in md
    assert "## Warnings" in md
    assert "## Materials" in md
    assert "## Steps" in md


def test_format_protocol_detail_materials():
    """Materials section lists reagents with vendors and SKUs."""
    data = _load_demo_json("demo_protocol.json")
    md = format_protocol_detail(data)
    assert "TRIzol Reagent" in md
    assert "Thermo Fisher" in md
    assert "SKU: 15596026" in md
    assert "Chloroform" in md
    assert "Sigma-Aldrich" in md


def test_format_protocol_detail_steps():
    """Steps section includes all 8 protocol steps."""
    data = _load_demo_json("demo_protocol.json")
    md = format_protocol_detail(data)
    assert "Step 1." in md
    assert "Step 8." in md
    assert "TRIzol Reagent" in md
    assert "NanoDrop" in md


def test_format_protocol_detail_stats():
    """Stats line includes views, steps, and exports."""
    data = _load_demo_json("demo_protocol.json")
    md = format_protocol_detail(data)
    assert "4520" in md
    assert "312" in md


def test_format_protocol_detail_disclaimer():
    """Protocol report includes the safety disclaimer."""
    data = _load_demo_json("demo_protocol.json")
    md = format_protocol_detail(data)
    assert DISCLAIMER in md


def test_format_protocol_detail_minimal():
    """Handles a minimal protocol object without crashing."""
    data = {"payload": {"title": "Minimal", "steps": []}}
    md = format_protocol_detail(data)
    assert "# Minimal" in md
    assert DISCLAIMER in md


# ---------------------------------------------------------------------------
# Formatters — steps
# ---------------------------------------------------------------------------


def test_format_steps_count():
    """Steps formatter reports correct step count."""
    data = {"steps": [{"step": "Do A"}, {"step": "Do B"}, {"step": "Do C"}]}
    md = format_steps(data, "test-123")
    assert "3 steps" in md
    assert "Step 1." in md
    assert "Step 3." in md


def test_format_steps_payload_list():
    """Steps formatter handles v4 API response where steps are in 'payload' list."""
    data = {"payload": [{"step": "Do A"}, {"step": "Do B"}], "status_code": 0, "items_count": 2}
    md = format_steps(data, "test-456")
    assert "2 steps" in md
    assert "Step 1." in md
    assert "Step 2." in md


def test_format_steps_json_draft_parsing():
    """Steps with draft JSON content are parsed to plain text."""
    draft = json.dumps({"blocks": [{"text": "Extracted text here"}]})
    data = {"steps": [{"step": draft}]}
    md = format_steps(data, "test")
    assert "Extracted text here" in md


def test_format_steps_empty():
    """Empty steps list renders gracefully."""
    data = {"steps": []}
    md = format_steps(data, "test")
    assert "0 steps" in md
    assert DISCLAIMER in md


# ---------------------------------------------------------------------------
# Token persistence
# ---------------------------------------------------------------------------


def test_save_and_load_tokens(tmp_path):
    """Tokens round-trip through save/load."""
    token_file = tmp_path / "tokens.json"
    with patch("protocols_io.TOKEN_FILE", token_file), \
         patch("protocols_io.CONFIG_DIR", tmp_path):
        save_tokens({"access_token": "abc123", "token_type": "bearer"})
        result = load_tokens()
    assert result["access_token"] == "abc123"
    assert result["token_type"] == "bearer"
    assert "saved_at" in result


def test_load_tokens_missing(tmp_path):
    """Missing token file returns None."""
    with patch("protocols_io.TOKEN_FILE", tmp_path / "nonexistent.json"):
        result = load_tokens()
    assert result is None


# ---------------------------------------------------------------------------
# URL / ID parsing
# ---------------------------------------------------------------------------


def test_parse_full_url():
    """Full protocols.io URL is parsed to URI slug."""
    result = _parse_protocol_id("https://www.protocols.io/view/one-blot-western-grhbv36")
    assert result == "one-blot-western-grhbv36"


def test_parse_full_url_trailing_slash():
    """Trailing slash is stripped."""
    result = _parse_protocol_id("https://www.protocols.io/view/some-protocol-abc123/")
    assert result == "some-protocol-abc123"


def test_parse_full_url_with_query_params():
    """Query params after the slug are stripped."""
    result = _parse_protocol_id("https://www.protocols.io/view/my-proto-xyz?version=2&foo=bar")
    assert result == "my-proto-xyz"


def test_parse_uri_slug():
    """Plain URI slug passes through unchanged."""
    result = _parse_protocol_id("lysis-buffer-20-ml-c4gytv")
    assert result == "lysis-buffer-20-ml-c4gytv"


def test_parse_numeric_id():
    """Numeric ID passes through as string."""
    result = _parse_protocol_id("30756")
    assert result == "30756"


def test_parse_doi():
    """DOI passes through unchanged."""
    result = _parse_protocol_id("10.17504/protocols.io.baaciaaw")
    assert result == "10.17504/protocols.io.baaciaaw"


def test_parse_whitespace():
    """Leading/trailing whitespace is stripped."""
    result = _parse_protocol_id("  some-protocol-abc  ")
    assert result == "some-protocol-abc"


# ---------------------------------------------------------------------------
# API functions (mocked)
# ---------------------------------------------------------------------------


def _mock_response(json_data, status_code=200, headers=None):
    """Build a mock requests.Response."""
    resp = MagicMock()
    resp.status_code = status_code
    resp.json.return_value = json_data
    resp.text = json.dumps(json_data)
    resp.headers = headers if headers is not None else {}
    return resp


def _unlimited_limiter() -> _RateLimiter:
    """Avoid client-side throttle interfering with API unit tests."""
    return _RateLimiter(max_calls=10_000, window=1.0)


@patch("protocols_io.requests.get")
def test_search_protocols(mock_get):
    """search_protocols calls the correct API endpoint."""
    demo_data = _load_demo_json("demo_search_results.json")
    mock_get.return_value = _mock_response(demo_data)

    with patch("protocols_io.get_access_token", return_value="fake_token"), \
         patch("protocols_io._rate_limiter", _unlimited_limiter()):
        result = search_protocols("RNA extraction")

    assert result is not None
    assert len(result["items"]) == 5
    call_url = mock_get.call_args[0][0]
    assert "/api/v3/protocols" in call_url


@patch("protocols_io.requests.get")
def test_get_protocol(mock_get):
    """get_protocol calls the v4 endpoint."""
    demo_data = _load_demo_json("demo_protocol.json")
    mock_get.return_value = _mock_response(demo_data)

    with patch("protocols_io.get_access_token", return_value="fake_token"), \
         patch("protocols_io._rate_limiter", _unlimited_limiter()):
        result = get_protocol(12001)

    assert result is not None
    assert result["payload"]["id"] == 12001
    call_url = mock_get.call_args[0][0]
    assert "/api/v4/protocols/12001" in call_url


@patch("protocols_io.requests.get")
def test_get_protocol_steps(mock_get):
    """get_protocol_steps calls the steps endpoint."""
    steps_data = {"steps": [{"step": "Do something"}], "status_code": 0}
    mock_get.return_value = _mock_response(steps_data)

    with patch("protocols_io.get_access_token", return_value="fake_token"), \
         patch("protocols_io._rate_limiter", _unlimited_limiter()):
        result = get_protocol_steps(12001)

    assert result is not None
    assert len(result["steps"]) == 1
    call_url = mock_get.call_args[0][0]
    assert "/api/v4/protocols/12001/steps" in call_url


@patch("protocols_io.requests.get")
def test_api_expired_token_returns_none(mock_get):
    """Expired token (status 1219) returns None with an error message."""
    expired_resp = _mock_response({"status_code": 1219}, status_code=200)
    mock_get.return_value = expired_resp

    with patch("protocols_io.get_access_token", return_value="old_token"), \
         patch("protocols_io._rate_limiter", _unlimited_limiter()):
        result = search_protocols("test")

    assert result is None
    assert mock_get.call_count == 1


@patch("protocols_io.time.sleep")
@patch("protocols_io.requests.get")
def test_api_429_retries_then_success(mock_get, mock_sleep):
    r429 = _mock_response({}, status_code=429, headers={"Retry-After": "1"})
    ok = _mock_response({"items": [], "status_code": 0}, status_code=200)
    mock_get.side_effect = [r429, ok]

    with patch("protocols_io.get_access_token", return_value="t"), \
         patch("protocols_io._rate_limiter", _unlimited_limiter()):
        result = search_protocols("x")

    assert result is not None
    assert mock_get.call_count == 2
    mock_sleep.assert_called_once_with(1)


@patch("protocols_io.time.sleep")
@patch("protocols_io.requests.get")
def test_api_429_exhausts_retries(mock_get, mock_sleep):
    r429 = _mock_response({}, status_code=429, headers={"Retry-After": "1"})
    mock_get.return_value = r429

    with patch("protocols_io.get_access_token", return_value="t"), \
         patch("protocols_io._rate_limiter", _unlimited_limiter()):
        result = search_protocols("x")

    assert result is None
    assert mock_get.call_count == 3


@patch("protocols_io.time.sleep")
@patch("protocols_io.requests.get")
def test_api_429_retry_after_capped(mock_get, mock_sleep):
    r429 = _mock_response({}, status_code=429, headers={"Retry-After": "9999"})
    ok = _mock_response({"items": [], "status_code": 0}, status_code=200)
    mock_get.side_effect = [r429, ok]

    with patch("protocols_io.get_access_token", return_value="t"), \
         patch("protocols_io._rate_limiter", _unlimited_limiter()):
        search_protocols("x")

    mock_sleep.assert_called_once_with(120)


@patch("protocols_io.time.sleep")
def test_rate_limiter_waits_when_full(mock_sleep):
    lim = _RateLimiter(max_calls=2, window=60.0)
    lim.wait()
    lim.wait()
    lim.wait()
    assert mock_sleep.called


# ---------------------------------------------------------------------------
# Formatters — peer-reviewed badge
# ---------------------------------------------------------------------------


def test_format_search_results_peer_reviewed_badge():
    """Peer-reviewed protocols show the ✅ badge; non-peer-reviewed do not."""
    data = _load_demo_json("demo_search_results.json")
    md = format_search_results(data, "RNA extraction")
    # "Gene calling with Prodigal" has peer_reviewed=1 in demo data
    assert "✅ Peer-reviewed method" in md


def test_format_search_results_url_is_plain_text():
    """URLs are rendered as plain text, not as markdown hyperlinks."""
    data = _load_demo_json("demo_search_results.json")
    md = format_search_results(data, "test")
    # A markdown link would look like [url](url); plain text has no leading '['
    assert "](https://www.protocols.io/view/" not in md
    assert "https://www.protocols.io/view/" in md


# ---------------------------------------------------------------------------
# search_protocols — filter parameters forwarded to API
# ---------------------------------------------------------------------------


@patch("protocols_io.requests.get")
def test_search_protocols_peer_reviewed_param(mock_get):
    """peer_reviewed=1 is forwarded as a query param to the API."""
    demo_data = _load_demo_json("demo_search_results.json")
    mock_get.return_value = _mock_response(demo_data)

    with patch("protocols_io.get_access_token", return_value="fake_token"), \
         patch("protocols_io._rate_limiter", _unlimited_limiter()):
        search_protocols("RNA extraction", peer_reviewed=1)

    call_kwargs = mock_get.call_args[1]
    assert call_kwargs["params"]["peer_reviewed"] == 1


@patch("protocols_io.requests.get")
def test_search_protocols_published_on_param(mock_get):
    """published_on timestamp is forwarded as a query param to the API."""
    demo_data = _load_demo_json("demo_search_results.json")
    mock_get.return_value = _mock_response(demo_data)

    ts = 1700000000
    with patch("protocols_io.get_access_token", return_value="fake_token"), \
         patch("protocols_io._rate_limiter", _unlimited_limiter()):
        search_protocols("RNA extraction", published_on=ts)

    call_kwargs = mock_get.call_args[1]
    assert call_kwargs["params"]["published_on"] == ts


@patch("protocols_io.requests.get")
def test_search_protocols_no_optional_params_by_default(mock_get):
    """peer_reviewed and published_on are absent from params when not supplied."""
    demo_data = _load_demo_json("demo_search_results.json")
    mock_get.return_value = _mock_response(demo_data)

    with patch("protocols_io.get_access_token", return_value="fake_token"), \
         patch("protocols_io._rate_limiter", _unlimited_limiter()):
        search_protocols("RNA extraction")

    call_kwargs = mock_get.call_args[1]
    assert "peer_reviewed" not in call_kwargs["params"]
    assert "published_on" not in call_kwargs["params"]


# ---------------------------------------------------------------------------
# _strip_html
# ---------------------------------------------------------------------------


def test_strip_html_removes_tags():
    assert _strip_html("<b>Bold text</b>") == "Bold text"


def test_strip_html_nested():
    assert _strip_html("<div><span>Hello</span></div>") == "Hello"


def test_strip_html_plain_passthrough():
    assert _strip_html("No tags here") == "No tags here"


def test_strip_html_empty():
    assert _strip_html("") == ""


# ---------------------------------------------------------------------------
# download_protocol_pdf
# ---------------------------------------------------------------------------


@patch("protocols_io.requests.get")
def test_download_protocol_pdf_success(mock_get, tmp_path):
    """PDF bytes are written to the specified output path."""
    fake_pdf = b"%PDF-1.4 fake content"
    resp = MagicMock()
    resp.status_code = 200
    resp.headers = {"Content-Type": "application/pdf"}
    resp.content = fake_pdf
    mock_get.return_value = resp

    out = tmp_path / "output.pdf"
    with patch("protocols_io.get_access_token", return_value=None):
        result = download_protocol_pdf("my-proto-abc123", output_path=out)

    assert result == out
    assert out.read_bytes() == fake_pdf


@patch("protocols_io.requests.get")
def test_download_protocol_pdf_http_error(mock_get, tmp_path):
    """Non-200 HTTP response returns None."""
    resp = MagicMock()
    resp.status_code = 404
    resp.headers = {}
    mock_get.return_value = resp

    with patch("protocols_io.get_access_token", return_value=None):
        result = download_protocol_pdf("missing-proto", output_path=tmp_path / "out.pdf")

    assert result is None


@patch("protocols_io.requests.get")
def test_download_protocol_pdf_not_pdf_content(mock_get, tmp_path):
    """Response with non-PDF content type and non-PDF magic bytes returns None."""
    resp = MagicMock()
    resp.status_code = 200
    resp.headers = {"Content-Type": "text/html"}
    resp.content = b"<html>Not a PDF</html>"
    mock_get.return_value = resp

    with patch("protocols_io.get_access_token", return_value=None):
        result = download_protocol_pdf("some-proto", output_path=tmp_path / "out.pdf")

    assert result is None


@patch("protocols_io.requests.get")
def test_download_protocol_pdf_default_path(mock_get, tmp_path):
    """When output_path is None, a slugified filename is used in cwd."""
    fake_pdf = b"%PDF-1.4 data"
    resp = MagicMock()
    resp.status_code = 200
    resp.headers = {"Content-Type": "application/pdf"}
    resp.content = fake_pdf
    mock_get.return_value = resp

    import os
    original_cwd = os.getcwd()
    os.chdir(tmp_path)
    try:
        with patch("protocols_io.get_access_token", return_value=None):
            result = download_protocol_pdf("my-cool-protocol-xyz")
        assert result is not None
        assert result.suffix == ".pdf"
        assert result.read_bytes() == fake_pdf
    finally:
        os.chdir(original_cwd)


@patch("protocols_io.requests.get")
def test_download_protocol_pdf_network_error(mock_get, tmp_path):
    """Network exception returns None."""
    import requests as req_lib
    mock_get.side_effect = req_lib.RequestException("connection refused")

    with patch("protocols_io.get_access_token", return_value=None):
        result = download_protocol_pdf("some-proto", output_path=tmp_path / "out.pdf")

    assert result is None


# ---------------------------------------------------------------------------
# _parse_protocol_id — dx.doi.org prefix
# ---------------------------------------------------------------------------


def test_parse_dx_doi_prefix():
    """dx.doi.org/ prefix is stripped to leave the bare DOI."""
    result = _parse_protocol_id("dx.doi.org/10.17504/protocols.io.abc123")
    assert result == "10.17504/protocols.io.abc123"


def test_parse_dx_doi_prefix_with_version():
    """dx.doi.org/ DOI with version suffix is parsed correctly."""
    result = _parse_protocol_id("dx.doi.org/10.17504/protocols.io.e6nvwwq2zvmk/v1")
    assert result == "10.17504/protocols.io.e6nvwwq2zvmk/v1"


# ---------------------------------------------------------------------------
# _sha256
# ---------------------------------------------------------------------------


def test_sha256_matches_hashlib(tmp_path):
    """_sha256 produces the same digest as hashlib.sha256 directly."""
    import hashlib
    content = b"test content for hashing"
    f = tmp_path / "test.bin"
    f.write_bytes(content)
    expected = hashlib.sha256(content).hexdigest()
    assert _sha256(f) == expected


# ---------------------------------------------------------------------------
# _write_reproducibility
# ---------------------------------------------------------------------------


def test_write_reproducibility_creates_files(tmp_path):
    """Reproducibility bundle creates commands.sh, checksums.sha256, environment.yml."""
    import argparse
    output_file = tmp_path / "report.md"
    output_file.write_text("# Report")
    args = argparse.Namespace(
        search="RNA extraction", protocol=None, steps=None, demo=False,
        filter="public", peer_reviewed=None, published_on=None,
        page_size=10, page=1, output=str(tmp_path),
    )
    _write_reproducibility(args, tmp_path, [output_file])
    repro = tmp_path / "reproducibility"
    assert (repro / "commands.sh").exists()
    assert (repro / "checksums.sha256").exists()
    assert (repro / "environment.yml").exists()


def test_write_reproducibility_checksum_is_correct(tmp_path):
    """Checksum in checksums.sha256 matches the actual file digest."""
    import argparse, hashlib
    output_file = tmp_path / "report.md"
    content = "# Reproducibility test"
    output_file.write_text(content)
    args = argparse.Namespace(
        search="test", protocol=None, steps=None, demo=False,
        filter="public", peer_reviewed=None, published_on=None,
        page_size=10, page=1, output=str(tmp_path),
    )
    _write_reproducibility(args, tmp_path, [output_file])
    checksum_text = (tmp_path / "reproducibility" / "checksums.sha256").read_text()
    expected = hashlib.sha256(content.encode()).hexdigest()
    assert expected in checksum_text
    assert "report.md" in checksum_text


def test_write_reproducibility_commands_contains_args(tmp_path):
    """commands.sh includes the search query and output path."""
    import argparse
    output_file = tmp_path / "search_results.md"
    output_file.write_text("# Results")
    args = argparse.Namespace(
        search="CRISPR knockout", protocol=None, steps=None, demo=False,
        filter="public", peer_reviewed=None, published_on=None,
        page_size=10, page=1, output=str(tmp_path),
    )
    _write_reproducibility(args, tmp_path, [output_file])
    cmd = (tmp_path / "reproducibility" / "commands.sh").read_text()
    assert "CRISPR knockout" in cmd
    assert str(tmp_path) in cmd


def test_write_reproducibility_skips_missing_files(tmp_path):
    """Files that don't exist are silently excluded from checksums."""
    import argparse
    missing = tmp_path / "nonexistent.md"
    args = argparse.Namespace(
        search="test", protocol=None, steps=None, demo=False,
        filter="public", peer_reviewed=None, published_on=None,
        page_size=10, page=1, output=str(tmp_path),
    )
    _write_reproducibility(args, tmp_path, [missing])
    checksum_text = (tmp_path / "reproducibility" / "checksums.sha256").read_text()
    assert checksum_text.strip() == ""


def test_write_reproducibility_environment_yml(tmp_path):
    """environment.yml contains the conda spec with requests dependency."""
    import argparse
    args = argparse.Namespace(
        search="test", protocol=None, steps=None, demo=False,
        filter="public", peer_reviewed=None, published_on=None,
        page_size=10, page=1, output=str(tmp_path),
    )
    _write_reproducibility(args, tmp_path, [])
    env = (tmp_path / "reproducibility" / "environment.yml").read_text()
    assert "clawbio-protocols-io" in env
    assert "requests" in env


# ---------------------------------------------------------------------------
# main() — --output saves files for --search and --steps
# ---------------------------------------------------------------------------


@patch("protocols_io.requests.get")
def test_main_search_output_saves_markdown(mock_get, tmp_path):
    """--search --output saves search_results.md and reproducibility bundle."""
    demo_data = _load_demo_json("demo_search_results.json")
    mock_get.return_value = _mock_response(demo_data)

    with patch("protocols_io.get_access_token", return_value="fake_token"), \
         patch("protocols_io._rate_limiter", _unlimited_limiter()), \
         patch("sys.argv", ["protocols_io.py", "--search", "RNA extraction",
                            "--output", str(tmp_path)]):
        main()

    assert (tmp_path / "report.md").exists()
    assert (tmp_path / "reproducibility" / "commands.sh").exists()
    assert (tmp_path / "reproducibility" / "checksums.sha256").exists()


@patch("protocols_io.requests.get")
def test_main_steps_output_saves_markdown(mock_get, tmp_path):
    """--steps --output saves steps.md and reproducibility bundle."""
    steps_data = {"steps": [{"step": "Do A"}, {"step": "Do B"}], "status_code": 0}
    mock_get.return_value = _mock_response(steps_data)

    with patch("protocols_io.get_access_token", return_value="fake_token"), \
         patch("protocols_io._rate_limiter", _unlimited_limiter()), \
         patch("sys.argv", ["protocols_io.py", "--steps", "12001",
                            "--output", str(tmp_path)]):
        main()

    assert (tmp_path / "report.md").exists()
    assert (tmp_path / "reproducibility" / "commands.sh").exists()
    assert (tmp_path / "reproducibility" / "checksums.sha256").exists()


@patch("protocols_io.requests.get")
def test_main_protocol_output_saves_pdf_and_repro(mock_get, tmp_path):
    """--protocol --output downloads PDF and writes reproducibility bundle."""
    demo_data = _load_demo_json("demo_protocol.json")
    fake_pdf = b"%PDF-1.4 fake"
    pdf_resp = MagicMock()
    pdf_resp.status_code = 200
    pdf_resp.headers = {"Content-Type": "application/pdf"}
    pdf_resp.content = fake_pdf
    mock_get.side_effect = [
        _mock_response(demo_data),  # get_protocol call
        pdf_resp,                   # download_protocol_pdf call
    ]

    with patch("protocols_io.get_access_token", return_value="fake_token"), \
         patch("protocols_io._rate_limiter", _unlimited_limiter()), \
         patch("sys.argv", ["protocols_io.py", "--protocol", "12001",
                            "--output", str(tmp_path)]):
        main()

    assert (tmp_path / "report.md").exists()
    pdfs = list(tmp_path.glob("*.pdf"))
    assert len(pdfs) == 1
    assert pdfs[0].read_bytes() == fake_pdf
    assert (tmp_path / "reproducibility" / "commands.sh").exists()
    assert (tmp_path / "reproducibility" / "checksums.sha256").exists()
