"""
test_base_client.py — BaseClient error reporting and 429 retry.

Runs against a local stdlib HTTP server, so no network is needed and the
real requests/urllib3 retry stack is exercised.
"""

import sys
import threading
import time
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

import pytest
import requests

SKILL_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(SKILL_DIR))

from gwas_lookup_api.base_client import BaseClient


class _Handler(BaseHTTPRequestHandler):
    def _respond(self):
        self.server.hits += 1
        status, headers, body = self.server.script.pop(0) if self.server.script else (200, {}, b'{"ok": true}')
        self.send_response(status)
        self.send_header("Content-Length", str(len(body)))
        for key, value in headers.items():
            self.send_header(key, value)
        self.end_headers()
        self.wfile.write(body)

    do_GET = do_POST = _respond

    def log_message(self, *args):
        pass


@pytest.fixture
def server():
    srv = ThreadingHTTPServer(("127.0.0.1", 0), _Handler)
    srv.hits = 0
    srv.script = []  # (status, headers, body) per request, then 200s
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    yield srv
    srv.shutdown()
    srv.server_close()


def _client(server):
    return BaseClient(f"http://127.0.0.1:{server.server_port}", rate_interval=0.01)


def test_http_error_keeps_the_server_explanation(server):
    server.script = [(400, {}, b'{"errors": [{"message": "Cannot query field \'rsId\'"}]}')]
    with pytest.raises(requests.HTTPError, match="Cannot query field 'rsId'"):
        _client(server).post("graphql", {"query": "{}"})


def test_retries_after_429(server):
    server.script = [(429, {"Retry-After": "1"}, b"")]
    assert _client(server).get("x") == {"ok": True}
    assert server.hits == 2


def test_honours_retry_after_instead_of_fixed_sleep(server):
    server.script = [(429, {"Retry-After": "1"}, b"")]
    start = time.monotonic()
    _client(server).get("x")
    assert 0.9 <= time.monotonic() - start < 1.8


def test_persistent_429_raises_http_error(server, monkeypatch):
    from urllib3.util import Retry

    monkeypatch.setattr(Retry, "sleep", lambda self, response=None: None)
    server.script = [(429, {}, b"slow down")] * 10
    with pytest.raises(requests.HTTPError, match="429"):
        _client(server).get("x")


def test_long_retry_after_is_capped(server, monkeypatch):
    """A quota-style Retry-After must not stall a lookup for an hour."""
    import urllib3.util.retry as retry_mod

    slept = []
    monkeypatch.setattr(retry_mod.time, "sleep", slept.append)
    server.script = [(429, {"Retry-After": "3600"}, b"")]
    assert _client(server).get("x") == {"ok": True}
    assert slept and max(slept) <= 60
