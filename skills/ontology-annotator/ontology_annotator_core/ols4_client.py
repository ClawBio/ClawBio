"""ols4_client.py — Caching client for the EBI OLS4 search API, with an
offline fixture mode used by --demo and by the test suite.

OLS4 docs: https://www.ebi.ac.uk/ols4/help
Endpoint used: GET https://www.ebi.ac.uk/ols4/api/search
"""

from __future__ import annotations

import hashlib
import json
import time
from pathlib import Path
from typing import Optional

import requests

from .scoring import normalise

BASE_URL = "https://www.ebi.ac.uk/ols4/api/search"
USER_AGENT = "ClawBio-OntologyAnnotator/0.1.0"
RATE_LIMIT_INTERVAL = 0.34  # ~3 req/sec — polite to a shared public EBI endpoint
CACHE_TTL = 86400 * 30  # ontology term labels change rarely; 30 days is generous
DEFAULT_ROWS = 5


def fixture_key(ontology: str, value: str) -> str:
    """Cache/fixture key: normalised value is what gets looked up, so
    'Lung' and 'lung' share one cache entry and one API call."""
    return f"{ontology.strip().lower()}:{normalise(value)}"


def load_fixtures(path: Path) -> dict:
    """Load the recorded-OLS4-responses fixture file.

    Accepts either the on-disk shape (`{"_meta": ..., "responses": {...}}`)
    or a bare `{key: response}` map, so tests can hand it either.
    """
    data = json.loads(Path(path).read_text())
    return data.get("responses", data) if isinstance(data, dict) and "responses" in data else data


class OLS4Client:
    """Rate-limited, disk-cached client for OLS4 `/api/search`.

    Two lookup modes, chosen at construction time:

    - Live (default): HTTP GET against OLS4, results cached to disk at
      `cache_dir` keyed by (ontology, normalised value).
    - Offline (`fixtures` given): every lookup is answered from the
      pre-recorded fixtures dict instead of the network. This is what
      `--demo` and the test suite use, so neither ever needs network access.
    """

    def __init__(
        self,
        cache_dir: Optional[Path] = None,
        use_cache: bool = True,
        fixtures: Optional[dict] = None,
        rows: int = DEFAULT_ROWS,
    ):
        self.cache_dir = cache_dir
        self.use_cache = use_cache and cache_dir is not None
        self.fixtures = fixtures
        self.rows = rows
        self._last_request_time = 0.0
        self._session: Optional[requests.Session] = None
        if self.use_cache and self.cache_dir:
            self.cache_dir.mkdir(parents=True, exist_ok=True)

    @property
    def session(self) -> requests.Session:
        if self._session is None:
            self._session = requests.Session()
            self._session.headers.update(
                {"Accept": "application/json", "User-Agent": USER_AGENT}
            )
        return self._session

    def _throttle(self) -> None:
        elapsed = time.time() - self._last_request_time
        if elapsed < RATE_LIMIT_INTERVAL:
            time.sleep(RATE_LIMIT_INTERVAL - elapsed)
        self._last_request_time = time.time()

    def _cache_path(self, ontology: str, value: str) -> Optional[Path]:
        if not self.cache_dir:
            return None
        digest = hashlib.sha256(fixture_key(ontology, value).encode()).hexdigest()[:16]
        return self.cache_dir / f"{ontology.strip().lower()}_{digest}.json"

    def _get_disk_cache(self, ontology: str, value: str) -> Optional[dict]:
        path = self._cache_path(ontology, value)
        if not path or not path.exists():
            return None
        try:
            data = json.loads(path.read_text())
            if time.time() - data.get("_cached_at", 0) < CACHE_TTL:
                return data.get("response")
        except (json.JSONDecodeError, KeyError):
            pass
        return None

    def _set_disk_cache(self, ontology: str, value: str, response: dict) -> None:
        path = self._cache_path(ontology, value)
        if not path:
            return
        path.write_text(
            json.dumps({"_cached_at": time.time(), "response": response}, indent=2)
        )

    def search(self, value: str, ontology: str) -> dict:
        """Return `{"status", "source", "docs"}` for `value` restricted to `ontology`.

        `status` is one of "ok", "no_fixture", or "error" — this method never
        raises for a zero-hit search (that is a normal, valid OLS4 response)
        or for a network failure; both come back with an empty `docs` list so
        callers can flag the row instead of crashing the whole run.
        """
        ontology = ontology.strip().lower()

        if self.fixtures is not None:
            key = fixture_key(ontology, value)
            raw = self.fixtures.get(key)
            if raw is None:
                return {
                    "status": "no_fixture",
                    "source": "fixture",
                    "message": f"No recorded fixture for {key!r} (offline mode).",
                    "docs": [],
                }
            return {"status": "ok", "source": "fixture", "docs": _extract_docs(raw)}

        if self.use_cache:
            cached = self._get_disk_cache(ontology, value)
            if cached is not None:
                return {"status": "ok", "source": "cache", "docs": _extract_docs(cached)}

        try:
            self._throttle()
            params = {"q": value, "ontology": ontology, "rows": self.rows}
            resp = self.session.get(BASE_URL, params=params, timeout=30)
            if resp.status_code == 429:
                time.sleep(2.0)
                resp = self.session.get(BASE_URL, params=params, timeout=30)
            resp.raise_for_status()
            data = resp.json()
        except Exception as exc:  # noqa: BLE001 — network can fail in many ways; degrade gracefully
            return {"status": "error", "source": "live", "message": str(exc), "docs": []}

        if self.use_cache:
            self._set_disk_cache(ontology, value, data)
        return {"status": "ok", "source": "live", "docs": _extract_docs(data)}


def _extract_docs(raw_response: dict) -> list[dict]:
    return (raw_response or {}).get("response", {}).get("docs", []) or []
