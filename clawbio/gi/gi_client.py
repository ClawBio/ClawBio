"""Client for the Genomic Intelligence API.

Shared across `gi-promoter`, `gi-splice`, `gi-enhancer`, `gi-chromatin`,
`gi-annotation`, and `gi-expression` skills.

Auth resolution order:
1. Explicit ``api_key=`` constructor arg.
2. ``GI_API_KEY`` environment variable.

If neither is supplied, ``resolve_api_key`` raises ``RuntimeError`` with
instructions. A shared hackathon-tier key is documented in ``.env.example``
at the repo root — ``cp .env.example .env && source .env`` puts it on the
environment. Heavier / production use: request an individual key at
contact@genomicintelligence.ai and ``export GI_API_KEY=gi_…``.

Base URL: ``GI_BASE_URL`` env, default ``https://api.genomicintelligence.ai``.

Contract reference: https://docs.genomicintelligence.ai
"""

from __future__ import annotations

import os
import time
from typing import Any, Dict, Optional

import requests


DEFAULT_BASE_URL = "https://api.genomicintelligence.ai"

MISSING_KEY_MESSAGE = (
    "GI_API_KEY is not set. The gi-* skills call the hosted Genomic "
    "Intelligence API (https://api.genomicintelligence.ai) and require a "
    "partner bearer key.\n\n"
    "Quick start (ClawBio hackathon): copy the shared key from .env.example:\n"
    "    cp .env.example .env && set -a && source .env && set +a\n\n"
    "For heavier / production use, request an individual key at "
    "contact@genomicintelligence.ai, then:\n"
    "    export GI_API_KEY=gi_yourkeyhere"
)


class GIError(RuntimeError):
    """Non-2xx response from the API. Mirrors the ``{error}`` envelope."""

    def __init__(self, status: int, body: Dict[str, Any]):
        err = (body or {}).get("error", {}) if isinstance(body, dict) else {}
        self.status = status
        self.code = err.get("code", "http_error")
        self.message = err.get("message", "")
        self.request_id = err.get("request_id")
        self.details = err.get("details")
        super().__init__(
            f"[{status} {self.code}] {self.message} (request_id={self.request_id})"
        )


def resolve_api_key(explicit: Optional[str] = None) -> str:
    """Apply the auth resolution order documented at module top.

    Raises ``RuntimeError`` with onboarding instructions if no key is found.
    """
    if explicit:
        return explicit
    env = os.environ.get("GI_API_KEY")
    if env:
        return env
    raise RuntimeError(MISSING_KEY_MESSAGE)


class Client:
    """Thin synchronous client for /v1/tasks/{task}/predict."""

    def __init__(
        self,
        api_key: Optional[str] = None,
        base_url: Optional[str] = None,
        timeout: float = 180.0,
    ) -> None:
        self.api_key = resolve_api_key(api_key)
        self.base_url = (base_url or os.environ.get("GI_BASE_URL") or DEFAULT_BASE_URL).rstrip("/")
        self.timeout = timeout
        self._session = requests.Session()
        self._session.headers.update(
            {
                "Authorization": f"Bearer {self.api_key}",
                "Content-Type": "application/json",
                "Accept": "application/json",
                "User-Agent": "ClawBio-GI-Skill/0.1.0",
            }
        )

    def _check(self, resp: requests.Response) -> Dict[str, Any]:
        try:
            body = resp.json()
        except ValueError:
            body = {"error": {"code": "non_json", "message": resp.text[:200]}}
        if not resp.ok:
            raise GIError(resp.status_code, body)
        return body

    @staticmethod
    def _build_body(
        sequence: str,
        sequence_name: str,
        model: Optional[str],
        options: Optional[Dict[str, Any]],
        tss_index: Optional[int],
    ) -> Dict[str, Any]:
        """Assemble a predict body. ``tss_index`` is expression-only.

        The API rejects unknown top-level fields on the expression route
        (``additionalProperties: false``), so optional fields are only
        included when set.
        """
        body: Dict[str, Any] = {"sequence": sequence, "sequence_name": sequence_name}
        if model is not None:
            body["model"] = model
        if options is not None:
            body["options"] = options
        if tss_index is not None:
            body["tss_index"] = tss_index
        return body

    def health(self) -> Dict[str, Any]:
        r = self._session.get(f"{self.base_url}/health", timeout=self.timeout)
        return self._check(r)

    def predict(
        self,
        task: str,
        sequence: str,
        sequence_name: str = "sequence",
        model: Optional[str] = None,
        options: Optional[Dict[str, Any]] = None,
        tss_index: Optional[int] = None,
    ) -> Dict[str, Any]:
        body = self._build_body(sequence, sequence_name, model, options, tss_index)
        r = self._session.post(
            f"{self.base_url}/v1/tasks/{task}/predict",
            json=body,
            timeout=self.timeout,
        )
        return self._check(r)

    def submit_async(
        self,
        task: str,
        sequence: str,
        sequence_name: str = "sequence",
        model: Optional[str] = None,
        options: Optional[Dict[str, Any]] = None,
        tss_index: Optional[int] = None,
    ) -> str:
        body = self._build_body(sequence, sequence_name, model, options, tss_index)
        r = self._session.post(
            f"{self.base_url}/v1/tasks/{task}/predict",
            headers={"Prefer": "respond-async"},
            json=body,
            timeout=self.timeout,
        )
        body = self._check(r)
        return body["data"]["job_id"]

    def get_job(self, job_id: str) -> requests.Response:
        return self._session.get(
            f"{self.base_url}/v1/tasks/jobs/{job_id}", timeout=self.timeout
        )

    def wait_for_job(
        self,
        job_id: str,
        poll_interval: float = 2.0,
        max_wait: float = 30 * 60,
        on_progress=None,
    ) -> Dict[str, Any]:
        deadline = time.monotonic() + max_wait
        while True:
            r = self.get_job(job_id)
            if r.status_code == 200:
                return r.json()
            if r.status_code == 202:
                if on_progress is not None:
                    try:
                        on_progress((r.json().get("data") or {}).get("progress") or {})
                    except Exception:
                        pass
                if time.monotonic() > deadline:
                    raise TimeoutError(f"job {job_id} did not finish within {max_wait}s")
                time.sleep(poll_interval)
                continue
            try:
                body = r.json()
            except ValueError:
                body = {"error": {"code": "non_json", "message": r.text[:200]}}
            raise GIError(r.status_code, body)


def read_fasta(path) -> tuple[str, str]:
    """Tiny FASTA parser (single record). Returns (sequence_name, sequence).

    Concatenates all non-header lines; uppercases; strips whitespace and
    non-ACGTN characters. Sufficient for the demo fixtures bundled in
    each gi-* skill; users with multi-record FASTA should pre-process.
    """
    from pathlib import Path
    name = None
    seq_parts: list[str] = []
    with open(Path(path)) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is None:
                    name = line[1:].split()[0] or "sequence"
                continue
            seq_parts.append("".join(c for c in line.upper() if c in "ACGTN"))
    seq = "".join(seq_parts)
    return name or "sequence", seq
