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
from typing import Any, Dict, Mapping, Optional

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


# IUPAC ambiguity codes. Listed so the parser can say *why* it is refusing:
# these are legitimate FASTA content the model cannot score, which is a
# different problem from a stray character and deserves a different hint.
_IUPAC_AMBIGUITY = "RYSWKMBDHV"


class FastaError(ValueError):
    """Malformed FASTA input, rejected rather than silently repaired.

    Subclasses ``ValueError`` so a caller doing broad input validation still
    catches it, while callers that want to distinguish input problems from
    API problems can catch this specifically.
    """


class GIError(RuntimeError):
    """Non-2xx response from the API. Mirrors the ``{error}`` envelope.

    Branch on ``code`` (a closed enum in the schema — treat an unlisted value
    as a generic failure, not a parse error), never on ``details``: the
    ``details`` payload is keyed on ``code`` and its shape varies by release, so
    it is carried through verbatim for display only. A 422 carries the declared
    ``{"errors": [...]}`` object; older responses sent a bare array instead,
    which is why nothing here parses it. In particular, the expression
    ``tss_index`` checks are a whole-body validator and report at
    ``loc: ["body"]``, never ``body.tss_index``.
    """

    def __init__(
        self,
        status: int,
        body: Dict[str, Any],
        headers: Optional[Mapping[str, str]] = None,
    ):
        err = (body or {}).get("error", {}) if isinstance(body, dict) else {}
        self.status = status
        self.code = err.get("code", "http_error")
        self.message = err.get("message", "")
        # Prefer the envelope's request_id; every error response carries it.
        # Fall back to the X-Request-Id header for robustness (e.g. a non-JSON
        # body from a proxy) — support tickets always need a correlation id.
        self.request_id = err.get("request_id") or (headers or {}).get("X-Request-Id")
        self.details = err.get("details")
        rid = self.request_id or "unset"
        super().__init__(f"[{status} {self.code}] {self.message} (request_id={rid})")


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
    """Thin synchronous client for the six per-task predict operations.

    The API publishes one operation per task — ``POST
    /v1/tasks/promoter/predict``, ``…/splice/…``, ``…/enhancer/…``,
    ``…/chromatin/…``, ``…/annotation/…``, ``…/expression/…`` — each with its
    own request schema (different ``minLength``, different ``options`` class).
    The paths differ only in that one segment, so the f-string below builds
    them all; it is not a claim that they share a contract. An unrecognised
    ``task`` segment is a ``404 not_found``, not a 422.
    """

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
            # http_error is a published enum value; the response arrived with a
            # status and body, it just was not JSON. Client-origin errors carry
            # no request_id, which distinguishes them from server codes.
            body = {"error": {"code": "http_error", "message": resp.text[:200]}}
        if not resp.ok:
            raise GIError(resp.status_code, body, resp.headers)
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

        Every per-task request model is ``additionalProperties: false``, and
        so is every per-task ``options`` model, so an unknown key is a hard
        ``422 validation_failed`` (``type: extra_forbidden``) rather than a
        silent ignore. Optional fields are therefore only included when set,
        and ``tss_index`` — declared on ``ExpressionPredictRequest`` alone —
        must stay unset for the other five tasks.
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
        # An async submit reads data.job_id, so a malformed 2xx has to fail here
        # rather than downstream: a missing data key raised KeyError and a
        # non-object data raised TypeError, neither of which reads as a bad
        # response to the caller.
        data = body.get("data") if isinstance(body, dict) else None
        job_id = data.get("job_id") if isinstance(data, dict) else None
        if not isinstance(job_id, str) or not job_id:
            raise GIError(
                r.status_code,
                {
                    "error": {
                        "code": "http_error",
                        "message": "async submit returned no job_id in data",
                    }
                },
                r.headers,
            )
        return job_id

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
                body = {"error": {"code": "http_error", "message": r.text[:200]}}
            raise GIError(r.status_code, body, r.headers)


def read_fasta(path) -> tuple[str, str]:
    """Parse a single-record FASTA. Returns (sequence_name, sequence).

    Rejects malformed input rather than repairing it. Earlier versions
    silently deleted every character outside ``ACGTN`` and concatenated a
    multi-record file into one chimeric sequence under the first record's
    name. Both are unrecoverable once they happen: deleting an IUPAC
    ambiguity code shifts every base after it, so the model scores a
    sequence the caller never supplied and returns a confident result with
    nothing to indicate the substitution.

    Whitespace, blank lines and lowercase input are still handled — those
    are formatting, not content.

    Raises:
        FastaError: more than one record, a base outside ``ACGTN``, or
            sequence appearing before the first header.
    """
    from pathlib import Path
    name = None
    record_names: list[str] = []
    seq_parts: list[str] = []
    offenders: dict[str, int] = {}
    with open(Path(path)) as fh:
        for lineno, line in enumerate(fh, start=1):
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                header = line[1:].split()[0] if line[1:].split() else "sequence"
                record_names.append(header)
                if name is None:
                    name = header
                continue
            if name is None:
                raise FastaError(
                    f"{path}: sequence on line {lineno} before any '>' header. "
                    f"Those bases would be scored under the first record's name, "
                    f"and the coordinates returned would not describe what you "
                    f"submitted. Add a header, or remove the stray lines."
                )
            # Whitespace anywhere in the line is formatting, not content: the
            # API strips newlines, spaces and tabs before measuring length, so
            # a space-grouped body (10-base blocks from a viewer or Sanger
            # output) must parse here too. Stripping it moves nothing in
            # coordinate space, which is what separates it from an ambiguity
            # code we refuse to guess at.
            upper = "".join(line.split()).upper()
            for char in upper:
                if char not in "ACGTN":
                    offenders.setdefault(char, lineno)
            seq_parts.append(upper)

    if len(record_names) > 1:
        shown = ", ".join(record_names[:3])
        more = f", … ({len(record_names)} total)" if len(record_names) > 3 else ""
        raise FastaError(
            f"{path}: expected a single FASTA record, found {len(record_names)} "
            f"({shown}{more}). Concatenating them would submit a chimeric "
            f"sequence under one name — split the file and submit one record "
            f"per request."
        )

    if offenders:
        detail = ", ".join(
            f"{char!r} (first at line {lineno})"
            for char, lineno in sorted(offenders.items(), key=lambda kv: kv[1])[:5]
        )
        ambiguity = sorted(c for c in offenders if c in _IUPAC_AMBIGUITY)
        hint = (
            " IUPAC ambiguity codes cannot be scored; resolve them to explicit "
            "bases or submit a different region."
            if ambiguity
            else " Remove or resolve them before submitting."
        )
        raise FastaError(
            f"{path}: sequence contains characters outside ACGTN: {detail}.{hint}"
        )

    return name or "sequence", "".join(seq_parts)
