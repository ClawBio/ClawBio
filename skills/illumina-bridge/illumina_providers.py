"""Metadata provider interfaces for Illumina Connected Analytics enrichment."""

from __future__ import annotations

import json
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any
from urllib.parse import urlparse

try:
    import requests
except ImportError:  # pragma: no cover - exercised in runtime envs without requests
    requests = None  # type: ignore[assignment]


DEFAULT_ICA_BASE_URL = "https://ica.illumina.com/ica/rest"
ICA_ACCEPT_HEADER = "application/vnd.illumina.v3+json"


@dataclass
class MetadataEnrichmentResult:
    """Structured result returned by metadata providers."""

    provider: str
    status: str
    warnings: list[str] = field(default_factory=list)
    project: dict[str, Any] = field(default_factory=dict)
    run: dict[str, Any] = field(default_factory=dict)
    samples: list[dict[str, Any]] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return {
            "provider": self.provider,
            "status": self.status,
            "warnings": list(self.warnings),
            "project": dict(self.project),
            "run": dict(self.run),
            "samples": [dict(sample) for sample in self.samples],
        }


class BaseMetadataProvider:
    """Provider interface for optional metadata enrichment."""

    provider_name = "none"

    def enrich(
        self,
        *,
        bundle_dir: Path,
        project_id: str | None = None,
        run_id: str | None = None,
        allow_mock: bool = False,
    ) -> MetadataEnrichmentResult:
        raise NotImplementedError


class NullMetadataProvider(BaseMetadataProvider):
    """Default provider that performs no enrichment."""

    provider_name = "none"

    def enrich(
        self,
        *,
        bundle_dir: Path,
        project_id: str | None = None,
        run_id: str | None = None,
        allow_mock: bool = False,
    ) -> MetadataEnrichmentResult:
        return MetadataEnrichmentResult(provider="none", status="disabled")


class ICAMetadataProvider(BaseMetadataProvider):
    """Metadata-only ICA connector."""

    provider_name = "ica"

    def __init__(
        self,
        *,
        api_key: str | None = None,
        base_url: str | None = None,
        session: Any = None,
    ):
        self._initialization_warnings: list[str] = []
        self.api_key = self._normalize_api_key(
            api_key if api_key is not None else os.environ.get("ILLUMINA_ICA_API_KEY")
        )
        resolved_base_url = base_url if base_url is not None else os.environ.get("ILLUMINA_ICA_BASE_URL")
        self.base_url = self._resolve_base_url(resolved_base_url)
        self.session = session
        if self.session is None and requests is not None:
            self.session = requests.Session()
        if self.session is not None and hasattr(self.session, "headers"):
            self.session.headers.update(
                {
                    "Accept": ICA_ACCEPT_HEADER,
                }
            )

    def __repr__(self) -> str:
        return (
            f"{self.__class__.__name__}("
            f"base_url={self.base_url!r}, "
            f"has_api_key={self.api_key is not None})"
        )

    @staticmethod
    def _normalize_api_key(value: str | None) -> str | None:
        if value is None:
            return None
        stripped = value.strip()
        return stripped or None

    @staticmethod
    def _is_allowed_base_url(candidate: str) -> bool:
        parsed = urlparse(candidate)
        host = parsed.hostname or ""
        return parsed.scheme == "https" and (
            host == "ica.illumina.com" or host.endswith(".illumina.com")
        )

    def _resolve_base_url(self, candidate: str | None) -> str:
        normalized = (candidate or DEFAULT_ICA_BASE_URL).strip().rstrip("/")
        if not normalized:
            normalized = DEFAULT_ICA_BASE_URL
        if self._is_allowed_base_url(normalized):
            return normalized
        self._initialization_warnings.append(
            "ILLUMINA_ICA_BASE_URL is not a trusted Illumina HTTPS endpoint; "
            f"falling back to {DEFAULT_ICA_BASE_URL}."
        )
        return DEFAULT_ICA_BASE_URL

    def _result_warnings(self, warnings: list[str] | None = None) -> list[str]:
        combined = list(self._initialization_warnings)
        if warnings:
            combined.extend(warnings)
        return combined

    def _fetch_json(self, endpoint: str) -> dict[str, Any]:
        if self.session is None:
            raise RuntimeError("requests is not installed; ICA metadata lookup is unavailable.")
        headers = {"X-API-Key": self.api_key} if self.api_key else {}
        response = self.session.get(f"{self.base_url}{endpoint}", timeout=30, headers=headers)
        response.raise_for_status()
        payload = response.json()
        if not isinstance(payload, dict):
            raise ValueError("ICA metadata response must be a JSON object.")
        return payload

    @staticmethod
    def _lookup_warning(resource: str, exc: Exception) -> str:
        """Describe the failed lookup without copying raw response or URL text."""
        if requests is not None and isinstance(exc, requests.HTTPError):
            status = exc.response.status_code if exc.response is not None else None
            code = f"HTTP {status}" if status is not None else "HTTP error"
            if status == 401:
                action = "Check ILLUMINA_ICA_API_KEY."
            elif status == 403:
                action = f"Check the API key's permissions for this {resource}."
            elif status == 404:
                flag = "--ica-project-id" if resource == "project" else "--ica-run-id"
                action = f"Check {flag} and access to the selected project."
            elif status == 429 or (status is not None and status >= 500):
                action = "Retry later and check ICA service availability."
            else:
                action = "Check ICA access and the endpoint configuration."
            detail = f"failed ({code}). {action}"
        elif isinstance(exc, ValueError):
            detail = "returned an invalid JSON response. Check ILLUMINA_ICA_BASE_URL and retry."
        elif requests is not None and isinstance(exc, requests.Timeout):
            detail = "timed out. Check network connectivity and retry."
        elif requests is not None and isinstance(exc, requests.ConnectionError):
            detail = "could not connect. Check the network and ILLUMINA_ICA_BASE_URL."
        else:
            detail = "request failed. Check ICA access and the endpoint configuration."
        return f"ICA {resource} lookup {detail} Continuing without ICA metadata enrichment."

    def _normalize_project(self, payload: dict[str, Any], fallback_id: str | None) -> dict[str, Any]:
        return {
            "id": payload.get("id", fallback_id or ""),
            "name": payload.get("name", ""),
            "active": payload.get("active"),
            # Kept for output compatibility; Project in ICA v3 has no status field.
            "status": payload.get("status", ""),
        }

    def _normalize_run(self, payload: dict[str, Any], fallback_id: str | None) -> dict[str, Any]:
        pipeline = payload.get("pipeline", {})
        pipeline_name = ""
        if isinstance(pipeline, dict):
            pipeline_name = pipeline.get("name", "") or pipeline.get("code", "")
        elif isinstance(pipeline, str):
            pipeline_name = pipeline

        return {
            "id": payload.get("id", fallback_id or ""),
            "name": payload.get("userReference") or payload.get("reference", ""),
            "status": payload.get("status", ""),
            "pipeline": pipeline_name,
        }

    def _load_mock_payload(self, bundle_dir: Path) -> dict[str, Any]:
        mock_path = bundle_dir / "mock_ica_metadata.json"
        if not mock_path.exists():
            raise FileNotFoundError(f"No mock ICA metadata found in bundle: {mock_path}")
        return json.loads(mock_path.read_text(encoding="utf-8"))

    def _result_from_payload(
        self,
        payload: dict[str, Any],
        *,
        status: str,
        warnings: list[str] | None = None,
        project_id: str | None = None,
        run_id: str | None = None,
    ) -> MetadataEnrichmentResult:
        warnings = self._result_warnings(warnings)
        run = self._normalize_run(payload.get("run", {}), run_id)
        analysis_status = run["status"] or "unknown"
        if analysis_status != "SUCCEEDED":
            warnings.append(
                f"ICA analysis status is {analysis_status}, not SUCCEEDED. "
                "Metadata enrichment does not confirm a successful analysis."
            )
        warnings.append(
            "ICA v3 analysis responses do not define sample-level metadata; "
            "sample-level enrichment is unavailable and ICA sample fields remain empty."
        )
        return MetadataEnrichmentResult(
            provider="ica",
            status=status,
            warnings=warnings,
            project=self._normalize_project(payload.get("project", {}), project_id),
            run=run,
            samples=[],
        )

    def enrich(
        self,
        *,
        bundle_dir: Path,
        project_id: str | None = None,
        run_id: str | None = None,
        allow_mock: bool = False,
    ) -> MetadataEnrichmentResult:
        if not project_id or not run_id:
            return MetadataEnrichmentResult(
                provider="ica",
                status="skipped",
                warnings=self._result_warnings(
                    ["ICA metadata requires both --ica-project-id and --ica-run-id (analysis ID); "
                     "continuing without ICA metadata enrichment."]
                ),
            )

        if not self.api_key:
            if allow_mock:
                try:
                    payload = self._load_mock_payload(bundle_dir)
                    return self._result_from_payload(
                        payload,
                        status="mocked-demo",
                        warnings=["Using bundled mock ICA metadata because no API key was configured."],
                        project_id=project_id,
                        run_id=run_id,
                    )
                except (FileNotFoundError, json.JSONDecodeError):
                    pass
            return MetadataEnrichmentResult(
                provider="ica",
                status="warning",
                warnings=self._result_warnings(
                    ["ILLUMINA_ICA_API_KEY is not set; continuing without ICA metadata enrichment."]
                ),
            )

        if self.session is None:
            return MetadataEnrichmentResult(
                provider="ica",
                status="warning",
                warnings=self._result_warnings(
                    ["The optional 'requests' dependency is not installed; continuing without ICA metadata enrichment."]
                ),
            )

        payloads = {}
        for resource, endpoint in (
            ("project", f"/api/projects/{project_id}"),
            ("analysis", f"/api/projects/{project_id}/analyses/{run_id}"),
        ):
            try:
                payloads[resource] = self._fetch_json(endpoint)
            except Exception as exc:
                return MetadataEnrichmentResult(
                    provider="ica",
                    status="warning",
                    warnings=self._result_warnings([self._lookup_warning(resource, exc)]),
                )

        combined = {
            "project": payloads["project"],
            "run": payloads["analysis"],
        }
        return self._result_from_payload(
            combined,
            status="enriched",
            project_id=project_id,
            run_id=run_id,
        )


def build_metadata_provider(name: str) -> BaseMetadataProvider:
    """Factory for provider selection."""

    if name == "ica":
        return ICAMetadataProvider()
    return NullMetadataProvider()
