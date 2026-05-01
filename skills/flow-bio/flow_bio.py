#!/usr/bin/env python3
"""
flow_bio.py — Flow.bio API bridge for ClawBio
===============================================
Authenticate, browse pipelines/samples/projects, search, upload data,
launch pipeline executions, and check run status on any Flow instance.

Usage:
    python flow_bio.py --demo
    python flow_bio.py --pipelines
    python flow_bio.py --samples
    python flow_bio.py --search "RNA-seq"
    python flow_bio.py --upload-sample --name S1 --sample-type RNA-Seq --reads1 R1.fq.gz
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

SKILL_DIR = Path(__file__).resolve().parent
_PROJECT_ROOT = SKILL_DIR.parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

try:
    from clawbio.common.report import write_result_json as _shared_write_result_json
except ImportError:
    _shared_write_result_json = None

# ---------------------------------------------------------------------------
# Flow API client
# ---------------------------------------------------------------------------

DEFAULT_BASE_URL = "https://app.flow.bio/api"


class FlowClient:
    """Lightweight REST client for the Flow.bio API.

    Authentication mirrors the flowbio Python library:
      1. POST /login returns an access token (5-min lifespan) in the body
         and sets an HTTP-only refresh token cookie (7-day lifespan).
      2. requests.Session stores the cookie automatically.
      3. Before each request, _ensure_token() checks the token age.
         If >4 minutes have elapsed, GET /token uses the refresh cookie
         to obtain a fresh access token — no re-login needed.
    """

    # Refresh the access token after this many seconds (token lifetime is
    # 5 minutes; we refresh a minute early to avoid mid-request expiry).
    _TOKEN_REFRESH_AFTER = 4 * 60  # 240 seconds

    def __init__(self, base_url: str | None = None, token: str | None = None):
        try:
            import requests as _requests
        except ImportError:
            print("ERROR: requests not installed. Run: pip install requests", file=sys.stderr)
            sys.exit(1)
        self._requests = _requests
        self.base_url = (base_url or os.environ.get("FLOW_URL", DEFAULT_BASE_URL)).rstrip("/")
        self._token: str | None = token or os.environ.get("FLOW_TOKEN")
        self._last_token_refresh: float | None = None
        self._session = _requests.Session()
        if self._token:
            self._session.headers["Authorization"] = f"Bearer {self._token}"
            self._last_token_refresh = time.time()

    # -- auth ---------------------------------------------------------------

    def login(self, username: str, password: str) -> dict:
        """Authenticate with username/password and store the JWT token.

        The server response also sets an HTTP-only refresh-token cookie
        on the session, which is used by refresh_token() to obtain fresh
        access tokens without re-sending credentials.
        """
        resp = self._session.post(f"{self.base_url}/login", json={"username": username, "password": password})
        resp.raise_for_status()
        data = resp.json()
        self._token = data.get("token")
        if self._token:
            self._session.headers["Authorization"] = f"Bearer {self._token}"
            self._last_token_refresh = time.time()
        return data

    def refresh_token(self) -> None:
        """Use the refresh-token cookie to obtain a fresh access token.

        Called automatically by _ensure_token() when the current access
        token is near expiry.  Can also be called manually.
        """
        resp = self._session.get(f"{self.base_url}/token")
        resp.raise_for_status()
        data = resp.json()
        new_token = data.get("token")
        if new_token:
            self._token = new_token
            self._session.headers["Authorization"] = f"Bearer {new_token}"
            self._last_token_refresh = time.time()

    def set_token(self, token: str) -> None:
        """Set an existing JWT token (e.g. from environment)."""
        self._token = token
        self._session.headers["Authorization"] = f"Bearer {token}"
        self._last_token_refresh = time.time()

    @property
    def is_authenticated(self) -> bool:
        return self._token is not None

    def _ensure_token(self) -> None:
        """Refresh the access token if it is near expiry.

        Mirrors the flowbio v1 client pattern: check token age before
        every authenticated request, and silently refresh if stale.
        """
        if (
            self._last_token_refresh is not None
            and (time.time() - self._last_token_refresh) > self._TOKEN_REFRESH_AFTER
        ):
            try:
                self.refresh_token()
            except Exception:
                pass  # best-effort; the next request will fail with 401 if truly expired

    # -- discovery ----------------------------------------------------------

    def get_pipelines(self) -> list[dict]:
        """List available pipelines."""
        return self._get("/pipelines")

    def get_pipeline(self, pipeline_id: str) -> dict:
        """Get details for a specific pipeline."""
        return self._get(f"/pipelines/{pipeline_id}")

    def get_samples_owned(self) -> list[dict]:
        """List samples owned by the authenticated user."""
        return self._get_paginated("/samples/owned", "samples")

    def get_samples_shared(self) -> list[dict]:
        """List samples shared with the authenticated user."""
        return self._get_paginated("/samples/shared", "samples")

    def get_sample(self, sample_id: str) -> dict:
        """Get details for a specific sample."""
        return self._get(f"/samples/{sample_id}")

    def get_projects_owned(self) -> list[dict]:
        """List projects owned by the authenticated user."""
        return self._get_paginated("/projects/owned", "projects")

    def get_project(self, project_id: str) -> dict:
        """Get details for a specific project."""
        return self._get(f"/projects/{project_id}")

    def get_organisms(self) -> list[dict]:
        """List available organisms."""
        return self._get("/organisms")

    def get_sample_types(self) -> list[dict]:
        """List available sample types."""
        return self._get("/samples/types")

    def get_metadata_attributes(self) -> list[dict]:
        """List metadata attributes for sample uploads."""
        return self._get("/samples/metadata")

    def get_executions_owned(self) -> list[dict]:
        """List executions owned by the authenticated user."""
        return self._get_paginated("/executions/owned", "executions")

    def get_execution(self, execution_id: str) -> dict:
        """Get details for a specific execution."""
        return self._get(f"/executions/{execution_id}")

    def get_data_owned(self) -> list[dict]:
        """List data owned by the authenticated user."""
        return self._get_paginated("/data/owned", "data")

    def get_data(self, data_id: str) -> dict:
        """Get details for a specific data item."""
        return self._get(f"/data/{data_id}")

    # -- search -------------------------------------------------------------

    def search(self, query: str) -> dict:
        """Full-text search across Flow resources."""
        return self._get("/search", params={"q": query})

    def search_samples(self, filters: dict[str, str]) -> dict:
        """Search samples with metadata and field filters.

        Supported filter keys (all optional, combinable):
          - name: substring match on sample name
          - sample_types: comma-separated type identifiers
          - organism: organism ID, or "true"/"false" for has/missing
          - owner: substring match on owner name/username
          - project: substring match on project name
          - Any metadata attribute identifier (e.g. purification_target,
            experimental_method): substring match on value, or "true"/"false"
            for has/missing
          - full_metadata: "true" to include all metadata (not just in_table)
          - sort: field to sort by (default: -created)

        Returns dict with "count", "page", and "samples" keys.
        """
        return self._get("/samples/search", params=filters)

    # -- actions ------------------------------------------------------------

    def upload_sample(
        self,
        name: str,
        sample_type: str,
        reads1: Path,
        reads2: Path | None = None,
        metadata: dict[str, str] | None = None,
        project_id: str | None = None,
        organism_id: str | None = None,
        chunk_size: int = 1_000_000,
    ) -> dict:
        """Upload a sample with chunked file transfer."""
        import uuid

        file_id = str(uuid.uuid4())
        files_to_upload: list[tuple[str, Path]] = [("reads1", reads1)]
        if reads2:
            files_to_upload.append(("reads2", reads2))

        last_response: dict = {}
        for field_name, file_path in files_to_upload:
            if not file_path.exists():
                raise FileNotFoundError(f"File not found: {file_path}")

            file_size = file_path.stat().st_size
            total_chunks = max(1, (file_size + chunk_size - 1) // chunk_size)

            try:
                from tqdm import tqdm
                progress = tqdm(total=total_chunks, desc=f"Uploading {file_path.name}", unit="chunk")
            except ImportError:
                progress = None

            with open(file_path, "rb") as fh:
                for chunk_idx in range(total_chunks):
                    blob = fh.read(chunk_size)
                    is_last = chunk_idx == total_chunks - 1

                    form_data: dict = {
                        "filename": file_path.name,
                        "chunk": str(chunk_idx),
                        "total_chunks": str(total_chunks),
                        "file_id": file_id,
                    }

                    if is_last and field_name == files_to_upload[-1][0]:
                        form_data["name"] = name
                        form_data["sample_type"] = sample_type
                        if metadata:
                            form_data["metadata"] = json.dumps(metadata)
                        if project_id:
                            form_data["project_id"] = project_id
                        if organism_id:
                            form_data["organism_id"] = organism_id

                    resp = self._session.post(
                        f"{self.base_url}/upload/sample",
                        data=form_data,
                        files={"file": (file_path.name, blob)},
                    )
                    resp.raise_for_status()
                    last_response = resp.json()

                    if progress:
                        progress.update(1)

            if progress:
                progress.close()

        return last_response

    def run_pipeline(
        self,
        pipeline_version_id: str,
        sample_ids: list[str] | None = None,
        data_ids: list[str] | None = None,
        params: dict | None = None,
        genome_id: str | None = None,
    ) -> dict:
        """Launch a pipeline execution."""
        body: dict = {}
        if sample_ids:
            body["sample_params"] = sample_ids
        if data_ids:
            body["data_params"] = data_ids
        if params:
            body["params"] = params
        if genome_id:
            body["genome"] = genome_id
        return self._post(f"/pipelines/versions/{pipeline_version_id}/run", json_body=body)

    # -- HTTP helpers -------------------------------------------------------

    def _get(self, path: str, params: dict | None = None) -> dict | list:
        self._ensure_token()
        resp = self._session.get(f"{self.base_url}{path}", params=params)
        resp.raise_for_status()
        return resp.json()

    def _get_paginated(self, path: str, items_key: str, max_pages: int = 10) -> list[dict]:
        """Fetch all pages of a paginated endpoint."""
        self._ensure_token()
        all_items: list[dict] = []
        for page in range(1, max_pages + 1):
            resp = self._session.get(f"{self.base_url}{path}", params={"page": page})
            resp.raise_for_status()
            data = resp.json()
            if isinstance(data, list):
                all_items.extend(data)
                break
            items = data.get(items_key, [])
            all_items.extend(items)
            if len(items) == 0 or page >= data.get("page_count", page):
                break
        return all_items

    def _post(self, path: str, json_body: dict | None = None, auth_required: bool = True) -> dict:
        if auth_required and not self.is_authenticated:
            print("ERROR: Not authenticated. Use --login first or set FLOW_TOKEN.", file=sys.stderr)
            sys.exit(1)
        self._ensure_token()
        resp = self._session.post(f"{self.base_url}{path}", json=json_body)
        resp.raise_for_status()
        return resp.json()


# ---------------------------------------------------------------------------
# Demo mode — live API with offline fallback
# ---------------------------------------------------------------------------

DEMO_CACHE = SKILL_DIR / "data" / "demo_cache.json"


def _load_demo_cache() -> dict:
    """Load pre-cached public Flow.bio data for offline demo mode."""
    if DEMO_CACHE.exists():
        return json.loads(DEMO_CACHE.read_text(encoding="utf-8"))
    return {}


def _flatten_pipelines(categories: list[dict]) -> list[dict]:
    """Flatten the nested categories → subcategories → pipelines response
    into a flat list of pipeline dicts, each tagged with its category."""
    flat: list[dict] = []
    for cat in categories:
        cat_name = cat.get("name", "")
        for subcat in cat.get("subcategories", []):
            subcat_name = subcat.get("name", "")
            for pipe in subcat.get("pipelines", []):
                flat.append({
                    **pipe,
                    "category": cat_name,
                    "subcategory": subcat_name,
                })
    return flat


def run_demo(
    output_dir: Path | None = None,
    base_url: str | None = None,
    username: str | None = None,
    password: str | None = None,
    token: str | None = None,
) -> dict:
    """Show a live overview of the Flow.bio instance.

    Public endpoints (pipelines, organisms, sample types) are always
    queried.  If credentials are provided (via arguments or env vars),
    owned samples, projects, and executions are included too.
    """
    if output_dir is None:
        output_dir = Path("/tmp/flow_demo")
    output_dir.mkdir(parents=True, exist_ok=True)

    url = base_url or os.environ.get("FLOW_URL", DEFAULT_BASE_URL)
    client = FlowClient(base_url=url, token=token)

    # Authenticate if credentials are available
    _username = username or os.environ.get("FLOW_USERNAME")
    _password = password or os.environ.get("FLOW_PASSWORD")
    authenticated = client.is_authenticated
    user_label = "(anonymous)"
    if not authenticated and _username and _password:
        try:
            resp = client.login(_username, _password)
            user_label = resp.get("user", {}).get("username", _username)
            authenticated = True
        except Exception as exc:
            print(f"  Warning: login failed ({exc}), continuing with public endpoints only", file=sys.stderr)
    elif authenticated:
        user_label = "(token)"

    ts = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")

    print()
    print("Flow Bio Bridge — Overview")
    print("=" * 50)
    print(f"Instance: {url}")
    print(f"User: {user_label}")
    print(f"Authenticated: {'yes' if authenticated else 'no (public endpoints only)'}")
    print()

    # -- Public endpoints (no auth required) --------------------------------
    # Try live API first; fall back to bundled demo cache for offline use.

    cache = _load_demo_cache()
    live = True

    # Pipelines (nested categories response)
    pipelines: list[dict] = []
    try:
        raw = client.get_pipelines()
        if isinstance(raw, list):
            pipelines = _flatten_pipelines(raw)
    except Exception:
        pipelines = cache.get("pipelines", [])
        live = False

    # Organisms
    organisms: list[dict] = []
    try:
        organisms = client.get_organisms()
        if not isinstance(organisms, list):
            organisms = []
    except Exception:
        organisms = cache.get("organisms", [])
        live = False

    # Sample types
    sample_types: list[dict] = []
    try:
        sample_types = client.get_sample_types()
        if not isinstance(sample_types, list):
            sample_types = []
    except Exception:
        sample_types = cache.get("sample_types", [])
        live = False

    if not live:
        print("(Using cached data — could not reach Flow API)")
        print()

    print(f"Pipelines ({len(pipelines)} available):")
    for p in pipelines:
        name = p.get("name", "?")
        desc = (p.get("description") or "")[:50]
        cat = p.get("subcategory", p.get("category", ""))
        print(f"  {name:30s} [{cat}]  {desc}")
    print()

    print(f"Organisms ({len(organisms)} available):")
    for o in organisms:
        print(f"  {o.get('id', '?'):>6s}  {o.get('name', o.get('latin_name', '?'))}")
    print()

    print(f"Sample Types ({len(sample_types)} available):")
    for t in sample_types:
        print(f"  {t.get('identifier', '?'):20s} {t.get('name', '?')}")
    print()

    # -- Protected endpoints (auth required) --------------------------------

    samples: list[dict] = []
    projects: list[dict] = []
    executions: list[dict] = []

    if authenticated:
        try:
            samples = client.get_samples_owned()
        except Exception as exc:
            print(f"  Warning: could not fetch samples ({exc})", file=sys.stderr)

        print(f"Samples ({len(samples)} owned):")
        for s in samples:
            sid = s.get("id", "?")
            name = s.get("name", "?")
            stype = s.get("sample_type", {})
            stype_name = stype.get("name", stype) if isinstance(stype, dict) else str(stype)
            org = s.get("organism", {})
            org_name = org.get("name", org) if isinstance(org, dict) else str(org or "—")
            proj = s.get("project_name", "—") or "—"
            print(f"  {sid:>8s}  {name:20s} {str(stype_name):12s} {str(org_name):25s} {proj}")
        print()

        try:
            projects = client.get_projects_owned()
        except Exception as exc:
            print(f"  Warning: could not fetch projects ({exc})", file=sys.stderr)

        print(f"Projects ({len(projects)} owned):")
        for p in projects:
            pid = p.get("id", "?")
            name = p.get("name", "?")
            desc = (p.get("description") or "")[:50]
            print(f"  {pid:>8s}  {name:25s} {desc}")
        print()

        try:
            executions = client.get_executions_owned()
        except Exception as exc:
            print(f"  Warning: could not fetch executions ({exc})", file=sys.stderr)

        print(f"Executions ({len(executions)} owned):")
        for e in executions:
            eid = e.get("id", "?")
            status = e.get("status", "?")
            pipe_name = e.get("pipeline_name", "?")
            pipe_ver = e.get("pipeline_version", "")
            created = str(e.get("created") or "?")[:19]
            ver_suffix = f" v{pipe_ver}" if pipe_ver else ""
            print(f"  {eid:>8s}  {pipe_name}{ver_suffix:30s} {status:12s} {created}")
        print()
    else:
        print("(Samples, projects, and executions require authentication.)")
        print("Set FLOW_USERNAME + FLOW_PASSWORD, or FLOW_TOKEN, to see your data.")
        print()

    # -- Write outputs ------------------------------------------------------

    result = {
        "mode": "live",
        "instance": url,
        "user": user_label,
        "authenticated": authenticated,
        "pipelines": pipelines,
        "organisms": organisms,
        "sample_types": sample_types,
        "samples": samples,
        "projects": projects,
        "executions": executions,
        "timestamp": ts,
    }
    (output_dir / "result.json").write_text(
        json.dumps(result, indent=2, default=str) + "\n", encoding="utf-8",
    )

    write_report(output_dir, "demo", result)

    repro = output_dir / "reproducibility"
    repro.mkdir(exist_ok=True)
    (repro / "commands.sh").write_text(
        "#!/usr/bin/env bash\n"
        f"# Flow overview — instance: {url}\n"
        f"# Date: {ts}\n\n"
        "python skills/flow-bio/flow_bio.py --demo --output /tmp/flow_demo\n",
        encoding="utf-8",
    )
    (repro / "environment.yml").write_text(
        f"flow_url: {url}\n"
        f"mode: live\n"
        f"date: {ts}\n",
        encoding="utf-8",
    )

    print(f"Report: {output_dir / 'report.md'}")
    print(f"Result: {output_dir / 'result.json'}")
    print(f"Reproducibility: {repro}/")
    print()
    return result


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _get_execution_pipeline_name(e: dict) -> str:
    """Extract pipeline name from an execution dict.

    The list endpoint returns flat ``pipeline_name`` / ``pipeline_version``
    fields, while the detail endpoint nests them under
    ``pipelineVersion.pipeline.name``.  This helper handles both.
    """
    # Flat fields from the list serializer
    name = e.get("pipeline_name")
    if name:
        ver = e.get("pipeline_version", "")
        return f"{name} v{ver}" if ver else name
    # Nested fields from the detail serializer
    pv = e.get("pipelineVersion", e.get("pipeline_version", {}))
    if isinstance(pv, dict):
        inner = pv.get("pipeline", {})
        name = inner.get("name", "?") if isinstance(inner, dict) else str(inner)
        ver = pv.get("version", pv.get("name", ""))
        return f"{name} v{ver}" if ver else name
    return str(pv) if pv else "?"


def _search_item_label(category: str, item: dict) -> str:
    """Pick the best display label for a search result item by category."""
    if category == "data":
        return item.get("filename", item.get("name", item.get("id", "?")))
    if category == "executions":
        return item.get("identifier", item.get("pipeline_name", item.get("id", "?")))
    return item.get("name", item.get("id", "?"))


def _is_auth_error(exc: Exception) -> bool:
    """Check if an exception is an HTTP 401 authentication error."""
    return hasattr(exc, "response") and getattr(getattr(exc, "response", None), "status_code", 0) == 401


# ---------------------------------------------------------------------------
# Report generator
# ---------------------------------------------------------------------------


def write_report(output_dir: Path, action: str, data: dict) -> Path:
    """Write a markdown report for a Flow.bio interaction."""
    report_path = output_dir / "report.md"
    ts = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")

    lines = [
        f"# Flow Bio Bridge Report",
        "",
        f"**Date**: {ts}",
        f"**Skill**: flow-bio",
        f"**Action**: {action}",
        f"**Instance**: {data.get('instance', data.get('flow_url', 'N/A'))}",
        "",
    ]

    if action == "demo":
        is_live = data.get("mode") == "live"
        heading = "## Flow.bio Overview" if is_live else "## Summary (Demo Mode)"
        lines.append(heading)
        lines.append("")
        if not is_live:
            lines.append("This report was generated with cached data. No Flow API calls were made.")
            lines.append("")
        lines.append(f"- **User**: {data.get('user', 'N/A')}")
        lines.append(f"- **Authenticated**: {'yes' if data.get('authenticated') else 'no'}")
        lines.append(f"- **Pipelines**: {len(data.get('pipelines', []))}")
        lines.append(f"- **Organisms**: {len(data.get('organisms', []))}")
        lines.append(f"- **Sample Types**: {len(data.get('sample_types', []))}")
        lines.append(f"- **Samples**: {len(data.get('samples', []))}")
        lines.append(f"- **Projects**: {len(data.get('projects', []))}")
        lines.append(f"- **Executions**: {len(data.get('executions', []))}")
        lines.append("")

        # Pipelines table
        if data.get("pipelines"):
            lines.append("### Pipelines")
            lines.append("")
            lines.append("| Name | Category | Description |")
            lines.append("|------|----------|-------------|")
            for p in data["pipelines"]:
                name = p.get("name", "?")
                cat = p.get("subcategory", p.get("category", ""))
                desc = (p.get("description") or "")[:60]
                lines.append(f"| {name} | {cat} | {desc} |")
            lines.append("")

        # Organisms table
        if data.get("organisms"):
            lines.append("### Organisms")
            lines.append("")
            lines.append("| ID | Name |")
            lines.append("|----|------|")
            for o in data["organisms"]:
                lines.append(f"| {o.get('id', '?')} | {o.get('name', o.get('latin_name', '?'))} |")
            lines.append("")

        # Sample types table
        if data.get("sample_types"):
            lines.append("### Sample Types")
            lines.append("")
            lines.append("| Identifier | Name |")
            lines.append("|------------|------|")
            for t in data["sample_types"]:
                lines.append(f"| {t.get('identifier', '?')} | {t.get('name', '?')} |")
            lines.append("")

        # Samples table (auth required)
        if data.get("samples"):
            lines.append("### Samples")
            lines.append("")
            lines.append("| ID | Name | Type | Organism | Project |")
            lines.append("|----|------|------|----------|---------|")
            for s in data["samples"]:
                sid = s.get("id", "?")
                name = s.get("name", "?")
                stype = s.get("sample_type", {})
                stype_name = stype.get("name", stype) if isinstance(stype, dict) else str(stype)
                org = s.get("organism", {})
                org_name = org.get("name", org) if isinstance(org, dict) else str(org or "—")
                proj = s.get("project_name", "—") or "—"
                lines.append(f"| {sid} | {name} | {stype_name} | {org_name} | {proj} |")
            lines.append("")

        # Executions table (auth required)
        if data.get("executions"):
            lines.append("### Executions")
            lines.append("")
            lines.append("| ID | Pipeline | Status | Created |")
            lines.append("|----|----------|--------|---------|")
            for e in data["executions"]:
                eid = e.get("id", "?")
                status = e.get("status", "?")
                pipe_name = _get_execution_pipeline_name(e)
                created = (str(e.get("created", "?")))[:19]
                lines.append(f"| {eid} | {pipe_name} | {status} | {created} |")

    elif action == "pipelines":
        lines.append("## Pipelines")
        lines.append("")
        lines.append("| Name | Category | Description |")
        lines.append("|------|----------|-------------|")
        for p in data.get("pipelines", []):
            name = p.get("name", "?")
            cat = p.get("subcategory", p.get("category", ""))
            desc = (p.get("description") or "")[:80]
            lines.append(f"| {name} | {cat} | {desc} |")

    elif action == "samples":
        lines.append("## Samples")
        lines.append("")
        lines.append("| ID | Name | Type | Organism | Project |")
        lines.append("|----|------|------|----------|---------|")
        for s in data.get("samples", []):
            sid = s.get("id", "?")
            name = s.get("name", "?")
            stype = s.get("sample_type", {})
            stype_name = stype.get("name", stype) if isinstance(stype, dict) else str(stype)
            org = s.get("organism", {})
            org_name = org.get("name", org) if isinstance(org, dict) else str(org or "—")
            proj = s.get("project_name", "—") or "—"
            lines.append(f"| {sid} | {name} | {stype_name} | {org_name} | {proj} |")

    elif action == "projects":
        lines.append("## Projects")
        lines.append("")
        lines.append("| ID | Name | Description |")
        lines.append("|----|------|-------------|")
        for p in data.get("projects", []):
            pid = p.get("id", "?")
            name = p.get("name", "?")
            desc = (p.get("description") or "")[:60]
            lines.append(f"| {pid} | {name} | {desc} |")

    elif action == "executions":
        lines.append("## Executions")
        lines.append("")
        lines.append("| ID | Pipeline | Status | Created |")
        lines.append("|----|----------|--------|---------|")
        for e in data.get("executions", []):
            eid = e.get("id", "?")
            status = e.get("status", "?")
            pipe_name = _get_execution_pipeline_name(e)
            created = str(e.get("created", "?"))[:19]
            lines.append(f"| {eid} | {pipe_name} | {status} | {created} |")

    elif action == "organisms":
        lines.append("## Organisms")
        lines.append("")
        lines.append("| ID | Name |")
        lines.append("|----|------|")
        for o in data.get("organisms", []):
            lines.append(f"| {o.get('id', '?')} | {o.get('name', o.get('latin_name', '?'))} |")

    elif action == "data":
        lines.append("## Data")
        lines.append("")
        lines.append("| ID | Filename | Type | Size |")
        lines.append("|----|----------|------|------|")
        for d in data.get("data_items", data.get("data", [])):
            lines.append(f"| {d.get('id', '?')} | {d.get('filename', d.get('name', '?'))} | {d.get('filetype', '?')} | {d.get('size', '?')} |")

    elif action == "metadata_attributes":
        lines.append("## Metadata Attributes")
        lines.append("")
        lines.append("| Name | Identifier | Required | Has Options | User Terms |")
        lines.append("|------|------------|----------|-------------|------------|")
        for attr in data.get("attributes", []):
            name = attr.get("name", "?")
            ident = attr.get("identifier", "?")
            req = "yes" if attr.get("required") else "no"
            opts = "yes" if attr.get("has_options") else "no"
            user = "yes" if attr.get("allow_user_terms") else "no"
            lines.append(f"| {name} | {ident} | {req} | {opts} | {user} |")

    elif action == "execution":
        lines.append("## Execution Details")
        lines.append("")
        e = data.get("execution", data)
        lines.append(f"- **ID**: {e.get('id', '?')}")
        lines.append(f"- **Status**: {e.get('status', '?')}")
        lines.append(f"- **Pipeline**: {_get_execution_pipeline_name(e)}")
        lines.append(f"- **Created**: {e.get('created', '?')}")
        lines.append(f"- **Started**: {e.get('started', 'N/A')}")
        lines.append(f"- **Finished**: {e.get('finished', 'N/A')}")

    elif action == "pipeline":
        lines.append("## Pipeline Details")
        lines.append("")
        lines.append(f"- **Name**: {data.get('name', '?')}")
        lines.append(f"- **ID**: {data.get('id', '?')}")
        lines.append(f"- **Description**: {data.get('description', 'N/A')}")

    elif action == "sample":
        lines.append("## Sample Details")
        lines.append("")
        lines.append(f"- **Name**: {data.get('name', '?')}")
        lines.append(f"- **ID**: {data.get('id', '?')}")
        stype = data.get("sample_type", data.get("sampleType", "?"))
        stype_name = stype.get("name", stype) if isinstance(stype, dict) else str(stype)
        lines.append(f"- **Type**: {stype_name}")
        org = data.get("organism", {})
        org_name = org.get("name", org) if isinstance(org, dict) else str(org or "N/A")
        lines.append(f"- **Organism**: {org_name}")
        metadata = data.get("metadata", {})
        if metadata:
            lines.append("")
            lines.append("### Metadata")
            lines.append("")
            lines.append("| Attribute | Value | Annotation |")
            lines.append("|-----------|-------|------------|")
            for attr_id, attr_data in metadata.items():
                if isinstance(attr_data, dict):
                    label = attr_data.get("attribute_name", attr_id)
                    value = attr_data.get("value") or "—"
                    annotation = attr_data.get("annotation") or ""
                    lines.append(f"| {label} | {value} | {annotation} |")
                else:
                    lines.append(f"| {attr_id} | {attr_data} | |")

    elif action == "search_samples":
        filters = data.get("filters", {})
        filter_desc = ", ".join(f"{k}={v}" for k, v in filters.items() if k != "full_metadata")
        lines.append("## Sample Search Results")
        lines.append("")
        lines.append(f"**Filters**: {filter_desc}")
        lines.append(f"**Found**: {data.get('count', 0)} samples")
        lines.append("")
        samples = data.get("samples", [])
        if samples:
            lines.append("| # | Name | ID | Type | Organism | Key Metadata |")
            lines.append("|---|------|----|------|----------|--------------|")
            for idx, s in enumerate(samples, 1):
                name = s.get("name", "?")
                sid = s.get("id", "?")
                stype = s.get("sample_type", {})
                stype_name = stype.get("name", stype) if isinstance(stype, dict) else str(stype)
                org = s.get("organism", {})
                org_name = org.get("name", org) if isinstance(org, dict) else str(org or "—")
                metadata = s.get("metadata", {})
                meta_parts = []
                for attr_id, attr_data in metadata.items():
                    if isinstance(attr_data, dict) and attr_data.get("value"):
                        meta_parts.append(f"{attr_data.get('attribute_name', attr_id)}: {attr_data['value']}")
                meta_str = "; ".join(meta_parts[:3])
                if len(meta_parts) > 3:
                    meta_str += f" (+{len(meta_parts) - 3} more)"
                lines.append(f"| {idx} | {name} | {sid} | {stype_name} | {org_name} | {meta_str} |")

    elif action == "search":
        lines.append("## Search Results")
        lines.append("")
        lines.append(f"Query: `{data.get('query', '')}`")
        lines.append("")
        results = data.get("results", data)
        if isinstance(results, dict):
            total = sum(len(v) for v in results.values() if isinstance(v, list))
            lines.append(f"**Total**: {total} results")
            lines.append("")
            for key, items in results.items():
                if isinstance(items, list) and items:
                    lines.append(f"### {key.title()} ({len(items)})")
                    lines.append("")
                    for idx, item in enumerate(items, 1):
                        if isinstance(item, dict):
                            label = _search_item_label(key, item)
                            item_id = item.get("id", "")
                            id_suffix = f" (`{item_id}`)" if item_id else ""
                            lines.append(f"{idx}. {label}{id_suffix}")
                        else:
                            lines.append(f"{idx}. {item}")
                    lines.append("")

    else:
        lines.append(f"## Result")
        lines.append("")
        lines.append(f"```json\n{json.dumps(data, indent=2, default=str)}\n```")

    lines.extend([
        "",
        "## Reproducibility",
        "",
        "See `reproducibility/commands.sh` to re-run this interaction.",
        "",
        "## Disclaimer",
        "",
        "ClawBio is a research and educational tool. It is not a medical device "
        "and does not provide clinical diagnoses. Consult a healthcare "
        "professional before making any medical decisions.",
    ])

    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return report_path


def write_result_json(output_dir: Path, action: str, data: dict) -> Path:
    """Write machine-readable result.json using the shared ClawBio envelope."""
    if _shared_write_result_json is not None:
        return _shared_write_result_json(
            output_dir=output_dir,
            skill="flow-bio",
            version="0.1.0",
            summary={"action": action},
            data=data,
        )
    # Fallback if clawbio.common is not importable (e.g. missing pandas)
    output_dir.mkdir(parents=True, exist_ok=True)
    envelope = {
        "skill": "flow-bio",
        "version": "0.1.0",
        "completed_at": datetime.now(timezone.utc).isoformat(),
        "summary": {"action": action},
        "data": data,
    }
    path = output_dir / "result.json"
    path.write_text(json.dumps(envelope, indent=2, default=str) + "\n", encoding="utf-8")
    return path


def write_reproducibility(output_dir: Path, command: str, flow_url: str) -> None:
    """Write reproducibility bundle."""
    repro = output_dir / "reproducibility"
    repro.mkdir(exist_ok=True)
    ts = datetime.now(timezone.utc).isoformat()
    (repro / "commands.sh").write_text(
        f"#!/usr/bin/env bash\n"
        f"# Reproduce this Flow.bio interaction\n"
        f"# Instance: {flow_url}\n"
        f"# Date: {ts}\n\n"
        f"{command}\n",
        encoding="utf-8",
    )
    (repro / "environment.yml").write_text(
        f"flow_url: {flow_url}\n"
        f"date: {ts}\n"
        f"requests_required: true\n",
        encoding="utf-8",
    )


# ---------------------------------------------------------------------------
# Display helpers
# ---------------------------------------------------------------------------


def _print_table(headers: list[str], rows: list[list[str]]) -> None:
    """Print a simple aligned table."""
    widths = [len(h) for h in headers]
    for row in rows:
        for i, cell in enumerate(row):
            if i < len(widths):
                widths[i] = max(widths[i], len(str(cell)))

    fmt = "  ".join(f"{{:<{w}}}" for w in widths)
    print(fmt.format(*headers))
    print(fmt.format(*("-" * w for w in widths)))
    for row in rows:
        padded = [str(row[i]) if i < len(row) else "" for i in range(len(headers))]
        print(fmt.format(*padded))


def _print_json(data: dict | list, indent: int = 2) -> None:
    """Pretty-print JSON data."""
    print(json.dumps(data, indent=indent, default=str))


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(
        description="Flow Bio Bridge — interact with the Flow.bio API",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # Auth
    auth_group = parser.add_argument_group("authentication")
    auth_group.add_argument("--login", action="store_true", help="Authenticate with Flow.bio")
    auth_group.add_argument("--username", metavar="USER", help="Flow username (or set FLOW_USERNAME)")
    auth_group.add_argument("--password", metavar="PASS", help="Flow password (or set FLOW_PASSWORD)")
    auth_group.add_argument("--token", metavar="TOKEN", help="Existing JWT token (or set FLOW_TOKEN)")
    auth_group.add_argument("--url", metavar="URL", help="Flow API base URL (default: https://app.flow.bio/api)")

    # Discovery (list)
    list_group = parser.add_argument_group("discovery")
    list_group.add_argument("--pipelines", action="store_true", help="List available pipelines")
    list_group.add_argument("--samples", action="store_true", help="List owned samples")
    list_group.add_argument("--projects", action="store_true", help="List owned projects")
    list_group.add_argument("--executions", action="store_true", help="List owned executions")
    list_group.add_argument("--organisms", action="store_true", help="List available organisms")
    list_group.add_argument("--sample-types", action="store_true", help="List available sample types")
    list_group.add_argument("--data", action="store_true", help="List owned data")
    list_group.add_argument("--metadata-attributes", action="store_true", help="List metadata attributes schema")

    # Detail (single resource)
    detail_group = parser.add_argument_group("details")
    detail_group.add_argument("--pipeline", "--pipeline-detail", metavar="ID", help="Get pipeline details")
    detail_group.add_argument("--sample", "--sample-detail", metavar="ID", help="Get sample details")
    detail_group.add_argument("--execution", "--execution-detail", metavar="ID", help="Get execution details")

    # Search
    parser.add_argument("--search", metavar="QUERY", help="Search Flow resources")
    parser.add_argument(
        "--search-samples", nargs="+", metavar="KEY=VALUE",
        help="Search samples by metadata/fields (e.g. purification_target=SNRPB organism=Hs name=iCLIP)",
    )

    # Upload
    upload_group = parser.add_argument_group("upload")
    upload_group.add_argument("--upload-sample", action="store_true", help="Upload a sample")
    upload_group.add_argument("--name", metavar="NAME", help="Sample name for upload")
    upload_group.add_argument("--sample-type", metavar="TYPE", help="Sample type (e.g. RNA-Seq)")
    upload_group.add_argument("--reads1", metavar="FILE", help="First reads file")
    upload_group.add_argument("--reads2", metavar="FILE", help="Second reads file (paired-end)")
    upload_group.add_argument("--organism", metavar="NAME", help="Organism name or ID")
    upload_group.add_argument("--project", metavar="ID", help="Project ID to add sample to")

    # Run
    run_group = parser.add_argument_group("run pipeline")
    run_group.add_argument("--run-pipeline", metavar="VERSION_ID", help="Run a pipeline version")
    run_group.add_argument("--run-samples", metavar="IDS", help="Comma-separated sample IDs for pipeline")
    run_group.add_argument("--run-data", metavar="IDS", help="Comma-separated data IDs for pipeline")
    run_group.add_argument("--run-params", metavar="JSON", help="Pipeline parameters as JSON string")
    run_group.add_argument("--genome", metavar="ID", help="Genome ID for pipeline run")

    # Output
    parser.add_argument("--output", metavar="DIR", help="Output directory for reports")
    parser.add_argument("--json", action="store_true", help="Output raw JSON instead of formatted tables")
    parser.add_argument("--demo", action="store_true", help="Run with cached data (offline, no credentials)")

    args = parser.parse_args()

    # No action specified
    has_action = any([
        args.demo, args.login, args.pipelines, args.samples, args.projects,
        args.executions, args.organisms, args.sample_types, args.data,
        args.metadata_attributes,
        args.pipeline, args.sample, args.execution, args.search,
        args.search_samples, args.upload_sample, args.run_pipeline,
    ])
    if not has_action:
        parser.print_help()
        sys.exit(0)

    # Demo mode
    if args.demo:
        out = Path(args.output) if args.output else None
        run_demo(
            output_dir=out,
            base_url=args.url,
            username=args.username,
            password=args.password,
            token=args.token,
        )
        return

    # Build client
    base_url = args.url or os.environ.get("FLOW_URL", DEFAULT_BASE_URL)
    token = args.token or os.environ.get("FLOW_TOKEN")
    client = FlowClient(base_url=base_url, token=token)

    # Login if requested or if we have credentials but no token
    if args.login or (not client.is_authenticated and (args.username or os.environ.get("FLOW_USERNAME"))):
        username = args.username or os.environ.get("FLOW_USERNAME", "")
        password = args.password or os.environ.get("FLOW_PASSWORD", "")
        if not username or not password:
            print("ERROR: --username and --password required (or set FLOW_USERNAME / FLOW_PASSWORD)", file=sys.stderr)
            sys.exit(1)
        print(f"Logging in to {base_url} as {username}...")
        resp = client.login(username, password)
        user_info = resp.get("user", {})
        print(f"Authenticated as: {user_info.get('username', username)}")
        print()

    output_dir = Path(args.output) if args.output else None

    # --- Discovery commands ---

    if args.pipelines:
        raw = client.get_pipelines()
        pipelines = _flatten_pipelines(raw) if isinstance(raw, list) else raw.get("pipelines", raw.get("results", [raw]))
        if args.json:
            _print_json(pipelines)
        else:
            print(f"\nFlow Pipelines ({len(pipelines)} available):\n")
            rows = []
            for p in pipelines:
                rows.append([
                    p.get("id", "?"),
                    p.get("name", "?"),
                    p.get("subcategory", p.get("category", "")),
                    (p.get("description") or "")[:60],
                ])
            _print_table(["ID", "Name", "Category", "Description"], rows)
            print()
        if output_dir:
            output_dir.mkdir(parents=True, exist_ok=True)
            write_result_json(output_dir, "pipelines", {"pipelines": pipelines, "flow_url": base_url})
            write_report(output_dir, "pipelines", {"pipelines": pipelines, "flow_url": base_url})
            write_reproducibility(output_dir, "python skills/flow-bio/flow_bio.py --pipelines", base_url)
        return

    if args.samples:
        try:
            samples = client.get_samples_owned()
        except Exception as exc:
            if _is_auth_error(exc):
                print("Authentication required to list owned samples.", file=sys.stderr)
                print("Set FLOW_USERNAME and FLOW_PASSWORD, or FLOW_TOKEN.", file=sys.stderr)
                sys.exit(1)
            raise
        if args.json:
            _print_json(samples)
        else:
            print(f"\nOwned Samples ({len(samples)}):\n")
            rows = []
            for s in samples:
                stype = s.get("sample_type", s.get("sampleType", {}))
                stype_name = stype.get("name", stype) if isinstance(stype, dict) else str(stype)
                org = s.get("organism", {})
                org_name = org.get("name", org) if isinstance(org, dict) else str(org or "—")
                rows.append([
                    s.get("id", "?"),
                    s.get("name", "?"),
                    stype_name,
                    org_name,
                ])
            _print_table(["ID", "Name", "Type", "Organism"], rows)
            print()
        if output_dir:
            output_dir.mkdir(parents=True, exist_ok=True)
            write_result_json(output_dir, "samples", {"samples": samples, "flow_url": base_url})
            write_report(output_dir, "samples", {"samples": samples, "flow_url": base_url})
            write_reproducibility(output_dir, "python skills/flow-bio/flow_bio.py --samples", base_url)
        return

    if args.projects:
        try:
            projects = client.get_projects_owned()
        except Exception as exc:
            if _is_auth_error(exc):
                print("Authentication required to list owned projects.", file=sys.stderr)
                print("Set FLOW_USERNAME and FLOW_PASSWORD, or FLOW_TOKEN.", file=sys.stderr)
                sys.exit(1)
            raise
        if args.json:
            _print_json(projects)
        else:
            print(f"\nOwned Projects ({len(projects)}):\n")
            rows = []
            for p in projects:
                rows.append([
                    p.get("id", "?"),
                    p.get("name", "?"),
                    (p.get("description") or "")[:50],
                ])
            _print_table(["ID", "Name", "Description"], rows)
            print()
        if output_dir:
            output_dir.mkdir(parents=True, exist_ok=True)
            write_result_json(output_dir, "projects", {"projects": projects, "flow_url": base_url})
            write_report(output_dir, "projects", {"projects": projects, "flow_url": base_url})
            write_reproducibility(output_dir, "python skills/flow-bio/flow_bio.py --projects", base_url)
        return

    if args.executions:
        try:
            executions = client.get_executions_owned()
        except Exception as exc:
            if _is_auth_error(exc):
                print("Authentication required to list executions.", file=sys.stderr)
                print("Set FLOW_USERNAME and FLOW_PASSWORD, or FLOW_TOKEN.", file=sys.stderr)
                sys.exit(1)
            raise
        if args.json:
            _print_json(executions)
        else:
            print(f"\nOwned Executions ({len(executions)}):\n")
            rows = []
            for e in executions:
                rows.append([
                    e.get("id", "?"),
                    _get_execution_pipeline_name(e),
                    e.get("status", "?"),
                    str(e.get("created", "?"))[:19],
                ])
            _print_table(["ID", "Pipeline", "Status", "Created"], rows)
            print()
        if output_dir:
            output_dir.mkdir(parents=True, exist_ok=True)
            write_result_json(output_dir, "executions", {"executions": executions, "flow_url": base_url})
            write_report(output_dir, "executions", {"executions": executions, "flow_url": base_url})
            write_reproducibility(output_dir, "python skills/flow-bio/flow_bio.py --executions", base_url)
        return

    if args.organisms:
        organisms = client.get_organisms()
        if args.json:
            _print_json(organisms)
        else:
            print(f"\nOrganisms ({len(organisms)} available):\n")
            rows = []
            for o in organisms:
                rows.append([
                    o.get("id", "?"),
                    o.get("name", o.get("latin_name", "?")),
                ])
            _print_table(["ID", "Name"], rows)
            print()
        if output_dir:
            output_dir.mkdir(parents=True, exist_ok=True)
            write_result_json(output_dir, "organisms", {"organisms": organisms, "flow_url": base_url})
            write_report(output_dir, "organisms", {"organisms": organisms, "flow_url": base_url})
            write_reproducibility(output_dir, "python skills/flow-bio/flow_bio.py --organisms", base_url)
        return

    if args.sample_types:
        types = client.get_sample_types()
        if args.json:
            _print_json(types)
        else:
            print(f"\nSample Types ({len(types)} available):\n")
            rows = []
            for t in types:
                rows.append([
                    t.get("identifier", "?"),
                    t.get("name", "?"),
                    (t.get("description") or "")[:50],
                ])
            _print_table(["Identifier", "Name", "Description"], rows)
            print()
        if output_dir:
            output_dir.mkdir(parents=True, exist_ok=True)
            write_result_json(output_dir, "sample_types", {"sample_types": types, "flow_url": base_url})
            write_report(output_dir, "sample_types", {"sample_types": types, "flow_url": base_url})
            write_reproducibility(output_dir, "python skills/flow-bio/flow_bio.py --sample-types", base_url)
        return

    if args.metadata_attributes:
        attributes = client.get_metadata_attributes()
        if args.json:
            _print_json(attributes)
        else:
            print(f"\nMetadata Attributes ({len(attributes)} defined):\n")
            for attr in attributes:
                ident = attr.get("identifier", "?")
                name = attr.get("name", "?")
                desc = attr.get("description", "")
                required = attr.get("required", False)
                allow_user = attr.get("allow_user_terms", False)
                has_opts = attr.get("has_options", False)
                req_label = " [required]" if required else ""
                print(f"  {name}{req_label}")
                print(f"    identifier: {ident}")
                if desc:
                    print(f"    description: {desc[:80]}")
                print(f"    has_options: {has_opts}  allow_user_terms: {allow_user}")
                # Show sample type requirements if present
                links = attr.get("sample_type_links", [])
                if links:
                    req_types = [lk.get("sample_type_name", "?") for lk in links if lk.get("required")]
                    if req_types:
                        print(f"    required for: {', '.join(req_types)}")
                print()
        if output_dir:
            output_dir.mkdir(parents=True, exist_ok=True)
            write_result_json(output_dir, "metadata_attributes", {"attributes": attributes, "flow_url": base_url})
            write_report(output_dir, "metadata_attributes", {"attributes": attributes, "flow_url": base_url})
            write_reproducibility(output_dir, "python skills/flow-bio/flow_bio.py --metadata-attributes", base_url)
        return

    if args.data:
        try:
            data_items = client.get_data_owned()
        except Exception as exc:
            if _is_auth_error(exc):
                print("Authentication required to list owned data.", file=sys.stderr)
                print("Set FLOW_USERNAME and FLOW_PASSWORD, or FLOW_TOKEN.", file=sys.stderr)
                sys.exit(1)
            raise
        if args.json:
            _print_json(data_items)
        else:
            print(f"\nOwned Data ({len(data_items)}):\n")
            rows = []
            for d in data_items:
                rows.append([
                    d.get("id", "?"),
                    d.get("filename", d.get("name", "?")),
                    d.get("filetype", "?"),
                    str(d.get("size", "?")),
                ])
            _print_table(["ID", "Filename", "Type", "Size"], rows)
            print()
        if output_dir:
            output_dir.mkdir(parents=True, exist_ok=True)
            write_result_json(output_dir, "data", {"data_items": data_items, "flow_url": base_url})
            write_report(output_dir, "data", {"data_items": data_items, "flow_url": base_url})
            write_reproducibility(output_dir, "python skills/flow-bio/flow_bio.py --data", base_url)
        return

    # --- Detail commands ---

    if args.pipeline:
        data = client.get_pipeline(args.pipeline)
        if args.json:
            _print_json(data)
        else:
            print(f"\nPipeline: {data.get('name', '?')}")
            print(f"ID: {data.get('id', '?')}")
            print(f"Description: {data.get('description', 'N/A')}")
            versions = data.get("versions", [])
            if versions:
                print(f"\nVersions ({len(versions)}):")
                for v in versions:
                    print(f"  - {v.get('id', '?')}: v{v.get('version', '?')}")
            print()
        if output_dir:
            output_dir.mkdir(parents=True, exist_ok=True)
            write_result_json(output_dir, "pipeline", {"pipeline": data, "flow_url": base_url})
            write_report(output_dir, "pipeline", {**data, "flow_url": base_url})
            write_reproducibility(output_dir, f"python skills/flow-bio/flow_bio.py --pipeline {args.pipeline}", base_url)
        return

    if args.sample:
        data = client.get_sample(args.sample)
        if args.json:
            _print_json(data)
        else:
            print(f"\nSample: {data.get('name', '?')}")
            print(f"ID: {data.get('id', '?')}")
            stype = data.get("sample_type", data.get("sampleType", {}))
            stype_name = stype.get("name", stype) if isinstance(stype, dict) else str(stype)
            print(f"Type: {stype_name}")
            org = data.get("organism", {})
            print(f"Organism: {org.get('name', '?') if isinstance(org, dict) else org}")
            metadata = data.get("metadata", {})
            if metadata:
                print(f"\nMetadata ({len(metadata)} attributes):")
                for attr_id, attr_data in metadata.items():
                    if isinstance(attr_data, dict):
                        label = attr_data.get("attribute_name", attr_id)
                        value = attr_data.get("value") or "—"
                        annotation = attr_data.get("annotation")
                        suffix = f" ({annotation})" if annotation else ""
                        print(f"  {label:30s} {value}{suffix}")
                    else:
                        print(f"  {attr_id:30s} {attr_data}")
            print()
        if output_dir:
            output_dir.mkdir(parents=True, exist_ok=True)
            write_result_json(output_dir, "sample", {"sample": data, "flow_url": base_url})
            write_report(output_dir, "sample", {**data, "flow_url": base_url})
            write_reproducibility(output_dir, f"python skills/flow-bio/flow_bio.py --sample {args.sample}", base_url)
        return

    if args.execution:
        data = client.get_execution(args.execution)
        if args.json:
            _print_json(data)
        else:
            print(f"\nExecution: {data.get('id', '?')}")
            print(f"Status: {data.get('status', '?')}")
            print(f"Pipeline: {_get_execution_pipeline_name(data)}")
            print(f"Created: {data.get('created', '?')}")
            print(f"Started: {data.get('started', 'N/A')}")
            print(f"Finished: {data.get('finished', 'N/A')}")
            print()
        if output_dir:
            output_dir.mkdir(parents=True, exist_ok=True)
            write_result_json(output_dir, "execution", {"execution": data, "flow_url": base_url})
            write_report(output_dir, "execution", {"execution": data, "flow_url": base_url})
            write_reproducibility(output_dir, f"python skills/flow-bio/flow_bio.py --execution {args.execution}", base_url)
        return

    # --- Search ---

    if args.search:
        query = args.search.strip()
        results = client.search(query)
        if args.json:
            _print_json(results)
        else:
            print(f"\nSearch results for \"{query}\":\n")
            if isinstance(results, dict):
                total = sum(len(v) for v in results.values() if isinstance(v, list))
                print(f"  Total: {total} results\n")
                for key, items in results.items():
                    if isinstance(items, list) and items:
                        print(f"  {key.title()} ({len(items)}):")
                        for idx, item in enumerate(items, 1):
                            if isinstance(item, dict):
                                label = _search_item_label(key, item)
                                item_id = item.get("id", "")
                                id_suffix = f"  (id: {item_id})" if item_id else ""
                                print(f"    {idx:>2}. {label}{id_suffix}")
                            else:
                                print(f"    {idx:>2}. {item}")
                        print()
            else:
                _print_json(results)
        if output_dir:
            output_dir.mkdir(parents=True, exist_ok=True)
            write_result_json(output_dir, "search", {"query": args.search, "results": results, "flow_url": base_url})
            write_report(output_dir, "search", {"query": args.search, "results": results, "flow_url": base_url})
            write_reproducibility(output_dir, f'python skills/flow-bio/flow_bio.py --search "{args.search}"', base_url)
        return

    # --- Search samples by metadata ---

    if args.search_samples:
        # Parse key=value pairs
        filters: dict[str, str] = {"full_metadata": "true"}
        for pair in args.search_samples:
            if "=" in pair:
                k, v = pair.split("=", 1)
                filters[k] = v
            else:
                # Treat bare words as name search
                filters["name"] = pair

        results = client.search_samples(filters)
        samples = results.get("samples", []) if isinstance(results, dict) else results
        count = results.get("count", len(samples)) if isinstance(results, dict) else len(samples)

        if args.json:
            _print_json(results)
        else:
            filter_desc = " ".join(f"{k}={v}" for k, v in filters.items() if k != "full_metadata")
            print(f"\nSample search: {filter_desc}")
            print(f"Found: {count} samples\n")
            for idx, s in enumerate(samples, 1):
                sid = s.get("id", "?")
                name = s.get("name", "?")
                stype = s.get("sample_type", {})
                stype_name = stype.get("name", stype) if isinstance(stype, dict) else str(stype)
                org = s.get("organism", {})
                org_name = org.get("name", org) if isinstance(org, dict) else str(org or "—")
                owner = s.get("owner_name", "")
                print(f"  {idx:>3}. {name}")
                print(f"       id: {sid}  type: {stype_name}  organism: {org_name}")
                # Show metadata values
                metadata = s.get("metadata", {})
                meta_parts = []
                for attr_id, attr_data in metadata.items():
                    if isinstance(attr_data, dict):
                        val = attr_data.get("value")
                        if val:
                            label = attr_data.get("attribute_name", attr_id)
                            meta_parts.append(f"{label}: {val}")
                if meta_parts:
                    print(f"       {', '.join(meta_parts)}")
                if owner:
                    print(f"       owner: {owner}")
                print()
        if output_dir:
            output_dir.mkdir(parents=True, exist_ok=True)
            write_result_json(output_dir, "search_samples", {"filters": filters, "count": count, "samples": samples, "flow_url": base_url})
            write_report(output_dir, "search_samples", {"filters": filters, "count": count, "samples": samples, "flow_url": base_url})
            cmd_filters = " ".join(f"{k}={v}" for k, v in filters.items() if k != "full_metadata")
            write_reproducibility(output_dir, f'python skills/flow-bio/flow_bio.py --search-samples {cmd_filters}', base_url)
        return

    # --- Upload ---

    if args.upload_sample:
        if not args.name or not args.sample_type or not args.reads1:
            print("ERROR: --name, --sample-type, and --reads1 are required for upload", file=sys.stderr)
            sys.exit(1)
        reads1 = Path(args.reads1)
        reads2 = Path(args.reads2) if args.reads2 else None
        print(f"Uploading sample '{args.name}' ({args.sample_type})...")
        result = client.upload_sample(
            name=args.name,
            sample_type=args.sample_type,
            reads1=reads1,
            reads2=reads2,
            project_id=args.project,
            organism_id=args.organism,
        )
        print(f"Upload complete: {json.dumps(result, indent=2)}")
        if output_dir:
            output_dir.mkdir(parents=True, exist_ok=True)
            write_result_json(output_dir, "upload", {"upload_result": result, "flow_url": base_url})
            cmd = f'python skills/flow-bio/flow_bio.py --upload-sample --name "{args.name}" --sample-type "{args.sample_type}" --reads1 {args.reads1}'
            if args.reads2:
                cmd += f" --reads2 {args.reads2}"
            write_reproducibility(output_dir, cmd, base_url)
        return

    # --- Run pipeline ---

    if args.run_pipeline:
        sample_ids = args.run_samples.split(",") if args.run_samples else None
        data_ids = args.run_data.split(",") if args.run_data else None
        params = json.loads(args.run_params) if args.run_params else None
        print(f"Launching pipeline version {args.run_pipeline}...")
        result = client.run_pipeline(
            pipeline_version_id=args.run_pipeline,
            sample_ids=sample_ids,
            data_ids=data_ids,
            params=params,
            genome_id=args.genome,
        )
        print(f"Execution launched: {json.dumps(result, indent=2)}")
        if output_dir:
            output_dir.mkdir(parents=True, exist_ok=True)
            write_result_json(output_dir, "run", {"execution": result, "flow_url": base_url})
            write_report(output_dir, "execution", {"execution": result})
            write_reproducibility(output_dir, f"python skills/flow-bio/flow_bio.py --run-pipeline {args.run_pipeline}", base_url)
        return


if __name__ == "__main__":
    main()
