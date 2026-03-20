#!/usr/bin/env python3
"""
protocols_io.py — protocols.io bridge for ClawBio
===================================================
Search, browse, and retrieve scientific protocols from protocols.io
via REST API with client token authentication.

Usage:
    python protocols_io.py --login
    python protocols_io.py --search "RNA extraction"
    python protocols_io.py --protocol 30756
    python protocols_io.py --steps 30756
    python protocols_io.py --demo
"""

from __future__ import annotations

import argparse
import collections
import getpass
import itertools
import json
import os
import sys
import threading
import time
from datetime import datetime, timezone
from pathlib import Path

try:
    import requests
except ImportError:
    print("ERROR: requests not installed. Run: pip install requests", file=sys.stderr)
    sys.exit(1)


class Spinner:
    """Animated terminal spinner that runs in a background thread."""

    FRAMES = ["⠋", "⠙", "⠹", "⠸", "⠼", "⠴", "⠦", "⠧", "⠇", "⠏"]

    def __init__(self, message: str = ""):
        self._message = message
        self._stop = threading.Event()
        self._thread: threading.Thread | None = None

    def _spin(self):
        for frame in itertools.cycle(self.FRAMES):
            if self._stop.is_set():
                break
            sys.stderr.write(f"\r  {frame} {self._message}")
            sys.stderr.flush()
            time.sleep(0.08)
        sys.stderr.write(f"\r  ✓ {self._message}\n")
        sys.stderr.flush()

    def __enter__(self):
        self._thread = threading.Thread(target=self._spin, daemon=True)
        self._thread.start()
        return self

    def __exit__(self, *_):
        self._stop.set()
        if self._thread:
            self._thread.join()

SKILL_DIR = Path(__file__).resolve().parent
DEMO_DIR = SKILL_DIR / "demo"
CONFIG_DIR = Path.home() / ".clawbio"
TOKEN_FILE = CONFIG_DIR / "protocols_io_tokens.json"

API_V3 = "https://www.protocols.io/api/v3"
API_V4 = "https://www.protocols.io/api/v4"
RATE_LIMIT = 100  # requests per 60-second window (protocols.io API)
RATE_WINDOW = 60.0
MAX_RETRIES_429 = 3

DISCLAIMER = (
    "*ClawBio is a research and educational tool. It is not a medical device "
    "and does not provide clinical diagnoses. Consult a healthcare professional "
    "before making any medical decisions.*"
)


# ---------------------------------------------------------------------------
# Rate limiting (100 requests / minute; 429 retry with Retry-After)
# ---------------------------------------------------------------------------


class _RateLimiter:
    """Thread-safe sliding-window limiter for outbound API GETs."""

    def __init__(self, max_calls: int = RATE_LIMIT, window: float = RATE_WINDOW):
        self._max = max_calls
        self._window = window
        self._timestamps: collections.deque[float] = collections.deque()
        self._lock = threading.Lock()

    def wait(self) -> None:
        with self._lock:
            now = time.monotonic()
            while self._timestamps and self._timestamps[0] <= now - self._window:
                self._timestamps.popleft()
            if len(self._timestamps) >= self._max:
                sleep_for = self._timestamps[0] - (now - self._window)
                if sleep_for > 0:
                    print(f"  Rate limit (client) — waiting {sleep_for:.1f}s", file=sys.stderr)
                    time.sleep(sleep_for)
            self._timestamps.append(time.monotonic())


_rate_limiter = _RateLimiter()


# ---------------------------------------------------------------------------
# Token management
# ---------------------------------------------------------------------------


def load_tokens() -> dict | None:
    """Load saved access token from disk."""
    if TOKEN_FILE.exists():
        try:
            data = json.loads(TOKEN_FILE.read_text(encoding="utf-8"))
            if data.get("access_token"):
                return data
        except (json.JSONDecodeError, KeyError):
            pass
    return None


def save_tokens(tokens: dict) -> None:
    """Persist access token to disk."""
    CONFIG_DIR.mkdir(parents=True, exist_ok=True)
    tokens["saved_at"] = datetime.now(timezone.utc).isoformat()
    TOKEN_FILE.write_text(json.dumps(tokens, indent=2), encoding="utf-8")
    TOKEN_FILE.chmod(0o600)
    print(f"  Token saved to {TOKEN_FILE}")


def get_access_token() -> str | None:
    """Resolve access token from env var or saved file."""
    env_token = os.environ.get("PROTOCOLS_IO_ACCESS_TOKEN")
    if env_token:
        return env_token
    tokens = load_tokens()
    if tokens and tokens.get("access_token"):
        return tokens["access_token"]
    return None


def token_login() -> str | None:
    """
    Login by pasting a client access token from protocols.io/developers.
    Verifies the token against the API and saves it locally.
    """
    print("\n  Paste your access token from https://www.protocols.io/developers")
    print("  (Log in → Your Applications → copy the 'Access Token')\n")
    token = getpass.getpass("  Access Token: ").strip()

    if not token:
        print("ERROR: No token provided.", file=sys.stderr)
        return None

    print("  Verifying token...")
    try:
        resp = requests.get(
            f"{API_V3}/session/profile",
            headers={"Authorization": f"Bearer {token}", "Accept": "application/json"},
            timeout=30,
        )
        data = resp.json()
    except Exception as e:
        print(f"ERROR: Could not verify token: {e}", file=sys.stderr)
        return None

    if resp.status_code != 200 or data.get("status_code") != 0:
        print(f"ERROR: Token rejected by protocols.io: {data.get('error_message', resp.text[:200])}", file=sys.stderr)
        return None

    save_tokens({"access_token": token, "token_type": "bearer"})
    user = data.get("user", {})
    print(f"  Logged in as: {user.get('name', 'unknown')} (@{user.get('username', '?')})")
    return token


# ---------------------------------------------------------------------------
# API helpers
# ---------------------------------------------------------------------------


def _headers(token: str | None = None) -> dict:
    """Build request headers with Bearer auth."""
    t = token or get_access_token()
    h = {"Accept": "application/json"}
    if t:
        h["Authorization"] = f"Bearer {t}"
    return h


def _api_get(url: str, params: dict | None = None, token: str | None = None) -> dict | None:
    """GET with client-side rate limiting and 429 handling (Retry-After)."""
    hdrs = _headers(token)

    for attempt in range(1, MAX_RETRIES_429 + 1):
        _rate_limiter.wait()

        try:
            resp = requests.get(url, headers=hdrs, params=params, timeout=30)
        except requests.RequestException as e:
            print(f"ERROR: Request failed: {e}", file=sys.stderr)
            return None

        if resp.status_code == 429:
            try:
                retry_after = int(resp.headers.get("Retry-After", "10"))
            except (TypeError, ValueError):
                retry_after = 10
            retry_after = max(1, min(retry_after, 120))
            print(
                f"  HTTP 429 Too Many Requests — retry in {retry_after}s "
                f"({attempt}/{MAX_RETRIES_429})",
                file=sys.stderr,
            )
            time.sleep(retry_after)
            continue

        try:
            data = resp.json()
        except (ValueError, requests.exceptions.JSONDecodeError):
            print(f"ERROR: Non-JSON response (HTTP {resp.status_code}): {resp.text[:200]}", file=sys.stderr)
            return None

        if data.get("status_code") == 1219:
            print("ERROR: Token expired. Run --login again to paste a new token.", file=sys.stderr)
            return None

        if resp.status_code != 200:
            msg = data.get("error_message", resp.text[:200]) if isinstance(data, dict) else resp.text[:200]
            print(f"API error {resp.status_code}: {msg}", file=sys.stderr)
            return None

        return data

    print("ERROR: Still rate-limited after retries.", file=sys.stderr)
    return None


# ---------------------------------------------------------------------------
# Core operations
# ---------------------------------------------------------------------------


def search_protocols(
    query: str,
    filter_type: str = "public",
    page_size: int = 10,
    page_id: int = 1,
    order_field: str = "activity",
    peer_reviewed: int | None = None,
) -> dict | None:
    """Search protocols.io for protocols matching a keyword query."""
    params: dict = {
        "filter": filter_type,
        "key": query,
        "order_field": order_field,
        "page_size": page_size,
        "page_id": page_id,
    }
    if peer_reviewed is not None:
        params["peer_reviewed"] = peer_reviewed
    return _api_get(f"{API_V3}/protocols", params=params)


def _parse_protocol_id(raw: str) -> str:
    """
    Extract a usable protocol identifier from various input formats:
    - Full URL:  https://www.protocols.io/view/some-protocol-slug-abc123  → some-protocol-slug-abc123
    - URI slug:  some-protocol-slug-abc123  → some-protocol-slug-abc123
    - Numeric:   30756  → 30756
    - DOI:       10.17504/protocols.io.abc123  → 10.17504/protocols.io.abc123
    """
    s = raw.strip().rstrip("/")
    if "protocols.io/view/" in s:
        s = s.split("protocols.io/view/")[-1]
    elif "protocols.io/api/" in s:
        pass
    s = s.split("?")[0].split("#")[0]
    return s


def get_protocol(protocol_id: str | int, content_format: str = "markdown") -> dict | None:
    """Retrieve full protocol detail by ID, URI, URL, or DOI."""
    pid = _parse_protocol_id(str(protocol_id))
    return _api_get(
        f"{API_V4}/protocols/{pid}",
        params={"content_format": content_format},
    )


def get_protocol_steps(protocol_id: str | int, content_format: str = "markdown") -> dict | None:
    """Retrieve protocol steps by ID, URI, URL, or DOI."""
    pid = _parse_protocol_id(str(protocol_id))
    return _api_get(
        f"{API_V4}/protocols/{pid}/steps",
        params={"content_format": content_format},
    )




# ---------------------------------------------------------------------------
# Formatters
# ---------------------------------------------------------------------------


def format_search_results(data: dict, query: str) -> str:
    """Render search results as markdown."""
    items = data.get("items", [])
    pagination = data.get("pagination", {})
    total = pagination.get("total_results", len(items))

    lines = [
        f"# Protocols.io Search: \"{query}\"\n",
        f"**{total} results found** (showing {len(items)})\n",
    ]

    for i, p in enumerate(items, 1):
        title = p.get("title", "Untitled")
        uri = p.get("uri", "")
        doi = p.get("doi", "")
        creator = p.get("creator", {}).get("name", "Unknown")
        published = p.get("published_on")
        pub_str = datetime.fromtimestamp(published, tz=timezone.utc).strftime("%Y-%m-%d") if published else "Draft"
        n_steps = p.get("number_of_steps") or p.get("stats", {}).get("number_of_steps", "?")
        url = f"https://www.protocols.io/view/{uri}" if uri else ""

        peer_reviewed = p.get("peer_reviewed")
        lines.append(f"## {i}. {title}\n")
        if peer_reviewed == 1:
            lines.append(f"- ✅ Peer-reviewed method")
        lines.append(f"- **Creator**: {creator}")
        lines.append(f"- **Published**: {pub_str}")
        lines.append(f"- **Steps**: {n_steps}")
        if doi:
            lines.append(f"- **DOI**: {doi}")
        if url:
            lines.append(f"- **URL**: {url}")
        lines.append("")

    lines.append(f"\n---\n{DISCLAIMER}")
    return "\n".join(lines)


def format_protocol_detail(data: dict) -> str:
    """Render a full protocol as markdown."""
    p = data.get("payload", data.get("protocol", data))
    title = p.get("title", "Untitled Protocol")
    doi = p.get("doi", "")
    uri = p.get("uri", "")
    creator = p.get("creator", {}).get("name", "Unknown")
    description = p.get("description", "")
    guidelines = p.get("guidelines", "")
    before_start = p.get("before_start", "")
    warning = p.get("warning", "")
    published = p.get("published_on")
    pub_str = datetime.fromtimestamp(published, tz=timezone.utc).strftime("%Y-%m-%d") if published else "Draft"

    authors = p.get("authors", [])
    author_str = ", ".join(a.get("name", "?") for a in authors) if authors else creator

    url = f"https://www.protocols.io/view/{uri}" if uri else ""

    lines = [
        f"# {title}\n",
        f"- **Authors**: {author_str}",
        f"- **Published**: {pub_str}",
    ]
    if doi:
        lines.append(f"- **DOI**: {doi}")
    if url:
        lines.append(f"- **URL**: {url}")

    stats = p.get("stats", {})
    n_steps = stats.get("number_of_steps") or p.get("number_of_steps", "?")
    if stats or p.get("number_of_steps"):
        lines.append(f"- **Views**: {stats.get('number_of_views', '?')} | "
                      f"**Steps**: {n_steps} | "
                      f"**Exports**: {stats.get('number_of_exports', '?')}")

    lines.append("")

    if description:
        lines.extend(["## Description\n", description, ""])
    if guidelines:
        lines.extend(["## Guidelines\n", guidelines, ""])
    if before_start:
        lines.extend(["## Before You Start\n", before_start, ""])
    if warning:
        lines.extend(["## Warnings\n", warning, ""])

    materials = p.get("materials", [])
    if materials:
        lines.append("## Materials\n")
        for m in materials:
            name = m.get("name", "Unknown reagent")
            vendor = m.get("vendor", {}).get("name", "")
            sku = m.get("sku", "")
            parts = [f"- **{name}**"]
            if vendor:
                parts.append(f"({vendor})")
            if sku:
                parts.append(f"[SKU: {sku}]")
            lines.append(" ".join(parts))
        lines.append("")

    steps = p.get("steps", [])
    if steps:
        lines.append("## Steps\n")
        for j, s in enumerate(steps, 1):
            section = s.get("section")
            if section:
                lines.append(f"### {_strip_html(section)}\n")
            step_text = s.get("step", "")
            if isinstance(step_text, str) and step_text.startswith("{"):
                try:
                    draft = json.loads(step_text)
                    blocks = draft.get("blocks", [])
                    step_text = "\n".join(b.get("text", "") for b in blocks)
                except json.JSONDecodeError:
                    pass
            lines.append(f"**Step {j}.**  {step_text}\n")
        lines.append("")

    lines.append(f"---\n{DISCLAIMER}")
    return "\n".join(lines)


def _strip_html(text: str) -> str:
    """Remove HTML tags from a string."""
    import re
    return re.sub(r"<[^>]+>", "", text).strip()


def format_steps(data: dict, protocol_id: str) -> str:
    """Render protocol steps as markdown."""
    # v4 /steps endpoint returns steps as a list under "payload";
    # fall back to "steps" key for backward-compat with mocked/demo data.
    payload = data.get("payload")
    if isinstance(payload, list):
        steps = payload
    else:
        steps = data.get("steps", [])
    lines = [
        f"# Protocol Steps — {protocol_id}\n",
        f"**{len(steps)} steps**\n",
    ]
    for j, s in enumerate(steps, 1):
        section = s.get("section")
        if section:
            lines.append(f"### {_strip_html(section)}\n")
        components = s.get("components", [])
        step_text = ""
        for comp in components:
            body = comp.get("source", {})
            if isinstance(body, dict):
                step_text += body.get("description", "")
            elif isinstance(body, str):
                step_text += body
        if not step_text:
            step_text = s.get("step", "(no content)")
            if isinstance(step_text, str) and step_text.startswith("{"):
                try:
                    draft = json.loads(step_text)
                    blocks = draft.get("blocks", [])
                    step_text = "\n".join(b.get("text", "") for b in blocks)
                except json.JSONDecodeError:
                    pass
        lines.append(f"**Step {j}.** {step_text}\n")

    lines.append(f"\n---\n{DISCLAIMER}")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Demo mode
# ---------------------------------------------------------------------------

def _load_demo_json(filename: str) -> dict:
    """Load a pre-cached demo JSON file from the demo/ directory."""
    path = DEMO_DIR / filename
    if not path.exists():
        print(f"ERROR: Demo file not found: {path}", file=sys.stderr)
        sys.exit(1)
    return json.loads(path.read_text(encoding="utf-8"))


def run_demo() -> None:
    """Run offline demo with pre-cached data."""
    print("\nProtocols.io Bridge — Demo Mode (offline)")
    print("=" * 50)

    demo_search = _load_demo_json("demo_search_results.json")
    demo_protocol = _load_demo_json("demo_protocol.json")

    print("\n--- Search Demo: \"RNA extraction\" (all results) ---\n")
    search_md = format_search_results(demo_search, "RNA extraction")
    print(search_md)

    print("\n--- Search Demo: \"RNA extraction\" --peer-reviewed ---\n")
    peer_reviewed_data = {
        **demo_search,
        "items": [p for p in demo_search.get("items", []) if p.get("peer_reviewed") == 1],
        "pagination": {
            **demo_search.get("pagination", {}),
            "total_results": sum(1 for p in demo_search.get("items", []) if p.get("peer_reviewed") == 1),
        },
    }
    pr_md = format_search_results(peer_reviewed_data, "RNA extraction (peer-reviewed)")
    print(pr_md)

    print("\n--- Protocol Detail Demo ---\n")
    detail_md = format_protocol_detail(demo_protocol)
    print(detail_md)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def _slugify(text: str, max_len: int = 60) -> str:
    """Turn a title or query into a safe filename slug."""
    import re
    s = text.lower().strip()
    s = re.sub(r"[^a-z0-9]+", "-", s)
    s = s.strip("-")[:max_len].rstrip("-")
    return s or "output"


def _dump_markdown(md: str, label: str) -> None:
    """Write markdown to an auto-named file in the current directory."""
    slug = _slugify(label)
    filename = f"{slug}.md"
    Path(filename).write_text(md, encoding="utf-8")
    print(f"  Saved to {filename}")


def _prompt_for_token() -> str | None:
    """Inline token prompt for commands that need auth but have no saved token."""
    print("  Get your token at: https://www.protocols.io/developers")
    print("  (Log in → Your Applications → copy the 'Access Token')\n")
    token = getpass.getpass("  Access Token (or press Enter to skip): ").strip()
    if not token:
        return None
    save_tokens({"access_token": token, "token_type": "bearer"})
    return token


def main() -> None:
    parser = argparse.ArgumentParser(
        description="protocols.io bridge — search, browse, and retrieve scientific protocols"
    )
    parser.add_argument("--login", action="store_true", help="Authenticate with access token")
    parser.add_argument("--search", type=str, help="Search protocols by keyword")
    parser.add_argument("--protocol", type=str, help="Retrieve full protocol by ID, URI, or DOI")
    parser.add_argument("--steps", type=str, help="Retrieve protocol steps by ID, URI, or DOI")
    parser.add_argument("--demo", action="store_true", help="Run offline demo with pre-cached data")
    parser.add_argument("--dump", action="store_true", help="Save output as a markdown file in the current directory")
    parser.add_argument("--page-size", type=int, default=10, help="Results per page (1-100)")
    parser.add_argument("--page", type=int, default=1, help="Page number")
    parser.add_argument("--filter", type=str, default="public",
                        choices=["public", "user_public", "user_private", "shared_with_user"],
                        help="Protocol filter type")
    parser.add_argument("--peer-reviewed", action="store_const", const=1, default=None,
                        help="Filter to peer-reviewed protocols only (omit to show all)")

    args = parser.parse_args()

    if args.demo:
        run_demo()
        return

    if args.login:
        result = token_login()
        if result:
            print("\n  Authentication successful!")
        else:
            print("\n  Authentication failed.", file=sys.stderr)
            sys.exit(1)
        return

    if args.search:
        token = get_access_token()
        if not token:
            print("  No access token found. Run --login first, or paste a token now.\n")
            token = _prompt_for_token()
            if not token:
                print("ERROR: Cannot search without an access token.", file=sys.stderr)
                sys.exit(1)

        with Spinner(f"Searching protocols.io for \"{args.search}\""):
            data = search_protocols(
                args.search,
                filter_type=args.filter,
                page_size=args.page_size,
                page_id=args.page,
                peer_reviewed=args.peer_reviewed,
            )
        if not data:
            print("ERROR: Search failed.", file=sys.stderr)
            sys.exit(1)

        report = format_search_results(data, args.search)
        print(report)
        if args.dump:
            _dump_markdown(report, f"search-{args.search}")
        return

    if args.protocol:
        with Spinner(f"Retrieving protocol {args.protocol}"):
            data = get_protocol(args.protocol)
        if not data:
            print("ERROR: Could not retrieve protocol.", file=sys.stderr)
            sys.exit(1)

        report = format_protocol_detail(data)
        print(report)
        if args.dump:
            p = data.get("payload", data.get("protocol", data))
            _dump_markdown(report, p.get("title", args.protocol))
        return

    if args.steps:
        with Spinner(f"Retrieving steps for protocol {args.steps}"):
            data = get_protocol_steps(args.steps)
        if not data:
            print("ERROR: Could not retrieve steps.", file=sys.stderr)
            sys.exit(1)

        report = format_steps(data, args.steps)
        print(report)
        if args.dump:
            _dump_markdown(report, f"steps-{args.steps}")
        return

    parser.print_help()


if __name__ == "__main__":
    main()
