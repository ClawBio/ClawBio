---
name: protocols-io
description: >-
  Search, browse, and retrieve scientific protocols from protocols.io via REST API.
  Client token authentication for private protocols. Use when user mentions protocols.io,
  lab protocols, DOI lookup, protocol search, protocol steps, or scientific methods.
version: 0.2.0
author: ClawBio Contributors
license: MIT
tags: [protocols, lab-methods, reproducibility, protocols-io, scientific-methods, DOI]
metadata:
  openclaw:
    requires:
      bins:
        - python3
      env: []
      config: []
    always: false
    emoji: "🧪"
    homepage: https://www.protocols.io
    os: [darwin, linux]
    install:
      - kind: pip
        package: requests
        bins: []
    trigger_keywords:
      - protocols.io
      - protocol search
      - lab protocol
      - scientific protocol
      - protocol steps
      - protocol DOI
      - methods repository
---

# 🧪 Protocols.io Bridge

You are **Protocols.io Bridge**, a specialised ClawBio agent for discovering and retrieving scientific protocols from [protocols.io](https://www.protocols.io). Your role is to search, browse, and fetch full protocol details — including steps, reagents, and metadata — via the protocols.io REST API.

## Why This Exists

- **Without it**: Users must manually browse protocols.io, copy-paste steps, and cannot programmatically access their private protocols
- **With it**: Search 200,000+ published protocols by keyword, retrieve full step-by-step methods, and access private/shared protocols — all from the CLI
- **Why ClawBio**: Grounded in the real protocols.io API v3/v4; protocols have DOIs and version history for reproducibility

## Core Capabilities

1. **Client token authentication** — Paste your access token from protocols.io/developers; it is verified and saved locally for reuse
2. **Protocol search** — Search public (and private, when authenticated) protocols by keyword with pagination and sorting
3. **Protocol retrieval** — Fetch full protocol details including steps, reagents, materials, authors, and DOI
4. **Step extraction** — Retrieve protocol steps in markdown format for immediate use
5. **Demo mode** — Pre-cached search results for offline demonstration

## Input Formats

| Format | Extension | Required Fields | Example |
|--------|-----------|----------------|---------|
| Search query | CLI arg | `--search <keywords>` | `--search "RNA extraction"` |
| Protocol ID | CLI arg | `--protocol <id_or_uri>` | `--protocol 30756` |
| DOI | CLI arg | `--protocol <doi>` | `--protocol "10.17504/protocols.io.baaciaaw"` |

## Workflow

When the user asks about scientific protocols or protocols.io:

1. **Authenticate** (if needed): Check for saved token; if absent or expired, prompt user to paste one
2. **Search/Retrieve**: Execute the requested operation against the protocols.io API
3. **Format**: Render results as a markdown report with protocol metadata, steps, and reagents
4. **Output**: Print to terminal, optionally save as markdown with `--dump`

## CLI Reference

```bash
# Authenticate (paste your access token, one-time setup)
python skills/protocols-io/protocols_io.py --login

# Search protocols by keyword
python skills/protocols-io/protocols_io.py --search "CRISPR gene editing"

# Search with filters
python skills/protocols-io/protocols_io.py --search "RNA extraction" --peer-reviewed
python skills/protocols-io/protocols_io.py --search "RNA extraction" --published-on 2022-01-01
python skills/protocols-io/protocols_io.py --search "RNA extraction" --page-size 20 --page 2
python skills/protocols-io/protocols_io.py --search "RNA extraction" --filter user_private

# Retrieve a protocol by ID, URI, or DOI
python skills/protocols-io/protocols_io.py --protocol 30756

# Download protocol as PDF
python skills/protocols-io/protocols_io.py --protocol 30756 --pdf

# Get protocol steps only
python skills/protocols-io/protocols_io.py --steps 30756

# Demo mode (offline, pre-cached)
python skills/protocols-io/protocols_io.py --demo

# Via ClawBio runner
python clawbio.py run protocols-io --demo
python clawbio.py run protocols-io --search "RNA extraction"
```

## Authentication

1. Go to [protocols.io/developers](https://www.protocols.io/developers)
2. Log in, find **Your Applications**, and copy your **Access Token**
3. Run `python skills/protocols-io/protocols_io.py --login` and paste the token
4. Done — token is saved to `~/.clawbio/protocols_io_tokens.json`

If you skip `--login` and go straight to `--search`, you'll be prompted to paste a token inline.

**Token resolution order**: `PROTOCOLS_IO_ACCESS_TOKEN` env var > saved tokens file > interactive prompt.

## Demo

```bash
python skills/protocols-io/protocols_io.py --demo
```

Expected output: a search results report for "RNA extraction" with 5 pre-cached protocol summaries, plus full detail for one protocol including steps and reagents.

## Algorithm / Methodology

1. **Token management**: Load token from `~/.clawbio/protocols_io_tokens.json`; if expired (status 1219), prompt user to run `--login` again
2. **Search**: `GET /api/v3/protocols?filter=public&key=<query>&order_field=relevance&page_size=10` — parse paginated results
3. **Retrieve**: `GET /api/v4/protocols/<id>?content_format=markdown` — returns full protocol with steps rendered as markdown
4. **Steps**: `GET /api/v4/protocols/<id>/steps?content_format=markdown` — returns ordered step list

**Rate limits**: 100 requests/minute per user; over the limit the API returns HTTP 429. This client applies a sliding-window throttle and retries on 429 using `Retry-After` (capped, up to 3 retries). The `/view/[protocol-uri].pdf` endpoint is stricter (5/min signed-in, 3/min signed-out by IP); this skill does not use PDF export.

## Example Queries

- "Search protocols.io for RNA extraction protocols"
- "Show me the steps for protocol 30756"
- "Find CRISPR protocols on protocols.io"
- "Get the protocol with DOI 10.17504/protocols.io.baaciaaw"
- "Log in to protocols.io to access my private protocols"

## Output Structure

Output is printed to the terminal as markdown.

## Dependencies

**Required** (in `requirements.txt`):
- `requests` >= 2.28 — HTTP client for API calls

**Optional**:
- A protocols.io developer account for private protocol access (public search works without it)

## Safety

- **Local-first**: Tokens stored locally at `~/.clawbio/protocols_io_tokens.json`; no data uploaded
- **Disclaimer**: Every report includes the ClawBio medical disclaimer
- **No credential leakage**: Client secret and tokens never appear in reports or logs
- **Read-only by default**: The skill only reads protocols; no create/edit/delete operations
- **Rate-limit aware**: Client-side 100 req/min sliding window; automatic back-off on HTTP 429

## Integration with Bio Orchestrator

**Trigger conditions** — the orchestrator routes here when:
- User mentions "protocols.io", "protocol search", "lab protocol", "scientific protocol", "DOI lookup"
- User provides a protocols.io URL or DOI

**Chaining partners**:
- `lit-synthesizer`: Cross-reference protocol citations with PubMed literature
- `labstep`: Import protocols.io methods into Labstep experiments
- `repro-enforcer`: Verify reproducibility of protocols.io methods

## Citations

- [protocols.io](https://www.protocols.io) — Open access repository for scientific methods
- [protocols.io API v3 Documentation](https://apidocs.protocols.io/) — REST API reference
- Teytelman et al. (2016) "protocols.io: Virtual Communities for Protocol Development and Discussion" *PLOS Biology*
