# Knowledge Guide Skill — Design Spec

**Date:** 2026-03-13
**Branch:** `feature/knowledge-guide` (off `feature/agent-reporting`)
**Status:** Approved design, pending implementation

## Purpose

After users receive ClawBio reports (PGx, PRS, equity, etc.), they often want to understand the underlying analysis — "what is variant calling?", "how does polygenic risk scoring work?". The knowledge-guide skill grounds these explanations in Galaxy Training Network (GTN) tutorials rather than generating answers from scratch, ensuring every explanation traces back to authoritative training material.

## Requirements

- **Standalone query mode**: free-text questions AND structured flags (`--topic`, `--tool`, `--concept`)
- **Report integration**: other skills embed a "Learn More" section by calling `api.py`
- **Deep content**: inline explanations pulled from GTN tutorial content, not just links
- **GTN as sole knowledge source**: no PubMed, no other APIs
- **Pre-cached mappings**: report "Learn More" sections use cached data (zero API latency)
- **Live fetch for standalone queries**: deep tutorial content fetched on demand

## Architecture

### File Structure

```
skills/knowledge-guide/
├── SKILL.md                        # Methodology doc
├── knowledge_guide.py              # CLI entry point (query engine)
├── api.py                          # Importable run() for other skills
├── gtn_client.py                   # GTN API fetcher (topics, tutorials, tool map)
├── gtn_cache_builder.py            # Builds gtn_cache.json + skill_recommendations.json
├── query_engine.py                 # Free-text → structured query matching logic
├── gtn_cache.json                  # Pre-built cache (committed to repo)
├── skill_recommendations.json      # Skill → tutorial mappings (committed to repo)
├── demo/
│   └── demo_queries.json           # Sample queries + expected output
└── tests/
    └── test_knowledge_guide.py
```

### Data Flow

```
GTN API (topics.json, top-tools.json, topic/*.json)
       ↓
gtn_cache_builder.py  ──→  gtn_cache.json (metadata, learning objectives, tool mappings)
                       ──→  skill_recommendations.json (clawbio skill → top 5 tutorials)
       ↓
knowledge_guide.py reads cache for:
  - Standalone queries (--query, --topic, --tool, --concept)
  - Deep content: live fetch via gtn_client.py when user wants full tutorial detail
       ↓
api.py exposes:
  - run(query=...) for standalone use
  - get_learn_more(skill_name=...) for report integration
```

## Cache Strategy

### `gtn_cache.json` — Main Knowledge Index

Built by `gtn_cache_builder.py` fetching three GTN API endpoints:

| Source | What We Store |
|---|---|
| `/api/topics.json` | 43 topics with names, summaries |
| `/api/topics/{id}.json` | Per-topic tutorial list: title, time estimate, difficulty, learning objectives, tools used, URL |
| `/api/top-tools.json` | Tool ID → tutorial mappings (reverse index) |

Stores ~400-500 tutorials with metadata. Full tutorial content is NOT cached (too large, goes stale) — fetched live on demand via `--deep`.

### `skill_recommendations.json` — Skill-to-Tutorial Mappings

Hand-seeded concepts + auto-enriched tutorial matches per ClawBio skill:

```json
{
  "pharmgx-reporter": {
    "concepts": ["pharmacogenomics", "CYP2D6", "drug metabolism", "CPIC guidelines"],
    "gtn_topics": ["variant-analysis"],
    "tutorials": [
      {
        "id": "variant-analysis/pharmacogenomics",
        "title": "Pharmacogenomics analysis",
        "time": "2h",
        "relevance": "Explains the variant calling pipeline behind PGx reports"
      }
    ]
  },
  "gwas-prs": {
    "concepts": ["polygenic risk", "GWAS", "effect sizes", "linkage disequilibrium"],
    "gtn_topics": ["variant-analysis", "statistics"],
    "tutorials": []
  }
}
```

The builder auto-matches by scanning each ClawBio skill's `SKILL.md` for keywords, then matching against GTN tutorial titles/descriptions. The `concepts` list is hand-seeded per skill to anchor matching.

### Refresh

`python gtn_cache_builder.py --refresh` re-fetches from GTN and rebuilds both files. GTN updates monthly-ish. Committed cache means the skill works fully offline.

## Query Engine

### Input Modes

| Flag | Example | Resolution |
|---|---|---|
| `--query "free text"` | `--query "what is variant calling?"` | Keyword extraction → match against tutorial titles, descriptions, learning objectives, topic names |
| `--topic <id>` | `--topic variant-analysis` | Direct lookup by topic ID |
| `--tool <name>` | `--tool fastqc` | Reverse lookup via tool→tutorial index |
| `--concept <term>` | `--concept "polygenic risk"` | Fuzzy match against concepts + tutorial metadata |

### Free-Text Matching (`query_engine.py`)

1. Tokenize and normalize (lowercase, strip stopwords)
2. Score each tutorial: title match, description match, learning objectives match, topic name match, tool name match
3. Weight title and learning objectives highest
4. Return top 5 ranked tutorials

Same scoring pattern as `galaxy_bridge/tool_recommender.py` — keyword matching with weighted signals. No ML, no embeddings.

### Output per Match

- Title, topic, URL, time estimate, difficulty
- Learning objectives (from cache)
- "Why this is relevant" line explaining the match

### Deep Content Pull

When `--deep` is passed, `gtn_client.py` fetches `/api/topics/{topic}/tutorials/{id}/tutorial.json` live and extracts key sections (introduction, methodology steps, key concepts). Only for standalone queries, never for report "Learn More" sections.

## Report Integration ("Learn More")

Other skills call `api.py`'s `get_learn_more(skill_name)`:

```python
from skills.knowledge_guide.api import get_learn_more

learn_more = get_learn_more("gwas-prs")
# Returns:
# {
#   "section_title": "Learn More",
#   "concepts": ["polygenic risk", "GWAS", "effect sizes"],
#   "tutorials": [
#     {"title": "...", "url": "...", "time": "2h", "relevance": "..."}
#   ],
#   "html": "<div class='learn-more'>...</div>"
# }
```

Key decisions:
- Reads only from `skill_recommendations.json` — zero network calls, zero latency
- Returns structured data (for `result.json`) and pre-rendered HTML (for HtmlReportBuilder)
- Unknown skill names return empty section gracefully
- Skills opt-in with ~3 lines of code (import, call, append)

HTML block: collapsible "Learn More" details section with tutorial cards (title, time, difficulty, relevance note). Links open the full GTN tutorial.

## CLI Interface

```bash
# Free-text query
python skills/knowledge-guide/knowledge_guide.py \
  --query "what is variant calling and why does it matter?" --output /tmp/kg_out

# Structured lookups
python skills/knowledge-guide/knowledge_guide.py --topic variant-analysis --output /tmp/kg_out
python skills/knowledge-guide/knowledge_guide.py --tool fastqc --output /tmp/kg_out
python skills/knowledge-guide/knowledge_guide.py --concept "polygenic risk" --output /tmp/kg_out

# Deep pull (fetches full tutorial content live)
python skills/knowledge-guide/knowledge_guide.py --query "RNA-seq DE" --deep --output /tmp/kg_out

# Skill "Learn More" lookup
python skills/knowledge-guide/knowledge_guide.py --skill gwas-prs --output /tmp/kg_out

# Refresh cache from GTN API
python skills/knowledge-guide/gtn_cache_builder.py --refresh

# Demo (offline, no API calls)
python skills/knowledge-guide/knowledge_guide.py --demo --output /tmp/kg_demo
```

## Output Structure

```
/tmp/kg_out/
├── report.md          # Markdown report with explanations + tutorial links
├── report.html        # Styled HTML report (HtmlReportBuilder)
├── result.json        # Standardized envelope (skill, version, summary, data)
└── tutorials/         # Deep-pulled tutorial content (only with --deep)
    └── variant-calling.md
```

### `result.json` Data Section

```json
{
  "query": "what is variant calling?",
  "mode": "free-text",
  "matches": 5,
  "tutorials": [],
  "deep_content_fetched": false
}
```

## Demo Mode

Bundled subset of `gtn_cache.json` with 3-4 pre-cached tutorials across variant analysis, transcriptomics, and single-cell. Runs a sample query and produces full output. No network required.

## What This Is NOT

- Not a chatbot hallucinating biology — every explanation traces to a GTN tutorial
- Not a replacement for tutorials — a signpost that pulls the right context
- Not a PubMed search tool — GTN is the sole knowledge source
- Not modifying existing skills — they opt in when ready
