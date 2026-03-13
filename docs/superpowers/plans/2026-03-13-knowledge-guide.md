# Knowledge Guide Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a GTN-backed educational skill that explains bioinformatics concepts grounded in Galaxy Training Network tutorials, with both standalone query mode and report-embedded "Learn More" integration.

**Architecture:** Two-layer design — a cache builder fetches GTN API data into committed JSON files, and a query engine reads the cache for fast lookups. Deep tutorial content is fetched live on demand. Other skills call `api.py` to embed "Learn More" sections using pre-computed skill→tutorial mappings.

**Tech Stack:** Python 3, requests (HTTP), JSON (cache), ClawBio common library (report.py, html_report.py)

---

## Chunk 1: Foundation — GTN Client, Cache Builder, and Cache Files

### Task 1: Create branch and skill directory

**Files:**
- Create: `skills/knowledge-guide/` (directory)

- [ ] **Step 1: Create feature branch off current branch**

```bash
git checkout -b feature/knowledge-guide
```

- [ ] **Step 2: Create skill directory structure**

```bash
mkdir -p skills/knowledge-guide/demo
mkdir -p skills/knowledge-guide/tests
```

Note: No commit here — git cannot track empty directories. The scaffold will be committed with Task 2's first files.

---

### Task 2: Build GTN client (`gtn_client.py`)

**Files:**
- Create: `skills/knowledge-guide/gtn_client.py`
- Test: `skills/knowledge-guide/tests/test_gtn_client.py`

- [ ] **Step 1: Write failing tests for GTN client**

```python
# skills/knowledge-guide/tests/test_gtn_client.py
"""Tests for the GTN API client."""
from __future__ import annotations
import json
import pytest
from pathlib import Path
import sys

SKILL_DIR = Path(__file__).resolve().parent.parent
if str(SKILL_DIR) not in sys.path:
    sys.path.insert(0, str(SKILL_DIR))


class TestFetchTopics:
    """Test fetch_topics returns structured topic data."""

    def test_returns_dict(self, monkeypatch):
        from gtn_client import fetch_topics

        sample = [
            {"name": "variant-analysis", "title": "Variant Analysis",
             "summary": "Detecting genetic variants", "material": []}
        ]
        monkeypatch.setattr(
            "gtn_client.requests.get",
            lambda url, **kw: type("R", (), {"json": lambda s: sample, "raise_for_status": lambda s: None})(),
        )
        result = fetch_topics()
        assert isinstance(result, list)
        assert result[0]["name"] == "variant-analysis"

    def test_handles_api_error(self, monkeypatch):
        from gtn_client import fetch_topics
        import requests as req

        def raise_error(url, **kw):
            raise req.ConnectionError("offline")

        monkeypatch.setattr("gtn_client.requests.get", raise_error)
        with pytest.raises(Exception):
            fetch_topics()


class TestFetchTopicDetail:
    """Test fetch_topic_detail returns tutorials for a topic."""

    def test_returns_tutorials(self, monkeypatch):
        from gtn_client import fetch_topic_detail

        sample = {
            "name": "variant-analysis",
            "title": "Variant Analysis",
            "materials": [
                {
                    "title": "Calling very rare variants",
                    "name": "rare-variants",
                    "hands_on": True,
                    "type": "tutorial",
                    "time_estimation": "3H",
                    "level": "Intermediate",
                    "objectives": ["Process duplex sequencing data"],
                    "tools": ["fastqc", "bwa_mem"],
                    "url": "/topics/variant-analysis/tutorials/rare-variants/tutorial.html",
                }
            ],
        }
        monkeypatch.setattr(
            "gtn_client.requests.get",
            lambda url, **kw: type("R", (), {"json": lambda s: sample, "raise_for_status": lambda s: None})(),
        )
        result = fetch_topic_detail("variant-analysis")
        assert result["name"] == "variant-analysis"
        assert len(result["materials"]) == 1
        assert result["materials"][0]["title"] == "Calling very rare variants"


class TestFetchToolMap:
    """Test fetch_tool_tutorial_map returns tool→tutorial index."""

    def test_returns_mapping(self, monkeypatch):
        from gtn_client import fetch_tool_tutorial_map

        sample = {
            "devteam/fastqc/fastqc": {
                "tool_id": [["toolshed.g2.bx.psu.edu/repos/devteam/fastqc/fastqc/0.74+galaxy1", "0.74+galaxy1"]],
                "tutorials": [
                    ["assembly/general-introduction", "An Introduction to Genome Assembly", "Assembly",
                     "/topics/assembly/tutorials/general-introduction/tutorial.html"]
                ],
            }
        }
        monkeypatch.setattr(
            "gtn_client.requests.get",
            lambda url, **kw: type("R", (), {"json": lambda s: sample, "raise_for_status": lambda s: None})(),
        )
        result = fetch_tool_tutorial_map()
        assert "devteam/fastqc/fastqc" in result
        assert len(result["devteam/fastqc/fastqc"]["tutorials"]) == 1


class TestFetchTutorialContent:
    """Test fetch_tutorial_content returns full tutorial body."""

    def test_returns_content(self, monkeypatch):
        from gtn_client import fetch_tutorial_content

        sample = {
            "name": "rare-variants",
            "title": "Calling very rare variants",
            "content": "# Introduction\nThis tutorial covers...",
        }
        monkeypatch.setattr(
            "gtn_client.requests.get",
            lambda url, **kw: type("R", (), {"json": lambda s: sample, "raise_for_status": lambda s: None})(),
        )
        result = fetch_tutorial_content("variant-analysis", "rare-variants")
        assert "content" in result
        assert "Introduction" in result["content"]

    def test_returns_none_on_404(self, monkeypatch):
        from gtn_client import fetch_tutorial_content
        import requests as req

        def raise_404(url, **kw):
            resp = req.Response()
            resp.status_code = 404
            raise req.HTTPError(response=resp)

        monkeypatch.setattr("gtn_client.requests.get", raise_404)
        result = fetch_tutorial_content("nonexistent", "nonexistent")
        assert result is None
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `python -m pytest skills/knowledge-guide/tests/test_gtn_client.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'gtn_client'`

- [ ] **Step 3: Implement GTN client**

```python
# skills/knowledge-guide/gtn_client.py
"""GTN API client — fetches topics, tutorials, and tool mappings from the
Galaxy Training Network (https://training.galaxyproject.org).

All functions return parsed JSON dicts. No caching here — that's handled
by gtn_cache_builder.py.
"""
from __future__ import annotations

import requests

GTN_BASE = "https://training.galaxyproject.org/training-material/api"


def fetch_topics() -> list[dict]:
    """Fetch all GTN topics (/api/topics.json)."""
    resp = requests.get(f"{GTN_BASE}/topics.json", timeout=30)
    resp.raise_for_status()
    return resp.json()


def fetch_topic_detail(topic_id: str) -> dict:
    """Fetch full detail for a single topic (/api/topics/{id}.json).

    Returns dict with 'name', 'title', 'summary', 'materials' (list of
    tutorial metadata including title, time_estimation, objectives, tools).
    """
    resp = requests.get(f"{GTN_BASE}/topics/{topic_id}.json", timeout=30)
    resp.raise_for_status()
    return resp.json()


def fetch_tool_tutorial_map() -> dict:
    """Fetch the tool→tutorial reverse index (/api/top-tools.json).

    Returns dict keyed by tool repo path (e.g. 'devteam/fastqc/fastqc'),
    each value has 'tool_id' (list of [full_id, version]) and 'tutorials'
    (list of [path, title, topic, url]).
    """
    resp = requests.get(f"{GTN_BASE}/top-tools.json", timeout=60)
    resp.raise_for_status()
    return resp.json()


def fetch_tutorial_content(topic_id: str, tutorial_id: str) -> dict | None:
    """Fetch full tutorial content for deep pull.

    Returns dict with 'name', 'title', 'content' (HTML/markdown body),
    or None if the tutorial is not found.
    """
    url = f"{GTN_BASE}/topics/{topic_id}/tutorials/{tutorial_id}/tutorial.json"
    try:
        resp = requests.get(url, timeout=30)
        resp.raise_for_status()
        return resp.json()
    except requests.HTTPError:
        return None
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `python -m pytest skills/knowledge-guide/tests/test_gtn_client.py -v`
Expected: All 5 tests PASS

- [ ] **Step 5: Commit**

```bash
git add skills/knowledge-guide/gtn_client.py skills/knowledge-guide/tests/test_gtn_client.py
git commit -m "feat(knowledge-guide): add GTN API client with tests"
```

---

### Task 3: Build cache builder (`gtn_cache_builder.py`)

**Files:**
- Create: `skills/knowledge-guide/gtn_cache_builder.py`
- Test: `skills/knowledge-guide/tests/test_cache_builder.py`

- [ ] **Step 1: Write failing tests for cache builder**

```python
# skills/knowledge-guide/tests/test_cache_builder.py
"""Tests for the GTN cache builder."""
from __future__ import annotations
import json
import sys
import tempfile
from pathlib import Path

import pytest

SKILL_DIR = Path(__file__).resolve().parent.parent
if str(SKILL_DIR) not in sys.path:
    sys.path.insert(0, str(SKILL_DIR))

# --- Sample GTN API responses for mocking ---

SAMPLE_TOPICS = [
    {"name": "variant-analysis", "title": "Variant Analysis",
     "summary": "Genetic variant detection", "material": []},
    {"name": "transcriptomics", "title": "Transcriptomics",
     "summary": "Gene expression analysis", "material": []},
]

SAMPLE_TOPIC_DETAIL = {
    "variant-analysis": {
        "name": "variant-analysis",
        "title": "Variant Analysis",
        "summary": "Genetic variant detection",
        "materials": [
            {
                "title": "Calling very rare variants",
                "name": "rare-variants",
                "hands_on": True,
                "type": "tutorial",
                "time_estimation": "3H",
                "level": "Intermediate",
                "objectives": ["Process duplex sequencing data", "Identify rare variants"],
                "tools": ["fastqc", "bwa_mem"],
            }
        ],
    },
    "transcriptomics": {
        "name": "transcriptomics",
        "title": "Transcriptomics",
        "summary": "Gene expression analysis",
        "materials": [
            {
                "title": "Reference-based RNA-Seq data analysis",
                "name": "ref-based",
                "hands_on": True,
                "type": "tutorial",
                "time_estimation": "4H",
                "level": "Intermediate",
                "objectives": ["Perform RNA-seq alignment", "Identify differentially expressed genes"],
                "tools": ["hisat2", "featurecounts", "deseq2"],
            }
        ],
    },
}

SAMPLE_TOOL_MAP = {
    "devteam/fastqc/fastqc": {
        "tool_id": [["toolshed.g2.bx.psu.edu/repos/devteam/fastqc/fastqc/0.74", "0.74"]],
        "tutorials": [
            ["variant-analysis/rare-variants", "Calling very rare variants",
             "Variant Analysis", "/topics/variant-analysis/tutorials/rare-variants/tutorial.html"],
        ],
    },
}


class TestBuildGtnCache:
    """Test build_gtn_cache produces a valid cache file."""

    def test_writes_cache_file(self, monkeypatch, tmp_path):
        from gtn_cache_builder import build_gtn_cache
        import gtn_client

        monkeypatch.setattr(gtn_client, "fetch_topics", lambda: SAMPLE_TOPICS)
        monkeypatch.setattr(gtn_client, "fetch_topic_detail",
                            lambda tid: SAMPLE_TOPIC_DETAIL[tid])
        monkeypatch.setattr(gtn_client, "fetch_tool_tutorial_map", lambda: SAMPLE_TOOL_MAP)

        cache_path = tmp_path / "gtn_cache.json"
        build_gtn_cache(output_path=cache_path)

        assert cache_path.exists()
        data = json.loads(cache_path.read_text())
        assert "topics" in data
        assert "tool_index" in data
        assert len(data["topics"]) == 2
        # Check tutorials are nested
        va = next(t for t in data["topics"] if t["name"] == "variant-analysis")
        assert len(va["tutorials"]) == 1
        assert va["tutorials"][0]["title"] == "Calling very rare variants"

    def test_atomic_write_preserves_old_on_failure(self, monkeypatch, tmp_path):
        """If fetching fails mid-build, old cache is preserved."""
        from gtn_cache_builder import build_gtn_cache
        import gtn_client

        # Write an existing cache
        cache_path = tmp_path / "gtn_cache.json"
        cache_path.write_text('{"old": true}')

        # Make fetch_topics fail
        monkeypatch.setattr(gtn_client, "fetch_topics",
                            lambda: (_ for _ in ()).throw(ConnectionError("offline")))

        with pytest.raises(ConnectionError):
            build_gtn_cache(output_path=cache_path)

        # Old cache should be preserved
        assert json.loads(cache_path.read_text()) == {"old": True}


class TestBuildSkillRecommendations:
    """Test build_skill_recommendations produces skill→tutorial mappings."""

    def test_writes_recommendations(self, monkeypatch, tmp_path):
        from gtn_cache_builder import build_skill_recommendations

        # Provide a pre-built cache
        cache_data = {
            "topics": [
                {
                    "name": "variant-analysis",
                    "title": "Variant Analysis",
                    "summary": "Genetic variant detection",
                    "tutorials": [
                        {
                            "title": "Calling very rare variants",
                            "name": "rare-variants",
                            "topic": "variant-analysis",
                            "time_estimation": "3H",
                            "level": "Intermediate",
                            "objectives": ["Process duplex sequencing data"],
                            "tools": ["fastqc", "bwa_mem"],
                            "url": "/topics/variant-analysis/tutorials/rare-variants/tutorial.html",
                        }
                    ],
                }
            ],
            "tool_index": {},
        }
        cache_path = tmp_path / "gtn_cache.json"
        cache_path.write_text(json.dumps(cache_data))

        recs_path = tmp_path / "skill_recommendations.json"
        build_skill_recommendations(cache_path=cache_path, output_path=recs_path)

        assert recs_path.exists()
        recs = json.loads(recs_path.read_text())
        # Should have entries for ClawBio skills
        assert isinstance(recs, dict)
        assert len(recs) > 0
        # Each entry should have concepts and tutorials
        for skill_name, entry in recs.items():
            assert "concepts" in entry
            assert "tutorials" in entry
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `python -m pytest skills/knowledge-guide/tests/test_cache_builder.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'gtn_cache_builder'`

- [ ] **Step 3: Implement cache builder**

```python
# skills/knowledge-guide/gtn_cache_builder.py
#!/usr/bin/env python3
"""
gtn_cache_builder.py — Build the GTN knowledge cache for ClawBio
================================================================
Fetches topics, tutorials, and tool mappings from the Galaxy Training
Network API and writes two committed JSON files:

  gtn_cache.json              — full topic + tutorial metadata
  skill_recommendations.json  — ClawBio skill → relevant tutorials

Usage:
    python gtn_cache_builder.py --refresh
"""
from __future__ import annotations

import argparse
import json
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

SKILL_DIR = Path(__file__).resolve().parent
if str(SKILL_DIR) not in sys.path:
    sys.path.insert(0, str(SKILL_DIR))

import gtn_client

GTN_TUTORIAL_BASE = "https://training.galaxyproject.org/training-material"

# Hand-seeded concept anchors per ClawBio skill.
# The builder matches these + SKILL.md keywords against GTN tutorials.
SKILL_CONCEPTS: dict[str, dict] = {
    "pharmgx-reporter": {
        "concepts": ["pharmacogenomics", "CYP2D6", "CYP2C19", "drug metabolism",
                      "CPIC guidelines", "star allele", "diplotype"],
        "gtn_topics": ["variant-analysis"],
    },
    "gwas-prs": {
        "concepts": ["polygenic risk score", "GWAS", "effect size",
                      "linkage disequilibrium", "PGS Catalog", "genome-wide association"],
        "gtn_topics": ["variant-analysis", "statistics"],
    },
    "equity-scorer": {
        "concepts": ["genomic diversity", "population genetics", "FST",
                      "heterozygosity", "PCA", "health equity"],
        "gtn_topics": ["variant-analysis", "sequence-analysis"],
    },
    "genome-compare": {
        "concepts": ["identity by state", "IBS", "genome comparison",
                      "ancestry estimation", "admixture"],
        "gtn_topics": ["variant-analysis"],
    },
    "scrna-orchestrator": {
        "concepts": ["single-cell RNA-seq", "scRNA-seq", "Scanpy", "clustering",
                      "UMAP", "marker genes", "doublet detection", "leiden"],
        "gtn_topics": ["single-cell", "transcriptomics"],
    },
    "clinpgx": {
        "concepts": ["pharmacogenomics", "gene-drug interaction", "PharmGKB",
                      "CPIC", "clinical annotation", "FDA drug label"],
        "gtn_topics": ["variant-analysis"],
    },
    "gwas-lookup": {
        "concepts": ["GWAS", "variant annotation", "rsID", "eQTL",
                      "PheWAS", "Open Targets", "GTEx"],
        "gtn_topics": ["variant-analysis"],
    },
    "nutrigx_advisor": {
        "concepts": ["nutrigenomics", "MTHFR", "folate metabolism",
                      "caffeine metabolism", "lactose tolerance", "vitamin D"],
        "gtn_topics": ["variant-analysis"],
    },
    "rnaseq-de": {
        "concepts": ["RNA-seq", "differential expression", "DESeq2",
                      "volcano plot", "pseudo-bulk", "count matrix"],
        "gtn_topics": ["transcriptomics"],
    },
    "galaxy-bridge": {
        "concepts": ["Galaxy", "bioinformatics tools", "workflows",
                      "tool discovery", "BioBlend"],
        "gtn_topics": ["galaxy-interface", "introduction"],
    },
    "claw-ancestry-pca": {
        "concepts": ["ancestry", "PCA", "population structure",
                      "admixture", "SGDP"],
        "gtn_topics": ["variant-analysis", "evolution"],
    },
}


def build_gtn_cache(output_path: Path | None = None) -> dict:
    """Fetch GTN data and build the cache. Atomic write — old cache
    preserved on failure.

    Args:
        output_path: Where to write gtn_cache.json. Defaults to SKILL_DIR/gtn_cache.json.

    Returns:
        The cache dict.
    """
    output_path = output_path or (SKILL_DIR / "gtn_cache.json")

    # Fetch all data (may raise on network error — old cache stays intact)
    print("Fetching GTN topics...")
    raw_topics = gtn_client.fetch_topics()

    print(f"Fetching detail for {len(raw_topics)} topics...")
    topics = []
    for t in raw_topics:
        topic_id = t.get("name", "")
        if not topic_id:
            continue
        try:
            detail = gtn_client.fetch_topic_detail(topic_id)
        except Exception as e:
            print(f"  WARN: skipping {topic_id}: {e}")
            continue

        tutorials = []
        for m in detail.get("materials", []):
            if m.get("type") != "tutorial":
                continue
            tutorials.append({
                "title": m.get("title", ""),
                "name": m.get("name", ""),
                "topic": topic_id,
                "time_estimation": m.get("time_estimation", ""),
                "level": m.get("level", ""),
                "objectives": m.get("objectives", []),
                "tools": m.get("tools", []),
                "url": f"{GTN_TUTORIAL_BASE}/topics/{topic_id}/tutorials/{m.get('name', '')}/tutorial.html",
            })

        topics.append({
            "name": topic_id,
            "title": detail.get("title", t.get("title", "")),
            "summary": detail.get("summary", t.get("summary", "")),
            "tutorials": tutorials,
        })

    print("Fetching tool→tutorial index...")
    tool_map = gtn_client.fetch_tool_tutorial_map()

    # Flatten tool map to simpler structure
    tool_index = {}
    for tool_key, data in tool_map.items():
        tool_index[tool_key] = [
            {
                "path": tut[0],
                "title": tut[1],
                "topic": tut[2],
                "url": f"{GTN_TUTORIAL_BASE}{tut[3]}",
            }
            for tut in data.get("tutorials", [])
        ]

    cache = {
        "built_at": datetime.now(timezone.utc).isoformat(),
        "gtn_base": GTN_TUTORIAL_BASE,
        "topic_count": len(topics),
        "tutorial_count": sum(len(t["tutorials"]) for t in topics),
        "tool_count": len(tool_index),
        "topics": topics,
        "tool_index": tool_index,
    }

    # Atomic write: write to temp, rename on success
    tmp_fd = tempfile.NamedTemporaryFile(
        mode="w", dir=output_path.parent, suffix=".tmp", delete=False
    )
    try:
        json.dump(cache, tmp_fd, indent=2, ensure_ascii=False)
        tmp_fd.close()
        Path(tmp_fd.name).replace(output_path)
    except Exception:
        Path(tmp_fd.name).unlink(missing_ok=True)
        raise

    total_tutorials = cache["tutorial_count"]
    print(f"Cache built: {len(topics)} topics, {total_tutorials} tutorials, "
          f"{len(tool_index)} tools → {output_path}")
    return cache


def build_skill_recommendations(
    cache_path: Path | None = None,
    output_path: Path | None = None,
) -> dict:
    """Build skill→tutorial recommendations from the cache + hand-seeded concepts.

    Args:
        cache_path: Path to gtn_cache.json. Defaults to SKILL_DIR/gtn_cache.json.
        output_path: Where to write. Defaults to SKILL_DIR/skill_recommendations.json.

    Returns:
        The recommendations dict.
    """
    cache_path = cache_path or (SKILL_DIR / "gtn_cache.json")
    output_path = output_path or (SKILL_DIR / "skill_recommendations.json")

    cache = json.loads(cache_path.read_text())

    # Build a flat list of all tutorials for scoring
    all_tutorials = []
    for topic in cache["topics"]:
        for tut in topic["tutorials"]:
            all_tutorials.append(tut)

    recommendations = {}
    for skill_name, seed in SKILL_CONCEPTS.items():
        concepts = seed["concepts"]
        preferred_topics = seed["gtn_topics"]

        # Score each tutorial against this skill's concepts
        scored = []
        for tut in all_tutorials:
            score = 0.0
            searchable = " ".join([
                tut.get("title", ""),
                " ".join(tut.get("objectives", [])),
                tut.get("topic", ""),
            ]).lower()

            # Concept matches
            for concept in concepts:
                if concept.lower() in searchable:
                    score += 5.0

            # Topic preference
            if tut.get("topic", "") in preferred_topics:
                score += 3.0

            # Objectives match (bonus for specific terms)
            for obj in tut.get("objectives", []):
                for concept in concepts:
                    if concept.lower() in obj.lower():
                        score += 2.0

            if score > 0:
                scored.append((score, tut))

        # Sort by score descending, take top 5
        scored.sort(key=lambda x: x[0], reverse=True)
        top_tutorials = []
        for score, tut in scored[:5]:
            top_tutorials.append({
                "id": f"{tut['topic']}/{tut['name']}",
                "title": tut["title"],
                "topic": tut["topic"],
                "time": tut.get("time_estimation", ""),
                "level": tut.get("level", ""),
                "url": tut.get("url", ""),
                "relevance_score": round(score, 1),
            })

        recommendations[skill_name] = {
            "concepts": concepts,
            "gtn_topics": preferred_topics,
            "tutorials": top_tutorials,
        }

    # Atomic write (same pattern as build_gtn_cache)
    tmp_fd = tempfile.NamedTemporaryFile(
        mode="w", dir=output_path.parent, suffix=".tmp", delete=False
    )
    try:
        json.dump(recommendations, tmp_fd, indent=2, ensure_ascii=False)
        tmp_fd.close()
        Path(tmp_fd.name).replace(output_path)
    except Exception:
        Path(tmp_fd.name).unlink(missing_ok=True)
        raise

    print(f"Skill recommendations built: {len(recommendations)} skills → {output_path}")
    return recommendations


def main():
    parser = argparse.ArgumentParser(description="Build GTN knowledge cache")
    parser.add_argument("--refresh", action="store_true",
                        help="Re-fetch from GTN API and rebuild cache files")
    parser.add_argument("--recommendations-only", action="store_true",
                        help="Rebuild skill_recommendations.json from existing cache")
    args = parser.parse_args()

    if args.recommendations_only:
        build_skill_recommendations()
        return

    if args.refresh or not (SKILL_DIR / "gtn_cache.json").exists():
        build_gtn_cache()
        build_skill_recommendations()
    else:
        print("Cache already exists. Use --refresh to rebuild.")


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `python -m pytest skills/knowledge-guide/tests/test_cache_builder.py -v`
Expected: All 3 tests PASS

- [ ] **Step 5: Commit**

```bash
git add skills/knowledge-guide/gtn_cache_builder.py skills/knowledge-guide/tests/test_cache_builder.py
git commit -m "feat(knowledge-guide): add GTN cache builder with atomic writes"
```

---

### Task 4: Build the initial cache (live GTN fetch)

**Files:**
- Create: `skills/knowledge-guide/gtn_cache.json` (auto-generated)
- Create: `skills/knowledge-guide/skill_recommendations.json` (auto-generated)

- [ ] **Step 1: Run the cache builder against live GTN API**

Run: `python skills/knowledge-guide/gtn_cache_builder.py --refresh`
Expected: Prints topic/tutorial/tool counts, writes both JSON files

- [ ] **Step 2: Verify cache structure**

Run: `python -c "import json; d=json.load(open('skills/knowledge-guide/gtn_cache.json')); print(f'Topics: {d[\"topic_count\"]}, Tutorials: {d[\"tutorial_count\"]}, Tools: {d[\"tool_count\"]}')"`
Expected: Topics: ~43, Tutorials: ~400-500, Tools: ~200+

- [ ] **Step 3: Verify recommendations structure**

Run: `python -c "import json; d=json.load(open('skills/knowledge-guide/skill_recommendations.json')); print(f'Skills: {len(d)}'); [print(f'  {k}: {len(v[\"tutorials\"])} tutorials') for k,v in d.items()]"`
Expected: 11 skills, each with 0-5 tutorials

- [ ] **Step 4: Commit cache files**

```bash
git add skills/knowledge-guide/gtn_cache.json skills/knowledge-guide/skill_recommendations.json
git commit -m "data(knowledge-guide): add initial GTN cache and skill recommendations"
```

---

## Chunk 2: Query Engine, CLI, and Demo

### Task 5: Build query engine (`query_engine.py`)

**Files:**
- Create: `skills/knowledge-guide/query_engine.py`
- Test: `skills/knowledge-guide/tests/test_query_engine.py`

- [ ] **Step 1: Write failing tests for query engine**

```python
# skills/knowledge-guide/tests/test_query_engine.py
"""Tests for the knowledge-guide query engine."""
from __future__ import annotations
import json
import sys
from pathlib import Path

import pytest

SKILL_DIR = Path(__file__).resolve().parent.parent
if str(SKILL_DIR) not in sys.path:
    sys.path.insert(0, str(SKILL_DIR))

# Minimal cache for testing
MINI_CACHE = {
    "topics": [
        {
            "name": "variant-analysis",
            "title": "Variant Analysis",
            "summary": "Genetic variant detection",
            "tutorials": [
                {
                    "title": "Calling very rare variants",
                    "name": "rare-variants",
                    "topic": "variant-analysis",
                    "time_estimation": "3H",
                    "level": "Intermediate",
                    "objectives": ["Process duplex sequencing data",
                                   "Identify rare variants without diploid assumptions"],
                    "tools": ["fastqc", "bwa_mem"],
                    "url": "https://training.galaxyproject.org/.../rare-variants/tutorial.html",
                },
                {
                    "title": "Exome sequencing variant calling",
                    "name": "exome-seq",
                    "topic": "variant-analysis",
                    "time_estimation": "4H",
                    "level": "Intermediate",
                    "objectives": ["Perform variant calling on exome data",
                                   "Annotate variants with functional information"],
                    "tools": ["bwa_mem", "freebayes", "snpeff"],
                    "url": "https://training.galaxyproject.org/.../exome-seq/tutorial.html",
                },
            ],
        },
        {
            "name": "transcriptomics",
            "title": "Transcriptomics",
            "summary": "Gene expression analysis",
            "tutorials": [
                {
                    "title": "Reference-based RNA-Seq data analysis",
                    "name": "ref-based",
                    "topic": "transcriptomics",
                    "time_estimation": "4H",
                    "level": "Intermediate",
                    "objectives": ["Perform RNA-seq alignment",
                                   "Identify differentially expressed genes"],
                    "tools": ["hisat2", "featurecounts", "deseq2"],
                    "url": "https://training.galaxyproject.org/.../ref-based/tutorial.html",
                },
            ],
        },
    ],
    "tool_index": {
        "devteam/fastqc/fastqc": [
            {"path": "variant-analysis/rare-variants",
             "title": "Calling very rare variants",
             "topic": "Variant Analysis",
             "url": "https://training.galaxyproject.org/.../rare-variants/tutorial.html"},
        ],
    },
}


class TestFreeTextQuery:
    """Test free-text query matching."""

    def test_matches_by_title(self):
        from query_engine import search_tutorials
        results = search_tutorials("variant calling", MINI_CACHE)
        assert len(results) > 0
        assert any("variant" in r["title"].lower() for r in results)

    def test_matches_by_objective(self):
        from query_engine import search_tutorials
        results = search_tutorials("differentially expressed genes", MINI_CACHE)
        assert len(results) > 0
        assert results[0]["topic"] == "transcriptomics"

    def test_returns_max_5(self):
        from query_engine import search_tutorials
        results = search_tutorials("analysis", MINI_CACHE)
        assert len(results) <= 5

    def test_empty_query_returns_empty(self):
        from query_engine import search_tutorials
        results = search_tutorials("", MINI_CACHE)
        assert results == []

    def test_no_matches_returns_empty(self):
        from query_engine import search_tutorials
        results = search_tutorials("quantum computing blockchain", MINI_CACHE)
        assert results == []


class TestTopicLookup:
    """Test direct topic lookup."""

    def test_finds_topic(self):
        from query_engine import lookup_topic
        result = lookup_topic("variant-analysis", MINI_CACHE)
        assert result is not None
        assert result["name"] == "variant-analysis"
        assert len(result["tutorials"]) == 2

    def test_unknown_topic_returns_none(self):
        from query_engine import lookup_topic
        result = lookup_topic("nonexistent-topic", MINI_CACHE)
        assert result is None


class TestToolLookup:
    """Test tool→tutorial reverse lookup."""

    def test_finds_tool(self):
        from query_engine import lookup_tool
        results = lookup_tool("fastqc", MINI_CACHE)
        assert len(results) > 0

    def test_unknown_tool_returns_empty(self):
        from query_engine import lookup_tool
        results = lookup_tool("nonexistent_tool_xyz", MINI_CACHE)
        assert results == []


class TestConceptLookup:
    """Test concept→tutorial fuzzy matching."""

    def test_matches_concept(self):
        from query_engine import lookup_concept

        recs = {
            "pharmgx-reporter": {
                "concepts": ["pharmacogenomics", "CYP2D6"],
                "tutorials": [
                    {"id": "variant-analysis/pharmgx", "title": "Pharmacogenomics"}
                ],
            }
        }
        results = lookup_concept("pharmacogenomics", MINI_CACHE, recs)
        assert len(results) > 0
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `python -m pytest skills/knowledge-guide/tests/test_query_engine.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'query_engine'`

- [ ] **Step 3: Implement query engine**

```python
# skills/knowledge-guide/query_engine.py
"""Query engine for the knowledge-guide skill.

Resolves free-text queries, topic lookups, tool lookups, and concept
lookups against the GTN cache. Scoring follows the same weighted-keyword
pattern used in galaxy_bridge/tool_recommender.py.
"""
from __future__ import annotations

import re

# Stopwords to strip from free-text queries
_STOPWORDS = frozenset({
    "a", "an", "the", "is", "are", "was", "were", "be", "been", "being",
    "have", "has", "had", "do", "does", "did", "will", "would", "could",
    "should", "may", "might", "shall", "can", "need", "must",
    "in", "on", "at", "to", "for", "of", "with", "by", "from", "as",
    "into", "through", "during", "before", "after", "above", "below",
    "between", "out", "off", "over", "under", "again", "further", "then",
    "and", "but", "or", "nor", "not", "so", "yet", "both", "either",
    "neither", "each", "every", "all", "any", "few", "more", "most",
    "other", "some", "such", "no", "only", "own", "same",
    "what", "which", "who", "whom", "this", "that", "these", "those",
    "how", "why", "when", "where",
    "i", "me", "my", "we", "our", "you", "your", "he", "she", "it",
    "they", "them", "its", "his", "her", "their",
    "about", "does", "matter",
})


def _tokenize(text: str) -> list[str]:
    """Lowercase, strip punctuation, remove stopwords."""
    words = re.findall(r"[a-z0-9]+(?:-[a-z0-9]+)*", text.lower())
    return [w for w in words if w not in _STOPWORDS]


def _flatten_tutorials(cache: dict) -> list[dict]:
    """Extract all tutorials from the cache into a flat list."""
    tutorials = []
    for topic in cache.get("topics", []):
        for tut in topic.get("tutorials", []):
            tutorials.append(tut)
    return tutorials


def search_tutorials(
    query: str,
    cache: dict,
    max_results: int = 5,
) -> list[dict]:
    """Free-text search across all cached tutorials.

    Scores each tutorial by keyword matches in title (weight 4),
    objectives (weight 3), topic title (weight 2), and tool names (weight 1).

    Returns top-N results sorted by score descending.
    """
    tokens = _tokenize(query)
    if not tokens:
        return []

    tutorials = _flatten_tutorials(cache)
    scored: list[tuple[float, dict]] = []

    for tut in tutorials:
        score = 0.0
        title_lower = tut.get("title", "").lower()
        objectives_text = " ".join(tut.get("objectives", [])).lower()
        topic_lower = tut.get("topic", "").lower()
        tools_text = " ".join(tut.get("tools", [])).lower()

        # Full phrase match (highest signal)
        query_lower = query.lower().strip()
        if query_lower in title_lower:
            score += 10.0
        if query_lower in objectives_text:
            score += 6.0

        # Per-token scoring
        for token in tokens:
            if token in title_lower:
                score += 4.0
            if token in objectives_text:
                score += 3.0
            if token in topic_lower:
                score += 2.0
            if token in tools_text:
                score += 1.0

        if score > 0:
            scored.append((score, tut))

    scored.sort(key=lambda x: x[0], reverse=True)

    results = []
    for score, tut in scored[:max_results]:
        results.append({
            "title": tut["title"],
            "name": tut.get("name", ""),
            "topic": tut.get("topic", ""),
            "time_estimation": tut.get("time_estimation", ""),
            "level": tut.get("level", ""),
            "objectives": tut.get("objectives", []),
            "url": tut.get("url", ""),
            "relevance_score": round(score, 1),
        })
    return results


def lookup_topic(topic_id: str, cache: dict) -> dict | None:
    """Direct topic lookup by ID. Returns topic dict or None."""
    for topic in cache.get("topics", []):
        if topic["name"] == topic_id:
            return topic
    return None


def lookup_tool(tool_name: str, cache: dict) -> list[dict]:
    """Reverse lookup: tool name → tutorials that use it.

    Searches the tool_index keys (partial match) and also scans
    tutorial tool lists for direct name matches.
    """
    tool_lower = tool_name.lower().strip()
    results = []
    seen = set()

    # 1. Check tool_index (keyed by repo path like devteam/fastqc/fastqc)
    for tool_key, tutorials in cache.get("tool_index", {}).items():
        if tool_lower in tool_key.lower():
            for tut in tutorials:
                key = tut.get("path", tut.get("title", ""))
                if key not in seen:
                    seen.add(key)
                    results.append(tut)

    # 2. Scan tutorial tool lists for direct matches
    for topic in cache.get("topics", []):
        for tut in topic.get("tutorials", []):
            for tool in tut.get("tools", []):
                if tool_lower in tool.lower():
                    key = f"{tut.get('topic', '')}/{tut.get('name', '')}"
                    if key not in seen:
                        seen.add(key)
                        results.append({
                            "path": key,
                            "title": tut["title"],
                            "topic": tut.get("topic", ""),
                            "url": tut.get("url", ""),
                            "time_estimation": tut.get("time_estimation", ""),
                            "level": tut.get("level", ""),
                            "objectives": tut.get("objectives", []),
                        })

    return results


def lookup_concept(
    concept: str,
    cache: dict,
    recommendations: dict | None = None,
) -> list[dict]:
    """Fuzzy concept lookup — checks skill_recommendations.json first,
    then falls back to free-text search.

    Returns list of tutorial matches.
    """
    concept_lower = concept.lower().strip()
    results = []
    seen = set()

    # 1. Check recommendations for matching concepts
    if recommendations:
        for skill_name, entry in recommendations.items():
            for c in entry.get("concepts", []):
                if concept_lower in c.lower() or c.lower() in concept_lower:
                    for tut in entry.get("tutorials", []):
                        if tut.get("id", "") not in seen:
                            seen.add(tut["id"])
                            results.append(tut)

    # 2. Fall back to free-text search
    if not results:
        results = search_tutorials(concept, cache)

    return results
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `python -m pytest skills/knowledge-guide/tests/test_query_engine.py -v`
Expected: All 9 tests PASS

- [ ] **Step 5: Commit**

```bash
git add skills/knowledge-guide/query_engine.py skills/knowledge-guide/tests/test_query_engine.py
git commit -m "feat(knowledge-guide): add query engine with weighted scoring"
```

---

### Task 6: Build CLI entry point (`knowledge_guide.py`)

**Files:**
- Create: `skills/knowledge-guide/knowledge_guide.py`
- Test: `skills/knowledge-guide/tests/test_cli.py`

- [ ] **Step 1: Write failing tests for CLI**

```python
# skills/knowledge-guide/tests/test_cli.py
"""Tests for the knowledge-guide CLI."""
from __future__ import annotations
import json
import subprocess
import sys
from pathlib import Path

import pytest

SKILL_DIR = Path(__file__).resolve().parent.parent
SCRIPT = SKILL_DIR / "knowledge_guide.py"


class TestDemoMode:
    """Test --demo flag runs without network and produces output."""

    def test_demo_produces_output(self, tmp_path):
        result = subprocess.run(
            [sys.executable, str(SCRIPT), "--demo", "--output", str(tmp_path)],
            capture_output=True, text=True, timeout=30,
        )
        assert result.returncode == 0, f"STDERR: {result.stderr}"
        assert (tmp_path / "report.md").exists()
        assert (tmp_path / "report.html").exists()
        assert (tmp_path / "result.json").exists()

    def test_demo_result_json_envelope(self, tmp_path):
        subprocess.run(
            [sys.executable, str(SCRIPT), "--demo", "--output", str(tmp_path)],
            capture_output=True, text=True, timeout=30,
        )
        data = json.loads((tmp_path / "result.json").read_text())
        assert data["skill"] == "knowledge-guide"
        assert "summary" in data
        assert "data" in data


class TestQueryMode:
    """Test --query flag with cache."""

    def test_query_produces_output(self, tmp_path):
        result = subprocess.run(
            [sys.executable, str(SCRIPT), "--query", "variant calling",
             "--output", str(tmp_path)],
            capture_output=True, text=True, timeout=30,
        )
        assert result.returncode == 0, f"STDERR: {result.stderr}"
        assert (tmp_path / "report.md").exists()
        assert (tmp_path / "result.json").exists()


class TestSkillMode:
    """Test --skill flag for Learn More lookups."""

    def test_skill_lookup(self, tmp_path):
        result = subprocess.run(
            [sys.executable, str(SCRIPT), "--skill", "pharmgx-reporter",
             "--output", str(tmp_path)],
            capture_output=True, text=True, timeout=30,
        )
        assert result.returncode == 0, f"STDERR: {result.stderr}"
        assert (tmp_path / "result.json").exists()
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `python -m pytest skills/knowledge-guide/tests/test_cli.py -v`
Expected: FAIL — script does not exist or missing implementation

- [ ] **Step 3: Implement CLI**

```python
# skills/knowledge-guide/knowledge_guide.py
#!/usr/bin/env python3
"""
knowledge_guide.py — GTN-backed educational guide for ClawBio
=============================================================
Explains bioinformatics concepts grounded in Galaxy Training Network
tutorials. Supports free-text queries, structured lookups, and
report-embedded "Learn More" sections.

Usage:
    python knowledge_guide.py --query "what is variant calling?" --output /tmp/kg
    python knowledge_guide.py --topic variant-analysis --output /tmp/kg
    python knowledge_guide.py --tool fastqc --output /tmp/kg
    python knowledge_guide.py --concept "polygenic risk" --output /tmp/kg
    python knowledge_guide.py --skill gwas-prs --output /tmp/kg
    python knowledge_guide.py --query "RNA-seq DE" --deep --output /tmp/kg
    python knowledge_guide.py --demo --output /tmp/kg_demo
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

SKILL_DIR = Path(__file__).resolve().parent
_PROJECT_ROOT = SKILL_DIR.parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))
if str(SKILL_DIR) not in sys.path:
    sys.path.insert(0, str(SKILL_DIR))

from clawbio.common.report import generate_report_header, generate_report_footer, write_result_json
from clawbio.common.html_report import HtmlReportBuilder

CACHE_PATH = SKILL_DIR / "gtn_cache.json"
RECS_PATH = SKILL_DIR / "skill_recommendations.json"

# Embedded mini-cache for demo mode (no network needed)
DEMO_CACHE = {
    "topics": [
        {
            "name": "variant-analysis",
            "title": "Variant Analysis",
            "summary": "Detecting genetic variants in sequencing data",
            "tutorials": [
                {
                    "title": "Calling very rare variants",
                    "name": "rare-variants",
                    "topic": "variant-analysis",
                    "time_estimation": "3H",
                    "level": "Intermediate",
                    "objectives": [
                        "Process duplex sequencing data",
                        "Identify rare variants without diploid assumptions",
                    ],
                    "tools": ["fastqc", "bwa_mem", "naive_variant_caller"],
                    "url": "https://training.galaxyproject.org/training-material/topics/variant-analysis/tutorials/rare-variants/tutorial.html",
                },
                {
                    "title": "Exome sequencing data analysis",
                    "name": "exome-seq",
                    "topic": "variant-analysis",
                    "time_estimation": "4H",
                    "level": "Intermediate",
                    "objectives": [
                        "Perform variant calling on whole-exome data",
                        "Annotate variants with functional impact",
                    ],
                    "tools": ["bwa_mem", "freebayes", "snpeff", "gemini"],
                    "url": "https://training.galaxyproject.org/training-material/topics/variant-analysis/tutorials/exome-seq/tutorial.html",
                },
            ],
        },
        {
            "name": "transcriptomics",
            "title": "Transcriptomics",
            "summary": "Gene expression analysis approaches",
            "tutorials": [
                {
                    "title": "Reference-based RNA-Seq data analysis",
                    "name": "ref-based",
                    "topic": "transcriptomics",
                    "time_estimation": "4H",
                    "level": "Intermediate",
                    "objectives": [
                        "Perform RNA-seq alignment with HISAT2",
                        "Identify differentially expressed genes with DESeq2",
                    ],
                    "tools": ["hisat2", "featurecounts", "deseq2"],
                    "url": "https://training.galaxyproject.org/training-material/topics/transcriptomics/tutorials/ref-based/tutorial.html",
                },
            ],
        },
        {
            "name": "single-cell",
            "title": "Single Cell",
            "summary": "Single-cell RNA-seq and single-cell analysis",
            "tutorials": [
                {
                    "title": "Clustering 3K PBMCs with Scanpy",
                    "name": "scrna-scanpy-pbmc3k",
                    "topic": "single-cell",
                    "time_estimation": "3H",
                    "level": "Intermediate",
                    "objectives": [
                        "Perform quality control on scRNA-seq data",
                        "Cluster cells and identify marker genes",
                    ],
                    "tools": ["scanpy_filter_cells", "scanpy_normalise", "scanpy_find_cluster"],
                    "url": "https://training.galaxyproject.org/training-material/topics/single-cell/tutorials/scrna-scanpy-pbmc3k/tutorial.html",
                },
            ],
        },
    ],
    "tool_index": {
        "devteam/fastqc/fastqc": [
            {
                "path": "variant-analysis/rare-variants",
                "title": "Calling very rare variants",
                "topic": "Variant Analysis",
                "url": "https://training.galaxyproject.org/training-material/topics/variant-analysis/tutorials/rare-variants/tutorial.html",
            },
        ],
    },
}


def _load_cache(demo: bool = False) -> dict:
    """Load the GTN cache (or demo cache)."""
    if demo:
        return DEMO_CACHE
    if not CACHE_PATH.exists():
        print(f"ERROR: Cache not found at {CACHE_PATH}", file=sys.stderr)
        print("Run: python gtn_cache_builder.py --refresh", file=sys.stderr)
        sys.exit(1)
    return json.loads(CACHE_PATH.read_text())


def _load_recommendations() -> dict:
    """Load skill recommendations if available."""
    if RECS_PATH.exists():
        return json.loads(RECS_PATH.read_text())
    return {}


def _format_tutorial_md(tut: dict, index: int = 1) -> str:
    """Format a single tutorial as markdown."""
    lines = [f"### {index}. {tut.get('title', 'Untitled')}"]

    meta = []
    if tut.get("topic"):
        meta.append(f"**Topic:** {tut['topic']}")
    if tut.get("time_estimation") or tut.get("time"):
        meta.append(f"**Time:** {tut.get('time_estimation') or tut.get('time')}")
    if tut.get("level"):
        meta.append(f"**Level:** {tut['level']}")
    if tut.get("relevance_score"):
        meta.append(f"**Relevance:** {tut['relevance_score']}")
    if meta:
        lines.append(" | ".join(meta))

    objectives = tut.get("objectives", [])
    if objectives:
        lines.append("")
        lines.append("**Learning objectives:**")
        for obj in objectives:
            lines.append(f"- {obj}")

    url = tut.get("url", "")
    if url:
        lines.append("")
        lines.append(f"[Full tutorial]({url})")

    lines.append("")
    return "\n".join(lines)


def _format_deep_content(content: dict) -> str:
    """Format deep tutorial content as markdown."""
    lines = [f"## Deep Content: {content.get('title', '')}"]
    body = content.get("content", "")
    if not body:
        lines.append("*(No content available for this tutorial)*")
        return "\n".join(lines)

    # Extract first ~2000 chars of content as a preview
    # The GTN tutorial content is HTML — extract text portions
    text = re.sub(r"<[^>]+>", "", body)  # strip HTML tags
    text = re.sub(r"\n{3,}", "\n\n", text)  # collapse whitespace
    preview = text[:2000]
    if len(text) > 2000:
        preview += "\n\n*(Content truncated — see full tutorial link above)*"

    lines.append("")
    lines.append(preview)
    return "\n".join(lines)


def run_query(args: argparse.Namespace) -> dict:
    """Execute a knowledge-guide query and return results."""
    from query_engine import search_tutorials, lookup_topic, lookup_tool, lookup_concept

    cache = _load_cache(demo=args.demo)
    recs = _load_recommendations()

    mode = "unknown"
    query_text = ""
    results = []

    if args.demo:
        mode = "demo"
        query_text = "variant calling"
        results = search_tutorials("variant calling", cache)
    elif args.query:
        mode = "free-text"
        query_text = args.query
        results = search_tutorials(args.query, cache)
    elif args.topic:
        mode = "topic"
        query_text = args.topic
        topic = lookup_topic(args.topic, cache)
        if topic:
            results = [
                {
                    "title": t["title"],
                    "name": t.get("name", ""),
                    "topic": topic["name"],
                    "time_estimation": t.get("time_estimation", ""),
                    "level": t.get("level", ""),
                    "objectives": t.get("objectives", []),
                    "url": t.get("url", ""),
                }
                for t in topic["tutorials"]
            ]
    elif args.tool:
        mode = "tool"
        query_text = args.tool
        results = lookup_tool(args.tool, cache)
    elif args.concept:
        mode = "concept"
        query_text = args.concept
        results = lookup_concept(args.concept, cache, recs)
    elif args.skill:
        mode = "skill"
        query_text = args.skill
        if args.skill in recs:
            entry = recs[args.skill]
            results = entry.get("tutorials", [])

    # Deep content fetch (live, standalone queries only)
    deep_content = None
    if getattr(args, "deep", False) and results and not args.demo:
        top = results[0]
        topic_id = top.get("topic", "")
        tut_name = top.get("name", top.get("path", "").split("/")[-1] if "/" in top.get("path", "") else "")
        if topic_id and tut_name:
            from gtn_client import fetch_tutorial_content
            deep_content = fetch_tutorial_content(topic_id, tut_name)

    return {
        "mode": mode,
        "query": query_text,
        "matches": len(results),
        "tutorials": results,
        "deep_content": deep_content,
        "deep_content_fetched": deep_content is not None,
    }


def generate_report(result: dict, output_dir: Path) -> None:
    """Generate report.md, report.html, and result.json."""
    output_dir.mkdir(parents=True, exist_ok=True)

    # --- report.md ---
    header = generate_report_header(
        title="Knowledge Guide",
        skill_name="knowledge-guide",
        extra_metadata={
            "Query": result["query"],
            "Mode": result["mode"],
            "Matches": str(result["matches"]),
        },
    )
    body_lines = []
    if not result["tutorials"]:
        body_lines.append("No matching tutorials found for this query.\n")
    else:
        body_lines.append("## Matching Tutorials\n")
        for i, tut in enumerate(result["tutorials"], 1):
            body_lines.append(_format_tutorial_md(tut, i))

    if result.get("deep_content"):
        body_lines.append(_format_deep_content(result["deep_content"]))

    footer = generate_report_footer()
    report_md = header + "\n".join(body_lines) + footer
    (output_dir / "report.md").write_text(report_md)

    # --- report.html ---
    builder = HtmlReportBuilder(title="Knowledge Guide", skill="knowledge-guide")
    builder.add_header_block(
        title="Knowledge Guide",
        subtitle=f"Query: {result['query']} ({result['mode']} mode)",
    )
    if result["tutorials"]:
        rows = []
        for tut in result["tutorials"]:
            rows.append([
                tut.get("title", ""),
                tut.get("topic", ""),
                tut.get("time_estimation", tut.get("time", "")),
                tut.get("level", ""),
                f'<a href="{tut.get("url", "#")}">Open</a>',
            ])
        builder.add_section(
            title="Matching Tutorials",
            content=builder.table(
                headers=["Tutorial", "Topic", "Time", "Level", "Link"],
                rows=rows,
            ),
        )
    else:
        builder.add_section(title="Results", content="No matching tutorials found.")
    html_path = output_dir / "report.html"
    html_path.write_text(builder.render())

    # --- result.json ---
    # Strip deep_content from result.json (too large)
    data = {k: v for k, v in result.items() if k != "deep_content"}
    write_result_json(
        output_dir=str(output_dir),
        skill="knowledge-guide",
        version="0.1.0",
        summary={
            "query": result["query"],
            "mode": result["mode"],
            "matches": result["matches"],
            "deep_content_fetched": result["deep_content_fetched"],
        },
        data=data,
    )

    # --- Deep content file (if fetched) ---
    if result.get("deep_content"):
        deep_dir = output_dir / "tutorials"
        deep_dir.mkdir(exist_ok=True)
        content = result["deep_content"]
        filename = f"{content.get('name', 'tutorial')}.md"
        (deep_dir / filename).write_text(_format_deep_content(content))

    print(f"Report written to {output_dir}/")


def main():
    parser = argparse.ArgumentParser(
        description="Knowledge Guide — GTN-backed educational guide for ClawBio"
    )
    parser.add_argument("--query", type=str, help="Free-text query")
    parser.add_argument("--topic", type=str, help="GTN topic ID (e.g. variant-analysis)")
    parser.add_argument("--tool", type=str, help="Galaxy tool name (e.g. fastqc)")
    parser.add_argument("--concept", type=str, help="Concept term (e.g. 'polygenic risk')")
    parser.add_argument("--skill", type=str, help="ClawBio skill name for Learn More lookup")
    parser.add_argument("--deep", action="store_true", help="Fetch full tutorial content (live)")
    parser.add_argument("--demo", action="store_true", help="Run demo with embedded data")
    parser.add_argument("--output", type=str, required=True, help="Output directory")

    args = parser.parse_args()

    # Require at least one query mode
    if not any([args.query, args.topic, args.tool, args.concept, args.skill, args.demo]):
        parser.error("Provide one of: --query, --topic, --tool, --concept, --skill, or --demo")

    result = run_query(args)
    generate_report(result, Path(args.output))


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `python -m pytest skills/knowledge-guide/tests/test_cli.py -v`
Expected: All 4 tests PASS

- [ ] **Step 5: Commit**

```bash
git add skills/knowledge-guide/knowledge_guide.py skills/knowledge-guide/tests/test_cli.py
git commit -m "feat(knowledge-guide): add CLI with query, topic, tool, concept, skill, and demo modes"
```

---

### Task 7: Create demo data and SKILL.md

**Files:**
- Create: `skills/knowledge-guide/demo/demo_queries.json`
- Create: `skills/knowledge-guide/SKILL.md`

- [ ] **Step 1: Create demo queries file**

```json
[
  {"query": "variant calling", "expected_topic": "variant-analysis"},
  {"query": "RNA-seq differential expression", "expected_topic": "transcriptomics"},
  {"query": "single cell clustering", "expected_topic": "single-cell"},
  {"tool": "fastqc", "expected_tutorials_min": 1},
  {"skill": "pharmgx-reporter", "expected_concepts": ["pharmacogenomics"]}
]
```

- [ ] **Step 2: Create SKILL.md**

```markdown
---
name: knowledge-guide
description: >-
  GTN-backed educational guide — explains bioinformatics concepts grounded
  in Galaxy Training Network tutorials. Supports free-text queries, structured
  lookups by topic/tool/concept, and report-embedded "Learn More" sections.
version: 0.1.0
metadata:
  openclaw:
    category: education
    data_source: Galaxy Training Network (training.galaxyproject.org)
---

# Knowledge Guide

Explains bioinformatics concepts by grounding every answer in Galaxy Training
Network (GTN) tutorials. No hallucinated biology — every explanation traces
back to an authoritative training resource.

## Core Capabilities

1. **Free-text search** — "what is variant calling?" → ranked tutorial matches
2. **Topic lookup** — browse all tutorials within a GTN topic
3. **Tool lookup** — find tutorials that use a specific Galaxy tool
4. **Concept lookup** — fuzzy match against skill-seeded concept lists
5. **Skill Learn More** — pre-computed tutorial recommendations per ClawBio skill
6. **Deep content pull** — live-fetch full tutorial content for inline explanation

## Input Modes

| Flag | Example | Description |
|------|---------|-------------|
| `--query` | `"what is variant calling?"` | Free-text keyword search |
| `--topic` | `variant-analysis` | Direct GTN topic lookup |
| `--tool` | `fastqc` | Reverse tool→tutorial index |
| `--concept` | `"polygenic risk"` | Fuzzy concept matching |
| `--skill` | `gwas-prs` | ClawBio skill Learn More |
| `--deep` | (flag) | Fetch full tutorial content live |

## Scoring Weights (Free-Text Mode)

| Signal | Weight | Bonus |
|--------|--------|-------|
| Title match | 4.0/token | +10.0 phrase |
| Objectives match | 3.0/token | +6.0 phrase |
| Topic match | 2.0/token | — |
| Tool match | 1.0/token | — |

## Output Structure

```text
output_dir/
├── report.md
├── report.html
├── result.json
└── tutorials/        (only with --deep)
    └── <name>.md
```

## Dependencies

- `requests` (for live GTN API calls and --deep mode)

## Disclaimer

*ClawBio is a research and educational tool. It is not a medical device
and does not provide clinical diagnoses. Consult a healthcare professional
before making any medical decisions.*
```

- [ ] **Step 3: Verify demo mode end-to-end**

Run: `python skills/knowledge-guide/knowledge_guide.py --demo --output /tmp/kg_demo`
Expected: Creates `/tmp/kg_demo/report.md` and `/tmp/kg_demo/result.json` with 2+ tutorial matches

- [ ] **Step 4: Commit**

```bash
git add skills/knowledge-guide/demo/ skills/knowledge-guide/SKILL.md
git commit -m "docs(knowledge-guide): add SKILL.md and demo queries"
```

---

## Chunk 3: API Integration and Registration

> **Prerequisite:** Chunks 1 and 2 must be complete before starting this chunk. Task 8's `api.py` imports from `knowledge_guide.py` (Chunk 2) and reads cache files (Chunk 1).

### Task 8: Build importable API (`api.py`)

**Files:**
- Create: `skills/knowledge-guide/api.py`
- Test: `skills/knowledge-guide/tests/test_api.py`

- [ ] **Step 1: Write failing tests for api.py**

```python
# skills/knowledge-guide/tests/test_api.py
"""Tests for the knowledge-guide importable API."""
from __future__ import annotations
import json
import sys
from pathlib import Path

import pytest

SKILL_DIR = Path(__file__).resolve().parent.parent
if str(SKILL_DIR) not in sys.path:
    sys.path.insert(0, str(SKILL_DIR))


class TestRun:
    """Test the run() entry point."""

    def test_run_with_query(self):
        from api import run
        result = run(options={"query": "variant calling", "demo": True})
        assert "summary" in result
        assert result["summary"]["matches"] > 0

    def test_run_with_output_dir(self, tmp_path):
        from api import run
        result = run(options={
            "query": "variant calling",
            "demo": True,
            "output_dir": str(tmp_path),
        })
        assert (tmp_path / "result.json").exists()


class TestGetLearnMore:
    """Test the get_learn_more() function for report integration."""

    def test_known_skill(self, tmp_path):
        """Should return tutorials for a registered skill."""
        from api import get_learn_more

        # Create a minimal recommendations file
        recs = {
            "pharmgx-reporter": {
                "concepts": ["pharmacogenomics", "CYP2D6"],
                "gtn_topics": ["variant-analysis"],
                "tutorials": [
                    {
                        "id": "variant-analysis/pharmgx",
                        "title": "Pharmacogenomics tutorial",
                        "topic": "variant-analysis",
                        "time": "2h",
                        "url": "https://example.com/tutorial",
                        "relevance_score": 10.0,
                    }
                ],
            }
        }
        recs_path = tmp_path / "skill_recommendations.json"
        recs_path.write_text(json.dumps(recs))

        result = get_learn_more("pharmgx-reporter", recommendations_path=recs_path)
        assert result["section_title"] == "Learn More"
        assert len(result["tutorials"]) == 1
        assert "html" in result

    def test_unknown_skill_returns_empty(self, tmp_path):
        """Should return empty section for unregistered skill."""
        from api import get_learn_more

        recs_path = tmp_path / "skill_recommendations.json"
        recs_path.write_text("{}")

        result = get_learn_more("nonexistent-skill", recommendations_path=recs_path)
        assert result["tutorials"] == []
        assert result["html"] == ""
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `python -m pytest skills/knowledge-guide/tests/test_api.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'api'`

- [ ] **Step 3: Implement api.py**

```python
# skills/knowledge-guide/api.py
"""Importable API for the knowledge-guide skill."""
from __future__ import annotations
import json
import sys
from pathlib import Path

_SKILL_DIR = Path(__file__).resolve().parent
_PROJECT_ROOT = _SKILL_DIR.parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))
if str(_SKILL_DIR) not in sys.path:
    sys.path.insert(0, str(_SKILL_DIR))

from clawbio.common.report import write_result_json

RECS_PATH = _SKILL_DIR / "skill_recommendations.json"


def run(genotypes: dict[str, str] | None = None, options: dict | None = None) -> dict:
    """Run a knowledge-guide query.

    Args:
        genotypes: Not used (present for API consistency).
        options: Dict with keys:
            - query: str — free-text query
            - topic: str — GTN topic ID
            - tool: str — Galaxy tool name
            - concept: str — concept term
            - skill: str — ClawBio skill name
            - deep: bool — fetch full tutorial content
            - demo: bool — use embedded demo data
            - output_dir: str — if provided, write report files

    Returns:
        Result dict with summary and data.
    """
    import argparse
    from knowledge_guide import run_query, generate_report

    options = options or {}

    # Build a Namespace to match CLI expectations
    args = argparse.Namespace(
        query=options.get("query"),
        topic=options.get("topic"),
        tool=options.get("tool"),
        concept=options.get("concept"),
        skill=options.get("skill"),
        deep=options.get("deep", False),
        demo=options.get("demo", False),
    )

    result = run_query(args)

    output_dir = options.get("output_dir")
    if output_dir:
        generate_report(result, Path(output_dir))

    return {
        "summary": {
            "query": result["query"],
            "mode": result["mode"],
            "matches": result["matches"],
        },
        "data": {
            "tutorials": result["tutorials"],
            "deep_content_fetched": result["deep_content_fetched"],
        },
    }


def get_learn_more(
    skill_name: str,
    recommendations_path: Path | None = None,
) -> dict:
    """Get a pre-built 'Learn More' section for a ClawBio skill.

    Reads from skill_recommendations.json (zero network calls).
    Returns structured data + pre-rendered HTML block.

    Args:
        skill_name: ClawBio skill name (e.g. 'gwas-prs').
        recommendations_path: Override path to recommendations JSON.

    Returns:
        Dict with:
            section_title: str
            concepts: list[str]
            tutorials: list[dict]
            html: str (pre-rendered HTML block, empty if no tutorials)
    """
    recs_path = recommendations_path or RECS_PATH

    if recs_path.exists():
        recs = json.loads(recs_path.read_text())
    else:
        recs = {}

    entry = recs.get(skill_name, {})
    concepts = entry.get("concepts", [])
    tutorials = entry.get("tutorials", [])

    # Pre-render HTML block
    html = ""
    if tutorials:
        items = []
        for tut in tutorials:
            url = tut.get("url", "#")
            title = tut.get("title", "Tutorial")
            time_est = tut.get("time", tut.get("time_estimation", ""))
            level = tut.get("level", "")
            meta = " | ".join(filter(None, [time_est, level]))
            items.append(
                f'<li><a href="{url}" target="_blank">{title}</a>'
                f'{f" <small>({meta})</small>" if meta else ""}</li>'
            )
        html = (
            '<details class="learn-more" open>'
            "<summary><strong>Learn More</strong></summary>"
            f'<ul>{"".join(items)}</ul>'
            "</details>"
        )

    return {
        "section_title": "Learn More",
        "concepts": concepts,
        "tutorials": tutorials,
        "html": html,
    }
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `python -m pytest skills/knowledge-guide/tests/test_api.py -v`
Expected: All 4 tests PASS

- [ ] **Step 5: Commit**

```bash
git add skills/knowledge-guide/api.py skills/knowledge-guide/tests/test_api.py
git commit -m "feat(knowledge-guide): add importable API with get_learn_more()"
```

---

### Task 9: Register skill in catalog.json and CLAUDE.md

**Files:**
- Modify: `skills/catalog.json`
- Modify: `CLAUDE.md`

- [ ] **Step 1: Add knowledge-guide to catalog.json**

Add this entry to the `skills` array in `skills/catalog.json` (alphabetical position, after `gwas-prs`):

```json
{
  "name": "knowledge-guide",
  "cli_alias": "guide",
  "description": "GTN-backed educational guide — explains bioinformatics concepts grounded in Galaxy Training Network tutorials, with free-text queries and report-embedded Learn More sections",
  "version": "0.1.0",
  "status": "mvp",
  "has_script": true,
  "has_tests": true,
  "has_demo": true,
  "demo_command": "python skills/knowledge-guide/knowledge_guide.py --demo --output /tmp/kg_demo",
  "dependencies": [],
  "tags": [
    "education",
    "tutorials",
    "GTN",
    "Galaxy-Training-Network",
    "knowledge-base"
  ],
  "trigger_keywords": [
    "what is",
    "how does",
    "explain",
    "learn more",
    "tutorial",
    "guide",
    "GTN"
  ],
  "chaining_partners": [
    "galaxy-bridge",
    "pharmgx-reporter",
    "gwas-prs",
    "scrna-orchestrator"
  ]
}
```

- [ ] **Step 2: Add routing entry to CLAUDE.md**

Add this row to the Skill Routing Table in `CLAUDE.md`:

```
| Educational guide, "what is X", "how does Y work", explain analysis, learn more, GTN tutorials | `skills/knowledge-guide/` | Run `knowledge_guide.py` |
```

Add to the CLI Reference section:

```bash
# Knowledge Guide — GTN-backed educational explainer
python skills/knowledge-guide/knowledge_guide.py \
  --query "what is variant calling?" --output <report_dir>
python skills/knowledge-guide/knowledge_guide.py \
  --topic variant-analysis --output <report_dir>
python skills/knowledge-guide/knowledge_guide.py \
  --tool fastqc --output <report_dir>
python skills/knowledge-guide/knowledge_guide.py \
  --concept "polygenic risk" --output <report_dir>
python skills/knowledge-guide/knowledge_guide.py \
  --skill gwas-prs --output <report_dir>
python skills/knowledge-guide/knowledge_guide.py \
  --query "RNA-seq DE" --deep --output <report_dir>
python skills/knowledge-guide/knowledge_guide.py --demo --output /tmp/kg_demo
```

Add to the Demo Data table:

```
| Knowledge Guide demo (3 topics, embedded) | `--demo` flag | knowledge-guide |
```

Add to the Demo Commands section:

```bash
# Knowledge Guide demo
python skills/knowledge-guide/knowledge_guide.py --demo --output /tmp/kg_demo
```

- [ ] **Step 3: Commit**

```bash
git add skills/catalog.json CLAUDE.md
git commit -m "docs: register knowledge-guide in catalog.json and CLAUDE.md routing table"
```

---

### Task 10: Run full test suite and verify

**Files:** (none — verification only)

- [ ] **Step 1: Run all knowledge-guide tests**

Run: `python -m pytest skills/knowledge-guide/tests/ -v`
Expected: All tests PASS (5 gtn_client + 3 cache_builder + 9 query_engine + 4 cli + 4 api = ~25 tests)

- [ ] **Step 2: Run demo end-to-end**

Run: `python skills/knowledge-guide/knowledge_guide.py --demo --output /tmp/kg_final_demo`
Expected: Creates report.md, result.json with tutorial matches

- [ ] **Step 3: Test query mode against live cache**

Run: `python skills/knowledge-guide/knowledge_guide.py --query "what is variant calling?" --output /tmp/kg_query_test`
Expected: Creates report with relevant variant analysis tutorials

- [ ] **Step 4: Test skill Learn More mode**

Run: `python skills/knowledge-guide/knowledge_guide.py --skill pharmgx-reporter --output /tmp/kg_skill_test`
Expected: Creates report with pharmacogenomics-related tutorials

- [ ] **Step 5: Final commit**

```bash
git add skills/knowledge-guide/SKILL.md skills/knowledge-guide/*.py \
      skills/knowledge-guide/tests/*.py skills/knowledge-guide/demo/ \
      skills/knowledge-guide/gtn_cache.json skills/knowledge-guide/skill_recommendations.json
git commit -m "feat(knowledge-guide): complete MVP — GTN-backed educational guide skill"
```
