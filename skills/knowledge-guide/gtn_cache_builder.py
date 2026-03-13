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

    # Atomic write
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
