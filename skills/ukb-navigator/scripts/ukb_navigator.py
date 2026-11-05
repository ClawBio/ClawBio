#!/usr/bin/env python3
"""
ukb_navigator.py — UK Biobank Schema Navigator (ClawBio Skill)
==============================================================
Semantic search across UK Biobank's 12,000+ data fields.
Embeds UKB schema into ChromaDB, then provides natural-language field discovery.

Usage:
    python ukb_navigator.py --query "blood pressure" --output /tmp/ukb_report
    python ukb_navigator.py --field 21001 --output /tmp/ukb_report
    python ukb_navigator.py --demo --output /tmp/ukb_demo

Dependencies: chromadb, voyageai (optional)
"""

import argparse
import csv
import json
import os
import sys
from datetime import datetime
from pathlib import Path

SKILL_DIR = Path(__file__).resolve().parent
DEMO_DATA_DIR = SKILL_DIR / "demo_data"
DEFAULT_DB_PATH = SKILL_DIR / "embeddings"
DEFAULT_COLLECTION = "ukb_schema"


# ---------------------------------------------------------------------------
# Embedding helpers
# ---------------------------------------------------------------------------

def _get_embed_fn():
    """Return a Voyage AI embedding function if available, else None (ChromaDB default)."""
    api_key = os.environ.get("VOYAGE_API_KEY")
    if not api_key:
        return None
    try:
        from chromadb.utils.embedding_functions import OpenAIEmbeddingFunction
        # Voyage AI is compatible with the OpenAI embedding interface
        import voyageai  # noqa: F401
        from chromadb.utils.embedding_functions import OpenAIEmbeddingFunction  # noqa: F811

        class VoyageEmbedder:
            def __init__(self):
                self.client = voyageai.Client(api_key=api_key)
                self.model = "voyage-3"

            def __call__(self, input):
                result = self.client.embed(input, model=self.model)
                return result.embeddings

        return VoyageEmbedder()
    except ImportError:
        return None


def _get_collection(db_path: Path, collection_name: str):
    """Get or create ChromaDB collection."""
    import chromadb

    db_path.mkdir(parents=True, exist_ok=True)
    client = chromadb.PersistentClient(path=str(db_path))
    embed_fn = _get_embed_fn()
    kwargs = {"name": collection_name}
    if embed_fn:
        kwargs["embedding_function"] = embed_fn
    return client.get_or_create_collection(**kwargs)


# ---------------------------------------------------------------------------
# Schema embedding
# ---------------------------------------------------------------------------

def embed_schema(schema_csv: Path, schema_txt: Path = None,
                 db_path: Path = DEFAULT_DB_PATH,
                 collection_name: str = DEFAULT_COLLECTION,
                 force: bool = False) -> str:
    """Embed UKB schema CSV (and optional schema_27.txt) into ChromaDB."""
    collection = _get_collection(db_path, collection_name)

    if not force and collection.count() > 0:
        return f"Collection already has {collection.count()} documents. Use --force to re-embed."

    if not schema_csv.exists():
        return f"Schema file not found: {schema_csv}"

    # Clear if forcing
    if force and collection.count() > 0:
        existing = collection.get()
        if existing["ids"]:
            collection.delete(ids=existing["ids"])

    docs, ids, metadatas = [], [], []

    # Parse schema CSV
    with open(schema_csv, "r", encoding="utf-8", errors="replace") as f:
        reader = csv.DictReader(f)
        for i, row in enumerate(reader):
            parts = [f"{k}: {v.strip()}" for k, v in row.items() if v and v.strip()]
            if parts:
                docs.append("\n".join(parts)[:2000])
                ids.append(f"ukb_schema_{i}")
                metadatas.append({"source": "ukb_schema", "row_index": i})

    # Parse schema_27.txt if present
    if schema_txt and schema_txt.exists():
        text = schema_txt.read_text(encoding="utf-8", errors="replace")
        lines = text.splitlines()
        chunk, chunk_chars, chunk_idx = [], 0, 0
        for line in lines:
            chunk.append(line)
            chunk_chars += len(line) + 1
            if chunk_chars >= 1000:
                docs.append("\n".join(chunk))
                ids.append(f"ukb_schema27_{chunk_idx}")
                metadatas.append({"source": "ukb_schema27", "chunk_index": chunk_idx})
                chunk, chunk_chars = [], 0
                chunk_idx += 1
        if chunk:
            docs.append("\n".join(chunk))
            ids.append(f"ukb_schema27_{chunk_idx}")
            metadatas.append({"source": "ukb_schema27", "chunk_index": chunk_idx})

    if not docs:
        return "No documents found to embed."

    # Batch add
    batch_size = 100
    for start in range(0, len(docs), batch_size):
        end = min(start + batch_size, len(docs))
        collection.add(
            documents=docs[start:end],
            ids=ids[start:end],
            metadatas=metadatas[start:end],
        )

    return f"Embedded {collection.count()} documents to '{collection_name}'."


# ---------------------------------------------------------------------------
# Query
# ---------------------------------------------------------------------------

def query_schema(question: str, n_results: int = 10,
                 db_path: Path = DEFAULT_DB_PATH,
                 collection_name: str = DEFAULT_COLLECTION) -> list[dict]:
    """Semantic search against embedded UKB schema. Returns list of matches."""
    collection = _get_collection(db_path, collection_name)

    if collection.count() == 0:
        return []

    results = collection.query(query_texts=[question], n_results=n_results)

    if not results["documents"] or not results["documents"][0]:
        return []

    matches = []
    for doc, meta, dist in zip(
        results["documents"][0],
        results["metadatas"][0],
        results["distances"][0],
    ):
        similarity = max(0, 1 - dist) if dist else 0
        matches.append({
            "text": doc[:500],
            "source": meta.get("source", "unknown"),
            "similarity": round(similarity, 3),
            "metadata": meta,
        })

    return matches


def field_lookup(field_id: str, **kwargs) -> list[dict]:
    """Direct lookup by UK Biobank field ID."""
    return query_schema(f"UK Biobank field {field_id}", **kwargs)


# ---------------------------------------------------------------------------
# Demo data
# ---------------------------------------------------------------------------

DEMO_RESULTS = [
    {
        "text": "Field ID: 4080\nField: Systolic blood pressure, automated reading\nCategory: Blood pressure\nUnits: mmHg\nParticipants: 475,014\nInstances: 4\nDescription: Automated reading of systolic blood pressure measured during assessment visit.",
        "source": "ukb_schema",
        "similarity": 0.912,
        "metadata": {"source": "ukb_schema", "row_index": 4080},
    },
    {
        "text": "Field ID: 4079\nField: Diastolic blood pressure, automated reading\nCategory: Blood pressure\nUnits: mmHg\nParticipants: 475,014\nInstances: 4\nDescription: Automated reading of diastolic blood pressure measured during assessment visit.",
        "source": "ukb_schema",
        "similarity": 0.897,
        "metadata": {"source": "ukb_schema", "row_index": 4079},
    },
    {
        "text": "Field ID: 93\nField: Systolic blood pressure, manual reading\nCategory: Blood pressure\nUnits: mmHg\nParticipants: 475,014\nInstances: 2\nDescription: Manual reading of systolic blood pressure with sphygmomanometer.",
        "source": "ukb_schema",
        "similarity": 0.881,
        "metadata": {"source": "ukb_schema", "row_index": 93},
    },
    {
        "text": "Field ID: 94\nField: Diastolic blood pressure, manual reading\nCategory: Blood pressure\nUnits: mmHg\nParticipants: 475,014\nInstances: 2\nDescription: Manual reading of diastolic blood pressure with sphygmomanometer.",
        "source": "ukb_schema",
        "similarity": 0.869,
        "metadata": {"source": "ukb_schema", "row_index": 94},
    },
    {
        "text": "Field ID: 6150\nField: Vascular/heart problems diagnosed by doctor\nCategory: Medical conditions\nDescription: Self-reported vascular or heart problems diagnosed by a doctor, including hypertension.",
        "source": "ukb_schema",
        "similarity": 0.743,
        "metadata": {"source": "ukb_schema", "row_index": 6150},
    },
    {
        "text": "Field ID: 6153\nField: Medication for cholesterol, blood pressure or diabetes\nCategory: Medications\nDescription: Self-reported medication for cholesterol, blood pressure, diabetes, or hormone replacement therapy.",
        "source": "ukb_schema",
        "similarity": 0.721,
        "metadata": {"source": "ukb_schema", "row_index": 6153},
    },
    {
        "text": "Field ID: 20002\nField: Non-cancer illness code, self-reported\nCategory: Medical conditions\nDescription: Self-reported non-cancer illnesses (including code 1065 = hypertension).",
        "source": "ukb_schema",
        "similarity": 0.698,
        "metadata": {"source": "ukb_schema", "row_index": 20002},
    },
    {
        "text": "Field ID: 41270\nField: Diagnoses - ICD10 (hospital inpatient)\nCategory: Hospital episodes\nDescription: ICD-10 diagnosis codes from Hospital Episode Statistics. Includes I10-I15 for hypertensive diseases.",
        "source": "ukb_schema",
        "similarity": 0.682,
        "metadata": {"source": "ukb_schema", "row_index": 41270},
    },
]


# ---------------------------------------------------------------------------
# Report generation
# ---------------------------------------------------------------------------

def generate_report(query: str, matches: list[dict], output_dir: Path,
                    is_demo: bool = False):
    """Generate markdown report and CSV from search results."""
    output_dir.mkdir(parents=True, exist_ok=True)

    # Markdown report
    lines = [
        f"# UKB Navigator Report",
        f"",
        f"**Query**: {query}",
        f"**Date**: {datetime.now().strftime('%Y-%m-%d %H:%M')}",
        f"**Results**: {len(matches)} matches",
    ]
    if is_demo:
        lines.append("**Mode**: Demo (pre-cached results)")
    lines.append("")
    lines.append("---")
    lines.append("")

    for i, m in enumerate(matches, 1):
        lines.append(f"## [{i}] Similarity: {m['similarity']:.3f} ({m['source']})")
        lines.append("")
        lines.append(m["text"])
        lines.append("")

    lines.append("---")
    lines.append("")
    lines.append("*ClawBio is a research and educational tool. It is not a medical device "
                 "and does not provide clinical diagnoses. Consult a healthcare professional "
                 "before making any medical decisions.*")

    report_path = output_dir / "report.md"
    report_path.write_text("\n".join(lines), encoding="utf-8")

    # CSV
    csv_path = output_dir / "matched_fields.csv"
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["rank", "similarity", "source", "text"])
        writer.writeheader()
        for i, m in enumerate(matches, 1):
            writer.writerow({
                "rank": i,
                "similarity": m["similarity"],
                "source": m["source"],
                "text": m["text"].replace("\n", " | "),
            })

    # Reproducibility
    repro_dir = output_dir / "reproducibility"
    repro_dir.mkdir(exist_ok=True)
    cmd = f'python {Path(__file__).name} --query "{query}" --output {output_dir}'
    if is_demo:
        cmd = f'python {Path(__file__).name} --demo --output {output_dir}'
    (repro_dir / "commands.sh").write_text(f"#!/bin/bash\n{cmd}\n", encoding="utf-8")

    return report_path


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="UKB Navigator — Semantic search across UK Biobank schema"
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--query", "-q", help="Natural language search query")
    group.add_argument("--field", "-f", help="Look up a specific UKB field ID")
    group.add_argument("--demo", action="store_true", help="Run demo with pre-cached results")
    group.add_argument("--embed", help="Path to ukb_schema.csv to embed")

    parser.add_argument("--output", "-o", help="Output directory for report")
    parser.add_argument("--n-results", type=int, default=10, help="Number of results (default: 10)")
    parser.add_argument("--db-path", type=str, default=str(DEFAULT_DB_PATH),
                        help="ChromaDB storage path")
    parser.add_argument("--schema-txt", type=str, help="Optional schema_27.txt path")
    parser.add_argument("--force", action="store_true", help="Force re-embedding")

    args = parser.parse_args()

    if args.embed:
        schema_csv = Path(args.embed)
        schema_txt = Path(args.schema_txt) if args.schema_txt else None
        result = embed_schema(schema_csv, schema_txt, Path(args.db_path), force=args.force)
        print(result)
        return

    if args.demo:
        query = "blood pressure and hypertension"
        matches = DEMO_RESULTS
        out_dir = Path(args.output) if args.output else Path("/tmp/ukb_demo")
        report_path = generate_report(query, matches, out_dir, is_demo=True)
        print(f"Demo report: {report_path}")
        # Also print summary to stdout
        print(f"\nUKB Navigator Demo: {query}")
        print("=" * 50)
        for i, m in enumerate(matches, 1):
            first_line = m["text"].split("\n")[0] if "\n" in m["text"] else m["text"][:80]
            print(f"  [{i}] ({m['similarity']:.3f}) {first_line}")
        return

    # Live query
    db_path = Path(args.db_path)
    if args.field:
        question = f"UK Biobank field {args.field}"
    else:
        question = args.query

    matches = query_schema(question, n_results=args.n_results, db_path=db_path)

    if not matches:
        print("No matches found. Run --embed first to index the UKB schema.")
        sys.exit(1)

    # Print summary
    print(f"UKB Schema Search: {question}")
    print("=" * 50)
    for i, m in enumerate(matches, 1):
        first_line = m["text"].split("\n")[0] if "\n" in m["text"] else m["text"][:80]
        print(f"  [{i}] ({m['similarity']:.3f}) {first_line}")

    # Generate report if output specified
    if args.output:
        out_dir = Path(args.output)
        report_path = generate_report(question, matches, out_dir)
        print(f"\nFull report: {report_path}")


if __name__ == "__main__":
    main()
