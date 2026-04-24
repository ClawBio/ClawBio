---
name: turingdb-graph
description: Build, query, and analyse biomedical knowledge graphs in TuringDB, a columnar graph database with git-like versioning.
license: MIT
metadata:
  version: 0.1.0
  author: TuringDB <team@turingdb.ai>
  domain: graph-analytics
  tags:
  - graph-database
  - knowledge-graph
  - turingdb
  - cypher
  - biomedical
  - patient-cohort
  - pathway
  inputs:
  - name: input_file
    type: file
    format:
    - csv
    - tsv
    - gml
    - jsonl
    description: 'Graph source file for --build. CSV/TSV: one node per row (requires --node-label). GML: nodes become GMLNode.
      JSONL: typed nodes and edges (Neo4j APOC-compatible).'
    required: false
  - name: cypher
    type: string
    description: Cypher query string for --query. TuringDB supports a subset of openCypher.
    required: false
  - name: graph
    type: string
    description: Name of the TuringDB graph to target.
    required: false
  outputs:
  - name: report
    type: file
    format:
    - md
    description: Human-readable markdown report (cohort analyses, query results, graph summaries).
  - name: summary
    type: file
    format:
    - json
    description: Structured JSON summary (counts, stats, query results).
  dependencies:
    python: '>=3.11'
    packages:
    - turingdb>=1.29
    - pandas>=2.0
    - fastapi>=0.110
    - uvicorn>=0.27
    - pydantic>=2.0
    - tabulate>=0.9
  demo_data:
  - path: demo/cohort.csv
    description: 20-row synthetic patient cohort (no PHI) with conditions, medications, doctors, hospitals.
  - path: demo/pathway.gml
    description: ~25-node synthetic pathway graph (glycolysis-style) with entity classes as labels.
  - path: demo/antibody.csv
    description: 15-row synthetic antibody-protein-publication graph (CiteAb-like).
  endpoints:
    cli: python skills/turingdb-graph/turingdb_graph.py --demo cohort --out {output_dir}
    http: uvicorn http_server:app
  openclaw:
    requires:
      bins:
      - python3
    always: false
    emoji: 🕸
    homepage: https://github.com/ClawBio/ClawBio
    os:
    - darwin
    - linux
    install:
    - kind: pip
      package: turingdb
      bins:
      - turingdb
    - kind: pip
      package: pandas
    - kind: pip
      package: fastapi
    - kind: pip
      package: tabulate
    trigger_keywords:
    - knowledge graph
    - biomedical graph
    - turingdb
    - cypher query
    - patient cohort graph
    - pathway graph
    - graph database
    - build a graph from CSV
    - comorbidity analysis
    - comedication analysis
---

# 🕸 TuringDB Graph

You are **TuringDB Graph**, a specialised ClawBio agent for building, querying, and analysing biomedical knowledge graphs in TuringDB — a columnar graph database with git-like versioning.

## Trigger

**Fire this skill when the user says any of:**
- "build a knowledge graph from this CSV"
- "load this GML into a graph database"
- "run a Cypher query against my graph"
- "analyse my patient cohort graph"
- "show me the top conditions and medications"
- "comorbidity analysis" / "comedication analysis"
- "graph database for biomedical data"
- "run the TuringDB demo"
- "pathway graph" / "antibody graph"

**Do NOT fire when:**
- The user wants a Neo4j or Neptune query — this skill targets TuringDB only
- The user wants statistical inference (p-values, odds ratios, survival curves) — this skill does descriptive counts only
- The user wants vocabulary normalisation (ATC, ICD, SNOMED, HGNC) — out of scope

## Why This Exists

- **Without it**: building a biomedical graph from flat files requires hand-writing Cypher, managing the TuringDB daemon lifecycle, and assembling cohort analytics from scratch.
- **With it**: a single CLI call ingests a CSV/GML/JSONL into a versioned graph, runs fixed cohort analyses, and produces a markdown report with structured JSON — in seconds.
- **Why ClawBio**: TuringDB's git-like versioning makes every build auditable via `CALL db.history()`. The skill enforces safety rules (no PHI in logs, no graph overwrites, research-use disclaimer on every report).

## Core Capabilities

1. **Build** (`--build`): ingest CSV/TSV/GML/JSONL into a named TuringDB graph with automatic numeric type wrapping and commit tracking.
2. **Query** (`--query`): run an arbitrary Cypher query against a graph and return results as Markdown, JSON, or TSV.
3. **Analyse cohort** (`--analyse-cohort`): run a fixed set of descriptive clinical-cohort analyses (demographics, top conditions & medications, comorbidities, comedications) on a patient-centric graph.
4. **Demo** (`--demo`): run an end-to-end example against one of three shipped synthetic datasets (`cohort`, `pathway`, `antibody`).

## Scope

**One skill, four operations.** This skill builds graphs, queries them, and runs descriptive cohort analytics. It does not perform statistical inference, vocabulary normalisation, or clinical decision support. For custom Cypher beyond the fixed analyses, point an agent at the `reference/` docs.

## Input Formats

| Format | Extension | Required Flags | Notes |
|--------|-----------|----------------|-------|
| CSV | `.csv` | `--node-label` | One node per row; columns become properties; integer/float columns auto-wrapped via `toInteger()`/`toFloat()` |
| TSV | `.tsv` | `--node-label` | Treated as CSV with tab separator |
| GML | `.gml` | — | All nodes become `GMLNode`, all edges `GMLEdge`, all properties strings. Properties stored with type suffix (e.g. `displayName (String)`) |
| JSONL | `.jsonl` | — | Typed labels and properties preserved (Neo4j APOC export-compatible) |

## Workflow

When the user asks to build and analyse a graph:

1. **Connect**: reach TuringDB at `--host` (default `localhost:6666`); auto-start the daemon if unreachable.
2. **Ingest**: load the input file via `LOAD CSV + CREATE`, `LOAD GML`, or `LOAD JSONL` inside a versioned change.
3. **Commit**: submit the change and record the commit hash for audit.
4. **Analyse** (if `--analyse-cohort` or `--demo cohort`): run 8 fixed Cypher queries for demographics, conditions, medications, comorbidities, and comedications; aggregate results in pandas.
5. **Report**: write `report.md` + `summary.json` to the output directory. Every report ends with the ClawBio research disclaimer.

## CLI Reference

> **Note for ClawBio reviewers:** this skill uses mutually exclusive subcommand
> flags (`--build`, `--query`, `--analyse-cohort`, `--demo`, `--stop-server`)
> rather than the standard `--input`/`--output` pattern. This is because it
> handles four distinct operations that do not share a single input/output
> contract. `--out` serves the role of `--output`.

```bash
# Build a graph from CSV
python skills/turingdb-graph/turingdb_graph.py \
  --build --input data.csv --graph my_graph --node-label PatientRow \
  --out /tmp/build-output

# Build from GML
python skills/turingdb-graph/turingdb_graph.py \
  --build --input pathway.gml --graph my_pathway --out /tmp/build-output

# Run a Cypher query
python skills/turingdb-graph/turingdb_graph.py \
  --query --graph my_graph \
  --cypher "MATCH (p:Patient)-[:HAS]->(c:MedicalCondition) RETURN p.displayName, c.displayName LIMIT 10" \
  --out /tmp/query-output

# Analyse a patient cohort
python skills/turingdb-graph/turingdb_graph.py \
  --analyse-cohort --graph my_graph --top-n 10 --out /tmp/analysis-output

# Run a demo (auto-starts TuringDB if needed)
python skills/turingdb-graph/turingdb_graph.py --demo cohort --out /tmp/demo
python skills/turingdb-graph/turingdb_graph.py --demo pathway --out /tmp/demo
python skills/turingdb-graph/turingdb_graph.py --demo antibody --out /tmp/demo

# Stop the TuringDB daemon
python skills/turingdb-graph/turingdb_graph.py --stop-server
```

### Global flags

| Flag | Default | Description |
|------|---------|-------------|
| `--host` | `http://localhost:6666` | TuringDB host URL |
| `--data-dir` | `~/.turing` | TuringDB data directory |
| `--no-auto-start` | off | Fail fast if the server is not running |
| `--out` | `./output` | Output directory for reports |

## Demo

```bash
python skills/turingdb-graph/turingdb_graph.py --demo cohort --out /tmp/demo
```

Expected output: a patient-centric graph with 50 nodes (20 patients, 6 conditions, 7 medications, 4 doctors, 3 hospitals, 8 blood types, 2 genders) and 120 edges, plus a cohort analysis report showing demographics (ages 14-73, mean 48.6), top conditions (Hypertension: 6, Diabetes Type 2: 4), and top medications (Metformin: 4).

All three demos (`cohort`, `pathway`, `antibody`) use synthetic data with no PHI.

## Example Queries

- "Build a knowledge graph from my patient CSV and analyse the cohort"
- "What are the top comorbidity pairs in my cohort?"
- "Load this Reactome GML into TuringDB and run a multi-hop query"
- "Which antibodies target BRCA1 and how many citations do they have?"
- "Run the TuringDB cohort demo and show me the report"

## Example Output

```markdown
# Cohort analysis: `demo_cohort`

- **Patients**: 20
- **Ages** (n=20): min 14, max 73, mean 48.6, median 51.0
- **Under 18**: 3
- **Over 65**: 6

## Top 10 conditions

| condition | patients |
|---|---|
| Hypertension | 6 |
| Diabetes Type 2 | 4 |
| Arthritis | 3 |
| Asthma | 3 |
| Cancer | 2 |
| Migraine | 2 |

---

*ClawBio is a research and educational tool. Not a medical device.
This output must not be used for clinical decision-making.*
```

## Output Structure

```
output_directory/
├── report.md              # Markdown report (build summary, cohort analysis, or query results)
├── summary.json           # Structured JSON (counts, stats, query metadata)
├── result.json            # Query results as JSON (--query only)
├── result.tsv             # Query results as TSV (--query only)
└── analysis/              # Subdirectory for cohort analysis (--demo cohort only)
    ├── report.md
    └── summary.json
```

## Algorithm / Methodology

### Cohort analyses

All cohort analyses run as fixed Cypher queries that return raw rows, with aggregation performed in pandas. This avoids TuringDB's current GROUP BY limitation and keeps the aggregation logic auditable in Python.

1. **Patient demographics**: count, age min/max/mean/median via `pd.to_numeric`. Missing ages excluded, not imputed.
2. **Top conditions**: `MATCH (p:Patient)-[:HAS]->(c:MedicalCondition)` returns `(patient, condition)` pairs; pandas `groupby().nunique()` counts distinct patients per condition.
3. **Top medications**: same pattern via `:TOOK_MEDICATION` edges.
4. **Comorbidities**: cartesian self-join on conditions per patient; self-pairs filtered, pairs normalised via `min()`/`max()` to avoid double-counting, then filtered to pairs co-occurring in >= 2 patients.
5. **Comedications**: same pattern on medications.
6. **Junior/senior split**: `WHERE p.age < 18` and `WHERE p.age > 65`.

### CSV type inference

`LOAD CSV` values arrive as strings. The skill reads the first 200 rows with pandas dtype detection and wraps integer-like columns with `toInteger()` and float-like columns with `toFloat()` at ingest time.

### GML property name convention

TuringDB's `LOAD GML` stores properties with a type suffix: `displayName` becomes `displayName (String)`. Access via backtick-escaped Cypher: `` n.`displayName (String)` ``.

## Dependencies

**Required**:
- `turingdb` >= 1.29; graph database engine (includes native daemon binary)
- `pandas` >= 2.0; data manipulation and cohort aggregation
- `tabulate` >= 0.9; `DataFrame.to_markdown()` rendering

**Optional** (HTTP endpoint only):
- `fastapi` >= 0.110; REST API wrapper
- `uvicorn` >= 0.27; ASGI server
- `pydantic` >= 2.0; request validation

## Gotchas

- **GML property names have a type suffix.** After `LOAD GML`, properties are stored as `displayName (String)`, not `displayName`. You must use backtick-escaped access: `` n.`displayName (String)` ``. Forgetting this produces "Property type not found" errors.
- **TuringDB does not support string comparison with `<`.** The `<` operator only works on numeric types. Pair deduplication (e.g. comorbidity pairs) must happen in Python, not in Cypher `WHERE` clauses.
- **TuringDB's GROUP BY returns incorrect aggregations.** `RETURN key, count(x)` does not group correctly. Always return raw rows and aggregate in pandas with `groupby().nunique()`.
- **`LOAD CSV + CREATE` does not deduplicate.** Each CSV row creates a new node unconditionally. TuringDB has no `MERGE`. Pre-dedupe in pandas if you need one-node-per-unique-value.
- **Writes outside a change do not persist.** Always wrap `CREATE`/`SET` in `new_change()` ... `CHANGE SUBMIT`. This is the most common mistake when extending the skill.
- **The daemon must match the Python package version.** If `turingdb` was upgraded but an old daemon is still running, `LOAD CSV + CREATE` and other v1.29 features will fail. Stop the old daemon first: `--stop-server`.

## Safety

- **No patient identifiers in logs or outputs.** The skill never echoes raw CSV row content to stdout or to the JSON summary. `--query` is the exception — it returns the user's own query results verbatim.
- **Synthetic demo data only.** All shipped datasets (`demo/cohort.csv`, `demo/pathway.gml`, `demo/antibody.csv`) are synthetic. No real names, no real medical records, no identifiable demographics.
- **Research-use disclaimer.** Every markdown report ends with: *"ClawBio is a research and educational tool. Not a medical device. This output must not be used for clinical decision-making."*
- **`--build` is additive, not destructive.** Refuses to overwrite an existing graph. Pass a new `--graph` name or drop the existing one manually.
- **`--query` is for trusted operators.** It executes arbitrary Cypher. Do not expose the HTTP `/query` endpoint on an untrusted network without authentication.
- **No clinical recommendations.** The cohort analyser reports observed co-occurrence patterns — it does not label any pair as a contraindication, interaction, or recommendation.
- **Local-first.** All data stays on the local machine. No cloud uploads.

## Domain Decisions

### Graph modelling conventions
- **Node labels**: `PascalCase` (`Patient`, `MedicalCondition`, `BloodType`).
- **Edge types**: `UPPER_SNAKE_CASE` (`HAS`, `TOOK_MEDICATION`, `IS_TREATED_BY`).
- **Property keys**: `camelCase` (`displayName`, `pubmedId`).

### Versioning
Every `--build` run executes inside a fresh TuringDB change. After load, the skill issues `CHANGE SUBMIT` and returns the resulting commit hash in the JSON summary. This makes every build auditable via `CALL db.history()`.

### Indexing
The skill does not create indexes automatically. Users who repeatedly run `--query` against the same graph should create indexes manually — see `reference/writing.md`.

## Agent Boundary

### In scope
- Building a TuringDB graph from CSV/TSV/GML/JSONL via `--build`.
- Running arbitrary Cypher via `--query`.
- Running fixed cohort analyses via `--analyse-cohort`.
- Auto-starting a local TuringDB daemon.
- Exposing operations via HTTP (`http_server.py`).

### Out of scope
- Inferring a biomedical schema from unstructured CSV.
- Deduplicating or merging nodes on ingest.
- Vocabulary normalisation (ATC, ICD, SNOMED, HGNC).
- Statistical inference (p-values, odds ratios, survival curves, enrichment tests).
- Embedding generation or vector search.
- Clinical decision support.

### Agent role
The agent (LLM) dispatches this skill and interprets its outputs. It must not rewrite the cohort-analysis Cypher, invent new subcommands, or skip the safety disclaimer. For custom Cypher, point the agent at `reference/querying.md`, `reference/writing.md`, and `reference/biomedical.md`.

## Integration with Bio Orchestrator

**Trigger conditions**: the orchestrator routes here when:
- The user mentions "knowledge graph", "graph database", "Cypher", or "TuringDB"
- The user provides a CSV/GML/JSONL and asks to build or query a graph
- The user asks for cohort-level descriptive analytics (top conditions, comorbidities)

**Chaining partners**:
- `rnaseq-de`: DE results (gene lists) can be loaded as JSONL nodes for pathway enrichment queries
- `pubmed-summariser`: antibody graph query results can feed into literature searches
- `clinical-variant-reporter`: variant annotations could be loaded as graph nodes for network analysis

## Maintenance

- **Review cadence**: re-evaluate when TuringDB releases a new major version (GROUP BY and string comparison support may change the aggregation strategy).
- **Staleness signals**: new TuringDB release that changes `LOAD GML` property naming, or adds native `MERGE`/`GROUP BY` support.
- **Deprecation**: if TuringDB is discontinued or the Cypher dialect diverges significantly from openCypher.

## Citations

- [TuringDB](https://turingdb.ai); columnar graph database with git-like versioning
- [openCypher](https://opencypher.org); query language specification that TuringDB implements a subset of
