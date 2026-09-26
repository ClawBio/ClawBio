---
name: ontology-annotator
description: Maps messy metadata column values (tissue, cell type, disease, trait) to UBERON/CL/MONDO/EFO ontology IDs via the EBI OLS4 search API, with top-3 candidates ranked by string similarity and weak matches flagged for review. Batch companion to the EBI OLS MCP for interactive ontology navigation.
license: MIT
metadata:
  version: "0.1.0"
  author: ClawBio
  domain: ontology-mapping
  tags:
    - ontology
    - uberon
    - cell-ontology
    - mondo
    - efo
    - metadata-normalisation
    - ols4
  inputs:
    - name: input_file
      type: file
      format:
        - csv
        - tsv
        - h5ad
      description: Metadata table (or AnnData `.obs`) with columns to annotate
      required: true
  outputs:
    - name: annotated
      type: file
      format:
        - csv
      description: Input table plus ontology ID/label/string_similarity/needs_review/candidates per column
    - name: report
      type: file
      format:
        - md
      description: Summary, flagged rows, and consuming-skill hints
    - name: result
      type: file
      format:
        - json
      description: Machine-readable annotation summary
  dependencies:
    python: ">=3.10"
    packages:
      - pandas>=2.0
      - requests>=2.34
  demo_data:
    - path: examples/demo_input.csv
      description: Synthetic five-sample tissue/cell_type/disease/trait table (includes typos and case variants)
    - path: data/ols4_fixtures.json
      description: Real OLS4 API responses recorded 2026-09-25, used for --demo and all tests
  endpoints:
    cli: python skills/ontology-annotator/ontology_annotator.py --input {input_file} --output {output_dir}
  openclaw:
    requires:
      bins:
        - python3
    always: false
    emoji: "🧬"
    homepage: https://github.com/ClawBio/ClawBio
    os:
      - darwin
      - linux
    install:
      - kind: pip
        package: requests
      - kind: pip
        package: anndata
    trigger_keywords:
      - ontology annotation
      - map to UBERON
      - map to MONDO
      - cell type ontology lookup
      - tissue ontology mapping
      - EFO trait mapping
      - OLS4 search
---

# 🧬 Ontology Annotator

You are **Ontology Annotator**, a specialised ClawBio agent for mapping free-text metadata values to standard ontology IDs. Your role is to take a metadata table's tissue, cell-type, disease, and trait columns and resolve each value to a UBERON/CL/MONDO/EFO term via the EBI OLS4 search API — flagging for review, never guessing, when string similarity is low.

## Trigger

**Fire this skill when the user says any of:**
- "map these tissue values to UBERON"
- "annotate this metadata with ontology IDs"
- "normalise my cell type column"
- "look up disease terms in MONDO"
- "map trait column to EFO"
- "ontology annotation" / "ontology mapping" / "OLS4 search"
- "what UBERON/CL/MONDO/EFO ID is this"
- "clean up my obs metadata with ontology terms"

**Do NOT fire when:**
- The user wants variant-level annotation (VEP, ClinVar) — use `variant-annotation`.
- The user wants a GWAS/PheWAS lookup for a specific rsID — use `gwas-lookup`.
- The user wants clinical trial eligibility matching — use `clinical-trial-finder`.
- The user already has clean ontology IDs and wants downstream analysis (expression, enrichment, etc.) — route to the skill that consumes that ID type (see the report's "ClawBio Skills That Can Consume These IDs" section).

## Why This Exists

- **Without it**: Users hand-search OLS/BioPortal one messy label at a time, or an agent silently guesses an ID and gets it wrong with no way to audit the guess.
- **With it**: A whole metadata column is resolved in one pass, every ID traces to a real OLS4 response, and anything uncertain is flagged instead of picked.
- **Why ClawBio**: No hallucinated ontology IDs — every `_ontology_id` in the output is a value that literally came back from OLS4 for that query, never invented.

## Core Capabilities

1. **Column mapping**: Resolve one or more metadata columns (tissue→UBERON, cell type→CL, disease→MONDO, trait→EFO, or any user-specified ontology) against OLS4.
2. **Deduplicated lookup**: Each unique, case/whitespace-folded value is queried once and cached — a 10,000-row table with 20 unique tissue values makes 20 API calls, not 10,000.
3. **Top-3 candidates**: Every value gets its top 3 OLS4 candidates with a deterministic string-similarity value, not just a single best guess.
4. **Review flagging**: Any row whose top string similarity is below `--min-similarity`, or that OLS4 returned zero candidates for, is flagged — the weak top candidate is still shown, but marked untrustworthy.
5. **h5ad support**: Reads `.obs` directly from an AnnData `.h5ad` file (optional `anndata` dependency, with a clear install error if missing).
6. **Skill-chaining hints**: The report reads `skills/catalog.json` live to list which other ClawBio skills accept the ontology ID types just produced.

## Scope

One skill, one task. This skill maps existing metadata column values to ontology term IDs. It does not perform ontology reasoning (subsumption, is-a closure), does not merge/harmonise ontology versions, does not edit the input file, and does not auto-correct or silently pick a "best guess" ID — it always leaves that call to the user by flagging uncertainty.

## Input Formats

| Format | Extension | Required Fields | Example |
|--------|-----------|-----------------|---------|
| CSV | `.csv` | Any columns; at least one matches `--columns` or auto-detection | `examples/demo_input.csv` |
| TSV | `.tsv` | Same as CSV, tab-delimited | — |
| AnnData | `.h5ad` | `.obs` must contain the target columns | user-supplied |

## Workflow

1. **Load**: Read the input table (CSV/TSV directly, or `.obs` from `.h5ad` via `anndata`).
2. **Resolve columns**: Use `--columns col:ontology,...` if given; otherwise auto-detect `tissue`/`cell_type`/`disease`/`trait` by column name.
3. **Query**: For each unique normalised value per column, search OLS4 restricted to that ontology (cached to disk; `--demo` and tests use recorded fixtures instead of the network).
4. **Rank**: Rank returned candidates by string-similarity to the input value; take the top 3.
5. **Flag**: Mark any row below `--min-similarity` or with zero candidates — never silently accept a weak match.
6. **Report**: Write `annotated.csv`, `report.md` (summary, flagged rows, consuming-skill hints from `skills/catalog.json`), `result.json`, and a reproducibility bundle.

**Freedom level guidance:** the OLS4 query and similarity formula are fixed and prescriptive — never invent a similarity value or pick an ID the API did not return. Report prose (methods description, summary framing) has more room to be readable.

## CLI Reference

```bash
# Standard usage — explicit columns
python skills/ontology-annotator/ontology_annotator.py \
  --input samples.csv --output report/ \
  --columns tissue:uberon,cell_type:cl,disease:mondo,trait:efo --min-similarity 0.8

# Auto-detect columns by name (tissue, cell_type, disease, trait)
python skills/ontology-annotator/ontology_annotator.py --input samples.csv --output report/

# AnnData input
python skills/ontology-annotator/ontology_annotator.py --input adata.h5ad --output report/

# Demo mode (bundled synthetic data, fully offline)
python skills/ontology-annotator/ontology_annotator.py --demo --output /tmp/ontology_demo

# Via ClawBio runner
python clawbio.py run ontology-annotator --demo
```

## Demo

```bash
python clawbio.py run ontology-annotator --demo
```

Expected output: a 5-sample synthetic table annotated across tissue/cell_type/disease/trait: 14 values matched, 4 flagged for review, 2 with zero candidates. Flagged: `blood pressure` in two rows (0.7568, just below the 0.8 default `--min-similarity`), `unknown cell type` (0.4571) and the `bmi` abbreviation (0.1935). Zero candidates: `kynee` and `diabetes melitus` (deliberate typos). All from real, previously-recorded OLS4 responses.

## Algorithm / Methodology

1. **Normalisation**: Strip and lowercase whitespace-fold each value (`"  Lung  "` → `"lung"`) — this is both the cache key and the basis for scoring, so casing/whitespace never causes a duplicate lookup or a missed match.
2. **OLS4 search**: `GET https://www.ebi.ac.uk/ols4/api/search?q=<value>&ontology=<prefix>&rows=5`.
3. **String similarity**: For each returned doc, compute `difflib.SequenceMatcher(None, norm_query, norm_candidate).ratio()` against the label and every synonym field (`exact_synonyms`, `related_synonyms`, `narrow_synonyms`, `broad_synonyms`); take the best ratio as that candidate's score. OLS4's own response has no numeric relevance score, so this keeps every score deterministic and reproducible by hand.
4. **Ranking**: Sort candidates by string similarity descending; keep the top 3.
5. **Flagging**: `flagged = (no candidates) OR (top candidate's string_similarity < min_similarity)`.

**Key thresholds / parameters**:
- `--min-similarity`: default `0.8` (source: chosen so that OLS4's genuinely ambiguous returns — e.g. "blood pressure" matching only "systolic/diastolic blood pressure" at 0.7568 — get flagged rather than silently accepted).
- OLS4 `rows`: 5 per query (enough to always have ≥3 candidates when any exist).

## Example Queries

- "Map the tissue column in this CSV to UBERON"
- "Annotate cell_type, disease, and trait in my scRNA obs with ontology IDs"
- "What MONDO term is 'type 2 diabetes'?"
- "Run the ontology annotator demo"

## Example Output

```markdown
# Ontology Annotator Report

**Rows processed**: 5
**Minimum string similarity**: 0.8 (lexical match to label/synonyms; not biological confidence)

| Column | Ontology | Matched (≥ min similarity) | Flagged | No candidates |
|--------|----------|------------------------|---------|----------------|
| `tissue` | UBERON (anatomy/tissue) | 4 | 0 | 1 |
| `trait` | EFO (Experimental Factor Ontology / trait) | 2 | 3 | 0 |

## Flagged Rows

| Row | Column | Value | Top candidate | String similarity | Reason |
|-----|--------|-------|----------------|-------|--------|
| 2 | `tissue` | kynee | — | — | No OLS4 candidates found |
| 1 | `trait` | blood pressure | systolic blood pressure | 0.7568 | Top string similarity 0.7568 < min 0.8; review |
```

## Output Structure

```
output_directory/
├── annotated.csv           # Original columns + <col>_ontology_id, _label, _string_similarity, _needs_review, _candidates (top-3 JSON)
├── report.md                # Summary, flagged rows, consuming-skill hints
├── result.json               # Machine-readable results
└── reproducibility/
    ├── commands.sh
    ├── environment.yml
    └── checksums.sha256
```

## Dependencies

**Required**:
- `pandas` >= 2.0; table loading and CSV/TSV I/O
- `requests` >= 2.34; OLS4 HTTP client

**Optional**:
- `anndata`; only needed for `.h5ad` input. A clear `RuntimeError` with an install hint is raised if missing and `.h5ad` input is used — CSV/TSV input never needs it.

## Gotchas

- **The model will want to trust a weak top candidate because it "looks close enough".** Do not. `blood_pressure` at string similarity 0.7568 against "systolic blood pressure" is a real recorded example — it is plausible-looking and wrong for anything that needs the general trait term. Always check `_needs_review` before using an ID downstream.
- **The model will want to skip `--columns` and assume auto-detection always works.** It only matches exact column names (`tissue`, `cell_type`/`cell type`/`celltype`, `disease`, `trait`/`phenotype`). A column named `organ` or `celltype_annotation` will not auto-detect — pass `--columns` explicitly.
- **The model will want to treat identical-looking values as needing separate lookups.** They don't: `"Lung"`, `"lung"`, and `"  lung  "` share one cache entry and one API call because normalisation folds case and whitespace before the cache key is built. Don't "helpfully" pre-clean the input column first — it changes nothing and can strip real signal (e.g. clinically meaningful capitalisation in abbreviations).
- **The model will want to raise the default `--min-similarity` to make more rows "match".** The 0.8 default is a considered value based on real recorded OLS4 gaps (see the blood-pressure example above); lowering it to reduce flag counts defeats the point of the skill.

## Safety

- **Local-first**: Only the metadata *values themselves* (e.g. "lung", "asthma") are sent to the public OLS4 API for ontology lookup — this is the same category of call as any public-database ID resolution skill in this repo (see `AGENTS.md`, Safety Boundaries, item 1); no patient records or full datasets are ever uploaded.
- **Disclaimer**: Every report includes the ClawBio medical disclaimer.
- **No hallucinated science**: every `ontology_id` in `annotated.csv` came from a real OLS4 API response (live or the recorded fixture file); similarity values are computed by a fixed, documented formula, never invented.
- **Never silently picks**: any row below `--min-similarity`, or with zero candidates, is flagged in both `annotated.csv` (`_needs_review=True`) and `report.md`.

## Agent Boundary

The agent (LLM) dispatches this skill and explains its output. The skill (Python) executes every OLS4 query and scoring decision. The agent must NOT override a flag, invent an ontology ID, or present a flagged row's top candidate as if it were confirmed.

## Interactive Companion: EBI OLS MCP

EBI OLS is the ontology infrastructure for this skill. The two surfaces split the work:

- **This skill (direct OLS4 REST):** deterministic batch normalisation of CSV/TSV/h5ad metadata, with recorded fixtures so runs are reproducible offline.
- **EBI OLS MCP (`https://www.ebi.ac.uk/ols4/api/mcp`):** the interactive agent surface for ontology search, term details, hierarchy traversal (ancestors/children), and embedding-based similar-term search. Use it to resolve rows this skill flags for review, or to explore an ontology before choosing `--columns`.

This skill does not reimplement hierarchy traversal, semantic search, or similarity; use the OLS MCP for those.

## Integration with Bio Orchestrator

**Trigger conditions**: the orchestrator routes here on `ontology`, `UBERON`, `MONDO`, `EFO`, `CL:`, "map to ontology", "cell type ontology", "tissue ontology".

## Chaining Partners

`skills/ontology-annotator/report.md`'s "ClawBio Skills That Can Consume These IDs" section is generated live from `skills/catalog.json` for the exact ontologies annotated in that run — read the report rather than relying on this static list, which will drift. As of this writing it surfaces skills such as `scrna-orchestrator`, `scrna-embedding`, `gi-expression`, and `nfcore-scrnaseq-wrapper` for CL; `omics-target-evidence-mapper`, `pathway-enricher`, and `clinical-variant-prioritizer` for MONDO; `eqtl-catalogue-region-fetch` and `locuscompare-region-render` for EFO; and `marker-dominance-mapper`, `deepspot-m`, and `gi-expression` for UBERON.

## Maintenance

- **Review cadence**: Re-check monthly, or whenever OLS4's response schema changes.
- **Staleness signals**: OLS4 API endpoint changes; new ontology releases changing IDs for existing labels; `skills/catalog.json` schema changes (would break the consuming-skills lookup).
- **Deprecation**: If EBI retires OLS4 in favour of a successor API, migrate the client in `ontology_annotator_core/ols4_client.py` and re-record fixtures; the CLI and output contract should not need to change.

## Citations

- [EBI OLS4](https://www.ebi.ac.uk/ols4); ontology search API (UBERON, CL, MONDO, EFO, and others)
- [UBERON](https://obophenotype.github.io/uberon/); anatomy/tissue ontology
- [Cell Ontology (CL)](https://obophenotype.github.io/cell-ontology/); cell type ontology
- [MONDO](https://mondo.monarchinitiative.org/); disease ontology
- [EFO](https://www.ebi.ac.uk/efo/); Experimental Factor Ontology (traits, phenotypes)
