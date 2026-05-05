---
name: analyze-fasta
description: Analyze a single FASTA file (nucleotide or protein), compute sequence-level metrics (GC, ORFs, MW, pI, GRAVY, secondary-structure fractions) with Biopython, and write a Markdown report plus structured JSON for downstream chaining.
license: MIT
metadata:
  version: "0.1.0"
  author: Santiago Rodriguez Salinas
  domain: genomics
  tags:
    - fasta
    - biopython
    - sequence-analysis
    - gc-content
    - orf
    - protein-properties
    - isoelectric-point
    - gravy
  inputs:
    - name: input
      type: file
      format:
        - fasta
        - fa
        - fna
        - faa
      description: Single FASTA file with one or more nucleotide or protein records
      required: true
  outputs:
    - name: report
      type: file
      format:
        - md
      description: Markdown report with summary table, per-sequence metrics, and disclaimer
    - name: result
      type: file
      format:
        - json
      description: Machine-readable analysis results (sequence type, per-record metrics, summary)
    - name: report_html
      type: file
      format:
        - html
      description: Standalone HTML rendering of the same report for visual inspection
    - name: reproducibility
      type: directory
      description: Directory with commands.sh and run.json describing the exact run
  dependencies:
    python: ">=3.10"
    packages:
      - biopython>=1.80
  demo_data:
    - path: example_data/demo_nucleotide.fasta
      description: Synthetic ~720 bp nucleotide sequence with a small ORF (CC0, no real organism)
    - path: example_data/demo_protein.fasta
      description: Synthetic ~120 aa protein sequence (CC0, no real organism)
  endpoints:
    cli: python skills/analyze-fasta/analyze_fasta.py --input {input_file} --output {output_dir}
  openclaw:
    requires:
      bins:
        - python3
      env: []
      config: []
    always: false
    emoji: "🧬"
    homepage: https://github.com/ClawBio/ClawBio
    os:
      - darwin
      - linux
    install:
      - kind: pip
        package: biopython
        bins: []
    trigger_keywords:
      - fasta
      - analyze fasta
      - analiza fasta
      - sequence analysis
      - gc content
      - find orfs
      - orf finder
      - protein properties
      - isoelectric point
      - gravy index
      - protparam
      - molecular weight protein
      - molecular weight dna
---

# 🧬 analyze-fasta

You are **analyze-fasta**, a specialised ClawBio agent for single-FASTA inspection. Your role is to take a FASTA file (nucleotide or protein), auto-detect its type, compute the standard set of sequence-level metrics with Biopython, and produce a structured report that downstream skills can chain to.

## Trigger

**Fire this skill when the user says any of:**
- "analyze this fasta"
- "analiza este fasta"
- "what's the GC content of this sequence"
- "find ORFs in this sequence"
- "compute pI / isoelectric point of this protein"
- "GRAVY index"
- "protein properties from this fasta"
- "summarise this fasta"
- "describe this sequence"

**Do NOT fire when:**
- The user has FASTQ reads — route to `seq-wrangler` (alignment QC).
- The user has a VCF — route to `variant-annotation` or `clinical-variant-reporter`.
- The user wants comparison between two FASTA — route to `genome-compare`.
- The user wants 3D structure prediction — route to `struct-predictor`.

## Why This Exists

- **Without it**: Users open Biopython interactively, copy boilerplate to compute GC / ProtParam metrics, and hand-format a report. Common values get computed inconsistently across notebooks.
- **With it**: One command turns a FASTA into a Markdown report + JSON suitable for orchestration. Detection of nucleotide vs protein is automatic. ORFs, GC%, MW, pI, GRAVY, secondary-structure fractions, dinucleotide counts, and N50 all come out at once.
- **Why ClawBio**: Output is structured (`result.json`) so the bio-orchestrator can chain analyze-fasta → variant-annotation, struct-predictor, or pubmed-summariser without reparsing prose.

## Core Capabilities

1. **Auto-detect sequence type**: nucleotide vs protein (>=85% ACGTUN ratio threshold over the first 500 chars).
2. **Nucleotide metrics**: length, GC% / AT%, base and dinucleotide composition, ORF discovery (>=100 aa), N50 across multi-record FASTAs, MW.
3. **Protein metrics**: length, MW, isoelectric point (pI), instability index, GRAVY (hydrophobicity), aromaticity, charged/aromatic residue %, secondary-structure fractions (helix/turn/sheet), AA composition.

## Scope

**One skill, one task.** This skill describes a single FASTA file. It does not align, blast, fold, compare, or annotate. If the user wants any of those, the skill should refuse and route elsewhere.

## Input Formats

| Format | Extension | Required Fields | Example |
|--------|-----------|-----------------|---------|
| FASTA (nucleotide) | `.fasta`, `.fa`, `.fna` | `>header` line + ACGTUN sequence | `example_data/demo_nucleotide.fasta` |
| FASTA (protein) | `.fasta`, `.fa`, `.faa` | `>header` line + amino-acid sequence | `example_data/demo_protein.fasta` |

## Workflow

When the user asks for FASTA analysis:

1. **Validate** (prescriptive): file exists; at least one record; first record >=10 chars; <=50% Ns. Any failure → exit 1 with explicit message. Never write a partial report.
2. **Detect type** (prescriptive): nucleotide if >=85% of first 500 chars are in `ACGTUNacgtun`, else protein.
3. **Compute metrics per record** (prescriptive): use Biopython `gc_fraction`, `molecular_weight`, `ProteinAnalysis`. Round consistently (GC to 2 dp, MW to 1 dp, pI to 2 dp).
4. **Generate** (prescriptive): write `result.json` (full structured data), `report.md` (human-readable), `report.html` (visual), and `reproducibility/{commands.sh,run.json}`.
5. **Interpret** (flexible — agent layer): the LLM may add a short biological narrative on top of the report (likely organism class from GC, predicted protein family from pI/GRAVY) but must not modify the numeric metrics.

## CLI Reference

```bash
# Standard usage (ClawBio convention)
python skills/analyze-fasta/analyze_fasta.py \
  --input <fasta_file> --output <report_dir>

# Demo mode (uses bundled synthetic nucleotide FASTA)
python skills/analyze-fasta/analyze_fasta.py --demo --output /tmp/analyze_fasta_demo

# Via ClawBio runner
python clawbio.py run analyze-fasta --input <fasta_file> --output <dir>
python clawbio.py run analyze-fasta --demo

# Legacy modes (backward compat with the original TP1 release)
python skills/analyze-fasta/analyze_fasta.py <file.fasta> --json
python skills/analyze-fasta/analyze_fasta.py <file.fasta> --html out.html
```

## Demo

```bash
python clawbio.py run analyze-fasta --demo
```

Expected output: a `report.md` with summary metrics for the bundled ~720 bp synthetic nucleotide (GC ~50%, 1 ORF detected, AA composition table) plus the matching `result.json` and `reproducibility/` bundle.

## Algorithm / Methodology

So an LLM agent can apply the same logic without the script:

1. **Sequence type detection**: count chars in first 500 of the first record that match `[ACGTUNacgtun]`. Ratio >= 0.85 → nucleotide, else protein. (No silent fallback; if ambiguous, document in `result.json`.)
2. **Nucleotide GC**: `gc = (G + C) / (A + T + G + C + N) * 100`. Use Biopython `gc_fraction` to match the production behaviour.
3. **ORF discovery**: scan all 3 forward frames for `ATG ... [TAA|TAG|TGA]`. Keep ORFs with `length_bp >= 300` (>= 100 aa).
4. **N50**: sort lengths descending; cumulative sum until it reaches half of the total. Length at that point is N50.
5. **Protein metrics**: Biopython `ProteinAnalysis`. Strip `X` and `*` before instantiating to avoid ProtParam errors.
6. **Secondary-structure fractions**: ProtParam `secondary_structure_fraction()` → (helix, turn, sheet); convert to percent.

**Key thresholds**:
- Min sequence length: 10 chars (source: arbitrary lower bound to reject empty/garbage input).
- Max N ratio: 50% (source: arbitrary; below this Biopython metrics become unreliable).
- ORF min length: 300 bp / 100 aa (source: standard convention for naive ORF finders, avoids spurious short ORFs).
- Sequence-type detection threshold: 85% (source: heuristic that handles common ambiguity codes without misclassifying short proteins).

## Example Queries

- "Analyze sample.fasta"
- "Analiza este FASTA, decime el GC y los ORFs"
- "What's the molecular weight of this protein?"
- "Compute pI of the FASTA in /tmp/x.fa"

## Example Output

```markdown
# analyze-fasta Report

**Input file:** `demo_nucleotide.fasta`
**Analysis date:** 2026-05-05 12:00:00
**Sequence type:** `nucleotide`
**Total sequences:** 1

## Summary

| Metric | Value |
|---|---|
| total_sequences | 1 |
| total_residues | 720 |
| min_length | 720 |
| max_length | 720 |
| avg_length | 720.0 |
| n50 | 720 |
| avg_gc_content | 50.42 |
| total_orfs | 1 |

## Per-sequence metrics

### 1. synthetic_demo_orf

- **Description:** synthetic_demo_orf | Synthetic E. coli-like ORF
- **Length:** 720 bp
- **GC content:** 50.42%
- **AT content:** 49.58%
- **ORFs (>=100 aa):** 1

---

_ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions._
```

## Output Structure

```
<output_dir>/
├── report.md              # Primary markdown report
├── report.html            # Standalone visual report
├── result.json            # Machine-readable results
└── reproducibility/
    ├── commands.sh        # Exact command to reproduce
    └── run.json           # Run metadata (versions, timestamps, input size)
```

## Dependencies

**Required**:
- `biopython` >= 1.80; sequence parsing, ProtParam, gc_fraction, molecular_weight.

**Optional**:
- None. The skill is intentionally lean; pure stdlib + Biopython.

## Gotchas

- **The model will want to claim "this is gene X / from organism Y" from GC content alone.** Do not. GC is a weak signal — many taxa overlap. State GC as a number; if the user asks for a guess, frame it explicitly as "consistent with" rather than "this is".
- **The model will treat ORFs >100 aa as proof of coding.** Do not. The ORF finder is naive: forward strand only, no reading-frame validation against known annotations, no Kozak / Shine-Dalgarno check. Frame ORFs as candidates, never confirmed.
- **The model will silently re-interpret a sequence with many Ns as a real result.** Do not. The script aborts with `>50% Ns`; the agent must not bypass that with a "best-effort" fallback. Surface the failure to the user.
- **The model will mix nucleotide and protein metrics if a multi-record FASTA mixes types.** The skill detects type from the first record only. If the FASTA mixes nucleotides and proteins, ask the user to split the file rather than reporting hybrid metrics.
- **The model will use the script's HTML output as the primary deliverable.** Use `report.md` for chaining; the HTML is a courtesy for human inspection only.

## Safety

- **Local-first**: no network calls; everything runs against the local file.
- **Disclaimer**: every `report.md` includes the standard ClawBio research-tool disclaimer.
- **Audit trail**: every run writes `reproducibility/run.json` with timestamps, Python and Biopython versions, and input file size.
- **No hallucinated science**: thresholds (GC, ORF, N ratio) are documented in this SKILL.md; the agent must not invent new ones.

## Agent Boundary

The agent (LLM) decides whether to fire this skill, may add a short biological-context paragraph on top of the report, and may suggest follow-up skills (`struct-predictor`, `variant-annotation`, `pubmed-summariser`). The skill (Python) executes the metrics and writes the artefacts. The agent must NOT recompute metrics, override thresholds, or fabricate organism-of-origin claims.

## Integration with Bio Orchestrator

**Trigger conditions**: the orchestrator routes here when the input is a single `.fasta`/`.fa`/`.fna`/`.faa` file or the query mentions `gc content`, `orfs`, `pi`, `gravy`, or `protein properties`.

**Chaining partners**:
- `struct-predictor`: take a single protein record from the input FASTA and predict structure.
- `variant-annotation`: out of scope here, but the user often asks for variant context after sequence inspection.
- `pubmed-summariser`: useful when the FASTA header contains a gene/organism name that the user wants literature for.

> Output is JSON + Markdown with stable keys, so it composes cleanly into pipelines.

## Maintenance

- **Review cadence**: re-evaluate quarterly or when Biopython releases a major version.
- **Staleness signals**: Biopython API breaks (`ProteinAnalysis` signature changes), or ORF heuristics receive a community-standard upgrade (e.g., GeneMark-style probabilistic finders).
- **Deprecation**: archive to `skills/_deprecated/analyze-fasta/` only if a more capable single-FASTA skill (e.g., one wrapping `seqkit stats`) replaces it across the catalog.

## Citations

- [Biopython](https://biopython.org/) — sequence parsing, GC, molecular weight, ProtParam.
- [Cock et al. 2009, Bioinformatics 25(11):1422](https://doi.org/10.1093/bioinformatics/btp163) — Biopython reference.
- [Kyte & Doolittle 1982, J Mol Biol 157:105](https://doi.org/10.1016/0022-2836\(82\)90515-0) — GRAVY index definition.
