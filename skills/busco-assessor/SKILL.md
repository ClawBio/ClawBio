---
name: busco-assessor
description: >-
  Genome, transcriptome, and protein completeness assessment via BUSCO v6.
  Agentic lineage routing from organism description, all three BUSCO modes,
  auto-lineage support, and full demo mode without the BUSCO binary.
license: MIT
metadata:
  version: "0.1.0"
  author: ClawBio Contributors
  domain: genomics
  tags:
    - busco
    - genome-completeness
    - assembly-qc
    - transcriptome
    - lineage
    - orthodb
    - hmmer
    - prokaryote
    - eukaryote
  inputs:
    - name: input
      type: file
      format: [fasta, fna, fa, faa]
      description: Assembly, transcriptome, or protein FASTA (required unless --demo)
      required: true
  outputs:
    - name: report
      type: file
      format: md
      description: Markdown completeness report with interpretation
    - name: result
      type: file
      format: json
      description: Machine-readable completeness scores (C/S/D/F/M/n)
    - name: busco_run
      type: directory
      description: Raw BUSCO outputs (short_summary.txt, full_table.tsv, short_summary.json)
    - name: reproducibility
      type: directory
      description: commands.sh, environment.yml, checksums.sha256
  dependencies:
    python: ">=3.10"
    packages: []
    external:
      - busco>=6.0 (runtime; not required for --demo)
      - hmmer>=3.1 (installed with BUSCO via conda)
      - sepp==4.5.5 (auto-lineage only — v4.5.6 is incompatible)
  demo_data:
    - path: "--demo flag"
      description: Synthetic 5-sequence FASTA with bacteria-like completeness C:95.2%[S:93.1%,D:2.1%],F:2.3%,M:2.5%,n:124
  endpoints:
    cli: python skills/busco-assessor/busco_assessor.py --input {input} --mode genome --output {output_dir}
  openclaw:
    requires:
      bins:
        - python3
      env: []
      config: []
    always: false
    emoji: "🧬"
    homepage: https://busco.ezlab.org/
    os:
      - darwin
      - linux
    install:
      - kind: conda
        package: busco=6.0.0
        channels: [bioconda, conda-forge]
      - kind: conda
        package: sepp=4.5.5
        channels: [bioconda, conda-forge]
    trigger_keywords:
      - "genome completeness"
      - "BUSCO score"
      - "BUSCO assessment"
      - "assembly quality"
      - "check my assembly"
      - "BUSCO genome mode"
      - "completeness metrics"
      - "assembly QC"
      - "transcriptome completeness"
      - "protein set completeness"
      - "auto-lineage"
      - "busco -m genome"
      - "how complete is my genome"
      - "BUSCO bacteria"
      - "run BUSCO"
---

# 🧬 BUSCO Assessor

You are the **busco-assessor**, a specialised ClawBio agent for genome, transcriptome, and protein-set completeness assessment. Your role is to run BUSCO v6 against the correct OrthoDB lineage dataset — inferred automatically from the user's organism description — and produce a reproducible, interpreted completeness report.

## Trigger

**Fire when the user says any of:**
- "genome completeness", "BUSCO score", "BUSCO assessment"
- "assembly quality", "check my assembly", "check assembly completeness"
- "BUSCO genome mode", "BUSCO transcriptome mode", "BUSCO proteins mode"
- "completeness metrics", "assembly QC", "how complete is my genome"
- "BUSCO bacteria", "run BUSCO", "busco -m genome"
- "auto-lineage", "transcriptome completeness", "protein set completeness"

**Do NOT fire when:**
- User wants to align reads → route to `seq-wrangler`
- User wants multi-tool QC aggregation across samples → route to `multiqc-reporter`
- User wants variant calling or annotation → route to `vcf-annotator`
- User wants protein structure prediction → route to `struct-predictor`
- User is asking about genome *assembly* (not quality assessment) → suggest external assemblers

## Why This Exists

- **Without it**: Users must manually browse ~100 OrthoDB lineage datasets, choose the correct `*_odb10/12` for their organism, construct the BUSCO command, and interpret C/S/D/F/M scores from raw text output.
- **With it**: A free-text organism description (e.g. "my E. coli assembly") is sufficient — the skill resolves the lineage, runs BUSCO, parses scores, and produces a structured report with interpretation.
- **Why ClawBio**: Completeness assessment is a prerequisite for downstream genomics (variant calling, annotation, pangenome analysis) and must be reproducible and interpretable without bioinformatics expertise.

## Core Capabilities

1. **Agentic lineage routing** — maps natural-language organism descriptions to the correct BUSCO lineage flag via a curated routing table (`LINEAGE_ROUTING`).
2. **Three assessment modes** — genome, transcriptome, proteins, each with appropriate tool dependencies.
3. **Auto-lineage support** — `--auto-lineage`, `--auto-lineage-euk`, `--auto-lineage-prok` with SEPP 4.5.5 compatibility enforcement.
4. **Score parsing and interpretation** — extracts C/S/D/F/M completeness from `short_summary.txt` and provides plain-language interpretation.
5. **Full demo without BUSCO binary** — synthetic FASTA and output files generated in Python; safe for CI/offline environments.
6. **Reproducibility bundle** — `commands.sh`, `environment.yml` (pinning busco=6.0.0 + sepp=4.5.5), `checksums.sha256`.

## Scope

One skill, one task: **BUSCO completeness assessment**. This skill does NOT assemble genomes, call variants, run read alignment, or annotate genes. For multi-sample QC aggregation of BUSCO results, chain to `multiqc-reporter` (BUSCO module).

## Input Formats

| Format | Extension | BUSCO Mode | Notes |
|--------|-----------|-----------|-------|
| Genome assembly | `.fna`, `.fa`, `.fasta` | `genome` | Scaffolds or contigs |
| Transcriptome | `.fna`, `.fa`, `.fasta` | `transcriptome` | Assembled transcripts |
| Protein sequences | `.faa`, `.fasta` | `proteins` | Amino-acid FASTA |

## Workflow

1. **Validate inputs** — check `--input` exists; check `busco` binary on PATH (skip in `--demo` mode).
2. **Resolve lineage** — apply this decision tree in order:
   - If `--lineage <dataset>` supplied → use it verbatim.
   - If `--auto-lineage*` flag supplied → use it verbatim.
   - If `--organism "<text>"` supplied → call `infer_lineage(text)` to map keywords to lineage flag.
   - If nothing supplied → default to `--auto-lineage` (requires SEPP 4.5.5).
3. **Build BUSCO command** — assemble CLI with `-i`, `-m`, `-c`, `--out-path`, `--out`, and resolved lineage flag.
4. **Execute BUSCO** — `subprocess.run` with 7200s timeout; raise `RuntimeError` on nonzero exit with last 10 stderr lines.
5. **Parse `short_summary.txt`** — regex extraction of C/S/D/F/M/n; glob both `short_summary.txt` and `short_summary.specific.*.txt` patterns.
6. **Parse `full_table.tsv`** — tab-separated rows (skip `#` comment lines); returns per-gene status table.
7. **Write `result.json`** — completeness scores + run parameters.
8. **Write `report.md`** — completeness table, score string, plain-language interpretation, top-10 gene results, disclaimer.
9. **Write reproducibility bundle** — `reproducibility/commands.sh`, `environment.yml`, `checksums.sha256`.

## CLI Reference

```bash
# Genome mode with explicit lineage
python skills/busco-assessor/busco_assessor.py \
  --input assembly.fna --mode genome --lineage bacteria_odb12 \
  --cpu 8 --output /tmp/busco_out

# Genome mode with auto-lineage (prokaryote)
python skills/busco-assessor/busco_assessor.py \
  --input assembly.fna --mode genome --auto-lineage-prok \
  --cpu 8 --output /tmp/busco_out

# Agentic: infer lineage from organism hint
python skills/busco-assessor/busco_assessor.py \
  --input assembly.fna --organism "fruit fly"--output /tmp/busco_out

# Transcriptome mode
python skills/busco-assessor/busco_assessor.py \
  --input transcriptome.fna --mode transcriptome --lineage insecta_odb10 \
  --output /tmp/busco_transcriptome

# Proteins mode
python skills/busco-assessor/busco_assessor.py \
  --input proteins.faa --mode proteins --lineage vertebrata_odb10 \
  --output /tmp/busco_proteins

# Offline demo (no BUSCO binary needed)
python skills/busco-assessor/busco_assessor.py --demo --output /tmp/busco_demo

# Live demo: downloads real S. cerevisiae Mito FASTA + NCBI taxonomy lineage lookup
python skills/busco-assessor/busco_assessor.py --demo-live --output /tmp/busco_live_demo
```

## Demo

### Offline demo (no internet, no BUSCO binary)
```bash
python skills/busco-assessor/busco_assessor.py --demo --output /tmp/busco_demo
```
**Expected:** bacteria-like completeness `C:95.2%[S:93.1%,D:2.1%],F:2.3%,M:2.5%,n:124` — fully synthetic, works in CI.

### Live demo (real data from Ensembl + NCBI Taxonomy)
```bash
python skills/busco-assessor/busco_assessor.py --demo-live --output /tmp/busco_live_demo
```
**What it does — 5 steps:**
1. Downloads *S. cerevisiae* mitochondrial chromosome (22 KB) from Ensembl Genomes release 62
2. Queries NCBI Taxonomy E-utilities API for `Saccharomyces cerevisiae` → resolves `saccharomycetes_odb10`
3. Runs BUSCO if installed, otherwise generates realistic synthetic output
4. Writes `report.md` with completeness table and mitochondrial-genome note
5. Writes reproducibility bundle (commands.sh pins `busco=6.0.0 sepp=4.5.5`)

**Expected output (no BUSCO binary):**
```markdown
Lineage: saccharomycetes_odb10   [NCBI Taxonomy API]
C:2.1%[S:2.1%,D:0.0%],F:0.9%,M:97.0%,n:2137
```
> The low completeness (2.1%) is correct and expected — the mito chromosome only encodes ~15–35 protein-coding genes; most of the 2137 BUSCO orthologs are nuclear genes. This is an educational feature, not a bug.

## NCBI Taxonomy Integration

When `--demo-live` is used (or `--organism` is passed with the `--ncbi` flag), the skill queries the NCBI E-utilities API to resolve the organism's taxonomic lineage and select the most specific BUSCO dataset automatically:

```
esearch  → https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&term={name}&retmode=json
            returns: {"esearchresult": {"idlist": ["4932"]}}

efetch   → https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id=4932&retmode=xml
            returns: XML with <LineageEx> containing {rank, ScientificName} pairs
```

The `NCBI_TO_BUSCO` table maps rank+name pairs (most-specific first) to BUSCO lineages. For *S. cerevisiae*:
- class `Saccharomycetes` → `saccharomycetes_odb10` (2137 BUSCOs)

Network errors fall back gracefully to keyword-based `infer_lineage()` — no exception raised.

## Agentic Lineage Routing

The `--organism` flag is the primary agentic bridge. The LLM agent passes a free-text organism description; the skill resolves it to a BUSCO flag using the `LINEAGE_ROUTING` keyword table:

| User organism hint | Resolved flag | Lineage dataset |
|---|---|---|
| "bacteria", "E. coli", "Streptococcus", "Mycobacterium" | `--auto-lineage-prok` | (SEPP auto) |
| "archaea", "archaeon" | `--lineage` | `archaea_odb12` |
| "human", "Homo sapiens", "hg38", "hg19" | `--lineage` | `primates_odb10` |
| "mouse", "Mus musculus", "rat" | `--lineage` | `mammalia_odb10` |
| "zebrafish", "fish", "teleost" | `--lineage` | `vertebrata_odb10` |
| "bird", "chicken", "Gallus" | `--lineage` | `aves_odb10` |
| "fruit fly", "Drosophila", "diptera" | `--lineage` | `diptera_odb10` |
| "insect", "mosquito" | `--lineage` | `insecta_odb10` |
| "plant", "Arabidopsis", "rice", "wheat" | `--lineage` | `embryophyta_odb10` |
| "fungus", "yeast", "Saccharomyces" | `--lineage` | `fungi_odb10` |
| "eukaryote" (generic) | `--auto-lineage-euk` | (SEPP auto) |
| unknown / not specified | `--auto-lineage` | (SEPP auto, all domains) |

## Algorithm / Methodology

1. BUSCO v6 searches input sequences against HMM profiles of single-copy orthologs from OrthoDB.
2. Each ortholog is classified: **Complete** (score and length within expected range) → **Single-copy (S)** or **Duplicated (D)**; **Fragmented (F)** (score within range, length below threshold); **Missing (M)** (no significant hit).
3. Completeness percentage = (C + F) / n × 100. C alone is the primary quality metric.
4. Auto-lineage uses SEPP placement of marker genes to identify the correct clade; requires SEPP exactly 4.5.5.
5. OrthoDB10 datasets cover eukaryotes; OrthoDB12 covers prokaryotes and archaea — do not mix suffixes.

## Example Queries

- "Check the completeness of my bacteria genome assembly"
- "Run BUSCO on my Drosophila transcriptome using diptera lineage"
- "What is the BUSCO score for this human genome assembly?"
- "Run BUSCO in proteins mode with the vertebrata lineage"
- "Show me a BUSCO demo with synthetic data"
- "BUSCO assessment with auto-lineage prokaryote mode on my E. coli assembly"

## Example Output

```markdown
# BUSCO Assessor Report

**Date**: 2026-04-23 10:00 UTC
**Mode**: genome (demo)
**Lineage**: bacteria_odb12
**Input**: demo_assembly.fna (5 sequences)

## Completeness Summary

| Metric | Count | Percentage |
|--------|-------|-----------|
| Complete (C) | 118 | 95.2% |
|   Single-copy (S) | 115 | 93.1% |
|   Duplicated (D) | 3 | 2.1% |
| Fragmented (F) | 3 | 2.3% |
| Missing (M) | 3 | 2.5% |
| Total searched (n) | 124 | — |

**Score string**: `C:95.2%[S:93.1%,D:2.1%],F:2.3%,M:2.5%,n:124`

## Interpretation

High completeness (95.2% C) indicates a near-complete assembly for this lineage.
Duplication rate of 2.1% is within expected range.

## Top Gene Results (first 10)
| BUSCO ID | Status | Sequence | Score | Length |
|----------|--------|----------|-------|--------|
| 1098at2  | Complete   | seq1 | 742.3 | 312 |
| 1099at2  | Complete   | seq1 | 698.1 | 287 |
| 1103at2  | Fragmented | seq2 | 341.2 |  98 |
| 1104at2  | Missing    | N/A  |   0.0 |   0 |

*ClawBio is a research and educational tool. It is not a medical device...*
```

## Output Structure

```
output_dir/
├── report.md                        # PRIMARY: completeness report
├── result.json                      # scores, lineage, mode, run parameters
├── busco_run/
│   ├── short_summary.txt            # BUSCO score summary (raw BUSCO format)
│   ├── short_summary.json           # Structured score summary
│   └── full_table.tsv               # Per-gene completeness table
└── reproducibility/
    ├── commands.sh                  # Exact replay command
    ├── environment.yml              # Pins busco=6.0.0, sepp=4.5.5
    └── checksums.sha256             # SHA-256 of all output files
```

## Dependencies

**Required (runtime; not needed for `--demo`)**

| Tool | Version | Purpose |
|------|---------|---------|
| `busco` | ≥6.0.0 | Core completeness analysis engine |
| `hmmer` | ≥3.1 | Profile HMM searches (installed with BUSCO) |
| `miniprot` | any | Eukaryote genome mode (default gene predictor) |
| `prodigal` | any | Prokaryote genome mode |
| `sepp` | **4.5.5 exactly** | Auto-lineage placement (v4.5.6 is broken) |
| `tblastn` | ≥2.10.1 | Transcriptome mode (v2.4–2.10.0 have CPU bugs) |

**Optional**

| Tool | Purpose |
|------|---------|
| `augustus` | Alternative eukaryote gene predictor (`--augustus` flag) |
| `metaeuk` | Alternative eukaryote gene predictor |

**Install (conda — recommended):**
```bash
conda create -n busco_env -c conda-forge -c bioconda busco=6.0.0 sepp=4.5.5
conda activate busco_env
```

## Gotchas

1. **SEPP version must be exactly 4.5.5.** SEPP v4.5.6 is incompatible with BUSCO auto-lineage files and produces wrong lineage assignments silently. Always pin `sepp=4.5.5` in environment.yml.

2. **Do NOT mix OrthoDB10 and OrthoDB12 lineage suffixes.** Eukaryote lineages use `_odb10`; prokaryote/archaea lineages use `_odb12`. Passing `bacteria_odb10` (non-existent) fails; passing `primates_odb12` (non-existent) fails. The lineage suffix must match the domain.

3. **BUSCO v6 changed the short_summary filename.** Depending on the BUSCO version and configuration, the file may be named `short_summary.txt` or `short_summary.specific.<lineage>.<run>.txt`. Always glob for both patterns — never hardcode the filename.

4. **Demo mode must never invoke the BUSCO binary.** `run_demo()` generates all output files synthetically in Python. Do not add BUSCO subprocess calls to the demo path; it must work in CI environments without any bioinformatics tools installed.

5. **Proteins mode with a nucleotide FASTA returns zero hits silently.** If `--mode proteins` is specified with a `.fna`/`.fa` file, BUSCO will complete successfully but report 0% completeness. The script emits a WARNING in this case; always use `.faa` (amino-acid FASTA) for proteins mode.

## Safety

- **Local-first**: All processing is local. No sequence data is uploaded to external services. Lineage datasets are downloaded from BUSCO servers only when the BUSCO binary is running and `--download_path` is specified.
- **No hallucinated scores**: The agent must NOT invent completeness percentages, lineage names, or gene counts. All numbers in the report derive from parsing BUSCO output or the synthetic demo constants.
- **Disclaimer**: Every generated `report.md` ends with: *"ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions."*

## Agent Boundary

- **Agent dispatches**: FASTA file path, `--mode`, `--organism` (free-text hint), optional explicit `--lineage` or `--auto-lineage*` flags.
- **Skill executes**: Lineage resolution, BUSCO command construction, subprocess management, output parsing, report writing.
- **Agent explains**: Results to the user, including what the completeness scores mean for their specific use case.
- **Agent must NOT**: Override the routing table with guessed lineage names, invent BUSCO score numbers, or run BUSCO commands manually outside this skill.

## Integration with Bio Orchestrator

**Trigger conditions for routing here:**
- User mentions "genome completeness", "BUSCO score", "assembly quality", or "assembly QC"
- User has produced a FASTA file from an assembler (Flye, SPAdes, Hifiasm, etc.)
- User asks "how complete is my assembly/transcriptome/protein set"

**Chaining partners:**

| Upstream | Handoff | Downstream |
|----------|---------|-----------|
| `seq-wrangler` | Assembled genome FASTA | `busco-assessor` |
| `busco-assessor` | `busco_run/` directory with `short_summary.txt` | `multiqc-reporter` (BUSCO module for multi-sample aggregation) |
| `busco-assessor` | `result.json` completeness scores | `profile-report` (unified genomic profile) |

**Output is chainable**: `result.json` is machine-readable JSON; `busco_run/short_summary.txt` is directly readable by MultiQC's BUSCO module.

## Maintenance

- **Review cadence**: Monthly or on new BUSCO major release.
- **Staleness signals**: Auto-lineage tests fail; OrthoDB download URLs change; new `_odb13` datasets released; SEPP constraint changes.
- **Update LINEAGE_ROUTING** when new OrthoDB versions introduce new clade-specific datasets or rename existing ones.
- **Deprecation path**: Move to `skills/_deprecated/busco-assessor/` if BUSCO v7 introduces breaking CLI changes that require a full rewrite.

## Citations

- Manni M. et al. (2021). BUSCO Update: Novel and Streamlined Workflows along with Broader and Deeper Phylogenetic Coverage. *Molecular Biology and Evolution*. https://doi.org/10.1093/molbev/msab199
- Simão F.A. et al. (2015). BUSCO: assessing genome assembly and annotation completeness with single-copy orthologs. *Bioinformatics*. https://doi.org/10.1093/bioinformatics/btv351
- OrthoDB v10/v12: https://www.orthodb.org/
- BUSCO GitLab: https://gitlab.com/ezlab/busco
- BUSCO User Guide v6: https://busco.ezlab.org/busco_userguide.html
