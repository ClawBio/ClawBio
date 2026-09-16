---
name: phylogenetics-builder
description: End-to-end ML phylogenetic tree inference — MSA, trimming, ModelFinder, IQ-TREE2/RAxML-NG.
license: MIT
metadata:
  openclaw:
    requires:
      bins:
      - python3
    always: false
    emoji: 🌳
    homepage: https://github.com/ClawBio/ClawBio
    os:
    - darwin
    - linux
    install:
    - kind: conda
      package: bioconda::iqtree
    - kind: conda
      package: bioconda::raxml-ng
    - kind: conda
      package: bioconda::mafft
    - kind: conda
      package: bioconda::muscle
    - kind: conda
      package: bioconda::trimal
    trigger_keywords:
    - phylogeny
    - phylogenetic tree
    - iqtree
    - raxml
    - maximum likelihood tree
    - mafft alignment
    - build tree from sequences
    - model selection
    - ModelFinder
    - bootstrap support
    - evolutionary tree
    - molecular phylogeny
  author: ClawBio
  version: 0.2.0
  demo_data:
  - path: demo_alignment.fasta
    description: "Synthetic 12-taxon primate alignment (500 bp, pre-aligned)"
  dependencies:
    python: '>=3.10'
    packages:
    - pandas>=2.0
    - biopython>=1.80
    - matplotlib>=3.5
    optional:
    - ete3>=3.1  # midpoint rooting
  domain: genomics
  endpoints:
    cli: python skills/phylogenetics-builder/phylogenetics_builder.py --input {input_file} --output {output_dir}
  inputs:
  - name: input_file
    type: file
    format:
    - fasta
    - fa
    - aln
    description: DNA or protein sequences (unaligned for full pipeline, aligned with --aligned flag)
    required: true
  outputs:
  - name: report
    type: file
    format:
    - md
    description: Full analysis report with pipeline summary and branch table
  - name: result
    type: file
    format:
    - json
    description: Machine-readable results (ClawBio output contract)
  - name: phylo_tree
    type: file
    format:
    - nwk
    description: Newick format tree with bootstrap support values
  - name: alignment
    type: file
    format:
    - fasta
    description: Alignments produced by this run, alignment/aligned.fasta after the MSA
      stage and alignment/trimmed.fasta after trimAl; omitted for a stage that did not run
---

# 🌳 Phylogenetics Builder

You are **Phylogenetics Builder**, a ClawBio agent for end-to-end maximum-likelihood phylogenetic tree inference. You run the full pipeline: MSA → trimming → model selection → tree inference → rooting → visualisation.

## Why This Exists

Maximum-likelihood phylogenetics requires correctly chaining at least five external tools (aligner → trimmer → model selector → tree engine → visualiser), each with non-obvious CLI quirks — conflicting flags between MUSCLE v3/v5, model-name format incompatibility between IQ-TREE and RAxML-NG, and different bootstrap confidence thresholds (UFBoot ≥ 95 vs standard ≥ 70). This skill encapsulates the correct invocation for all supported tools and handles their output differences automatically.

## Trigger

**Fire when the user says:**
- "build a phylogenetic tree from these sequences"
- "run phylogeny analysis" / "infer evolutionary tree"
- "run IQ-TREE on my FASTA" / "use RAxML"
- "what substitution model should I use?" (with a FASTA file present)
- "align and build a tree" / "MSA then tree"
- "bootstrap support values" / "UFBoot replicates"
- "midpoint root the tree" / "root with outgroup"
- "maximum likelihood tree from my sequences"

**Do NOT fire when:**
- The user wants k-mer/distance trees only → use `fastreer` instead
- The user wants variant-based trees from VCF → use `fastreer` instead
- The user wants protein structure prediction → use `struct-predictor`
- The user needs alignment only (no tree) → recommend mafft/muscle standalone

## Scope

**One skill, one task.** This skill infers a maximum-likelihood phylogenetic tree from DNA or protein sequences. It does not annotate variants, predict structures, or perform downstream comparative genomics. Each post-tree task chains to another skill.

Supported pipeline stages:

- **6 MSA algorithms**: mafft (default), muscle, clustalw, kalign, tcoffee, prank
- **Alignment trimming**: trimAl `-automated1` (removes gapped columns)
- **Automatic model selection**: IQ-TREE2 ModelFinder (`-m MFP`), BIC-selected
- **Two inference engines**: IQ-TREE2 (default) and RAxML-NG
- **Three bootstrap modes**: UFBoot (1 000 reps, threshold ≥ 95), standard Felsenstein (100 reps, threshold ≥ 70), triple support (UFBoot + aLRT + aBayes)
- **Post-inference rooting**: outgroup or midpoint (ETE3 primary, Bio.Phylo `root_at_midpoint` fallback)
- **Visualisation**: proportional phylogram via Bio.Phylo + matplotlib
- **Reproducibility bundle**: exact CLI command, Conda environment definition, SHA-256 checksums
- **Offline demo**: pre-computed 12-taxon primate tree — never refuses when binaries are absent

## Workflow

1. **Validate input** — parse FASTA, check ≥3 sequences, check alignment if `--aligned` is set.
2. **MSA** (skip if `--aligned`) — run the chosen aligner; default is `mafft --auto` for speed/quality balance. Alternative aligners: `muscle`, `clustalw`, `kalign`, `tcoffee`, `prank`.
3. **Trim** (skip with `--no-trim`) — run `trimal -automated1`. This removes gapped columns that add noise without information. Skip for protein alignments where all columns are informative.
4. **Model selection** (skip if `--model` provided) — run `iqtree2 -m MFP`. Parse `Best-fit model according to BIC:` from the `.iqtree` log. The selected model is passed directly to tree inference.
5. **Tree inference** — choose engine:
   - **IQ-TREE2** (default): `iqtree2 -s aln.fa -m MODEL --prefix prefix [-bb 1000] [-b 100] [-bb 1000 -alrt 1000 -abayes]`
   - **RAxML-NG**: first `raxml-ng --check` (validates input), then `raxml-ng --all --msa aln.fa --model MODEL --bs-trees N`
6. **Rooting** (optional):
   - Outgroup: pass `--outgroup TAXON` to tree builder (`-o` in IQ-TREE, `--outgroup` in RAxML-NG)
   - Midpoint: use ETE3 `t.get_midpoint_outgroup(); t.set_outgroup(midpoint)` post-inference
7. **Parse & render** — extract branch lengths and support values from Newick; draw proportional phylogram with Bio.Phylo + matplotlib.
8. **Report** — write `report.md`, `result.json` (ClawBio contract), `phylo_tree.nwk`, `alignment/aligned.fasta` and `alignment/trimmed.fasta` when those stages ran, `figures/phylogram.png`, `tables/branch_support.csv`, `reproducibility/`.

**Demo fallback (MANDATORY):** If no binaries are installed, skip to pre-computed tree from `examples/demo_tree.nwk`. Always show a result, never refuse.

**Freedom level per step:**
- Steps 1–6 are prescriptive (exact CLI flags matter).
- Report narrative can be written flexibly by the agent.

## CLI Reference

```bash
# Full pipeline: unaligned → MSA → trim → ModelFinder → IQ-TREE2
python skills/phylogenetics-builder/phylogenetics_builder.py \
  --input sequences.fasta --output /tmp/phylo

# Pre-aligned input (skip MSA)
python skills/phylogenetics-builder/phylogenetics_builder.py \
  --input aligned.fasta --output /tmp/phylo --aligned

# Choose MSA algorithm (mafft default)
python skills/phylogenetics-builder/phylogenetics_builder.py \
  --input sequences.fasta --output /tmp/phylo \
  --aligner muscle

# Standard bootstrap instead of UFBoot
python skills/phylogenetics-builder/phylogenetics_builder.py \
  --input aligned.fasta --output /tmp/phylo --aligned \
  --bootstrap standard

# Triple support: UFBoot + aLRT + aBayes
python skills/phylogenetics-builder/phylogenetics_builder.py \
  --input aligned.fasta --output /tmp/phylo --aligned \
  --bootstrap all

# Root by outgroup
python skills/phylogenetics-builder/phylogenetics_builder.py \
  --input aligned.fasta --output /tmp/phylo --aligned \
  --outgroup Mus_musculus,Rattus_norvegicus

# Midpoint rooting (requires ETE3)
python skills/phylogenetics-builder/phylogenetics_builder.py \
  --input aligned.fasta --output /tmp/phylo --aligned \
  --root midpoint

# Use RAxML-NG engine
python skills/phylogenetics-builder/phylogenetics_builder.py \
  --input aligned.fasta --output /tmp/phylo --aligned \
  --engine raxml-ng

# Skip trimming
python skills/phylogenetics-builder/phylogenetics_builder.py \
  --input aligned.fasta --output /tmp/phylo --aligned --no-trim

# Provide model explicitly (skip ModelFinder)
python skills/phylogenetics-builder/phylogenetics_builder.py \
  --input aligned.fasta --output /tmp/phylo --aligned \
  --model GTR+F+G4

# Demo mode (works offline, no binaries needed)
python skills/phylogenetics-builder/phylogenetics_builder.py \
  --demo --output /tmp/phylo_demo
```

### All flags

| Flag | Default | Description |
|------|---------|-------------|
| `--input FILE` | — | Input FASTA (unaligned or aligned) |
| `--output DIR` | — | Output directory |
| `--demo` | off | Run with built-in 12-taxon primate data |
| `--aligned` | off | Input is already aligned — skip MSA |
| `--aligner` | mafft | MSA algorithm: mafft / muscle / clustalw / kalign / tcoffee / prank |
| `--engine` | iqtree2 | Tree engine: iqtree2 / raxml-ng |
| `--model MODEL` | auto | Skip ModelFinder; use this substitution model |
| `--bootstrap` | ufboot | Bootstrap: ufboot / standard / all |
| `--outgroup TAXA` | — | Comma-separated outgroup taxon name(s) |
| `--root midpoint` | — | Midpoint rooting via ETE3 (post-inference) |
| `--no-trim` | off | Skip trimAl trimming |
| `--threads N` | 2 | CPU threads for tree inference |
| `--seed N` | 42 | Random seed for reproducibility |

## MSA Algorithm Guide

| Aligner | Speed (10 seqs) | Speed (250 seqs) | Recommendation |
|---------|-----------------|------------------|----------------|
| mafft | 4.4 s | 42 s | **Default** — best speed/quality balance |
| kalign | 0.5 s | 8 s | Fastest for large datasets (>100 seqs) |
| muscle | 5 s | 30 min | Good for protein alignments |
| clustalw | 5.6 s | 49 min | Legacy; avoid for large datasets |
| tcoffee | slow | very slow | Most accurate; use for ≤20 sequences |
| prank | slow | very slow | Codon-aware; use with `-codon` for coding DNA |

_Benchmarks on SUP35 gene dataset from NGS Handbook._

## Bootstrap Guide

| Mode | Flag | Speed | Use when |
|------|------|-------|----------|
| `ufboot` | `-bb 1000` | ~3 sec | Default; fast and reliable (threshold: ≥95) |
| `standard` | `-b 100` | ~3 min | Publication standard; slower Felsenstein bootstrap |
| `all` | `-bb 1000 -alrt 1000 -abayes` | ~5 sec | Need triple validation; parse with `/` delimiter |

Triple support labels format: `{alrt}/{abayes}/{ufb}` — thresholds: alrt > 70, abayes > 0.7, ufb > 95.

## Example Output

```markdown
# Phylogenetics Builder Report

### Pipeline Summary

| Parameter | Value |
|-----------|-------|
| Input | `sequences.fasta` |
| Taxa | 12 |
| Aligner | mafft |
| Trimming | trimAl -automated1 |
| Substitution model | `TIM3+F+G4` |
| Tree engine | iqtree2 |
| Bootstrap | UFBoot (1 000 replicates) |
| Rooting | unrooted |

### Pipeline Steps
- `msa:mafft`
- `trim:trimal`
- `modelfinder:TIM3+F+G4`
- `tree:iqtree2:ufboot`

### Branch Lengths & Support Values
| Node / Taxon | Branch Length | Support |
|:-------------|:-------------:|:-------:|
| Homo_sapiens | 0.01000 | 100 |
| Pan_troglodytes | 0.00800 | 98 |
...
```

## Output Structure

```
output_directory/
├── report.md                    # Primary markdown report with pipeline summary
├── result.json                  # Machine-readable ClawBio output contract
├── phylo_tree.nwk               # Newick format tree with bootstrap support
├── alignment/                  # aligned.fasta (MSA) and trimmed.fasta (trimAl), when produced
├── figures/
│   └── phylogram.png            # Proportional phylogram (matplotlib)
├── tables/
│   └── branch_support.csv       # Per-node branch lengths and support values
└── reproducibility/
    ├── commands.sh              # Exact CLI command used
    ├── environment.yml          # Conda environment definition
    └── checksums.sha256         # SHA-256 checksums of all outputs
```

## Gotchas

- **Prank writes to `{prefix}.best.fas`, not `{prefix}`**. If using prank as aligner, the skill auto-renames this file. If you call prank manually, remember to look for the `.best.fas` suffix.
- **ModelFinder adds `+F` to models; RAxML-NG rejects it**. `TIM3+F+G4` from IQ-TREE ModelFinder must be stripped to `TIM3+G4` for RAxML-NG. The skill handles this automatically via `adapt_model_for_engine()`. If you pass `--model` manually with `--engine raxml-ng`, omit the `+F`.
- **UFBoot threshold is 95, not 70**. A common mistake is applying the standard bootstrap threshold (70) to UFBoot values. UFBoot support of 70 is NOT reliable. Use ≥ 95 as the cutoff for UFBoot.
- **Triple support labels contain `/`**. With `--bootstrap all`, node labels encode `alrt/abayes/ufb` (e.g. `80.5/0.85/97`). Standard Newick readers interpret the whole string as a confidence value. Use `strsplit(label, "/")` in R or split by `/` in Python.
- **trimAl on protein alignments may over-trim**. For small protein alignments (<20 sequences, >200 aa), use `--no-trim` or use `-nogaps` strategy instead of `-automated1`, which can remove too many columns.
- **IQ-TREE `-T AUTO` can block tests**. Always specify explicit thread count (`-T 2`) in automated/test contexts to avoid IQ-TREE hanging on thread detection.
- **Pre-aligned input still needs equal-length sequences**. If you use `--aligned`, the skill validates that all sequences are the same length. Gaps (`-`) are allowed; just ensure no sequences were accidentally truncated.

## Safety

- **All processing is local** — sequences never leave the machine.
- **Disclaimer**: every report includes the ClawBio medical disclaimer.
- **No hallucinated models or distances** — model and bootstrap values come directly from tool output.
- **Warn before overwriting** — script warns if output directory is non-empty.

## Agent Boundary

The agent (LLM) dispatches to this skill and explains the results. The skill (Python script) executes all computation. The agent must NOT invent substitution model names, bootstrap values, or branch lengths.

## Integration with Bio Orchestrator

Route to this skill when the query matches any `trigger_keywords` or the intent is maximum-likelihood tree inference. The orchestrator passes the FASTA path and any user-specified flags; the skill owns all tool decisions internally.

After the run, read `result.json`:
- `chat_summary_lines` — surface to the user verbatim.
- `preferred_artifacts` — open the figure and tree file for the user.
- `run_mode == "demo-fallback"` — surface `contract_alerts[0]` to prompt IQ-TREE2 installation.
- `workflow_state == "completed"` — no retry needed.

Do **not** pass raw tool flags from the user directly to the CLI without validation; use the documented `--flag` surface only.

## Chaining Partners

| Skill | When to chain |
|-------|--------------|
| `fastreer` | User wants a fast k-mer distance tree without full MSA |
| `variant-annotation` | Annotate variants found in sequences before building tree |
| `genome-compare` | Compare multiple genomes before phylogenetic inference |
| `profile-report` | Add evolutionary context to a patient profile |
| `claw-ancestry-pca` | Population structure analysis complements phylogenetics |

## Dependencies

| Dependency | Version | Required | Purpose |
|------------|---------|----------|---------|
| python | ≥ 3.10 | yes | Runtime |
| biopython | ≥ 1.80 | yes | Newick I/O, `root_at_midpoint`, visualisation |
| matplotlib | ≥ 3.5 | yes | Phylogram rendering |
| pandas | ≥ 2.0 | yes | Branch support CSV export |
| iqtree2 | ≥ 2.0 | recommended | ModelFinder + default tree engine |
| raxml-ng | any | optional | Alternative tree engine |
| mafft | any | optional | Default MSA aligner |
| muscle | ≥ 5.0 | optional | Alternative MSA aligner (v5 `-align`/`-output` syntax) |
| trimal | any | optional | Alignment column trimming |
| ete3 | ≥ 3.1 | optional | Midpoint rooting (Bio.Phylo fallback if absent) |
| clustalw | any | optional | Legacy MSA aligner |
| kalign | ≥ 3 | optional | Fast MSA for large datasets |
| t_coffee | any | optional | High-accuracy MSA for ≤ 20 sequences |
| prank | any | optional | Codon-aware MSA |

Install all bioinformatics binaries:
```bash
conda install -c bioconda iqtree raxml-ng mafft muscle trimal clustalw kalign3 t_coffee prank
```

## Maintenance

- **Review cadence**: bi-annually, or when IQ-TREE2 / RAxML-NG release major versions.
- **Staleness signals**: IQ-TREE2 `-m MFP` syntax changes; RAxML-NG `--all` flag renamed.
- **Deprecation criteria**: superseded by a unified tree inference + annotation tool.

## Citations

- [IQ-TREE 2](https://doi.org/10.1093/molbev/msaa015) — Minh et al., 2020
- [RAxML-NG](https://doi.org/10.1093/bioinformatics/btz305) — Kozlov et al., 2019
- [ModelTest-NG](https://doi.org/10.1093/molbev/msz189) — Darriba et al., 2020
- [trimAl](https://doi.org/10.1093/bioinformatics/btp348) — Capella-Gutiérrez et al., 2009
- [MAFFT](https://doi.org/10.1093/molbev/mst010) — Katoh & Standley, 2013
- [MUSCLE v5](https://doi.org/10.1101/2021.06.20.449169) — Edgar, 2021 (bioRxiv); [MUSCLE v3](https://doi.org/10.1093/nar/gkh340) — Edgar, 2004
- [T-Coffee](https://doi.org/10.1006/jmbi.2000.4042) — Notredame et al., 2000
- [ClustalW](https://doi.org/10.1093/nar/22.22.4673) — Thompson et al., 1994
- [KAlign](https://doi.org/10.1186/1471-2105-6-298) — Lassmann & Sonnhammer, 2005
- [PRANK](https://doi.org/10.1126/science.1158395) — Löytynoja & Goldman, 2008
- [ETE3](https://doi.org/10.1093/molbev/msw046) — Huerta-Cepas et al., 2016
- [Bio.Phylo](https://doi.org/10.1186/1471-2105-13-209) — Talevich et al., 2012
