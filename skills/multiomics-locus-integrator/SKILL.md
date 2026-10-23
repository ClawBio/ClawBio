---
name: multiomics-locus-integrator
description: >-
  Integrate chromatin accessibility, histone marks, and 3D-contact (TAD/loop)
  evidence at a genomic locus into one scored regulatory-evidence table, with
  target genes assigned by proximity or loop contact.
license: MIT
metadata:
  version: "0.1.0"
  author: Your Name
  domain: genomics
  tags:
    - multiomics
    - chromatin
    - regulatory-genomics
    - enhancer
    - 3d-genome
    - locus-annotation
  inputs:
    - name: loci
      type: file
      format:
        - bed
      description: Candidate loci (chrom, start, end, id)
      required: true
    - name: gene_targets
      type: file
      format:
        - bed
      description: Candidate target gene coordinates with gene symbol and strand (chrom, start, end, gene_symbol, strand)
      required: true
    - name: accessibility_track
      type: file
      format:
        - bed
      description: ATAC/DNase peak calls (optional evidence layer)
      required: false
    - name: histone_track
      type: file
      format:
        - bed
      description: Histone ChIP-seq peak calls, e.g. H3K27ac/H3K4me1 (optional evidence layer)
      required: false
    - name: contact_evidence
      type: file
      format:
        - bedpe
      description: Loop/TAD anchor pairs, e.g. Hi-C or similar chromatin conformation capture loops (optional evidence layer)
      required: false
  outputs:
    - name: report
      type: file
      format:
        - md
      description: Per-locus evidence summary with confidence tiers
    - name: evidence_table
      type: file
      format:
        - tsv
      description: Machine-readable locus x layer evidence matrix with integrated score
    - name: result
      type: file
      format:
        - json
      description: Full machine-readable result, same content as the report
  dependencies:
    python: ">=3.11"
  demo_data:
    - path: demo_loci.bed
      description: 6 synthetic loci on a fictitious chr5 region, engineered to exercise every tier and gotcha
    - path: demo_atac.bed
      description: Synthetic ATAC peaks, limited genomic span (demonstrates not_assayed vs. absent)
    - path: demo_h3k27ac.bed
      description: Synthetic H3K27ac peaks
    - path: demo_loops.bedpe
      description: Synthetic loop anchors, including one that contacts a different gene than the nearest gene call
    - path: demo_targets.bed
      description: 4 synthetic candidate target genes
  endpoints:
    cli: python skills/multiomics-locus-integrator/multiomics_locus_integrator.py --loci {loci} --gene-targets {gene_targets} --output {output_dir}
  openclaw:
    requires:
      bins:
        - python3
    always: false
    emoji: "🧩"
    homepage: https://github.com/ClawBio/ClawBio
    os:
      - darwin
      - linux
    trigger_keywords:
      - multiomics integration
      - multi-omic
      - ATAC H3K27ac integration
      - chromatin loop target gene
      - enhancer gene evidence
      - integrate chromatin and variant evidence
      - regulatory evidence scoring
      - TAD loop locus integration
---

# 🧩 Multiomics Locus Integrator

You are **Multiomics Locus Integrator**, a specialised ClawBio agent for regulatory genomics. Your role is to join already-called chromatin and 3D-contact evidence at a set of candidate loci and produce one scored, auditable evidence table per locus.

## Trigger

**Fire this skill when the user says any of:**
- "integrate ATAC and H3K27ac at these loci"
- "which of these SNPs/peaks have chromatin support"
- "combine chromatin accessibility, histone marks, and loop data for this region"
- "score these loci for regulatory evidence"
- "does this SNP contact a different gene than its nearest gene"
- "multiomics integration" / "multi-omic evidence"

**Do NOT fire when:**
- The user wants peaks *called* from raw BAM/bigWig data — that's upstream of this skill (ATAC/ChIP-seq processing, chromatin conformation capture (Hi-C-style) loop calling). This skill only joins and scores already-called evidence.
- The user wants literature/database evidence aggregation across omics types — route to `omics-target-evidence-mapper` instead.
- The user wants statistical fine-mapping of a GWAS signal — route to `fine-mapping`.

## Why This Exists

- **Without it**: every lab re-derives the same join logic by hand (peak overlap → loop anchor overlap → nearest/contacted gene → score) for every locus list, in a spreadsheet, from scratch.
- **With it**: one deterministic, auditable scorer that takes whatever layers are available (even one of three) and produces a consistent per-locus confidence tier, with the evidence basis always shown.
- **Why ClawBio**: no existing skill composes across chromatin-layer evidence at the locus level — this is the missing glue between `gi-chromatin`/`gi-enhancer` (prediction) and `locuscompare-region-render`/`eqtl-catalogue-region-fetch` (region-level rendering and lookup).

## Core Capabilities

1. **Layer overlap scoring**: weighted integration of accessibility, histone, and contact evidence, using only the layers actually supplied.
2. **Not-assayed vs. absent distinction**: a layer that never covered a locus's region is reported as `not_assayed`, never silently conflated with a true negative.
3. **Dual target-gene assignment**: nearest gene by proximity, and loop-contacted gene if contact evidence is present, with disagreements between the two explicitly flagged rather than silently resolved.

## Scope

**One locus, one row.** This skill joins and scores pre-called evidence; it does not call peaks, call loops, or process raw sequencing data.

## Input Formats

| Format | Extension | Required Fields | Example |
|--------|-----------|-----------------|---------|
| Loci BED | `.bed` | chrom, start, end, id | `demo_loci.bed` |
| Gene targets BED | `.bed` | chrom, start, end, gene_symbol, strand (+/-) | `demo_targets.bed` |
| Accessibility/histone BED | `.bed` | chrom, start, end, peak_id | `demo_atac.bed` |
| Contact BEDPE | `.bedpe` | chrom1, start1, end1, chrom2, start2, end2, loop_id | `demo_loops.bedpe` |

## Workflow

1. **Validate**: require `--loci` and `--gene-targets`; require at least one evidence layer (accessibility, histone, or contact).
2. **Process**: for each locus, determine per-layer `checked` / `supported` / `not_assayed` status (a layer is `not_assayed` for a locus outside that layer's observed chromosomal span), compute the weighted score, assign a tier, and resolve nearest-gene and loop-contact target genes.
3. **Generate**: write `tables/evidence_table.tsv` (one row per locus) and `result.json` (full structured result).
4. **Report**: write `report.md` with the summary table and the reproducibility bundle.

## CLI Reference

```bash
# Standard usage
python skills/multiomics-locus-integrator/multiomics_locus_integrator.py \
  --loci loci.bed --gene-targets genes.bed \
  --accessibility-track atac.bed --histone-track h3k27ac.bed --contact-evidence loops.bedpe \
  --output report_dir

# Demo mode (synthetic data, no user files needed)
python skills/multiomics-locus-integrator/multiomics_locus_integrator.py --demo --output /tmp/demo
```

## Demo

```bash
python skills/multiomics-locus-integrator/multiomics_locus_integrator.py --demo --output /tmp/demo
```

Expected output: a 6-locus report on a fictitious chr5 region demonstrating all four tiers (`high`, `partial`, `unsupported`, `not_assayed`) and one nearest-gene-vs-contact-gene disagreement (locus `L4_disagree`).

All demo coordinates are synthetic (fictitious `chr5` positions, invented gene symbols `GENE_A`–`GENE_D`) — not derived from any real dataset.

## Algorithm / Methodology

1. For each locus, overlap against each supplied layer (BED/BEDPE interval intersection, half-open coordinates).
2. A layer is `not_assayed` for a locus if the locus falls outside that layer's observed min/max span on that chromosome (the layer file never covered this region) — reported separately from `unsupported` (checked, no overlap found).
3. Score = (sum of weights of *supported* layers) / (sum of weights of *checked* layers) — only layers actually supplied and spanning the locus count toward the denominator, so a locus is never penalised for a layer nobody ran.
4. Tier: `high` (score ≥ 0.75 **and** ≥ 2 distinct supported layer types — a single heavily-weighted layer alone cannot reach `high`), `partial` (any support), `unsupported` (checked, zero support), `not_assayed` (no layer covered this locus at all).
5. Target gene: nearest gene by proximity to a strand-aware TSS window (`TSS-1000bp` to `TSS+250bp` — standard promoter-proximal convention; `+`-strand TSS = BED start, `-`-strand TSS = BED end), **and separately** any gene contacted via a loop whose other anchor overlaps the locus. If they disagree, both are reported — contact evidence does not silently override proximity.

**Key thresholds**:
- Layer weights: accessibility 1.0, histone 1.0, contact 2.0 (contact evidence weighted higher — harder to obtain spuriously than a single peak call). Not currently user-overridable; a `--weights` flag is a natural v0.2 addition.
- High-tier threshold: score ≥ 0.75 and ≥ 2 distinct supported layer types.

## Example Output

```markdown
# Multiomics Locus Integrator Report

**Loci evaluated**: 6
**Layers supplied**: accessibility, contact, histone
**Tier counts**: {'high': 2, 'partial': 2, 'unsupported': 1, 'not_assayed': 1}
**Nearest-gene vs. contact-gene disagreements**: 1

| Locus | Coordinates | Score | Tier | Checked | Supported | Not assayed | Nearest gene | Contact gene | Disagreement |
|---|---|---:|---|---|---|---|---|---|---|
| L1_agree_high | chr5:1000100-1000101 | 1.0 | high | accessibility,contact,histone | accessibility,contact,histone | - | GENE_A | GENE_A | |
| L4_disagree | chr5:1030100-1030101 | 0.667 | partial | contact,histone | contact | accessibility | GENE_B | GENE_C | YES |
| L6_not_assayed | chr5:1050100-1050101 |  | not_assayed | - | - | accessibility,contact,histone | GENE_D |  | |

ClawBio is a research tool. Locus tiers are a starting point for manual review, not a substitute for experimental validation.
```

## Output Structure

```
output_directory/
├── report.md
├── result.json
├── tables/
│   └── evidence_table.tsv
└── reproducibility/
    ├── commands.sh
    ├── environment.yml
    └── checksums.sha256
```

## Dependencies

**Required**: none beyond the Python standard library and `clawbio.common.reproducibility`.

## Gotchas

- **A gene with no strand is a BED record with no usable TSS.** The upstream/downstream sense of "1kb upstream, 250bp downstream" only means something once strand picks which end of the interval is the start — `gene-targets` therefore requires a 5th strand column (`+`/`-`); earlier revisions silently measured distance to whichever gene-body edge was nearer, which is the wrong edge for every minus-strand gene (PR #678 review).
- **A `not_assayed` locus is not evidence of absence.** If a track was never generated over a region, every locus in it will otherwise read identically to a true negative unless this distinction is computed explicitly — this skill derives `not_assayed` from each layer's own observed chromosomal span, not from a hardcoded genome-wide assumption.
- **Enumerate every evidence file's actual coverage before trusting a `tier`.** A loop/TAD file can cover only part of a chromosome; the per-chromosome span check is what catches this, but it is only as good as the supplied file — a track with one stray far-off interval will inflate its apparent span and hide a real `not_assayed` region.
- **Nearest-gene and contact-gene disagreement must be surfaced, not silently resolved.** Early drafts defaulted to "contact wins" — this is wrong: both are real candidate answers with different evidence bases, and collapsing them hides exactly the kind of judgment call a reviewer needs to see (demonstrated by `L4_disagree` in the demo data).
- **A single heavily-weighted layer should not alone earn `high`.** Requiring ≥2 distinct supported layer types for the `high` tier (not just a score threshold) prevents one noisy, high-weight loop call from single-handedly producing a false high-confidence result.

## Safety

- **Local-first**: no data upload; all inputs are local BED/BEDPE files.
- **No hallucinated science**: layer weights and tier thresholds are declared, fixed parameters in this file, not inferred or tuned per-run.
- **Audit trail**: every run writes a reproducibility bundle (`commands.sh`, `environment.yml`, `checksums.sha256`).

## Agent Boundary

The agent (LLM) explains tiers and disagreements to the user. The skill (Python) computes every overlap, weight, and threshold deterministically; the agent must not override a tier or invent a target gene not present in the evidence table.

## Integration with Bio Orchestrator

**Chaining partners**: `gwas-catalog-region-fetch` / `eqtl-catalogue-region-fetch` (supply candidate loci), `gi-enhancer` / `gi-chromatin` (could supply a predicted layer for loci lacking real ChIP-seq/ATAC data), `locuscompare-region-render` (renders the resulting evidence table alongside LD/fine-mapping).

## Maintenance

- **v0.2 candidates**: variant/genotype evidence as a fifth layer; user-overridable `--weights`; bigWig signal-value support (currently BED/BEDPE presence-only) via optional `pyBigWig`.
- **Staleness signals**: none — this skill has no external database dependency.

## Citations

- Fulco, C.P. et al. (2019). Activity-by-contact model of enhancer-promoter regulation from thousands of CRISPR perturbations. *Nature Genetics*, 51, 1664–1669. (Conceptual basis for contact-weighted locus-gene assignment.)
