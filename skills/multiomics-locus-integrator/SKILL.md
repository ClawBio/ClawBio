---
name: multiomics-locus-integrator
description: >-
  Integrate chromatin accessibility, histone marks, and 3D-contact (TAD/loop)
  evidence at a genomic locus into one scored regulatory-evidence table, with
  target genes assigned by proximity or loop contact.
license: MIT
metadata:
  version: "0.1.0"
  author: Muhammad Adeel
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
    - name: accessibility_coverage_bed
      type: file
      format:
        - bed
      description: Real assayed footprint for accessibility_track (e.g. a capture panel's target regions); overrides the default chromosome-wide-span not_assayed heuristic for this layer with a tighter, per-interval gap check
      required: false
    - name: histone_coverage_bed
      type: file
      format:
        - bed
      description: Same as accessibility_coverage_bed, for histone_track
      required: false
    - name: contact_evidence
      type: file
      format:
        - bedpe
      description: Loop/TAD anchor pairs, e.g. Hi-C or similar chromatin conformation capture loops (optional evidence layer)
      required: false
    - name: contact_coverage_bed
      type: file
      format:
        - bed
      description: Real assayed footprint for contact_evidence (e.g. a Capture-C/HiChIP capture region); overrides the default per-loop-anchor-span not_assayed heuristic with a tighter, per-interval gap check
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
    cli: python skills/multiomics-locus-integrator/multiomics_locus_integrator.py --loci {loci} --gene-targets {gene_targets} --accessibility-track {accessibility_track} --output {output_dir}
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
| Gene targets BED | `.bed` | chrom, start, end, gene_symbol, strand (+/-) **or** standard BED6 (chrom, start, end, name, score, strand) | `demo_targets.bed` |
| Accessibility/histone BED | `.bed` | chrom, start, end, peak_id | `demo_atac.bed` |
| Coverage BED (optional) | `.bed` | chrom, start, end, region_id | n/a -- user-supplied capture/target regions |
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

# With a targeted/capture track, where a gap really does mean "never
# measured" rather than "no signal called" -- see Gotchas below
python skills/multiomics-locus-integrator/multiomics_locus_integrator.py \
  --loci loci.bed --gene-targets genes.bed \
  --accessibility-track atac.bed --accessibility-coverage-bed panel_targets.bed \
  --output report_dir

# Same idea for contact evidence: --contact-coverage-bed states the real
# capture/assayed footprint of the 3D-contact assay, so a gap inside it
# reads as checked-and-unsupported rather than not_assayed
python skills/multiomics-locus-integrator/multiomics_locus_integrator.py \
  --loci loci.bed --gene-targets genes.bed \
  --contact-evidence loops.bedpe --contact-coverage-bed capture_region.bed \
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
2. A layer is `not_assayed` for a locus if the locus falls outside that layer's observed min/max span on that chromosome (the layer file never covered this region) — reported separately from `unsupported` (checked, no overlap found). If an explicit `--accessibility-coverage-bed` / `--histone-coverage-bed` is supplied for that layer, its own merged, per-interval blocks are used instead of the track's chromosome-wide span (see Gotchas).
3. Score = (sum of weights of *supported* layers) / (sum of weights of *checked* layers) — only layers actually supplied and spanning the locus count toward the denominator, so a locus is never penalised for a layer nobody ran.
4. Tier: `high` (score ≥ 0.75 **and** ≥ 2 distinct supported layer types — a single heavily-weighted layer alone cannot reach `high`), `partial` (any support), `unsupported` (checked, zero support), `not_assayed` (no layer covered this locus at all).
5. Target gene: nearest gene by **true distance to its strand-aware TSS point** (`+`-strand TSS = BED start, `-`-strand TSS = BED end — not to the TSS window's edge, which could occasionally rank the wrong gene as nearest), **and separately** any gene contacted via a loop whose other anchor overlaps the locus (checked against the TSS window: `TSS-1000bp` to `TSS+250bp`, standard promoter-proximal convention). If they disagree, both are reported — contact evidence does not silently override proximity.
6. Chromosome names are compared with a `chr`-prefix naming difference normalized out (`"5"` and `"chr5"` are treated as the same chromosome in every overlap/span/nearest-gene comparison) — but never normalized in any output field, which always reports the chromosome exactly as its own source file spelled it. This does **not** detect or correct a genome-build mismatch (e.g. hg19 vs. hg38): that requires reference data this skill does not have, and two builds can share identical `chr5`-vs-`chr5` naming while meaning different sequences. Verify all inputs share one assembly yourself.

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

## Layer coverage (observed span per supplied track)

- **accessibility**: chr5:1000050-1021100
- **contact**: chr5:999950-1041100
- **histone**: chr5:1000050-1040200

| Locus | Coordinates | Score | Tier | Checked | Supported | Not assayed | Nearest gene | Distant | Contact gene | Disagreement |
|---|---|---:|---|---|---|---|---|---|---|---|
| L1_agree_high | chr5:1000100-1000101 | 1.0 | high | accessibility,contact,histone | accessibility,contact,histone | - | GENE_A | | GENE_A | |
| L4_disagree | chr5:1030100-1030101 | 0.667 | partial | contact,histone | contact | accessibility | GENE_B | | GENE_C | YES |
| L6_not_assayed | chr5:1050100-1050101 |  | not_assayed | - | - | accessibility,contact,histone | GENE_D | | | |

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

- **A gene with no strand is a BED record with no usable TSS.** The upstream/downstream sense of "1kb upstream, 250bp downstream" only means something once strand picks which end of the interval is the start — `gene-targets` therefore requires a 5th strand column (`+`/`-`); earlier revisions silently measured distance to whichever gene-body edge was nearer, which is the wrong edge for every minus-strand gene.
- **A `not_assayed` locus is not evidence of absence.** If a track was never generated over a region, every locus in it will otherwise read identically to a true negative unless this distinction is computed explicitly — this skill derives `not_assayed` from each layer's own observed chromosomal span, not from a hardcoded genome-wide assumption.
- **A nearest gene is only as meaningful as the gene list it came from.** With a sparse `--gene-targets` file, "nearest" can mean "the only one supplied, megabases away" — `target_gene_nearest_is_distant` flags any pick past `NEAREST_GENE_FAR_THRESHOLD_BP` (1Mb) rather than silently reporting a distant gene as if it were a confident local call. The gene and its real distance are still returned; nothing is hidden, only flagged.
- **Report the layer coverage, not just the per-locus `not_assayed` flag.** A track's observed span can be stretched by a single far-off interval, which would otherwise mask a real coverage gap elsewhere — the `## Layer coverage` section in `report.md` makes each track's actual span visible so a reader can judge this directly rather than inferring it from which loci happen to read `not_assayed`.
- **Enumerate every evidence file's actual coverage before trusting a `tier`.** A loop/TAD file can cover only part of a chromosome; the per-chromosome span check is what catches this, but it is only as good as the supplied file — a track with one stray far-off interval will inflate its apparent span and hide a real `not_assayed` region.
- **A chromosome-wide span is a deliberate choice, not an oversight -- and `--coverage-bed` is the opt-out.** `not_assayed` treats the gap between two peaks on the same chromosome as "assayed, no signal called" (a true negative) by default, matching a whole-genome assay where peaks are sparse by design -- not as missing data. For a targeted or capture-based track, where a gap really does mean "never measured," supply `--accessibility-coverage-bed` / `--histone-coverage-bed` / `--contact-coverage-bed` pointing at the actual assayed regions (e.g. the panel's target BED, or a Capture-C/HiChIP capture region): that layer's `not_assayed` is then decided from the coverage file's own merged, per-interval blocks instead of the track's span (or, for contact, instead of each loop's own anchor span). This is opt-in per layer and never changes the default for a layer that does not supply one; all three layers' overrides are independent of each other.
- **A `chr`-prefix naming difference is normalized; a genome-build difference is not, and cannot be.** `"5"` and `"chr5"` compare equal everywhere this skill compares chromosomes -- a locus on one naming convention against a track on the other used to silently read as no overlap / `not_assayed` instead of matching. Every output field still reports the chromosome exactly as its own source file spelled it; only the comparison is normalized. This skill has no reference genome to check against, so it cannot detect (and does not attempt to detect) a genome-build mismatch -- two different assemblies can use identical `chr5` labels for coordinates that mean different things. Confirm all your input files share one assembly before trusting a result.
- **Nearest-gene and contact-gene disagreement must be surfaced, not silently resolved.** Early drafts defaulted to "contact wins" — this is wrong: both are real candidate answers with different evidence bases, and collapsing them hides exactly the kind of judgment call a reviewer needs to see (demonstrated by `L4_disagree` in the demo data).
- **A locus's contact coverage comes from its own loop, not every loop on the chromosome.** Merging all loops' anchors into one chromosome-wide span would let a single distant, unrelated loop make a locus look "assayed for contact" when no loop is actually anywhere near it — this skill checks each loop's own span independently instead, so `not_assayed` stays as tight as the real evidence.
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
- **Resolved**: nearest-gene ranking now measures from the locus midpoint to each gene's true TSS point, carried alongside its regulatory window in a dedicated `GeneWindow` tuple type, rather than to the nearer edge of the window -- the window-edge method could occasionally pick the wrong gene as nearest.
- **Staleness signals**: none — this skill has no external database dependency.

## Citations

- Fulco, C.P. et al. (2019). Activity-by-contact model of enhancer-promoter regulation from thousands of CRISPR perturbations. *Nature Genetics*, 51, 1664–1669. (Conceptual basis for weighting contact evidence higher than a single peak call -- **not an implementation of the ABC score itself**. The real ABC model multiplies measured activity by measured contact frequency; this skill's score is a simpler weighted presence/absence sum over whichever evidence layers were supplied. A reader expecting an actual activity-times-contact calculation should look elsewhere; this is ABC-inspired, not ABC.)
