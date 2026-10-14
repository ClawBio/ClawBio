---
name: fastreer
description: >-
  Phylogenetic distance matrices and trees from VCF or FASTA data using the
  fastreeR hybrid Java/Python toolkit (VCF2TREE, VCF2DIST, DIST2TREE, FASTA2DIST).
license: GPL-3.0
metadata:
  version: "0.1.0"
  author: Anestis Gkanogiannis
  domain: phylogenetics
  tags:
    - phylogenetics
    - distance-matrix
    - tree-building
    - vcf
    - fasta
    - population-genomics
  inputs:
    - name: input_file
      type: file
      format:
        - vcf
        - vcf.gz
        - fasta
        - fasta.gz
        - fa
        - fa.gz
        - fas
        - fas.gz
        - dist
      description: >-
        VCF file (biallelic/multiallelic SNPs, compressed or plain),
        FASTA file (aligned or unaligned sequences, compressed or plain), or
        PHYLIP distance matrix (for DIST2TREE).
      required: true
  outputs:
    - name: tree
      type: file
      format:
        - nwk
      description: Newick phylogenetic tree (VCF2TREE / DIST2TREE)
    - name: distances
      type: file
      format:
        - dist
      description: PHYLIP distance matrix (VCF2DIST / FASTA2DIST)
    - name: report
      type: file
      format:
        - md
      description: Analysis summary with sample list and interpretation
    - name: result
      type: file
      format:
        - json
      description: Machine-readable metadata (samples, command, paths)
  dependencies:
    python: ">=3.10"
    packages:
      - fastreer>=2.2.0
    system:
      - java>=11
  demo_data:
    - path: examples/demo_samples.vcf
      description: Synthetic VCF with 5 samples and 20 biallelic SNPs on chr1
    - path: examples/demo_sequences.fasta
      description: Synthetic FASTA with 5 sequences of 60 bp each
  endpoints:
    cli: >-
      python skills/fastreer/fastreer.py --command {command} --input {input_file} --output {output_dir}
  openclaw:
    requires:
      bins:
        - python3
        - java
    always: false
    emoji: "🌳"
    homepage: https://github.com/ClawBio/ClawBio
    os:
      - darwin
      - linux
    install:
      - kind: pip
        package: fastreer
    trigger_keywords:
      - phylogenetic tree
      - distance matrix from VCF
      - VCF
      - VCF2TREE
      - VCF2DIST
      - fastreer
      - fastreeR
      - genomic distance
      - hierarchical clustering tree
      - k-mer distance
      - FASTA2DIST
      - variant distance matrix
      - population tree
      - DIST2TREE
      - cosine distance VCF
      - sample phylogeny
---

# fastreeR

You are **fastreeR**, a specialised ClawBio skill for computing phylogenetic
distance matrices and trees from genomic VCF or FASTA data using the
[fastreeR](https://github.com/gkanogiannis/fastreeR) hybrid Java/Python toolkit.

## Trigger

**Fire this skill when the user says any of:**
- "build a phylogenetic tree from my VCF"
- "compute a distance matrix from variants"
- "VCF2TREE", "VCF2DIST", "DIST2TREE", "FASTA2DIST"
- "fastreer" or "fastreeR"
- "how similar are my samples genetically"
- "genomic distance between samples"
- "population tree from VCF"
- "k-mer distance from FASTA"
- "hierarchical clustering of samples"
- "cosine distance from genotypes"
- "sample distance matrix"

**Do NOT fire when:**
- The user wants population genetics statistics (π, Tajima's D, Fst) → route to `dnasp`
- The user wants protein structure prediction → route to `struct-predictor`
- The user wants alignment (not tree building) → use `seq-wrangler`
- The user wants ancestry/PCA decomposition → route to `claw-ancestry-pca`
- The user wants variant annotation → route to `variant-annotation`

## Why This Exists

- **Without it**: Building phylogenetic trees from VCF requires awkward conversion steps
  (VCF → PLINK → distance matrix → external tree software) with no unified output.
- **With it**: One command converts a VCF or FASTA directly to a Newick tree or
  PHYLIP distance matrix, with optional bootstrap support and windowed analysis.
- **Why ClawBio**: fastreeR is purpose-built for large population VCFs; it streams
  data in O(n_samples²) RAM rather than loading everything into memory.

## Core Capabilities

1. **VCF2TREE**: Computes cosine dissimilarity between samples and builds a
   hierarchical clustering tree directly from a VCF, with optional bootstrap resampling.
2. **VCF2DIST / FASTA2DIST**: Exports the underlying PHYLIP distance matrix for use
   in downstream tools (R, Python, ape, BioPython).
3. **Windowed analysis**: Streams per-window trees or matrices across genomic regions
   via `--window-bp` or `--window-variants`.

## Scope

This skill computes pairwise genomic distances and
hierarchical trees from VCF or FASTA input. It does not perform alignment, variant
calling, variant annotation, or population genetics statistics.

## Input Formats

| Format | Extension | Required Fields | Example |
|--------|-----------|-----------------|---------|
| VCF | `.vcf`, `.vcf.gz` | GT genotype field; ≥2 samples | `samples.vcf.gz` |
| FASTA | `.fasta`, `.fasta.gz`, `.fa`, `.fa.gz`, `.fas`, `.fas.gz` | ≥2 sequences | `sequences.fasta` |
| PHYLIP dist | `.dist` | PHYLIP matrix header + rows | `distances.dist` |

## Workflow

When the user provides a VCF or FASTA:

1. **Validate**: Confirm input file exists; detect format from extension; check Java 11+ is installed
2. **Select command**:
   - VCF + want tree → `VCF2TREE`
   - VCF + want distances only → `VCF2DIST`
   - Distance matrix + want tree → `DIST2TREE`
   - FASTA + want k-mer distances → `FASTA2DIST`
3. **Run fastreeR**: Invoke via `fastreer.py` with appropriate flags (threads, mem, bootstrap)
4. **Generate outputs**: Write `tree.nwk` or `distances.dist`, `report.md`, `result.json`,
   and reproducibility bundle
5. **Explain**: Summarise the tree topology or distance range; note any bootstrap support

**Freedom levels:**
- Steps 1–3 (execution): prescriptive; exact flags must be used
- Step 5 (interpretation): flexible; reason from the Newick or distance values

## CLI Reference

```bash
# Newick tree from VCF (with bootstrap)
python skills/fastreer/fastreer.py \
  --command VCF2TREE --input samples.vcf.gz --bootstrap 100 \
  --threads 4 --output <report_dir>

# Distance matrix from VCF
python skills/fastreer/fastreer.py \
  --command VCF2DIST --input samples.vcf.gz --threads 4 \
  --output <report_dir>

# Tree from pre-computed distance matrix
python skills/fastreer/fastreer.py \
  --command DIST2TREE --input distances.dist --output <report_dir>

# K-mer distance from FASTA sequences
python skills/fastreer/fastreer.py \
  --command FASTA2DIST --input sequences.fasta --kmer 5 \
  --output <report_dir>

# Windowed analysis (100 kb windows)
python skills/fastreer/fastreer.py \
  --command VCF2DIST --input samples.vcf.gz --window-bp 100000 \
  --output <report_dir>

# Demo (no data needed)
python skills/fastreer/fastreer.py --demo --output /tmp/fastreer_demo

# Via ClawBio runner
python clawbio.py run fastreer --demo
python clawbio.py run fastreer --input samples.vcf.gz
```

## Demo

```bash
python clawbio.py run fastreer --demo
```

Expected output: VCF2TREE run on a synthetic 5-sample / 20-SNP VCF. Produces a
Newick tree (`tree.nwk`), `report.md` with sample list and interpretation, and a
reproducibility bundle. If Java / fastreeR is not installed, synthetic demo output
is generated to illustrate the expected format.

## Algorithm / Methodology

**VCF2TREE / VCF2DIST** (cosine dissimilarity from genotypes):

1. For each sample pair (i, j), compute the cosine dissimilarity over all biallelic
   variant sites as: `d(i,j) = 1 - cosine_similarity(gt_vector_i, gt_vector_j)`
   where genotypes are encoded as allele dosages (0/0→0, 0/1→1, 1/1→2).
2. Build an N×N PHYLIP distance matrix; emit to `.dist` file.
3. For tree building: apply average-linkage (UPGMA) hierarchical clustering to the
   distance matrix; emit Newick with optional bootstrap node labels.

**FASTA2DIST** (D2S k-mer distance):

1. For each sequence, compute k-mer frequency vectors (default k=4).
2. Apply the D2S statistic (Reinert et al. 2009) to compute pairwise distances.
3. Emit PHYLIP matrix.

**Key parameters**:
- `--threads`: parallelism for distance computation (default: 1)
- `--mem`: JVM heap in MB (default: 256; increase for >500 samples)
- `--bootstrap`: streaming bootstrap replicates from VCF (VCF2TREE only)
- `--kmer`: k-mer size for FASTA2DIST (default: 4; range 3–8 typical)

## Example Queries

- "Build a phylogenetic tree from my population VCF"
- "Compute a distance matrix between my 200 samples using VCF2DIST"
- "Run fastreer on sequences.fasta with k=5"
- "Show me how similar my samples are genetically"
- "Use DIST2TREE to convert my distance matrix to a Newick tree"

## Example Output

```markdown
# fastreeR Report

**Command**: `VCF2TREE`
**Input**: `demo_samples.vcf` (5 samples, 20 variants)
**Date**: 2026-05-11

## Samples (5)
  - SAMPLE1
  - SAMPLE2
  - SAMPLE3
  - SAMPLE4
  - SAMPLE5

## Phylogenetic Tree

**Output format**: Newick
**File**: `tree.nwk`

((SAMPLE1:0.120,SAMPLE2:0.098):0.045,
 (SAMPLE3:0.110,(SAMPLE4:0.087,SAMPLE5:0.132):0.062):0.038);

SAMPLE1 and SAMPLE2 cluster together (distance 0.12), suggesting greater
genomic similarity relative to SAMPLE3–5. SAMPLE4 and SAMPLE5 are the
second closest pair (distance 0.087).
```

## Output Structure

```
output_directory/
├── report.md              # Summary: samples, tree/matrix preview, interpretation
├── result.json            # Machine-readable: command, samples, paths, metadata
├── tree.nwk               # Newick tree (VCF2TREE / DIST2TREE)
├── distances.dist         # PHYLIP distance matrix (VCF2DIST / FASTA2DIST)
└── reproducibility/
    ├── commands.sh        # Portable replay command ($CLAWBIO_ROOT / $OUTPUT_DIR)
    ├── environment.yml    # Conda recipe (fastreer + openjdk)
    ├── environment.txt    # Java version + pip fastreer version
    └── checksums.sha256   # SHA-256 of every output file
```

## Dependencies

**Required**:
- `fastreer >= 2.2.0` (install with `pip install fastreer`)
- **Java 11+**: the Python package wraps a Java backend; install via
  `sudo apt install default-jre` (Linux) or `brew install openjdk@17` (macOS)

**Optional**:
- `matplotlib`, for tree/heatmap visualisation in future versions

## Gotchas

- **Java version check**: The model may assume Python alone is sufficient. It is not.
  fastreeR's core is a Java application. Always check Java 11+ is present before
  running; emit a clear error if missing, not a cryptic JVM crash.

- **VCF must have sample columns**: Variant-only VCFs (no FORMAT/GT fields, no sample
  columns) will silently fail or produce empty output. Validate that `#CHROM` line
  has columns beyond FORMAT (i.e., at least one sample name).

- **JVM heap for large datasets**: The default `--mem 256` is insufficient for >500
  samples. Rule of thumb: `4 × n_samples² × n_threads / 1e6` MB. For 1000 samples
  with 8 threads: ~32 GB. Document this prominently or auto-compute a suggested value.

- **Windowed output is multi-block**: `--window-bp` produces a single file containing
  multiple concatenated PHYLIP matrices or Newick trees separated by comment lines.
  Do not attempt to parse it as a single matrix.

- **VCF2EMB is not included**: The embedding command requires downloading a 500MB
  BioFM language model. It is intentionally excluded from v0.1.0. If the user asks
  for variant embeddings, explain the requirement and the manual install steps.

- **Compressed VCF via stdin**: Piping `zcat input.vcf.gz | fastreer VCF2TREE -i -`
  works but the `-` stdin mode requires fastreeR ≥ 2.1.0 and may not stream on
  Windows. Use `-i input.vcf.gz` directly for portability.

## Safety

- **Local-first**: All computation is local. No genomic data is sent externally.
- **Disclaimer**: Every report includes the ClawBio medical disclaimer.
- **Audit trail**: `reproducibility/commands.sh` records the exact command run.
- **No hallucinated science**: All distance formulas trace to fastreeR source and cited papers.

## Agent Boundary

The agent (LLM) dispatches and explains results. The Python script (`fastreer.py`)
executes fastreeR and writes outputs. The agent must NOT invent tree topologies,
distance values, or bootstrap support figures; all must come from fastreeR output.

## Integration with Bio Orchestrator

**Trigger conditions**: the orchestrator routes here when:
- User provides a VCF and asks for tree/distance/phylogenetics
- User mentions `fastreer`, `fastreeR`, `VCF2TREE`, `VCF2DIST`, `FASTA2DIST`
- User asks "how similar are my samples" with a VCF or FASTA

**Chaining partners**:
- `dnasp`: Run DnaSP population statistics on the same VCF, then fastreeR for tree
- `variant-annotation`: Annotate variants first, then build a tree to visualise population structure
- `claw-ancestry-pca`: use PCA for admixture and fastreeR for hierarchical clustering; the two provide complementary views of population structure
- `seq-wrangler`: Align sequences first (seq-wrangler), then compute FASTA2DIST tree

## Maintenance

- **Review cadence**: Check for new fastreeR releases quarterly (PyPI: `pip index versions fastreer`)
- **Staleness signals**: New fastreeR CLI flags not exposed; VCF2EMB added to scope
- **Deprecation**: Archive to `skills/_deprecated/` if fastreeR is superseded or unmaintained

## Citations

- [fastreeR GitHub](https://github.com/gkanogiannis/fastreeR): source code, documentation, and Docker image
- [fastreeR PyPI](https://pypi.org/project/fastreer/): Python package
- [fastreeR Bioconductor](https://bioconductor.org/packages/fastreeR/): R/Bioconductor package
- Gkanogiannis A (2016) A scalable assembly-free variable selection algorithm for biomarker discovery from metagenomes. *BMC Bioinformatics* 17, 311. https://doi.org/10.1186/s12859-016-1186-3
- Reinert G et al. (2009) Alignment-free sequence comparison (I): statistics and power. *J Comput Biol* 16(12):1615-34. D2S k-mer statistic used in FASTA2DIST.
