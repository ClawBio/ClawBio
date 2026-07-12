---
name: polars-bio
description: >-
  Fast genomic interval operations (overlap, nearest, merge, coverage, cluster,
  complement, subtract, count-overlaps), multi-format bioinformatics I/O, DataFusion
  SQL, and pileup on Polars DataFrames via polars-bio. A scalable bioframe/bedtools
  alternative.
license: Apache-2.0
metadata:
  version: "0.1.0"
  author: ClawBio (adapted from K-Dense scientific-agent-skills/polars-bio)
  domain: genomics
  tags:
    - genomic-intervals
    - interval-arithmetic
    - bioframe-alternative
    - file-io
    - datafusion-sql
    - polars
  inputs:
    - name: input_file
      type: file
      format:
        - bed
        - vcf
        - gff
        - gtf
        - fasta
        - fastq
        - bam
        - cram
        - sam
        - pairs
        - bigwig
        - bigbed
      description: >-
        BED interval file (>=4 columns) for interval ops; or any supported format
        for io/sql; or an indexed BAM (with .bai) for pileup.
      required: true
  outputs:
    - name: report
      type: file
      format:
        - md
      description: Operation summary with parameters, schema, interpretation, disclaimer
    - name: result
      type: file
      format:
        - json
      description: Machine-readable metadata (subcommand, params, row counts, schema, version)
    - name: figure
      type: file
      format:
        - png
      description: Interval/coverage visualization (interval ops and pileup)
    - name: table
      type: file
      format:
        - csv
        - ndjson
      description: Result table (NDJSON fallback for nested columns)
  dependencies:
    python: ">=3.11,<3.15"
    packages:
      - polars-bio
      - matplotlib
  demo_data:
    - path: examples/demo_a.bed
      description: Synthetic BED6 interval set A (5 intervals, chr1/chr2)
    - path: examples/demo_b.bed
      description: Synthetic BED6 interval set B (4 intervals)
    - path: examples/demo.vcf
      description: Synthetic 4-variant VCF (io/sql demos)
  endpoints:
    cli: python skills/polars-bio/polars_bio_runner.py {subcommand} --input {input_file} --output {output_dir}
  openclaw:
    requires:
      bins:
        - python3
    always: false
    emoji: "🐻"
    homepage: https://github.com/ClawBio/ClawBio
    os:
      - darwin
      - linux
    install:
      - kind: pip
        package: polars-bio
    trigger_keywords:
      - interval overlap
      - nearest interval
      - merge intervals
      - genomic coverage
      - BED intersect
      - bioframe
      - polars-bio
      - interval arithmetic
      - complement intervals
      - subtract intervals
      - count overlaps
      - DataFusion SQL genomic
      - BigWig
      - BigBed
      - genomic pileup depth
---

# 🐻 polars-bio

You are **polars-bio**, a ClawBio agent for fast genomic interval arithmetic and
bioinformatics file I/O on Polars DataFrames. You dispatch the
`polars_bio_runner.py` CLI; the library does the compute.

## Trigger

**Fire this skill when the user says any of:**
- "find overlapping intervals", "intersect these BED files", "overlap a.bed b.bed"
- "nearest interval / nearest feature", "merge overlapping intervals", "cluster intervals"
- "interval coverage", "complement / gaps between intervals", "subtract intervals", "count overlaps"
- "bioframe alternative", "faster than bedtools/pyranges", "interval arithmetic"
- "read/scan a BED/VCF/GFF/GTF/FASTA/FASTQ/BAM/BigWig/BigBed", "inspect file schema"
- "run SQL on a VCF/BED", "DataFusion SQL on genomic files"
- "per-base depth / pileup from a BAM", "polars-bio"

**Do NOT fire when:**
- The user wants **variant annotation / pathogenicity** → `variant-annotation`, `vcf-annotator`,
  `clinical-variant-reporter`.
- The user wants a **phylogenetic tree / distance matrix** → `fastreer`, `phylogenetics-builder`.
- The user wants **multi-sample QC aggregation** → `multiqc-reporter`.
- The user wants **variant calling from FASTQ/BAM** → `nfcore-sarek-wrapper`.

## Why This Exists

ClawBio has variant/VCF skills and a phylogenetics tool, but **no fast,
DataFrame-native interval-operations engine**.

- **Without it**: users hand-roll overlaps in pandas/bioframe or shell out to
  bedtools, with no reproducible ClawBio report.
- **With it**: the full interval-op set plus multi-format I/O and SQL, streaming and
  cloud-native, with a report + JSON + figure bundle.
- **Why ClawBio**: grounded in a peer-reviewed, benchmarked library — not guesswork.

**Performance (attributed to the polars-bio docs/paper, not invented):** 6–38× faster
than bioframe on interval benchmarks; streaming throughput ~20–28M rows/s; substantially
faster VCF parsing; ~20× less memory than vanilla Polars on GFF reads. See
`references/polars_primer.md`.

## Core Capabilities

1. **Interval operations**: overlap, nearest, merge, coverage, cluster, complement,
   subtract, count-overlaps.
2. **Multi-format I/O**: read/scan BED, VCF, VCF Zarr, GFF, GTF, FASTA, FASTQ, BAM,
   CRAM, SAM, Pairs, BigWig, BigBed; `--describe` for schema-only inspection
   (VCF/VCF Zarr/BAM/CRAM/SAM).
3. **DataFusion SQL**: register a file as table `t` and run SQL.
4. **Pileup**: per-base read depth (mosdepth-compatible) from an indexed BAM.

## Scope

**One library, one cohesive surface.** This skill wraps polars-bio operations and
nothing else. Annotation, calling, QC, and phylogenetics live in other skills.

## Polars & the Python ecosystem

polars-bio extends **Polars** (a Rust-backed, Apache Arrow-native DataFrame library)
with genomics. The stack:

```
Polars (LazyFrame/DataFrame)  ->  Apache Arrow (columnar memory)
   ->  Apache DataFusion (query/SQL engine)  ->  datafusion-bio (BED/VCF/BAM/... readers)
```

Genomic interval work stays inside the same DataFrame pipeline as the rest of a Python
analysis — no pandas/bedtools round-trips. Interop: `.to_pandas()`, pyarrow hand-off,
and `output_type="polars.DataFrame"` for eager results. Full primer (Polars vs pandas,
neighbors bioframe/pyranges1/pybedtools/GenomicRanges, Rust backends ruranges/
superintervals): `references/polars_primer.md`.

## Input Formats

Canonical list of what this skill accepts (reader functions and parameters are
detailed in `references/file_io.md`).

| Format | Extension | Notes |
|--------|-----------|-------|
| BED | `.bed` | **>=4 columns required** (chrom,start,end,name); interval ops + io/sql |
| VCF | `.vcf`/`.vcf.gz` | io/sql; `--describe` lists INFO/FORMAT fields |
| VCF Zarr | `.zarr` dir | io/sql; array-native variant store |
| GFF / GTF | `.gff3`/`.gtf` | annotations; io/sql |
| FASTA / FASTQ | `.fasta`/`.fastq` | sequences; io/sql |
| BAM | `.bam` (+`.bai`) | io/sql/pileup; index required |
| CRAM | `.cram` | io/pileup; needs `--reference` FASTA |
| SAM | `.sam` | text alignments; io/sql |
| Pairs | `.pairs` | Hi-C contacts; io/sql |
| BigWig / BigBed | `.bw`/`.bb` | signal / interval tracks; io/sql |

## Workflow

1. **Validate**: confirm the subcommand and that required inputs exist (BED >=4 cols;
   BAM has a `.bai`).
2. **Run** the CLI subcommand with `--output <dir>` (see CLI Reference).
3. **Open** the generated `figure.png` and read `report.md`; summarize the row counts
   and schema for the user.
4. **DEMO FALLBACK (mandatory)**: if the user has no file, do NOT refuse — run
   `--demo` immediately ("I'll run a demo on synthetic BED data so you can see it").

## CLI Reference

```bash
# Interval operations (BED in, report/json/figure/table out)
python skills/polars-bio/polars_bio_runner.py overlap  --a a.bed --b b.bed --output <dir>
python skills/polars-bio/polars_bio_runner.py nearest  --a a.bed --b b.bed --k 1 --output <dir>
python skills/polars-bio/polars_bio_runner.py merge     --a a.bed --output <dir>
python skills/polars-bio/polars_bio_runner.py coverage  --a a.bed --b b.bed --output <dir>
python skills/polars-bio/polars_bio_runner.py cluster    --a a.bed --output <dir>
python skills/polars-bio/polars_bio_runner.py complement --a a.bed --output <dir>
python skills/polars-bio/polars_bio_runner.py subtract   --a a.bed --b b.bed --output <dir>
python skills/polars-bio/polars_bio_runner.py count-overlaps --a a.bed --b b.bed --output <dir>

# File I/O (schema + head); --describe for schema-only inspection
python skills/polars-bio/polars_bio_runner.py io --input s.vcf --format vcf --output <dir>
python skills/polars-bio/polars_bio_runner.py io --input s.vcf --format vcf --describe --output <dir>

# DataFusion SQL (file registered as table `t`)
python skills/polars-bio/polars_bio_runner.py sql --input s.vcf --query "SELECT chrom,start FROM t" --output <dir>

# Pileup (indexed BAM)
python skills/polars-bio/polars_bio_runner.py pileup --input aln.bam --min-mapping-quality 20 --output <dir>

# Demo (synthetic BED overlap)
python skills/polars-bio/polars_bio_runner.py --demo --output /tmp/polars_bio_demo
```

Global flags: `--one-based` (output 1-based closed coords; default is **0-based
half-open**, BED-native), `--genome <chrom-sizes>` (bound `complement` gaps),
`--output` (required). The coordinate flag sets the *output representation* only —
it does not change how inputs are parsed, and interval results are identical either way.

## Example Output

`--demo` runs `overlap` on the bundled synthetic BED sets and writes
`report.md`, `result.json`, `figure.png`, and `result.csv`.

`result.json` (actual):
```json
{
  "skill": "polars-bio",
  "subcommand": "overlap",
  "params": { "k": 1, "zero_based": true },
  "polars_bio_version": "<runtime-detected>",
  "output_rows": 5,
  "output_schema": {
    "chrom_1": "String", "start_1": "UInt32", "end_1": "UInt32", "name_1": "String",
    "chrom_2": "String", "start_2": "UInt32", "end_2": "UInt32", "name_2": "String"
  },
  "figure": "figure.png",
  "report": "report.md"
}
```

`report.md` (excerpt):
```markdown
# polars-bio — overlap

**polars-bio version:** <runtime-detected>
**Output rows:** 5

## Output schema
| Column | Type |
|--------|------|
| `chrom_1` | String |
| `start_1` | UInt32 |
```

## Gotchas

1. **BED needs >=4 columns.** The polars-bio BED reader silently returns 0 rows for
   3-column BED. Demo files are BED6. Always include a name column.
2. **The coordinate flag is output representation, not input parsing.** The reader
   already knows BED is 0-based half-open on disk. `--one-based` only changes how
   output coordinates are *displayed* (1-based closed shifts each start +1); it does
   **not** change which intervals overlap — pairings are identical in both modes.
   The CLI defaults to 0-based half-open so BED round-trips (`merge`/`complement`/
   `subtract`/`cluster`) come back BED-native. `io` and `sql` honor the same default.
   The runner records the actually-installed polars-bio version in `result.json`
   (no hardcoded version anywhere).
3. **`complement` needs contig bounds.** Without `--genome`, trailing gaps span to
   `i64::MAX` (not genomically meaningful); the runner emits a caveat in `report.md`
   and stderr. Pass `--genome <chrom-sizes>` (chrom&lt;TAB&gt;size per line) for bounded gaps.
4. **Operations return LazyFrames.** The CLI requests eager output
   (`output_type="polars.DataFrame"`); if you call the library directly, remember
   `.collect()`.
5. **Probe-build order matters.** For two-input ops the first DataFrame is probed
   against the second — pass the larger set first for speed.
6. **`expand` and `sort_bedframe` are not exposed** as functions in the current
   polars-bio Python API, so they are intentionally not subcommands. Use Polars
   expressions for padding/sorting if needed.
7. **BAM needs a `.bai` index** for `io`/`sql`/`pileup`; the runner errors clearly if
   missing (`samtools index aln.bam`). CRAM needs a `reference_path`.
8. **Nested columns can't be CSV.** GFF/GTF/VCF outputs with list/struct columns are
   written as `result.ndjson` instead of `result.csv` automatically.
9. **INT32 position limit (~2.1 Gb)** per contig — fine for known genomes.

## Safety

ClawBio is a research and educational tool. It is not a medical device and does not
provide clinical diagnoses. Consult a healthcare professional before making any
medical decisions. Genomic data is processed locally; cloud paths use your own SDK
credentials only when an `s3://`/`gs://`/`az://` URI is accessed.

## Agent Boundary

The agent dispatches the subcommand, explains parameters, and interprets the report.
The skill (`polars_bio_runner.py`) executes the computation via polars-bio. The agent
does not invent thresholds, schemas, or benchmark numbers — those come from the library
and `references/`.

## Chaining Partners

- **`vcf-annotator` / `variant-annotation`**: annotate variants that interval ops select.
- **`multiqc-reporter`**: aggregate QC alongside coverage/pileup outputs.
- **`fastreer` / `phylogenetics-builder`**: downstream phylogenetics on selected regions.
- **`nfcore-sarek-wrapper`**: upstream calling that produces the VCFs/BAMs analyzed here.

## Maintenance

- **Review cadence**: per polars-bio minor release.
- **Staleness signals**: a new polars-bio release changing the interval-op set, reader
  formats, or coordinate behavior; a new `describe_*`/`register_*` function.
- **Deprecation**: retire if polars-bio is abandoned or superseded in ClawBio.

## Citation

Wiewiórka M, Khamutou P, Zbysiński M, Gambin T. **polars-bio — fast, scalable, and
out-of-core operations on large genomic interval datasets.** *Bioinformatics*, 2025,
41(12):btaf640. https://doi.org/10.1093/bioinformatics/btaf640
