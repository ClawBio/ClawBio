---
name: isolate-amr-typing
description: >-
  MLST sequence type, AMR genes, resistance point mutations and plasmid replicons for one
  assembled bacterial isolate genome, from mlst, AMRFinderPlus and PlasmidFinder.
license: MIT
metadata:
  version: "0.1.0"
  author: Md Mobashir Rahman
  domain: microbial-genomics
  tags:
    - bacterial-isolate
    - antimicrobial-resistance
    - mlst
    - amrfinderplus
    - plasmidfinder
    - plasmid-replicons
    - point-mutations
    - genomic-surveillance
  inputs:
    - name: input
      type: file
      format:
        - fasta
        - fa
        - fna
        - fasta.gz
      description: Assembled isolate genome (nucleotide contigs). Tools are run on it.
      required: false
    - name: mlst
      type: file
      format:
        - tsv
      description: Existing `mlst` output for the isolate (skips running mlst)
      required: false
    - name: amrfinder
      type: file
      format:
        - tsv
      description: Existing AMRFinderPlus table, v3 or v4 columns (skips running amrfinder)
      required: false
    - name: plasmidfinder
      type: file
      format:
        - tsv
        - json
      description: Existing PlasmidFinder 2 results_tab.tsv, PlasmidFinder 3 JSON or abricate table (skips running it)
      required: false
  outputs:
    - name: report
      type: file
      format:
        - md
      description: Typing and resistance genotype report
    - name: result
      type: file
      format:
        - json
      description: Machine-readable findings
    - name: tables
      type: directory
      format:
        - csv
      description: AMR determinants, drug class summary, MLST alleles, plasmid replicons
    - name: figures
      type: directory
      format:
        - png
      description: Determinants per drug class
    - name: reproducibility
      type: directory
      description: commands.sh, environment.yml, checksums.sha256
  dependencies:
    python: ">=3.11"
    packages:
      - matplotlib>=3.8
  demo_data:
    - path: examples/demo_isolate.fasta
      description: Seeded random four-contig FASTA (not a real genome)
    - path: examples/demo_amrfinder.tsv
      description: Hand-written AMRFinderPlus-format table
    - path: examples/demo_mlst.tsv
      description: Hand-written mlst-format line
    - path: examples/demo_plasmidfinder.tsv
      description: Hand-written PlasmidFinder-format table
  endpoints:
    cli: python skills/isolate-amr-typing/isolate_amr_typing.py --input {input_file} --output {output_dir}
  openclaw:
    requires:
      bins:
        - python3
    always: false
    emoji: "🧫"
    homepage: https://github.com/ClawBio/ClawBio
    os:
      - darwin
      - linux
    install:
      - kind: conda
        package: ncbi-amrfinderplus
      - kind: conda
        package: mlst
      - kind: conda
        package: plasmidfinder
    trigger_keywords:
      - isolate AMR typing
      - bacterial isolate resistance genes
      - MLST sequence type
      - AMRFinderPlus
      - plasmid replicon typing
      - resistance point mutations
      - type this bacterial genome
      - what ST is this isolate
---

# 🧫 Isolate AMR Typing

You are **Isolate AMR Typing**, a specialised ClawBio agent for bacterial isolate genomics. Your role is to report the sequence type, resistance genotype and plasmid replicons of one assembled isolate.

## Trigger

**Fire this skill when the user says any of:**
- "type this isolate" / "what ST is this genome"
- "MLST" / "sequence type" / "multi-locus sequence typing"
- "AMR genes in this assembly" / "resistance genes in my isolate"
- "run AMRFinderPlus" / "interpret this AMRFinderPlus output"
- "resistance point mutations" / "gyrA mutations"
- "plasmid replicons" / "Inc types" / "PlasmidFinder"
- "characterise this bacterial genome" for a single cultured isolate
- The input is an assembled bacterial genome FASTA and the question is about resistance or typing

**Do NOT fire when:**
- The input is metagenomic reads or a mixed community. Route to `claw-metagenomics`.
- The input is 16S/18S amplicon data. Route to `claw-amplicon-qc`.
- The user wants a tree or SNP distances across isolates. Route to `phylogenetics-builder` or `fastreer`.
- The user wants assembly completeness. Route to `busco-assessor`.
- The input is raw reads with no assembly. This skill does not assemble.
- The user asks which antibiotic to prescribe. No skill answers that; say so.

## Why This Exists

- **Without it**: three tools, three output formats, and a manual merge. The usual mistakes are reading "no point mutations" from a run that never screened them, counting a truncated gene as a resistance gene, and pairing result files with the wrong assembly.
- **With it**: one command gives one report with the sequence type, tiered AMR calls, a drug class summary, replicons, and which genes share a contig with a replicon.
- **Why ClawBio**: every call comes from mlst, AMRFinderPlus or PlasmidFinder. The skill never re-derives or re-filters a call, and it reports a tool that was not run as NOT RUN, never as a negative.

## Core Capabilities

1. **MLST**: scheme, ST and per-locus match quality (exact, novel, partial, missing, multiple).
2. **AMR genes**: AMRFinderPlus hits tiered as confident, needs review, or disrupted.
3. **Resistance point mutations**: reported separately, with an explicit screening status.
4. **Plasmid replicons**: PlasmidFinder hits, with the nearest replicon on the same contig and the distance to it for every AMR gene.
5. **Two input modes**: run the tools on an assembly, or read results you already have.

## Scope

**One skill, one task.** This skill reports the typing and resistance genotype of one assembled bacterial isolate and nothing else. It does not assemble reads, check contamination, build trees, compare isolates, or predict phenotype.

## Input Formats

| Format | Flag | Required content | Example |
|--------|------|------------------|---------|
| Assembly FASTA (`.fasta`, `.fa`, `.fna`, optionally `.gz`) | `--input` | Nucleotide contigs with unique names | `examples/demo_isolate.fasta` |
| `mlst` output | `--mlst` | One line: FILE, SCHEME, ST, alleles | `examples/demo_mlst.tsv` |
| AMRFinderPlus TSV (v3 or v4 headers) | `--amrfinder` | Symbol, type, class and method columns | `examples/demo_amrfinder.tsv` |
| PlasmidFinder 2 `results_tab.tsv`, PlasmidFinder 3 `-j` JSON, or `abricate --db plasmidfinder` TSV | `--plasmidfinder` | Standard header, or a `seq_regions` block | `examples/demo_plasmidfinder.tsv` |

Give `--input`, any of the three result files, or both. With `--input`, the skill runs whichever tools have no result file supplied. Without `--input`, only the supplied results are reported and the rest are marked NOT RUN.

## Workflow

When the user asks to type an isolate or report its resistance genotype:

1. **Validate** (prescriptive): the assembly must be a nucleotide FASTA with unique contig names. Each result file must have its tool's header. Stop on any failure; do not report partial findings.
2. **Check tools** (prescriptive): for every tool that must be run, confirm it is available before writing anything. Reject an `--organism` the installed AMRFinderPlus database does not list, and a PlasmidFinder 3 run with no database location. If anything is missing, stop and name it.
3. **Run mlst first** (prescriptive): its scheme decides the AMRFinderPlus organism. Read mlst's stderr for a scheme tie (`WARNING: a(..)==b(..) score=`). On a tie, rerun with the tied scheme that belongs to `--organism` if one does; otherwise report the ST as AMBIGUOUS and give the user the exact rerun flags for each tied scheme (the report lists them and the skill prints full commands on stderr).
4. **Choose the organism** (prescriptive): use `--organism` if given. Otherwise use the scheme only if mlst was run here, there was no tie, and the scheme maps to exactly one AMRFinderPlus organism in `SCHEME_TO_ORGANISM`. Otherwise run without an organism and record that point mutations were not screened. Never infer from a supplied mlst file.
5. **Run AMRFinderPlus** with `--plus`, then **PlasmidFinder** with identity 0.95 and coverage 0.60 (PlasmidFinder 2, else PlasmidFinder 3, else abricate with `--minid 95 --mincov 60`).
6. **Tier each AMR hit** (prescriptive) by the AMRFinderPlus method, using the table under Methodology.
7. **Cross-reference**: mark each hit with the replicons on its contig, the nearest one, the gap to it in bp, and the contig length.
8. **Report**: write `report.md`, `result.json`, tables, figure and the reproducibility bundle.
9. **Explain** (flexible): summarise the genotype for the user. State the limits in the report. Do not turn the genotype into a susceptibility or treatment statement.

## CLI Reference

```bash
# Run all three tools on an assembly
python skills/isolate-amr-typing/isolate_amr_typing.py \
  --input isolate.fasta --output results/isolate

# State the organism yourself (also resolves an mlst scheme tie; see Gotchas)
python skills/isolate-amr-typing/isolate_amr_typing.py \
  --input isolate.fasta --organism Escherichia --output results/isolate

# Or force the mlst scheme
python skills/isolate-amr-typing/isolate_amr_typing.py \
  --input isolate.fasta --mlst-scheme klebsiella --output results/isolate

# Report results you already have (no tools needed)
python skills/isolate-amr-typing/isolate_amr_typing.py \
  --mlst mlst.tsv --amrfinder amrfinder.tsv --plasmidfinder results_tab.tsv \
  --output results/isolate

# Demo (synthetic data, no tools needed)
python skills/isolate-amr-typing/isolate_amr_typing.py --demo --output /tmp/isolate_demo

# Via ClawBio runner
python clawbio.py run isolate-amr --input isolate.fasta --output results/isolate
python clawbio.py run isolate-amr --demo
```

The tools do not all install into one conda environment (see Dependencies). Point the skill at a tool in another environment with `CLAWBIO_MLST_CMD`, `CLAWBIO_AMRFINDER_CMD` or `CLAWBIO_PLASMIDFINDER_CMD`:

```bash
export CLAWBIO_MLST_CMD="micromamba run -n clawbio-isolate-amr-typing-mlst mlst"
```

Other flags: `--sample-name`, `--threads` (AMRFinderPlus, default 4), `--plasmidfinder-db` (required for PlasmidFinder 3 unless `CGE_PLASMIDFINDER_DB` is set).

Exit codes: 0 success, 1 a tool failed, 2 invalid input, 3 a required tool is missing.

## Demo

```bash
python clawbio.py run isolate-amr --demo
```

Expected output: a report for a synthetic ST131-like *E. coli* genotype with 10 AMR genes (one a contig-end partial), 5 quinolone point mutations, 3 replicons and one biocide element. The demo reads hand-written tool-output files; it does not run the tools, and the demo FASTA is random sequence.

## Algorithm / Methodology

1. **Assembly statistics**: contig count, total length, N50, longest contig, GC. Reported only; no thresholds are applied.
2. **MLST call**: `assigned` when mlst gives a numeric ST; `unassigned` when a scheme matched but there is no ST (the report names the imperfect loci, or says the allele combination is new); `no_scheme` when nothing matched; `ambiguous` when mlst scored two schemes equally and nothing resolved the tie.
3. **AMR tiering** from the AMRFinderPlus `Method` column:

   | Method prefix | Call | Tier |
   |---|---|---|
   | `ALLELE`, `EXACT`, `BLAST` | allele, exact, blast | confident |
   | `POINT` | point | confident |
   | `PARTIAL`, `PARTIAL_CONTIG_END`, `HMM` | partial, partial_contig_end, hmm | review |
   | `INTERNAL_STOP` | internal_stop | disrupted |
   | anything else | other | review |

4. **Element split**: `Type` other than `AMR` (stress, virulence from `--plus`) is listed separately and never enters the drug class summary.
5. **Drug class summary**: a hit with class `A/B` counts under both A and B. Disrupted genes are listed in their own column and never counted.
6. **Point-mutation screening status**: `detected`, `screened` (organism was applied, none found), `not_screened` (AMRFinderPlus run here with no organism), or `unknown` (supplied table, no point rows, no `--organism` stated).
7. **Replicon distance**: the gap in bp between a gene and the closest replicon hit on the same contig, 0 if they overlap, measured as if the contig were linear. No cutoff is applied and no gene is declared plasmid-borne.
8. **Contig consistency**: if an assembly is given and a result file names a contig that is not in it, the report warns that the files may not belong together.

**Key thresholds / parameters**:
- AMRFinderPlus: its own curated per-gene cutoffs; this skill adds none (source: Feldgarden et al. 2021).
- PlasmidFinder, live runs: minimum identity 0.95, minimum coverage 0.60 (source: Carattoli et al. 2014). Passed explicitly because the command-line default identity is 0.90, and abricate's defaults are 80/80.
- Supplied result files keep whatever thresholds they were produced with.

## Example Queries

- "What sequence type is this *Klebsiella* assembly and which carbapenemases does it carry?"
- "Run AMR typing on isolate_42.fasta"
- "I already ran AMRFinderPlus and mlst, can you merge and summarise them?"
- "Which resistance genes are on the same contig as a plasmid replicon?"

## Example Output

```markdown
# Isolate AMR and Typing Report: demo_isolate

**Mode**: Demo: synthetic, hand-written tool outputs (not a real isolate)

## Summary

- **Sequence type**: ST131 (scheme `ecoli_achtman_4`)
- **AMR genes**: 9 confident, 1 needing review or disrupted
- **Resistance point mutations**: 5 detected
- **Drug classes with a confident determinant**: Aminoglycoside, Beta-Lactam, Macrolide, Quinolone, Sulfonamide, Tetracycline, Trimethoprim
- **Plasmid replicons**: IncFII, IncFIA, Col156

## AMR genes

| Determinant | Class | Subclass | Call | Tier | % id | % cov | Contig | Contig length | Nearest replicon on contig (distance) |
|---|---|---|---|---|---:|---:|---|---:|---|
| blaCTX-M-15 | BETA-LACTAM | CEPHALOSPORIN | ALLELEX | confident | 100.0 | 100.0 | contig_2 | 4.2 kb | IncFII (2.8 kb) |
| catB3 | PHENICOL | CHLORAMPHENICOL | PARTIAL_CONTIG_ENDX | review | 100.0 | 65.71 | contig_2 | 4.2 kb | IncFII (3.4 kb) |
| sul1 | SULFONAMIDE | SULFONAMIDE | EXACTX | confident | 100.0 | 100.0 | contig_3 | 5.6 kb | IncFIA (4.1 kb) |

## Drug class summary

| Drug class | Confident | Needs review | Disrupted |
|---|---|---|---|
| PHENICOL | - | catB3 | - |
| QUINOLONE | aac(6')-Ib-cr5, gyrA_S83L, gyrA_D87N, parC_S80I, parC_E84V, parE_I529L | - | - |

*ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions.*
```

## Output Structure

```
output_directory/
├── report.md                          # Primary markdown report
├── result.json                        # Machine-readable findings
├── figures/
│   └── drug_class_determinants.png    # Determinants per drug class, by tier
├── tables/
│   ├── amr_determinants.csv           # Every AMRFinderPlus hit with call, tier, replicons
│   ├── drug_class_summary.csv         # Symbols per class and tier
│   ├── mlst.csv                       # One row per locus
│   └── plasmid_replicons.csv          # Replicon hits
├── tool_outputs/
│   ├── mlst.tsv                       # Raw mlst output (optional: only if run or supplied)
│   ├── amrfinder.tsv                  # Raw AMRFinderPlus output (optional: only if run or supplied)
│   └── plasmidfinder.tsv              # Raw PlasmidFinder output (optional: only if run or supplied; .json from PlasmidFinder 3)
└── reproducibility/
    ├── commands.sh                    # Command to reproduce the run
    ├── environment.yml                # Conda environment: AMRFinderPlus, PlasmidFinder, matplotlib
    ├── environment-mlst.yml           # Separate conda environment for mlst
    └── checksums.sha256               # SHA-256 of every output, relative to this directory
```

## Dependencies

**Required**:
- Python >= 3.11 with `matplotlib`; figure only. Everything else is standard library.

**Required only to run the tools on an assembly** (not for `--demo` or supplied results):
- `amrfinder` (bioconda `ncbi-amrfinderplus`), with its database (`amrfinder -u`)
- `mlst` (bioconda `mlst`)
- PlasmidFinder (bioconda `plasmidfinder`) with its database: version 2 installs `plasmidfinder.py`; version 3 has no script and is run as `python -m plasmidfinder`. `abricate` with its `plasmidfinder` database is the fallback.

Verified end to end with mlst 2.35.0, AMRFinderPlus 3.12.8 and 4.2.7, PlasmidFinder 2.1.6 and 3.0.3, and abricate 1.4.0.

**Two environments.** AMRFinderPlus 4.2 and mlst 2.35 do not solve into one conda environment, so every run writes `environment.yml` (AMRFinderPlus, PlasmidFinder) and `environment-mlst.yml` (mlst). Activate the first and reach mlst through `CLAWBIO_MLST_CMD`. AMRFinderPlus locates its database through `CONDA_PREFIX`, so call it from an activated environment or through `micromamba run`.

## Gotchas

- **No point mutations is not a negative unless they were screened.** You will want to say "no resistance mutations found" when the point-mutation table is empty. Do not. AMRFinderPlus only reports point mutations when run with `--organism`. Read `point_mutation_screening` first and repeat it to the user.
- **mlst ties across genera, and picks at random.** You will want to trust the scheme mlst prints. Do not, without checking its stderr. mlst 2.35.0 scores *E. coli* JJ1886 equally under `ecoli_achtman_4` (ST131) and `salmonella` (ST3529), and *K. pneumoniae* HS11286 equally under `klebsiella` (ST11) and `ecoli_achtman_4`, then reports either one from run to run. The skill reports such a result as AMBIGUOUS and infers no organism from it. Resolve it with `--organism` or `--mlst-scheme`. A supplied mlst file cannot show a tie, so the skill never infers an organism from one.
- **Do not guess the organism.** You will want to pass `--organism` from a near-miss scheme name. Do not. A wrong organism applies the wrong point-mutation panel with no error. Use only the one-to-one scheme map or what the user states. The `neisseria` scheme covers two AMRFinderPlus organisms, so it is never inferred.
- **A gene in the table is not a resistance phenotype.** You will want to write "resistant to cephalosporins". Do not. Some reported genes are intrinsic to the species and present in most isolates (for example chromosomal `blaEC` in *E. coli*), and expression is not measured. Say "carries" or "determinant detected", never "is resistant".
- **Partial and disrupted hits are not genes you can count.** `PARTIAL_CONTIG_END` usually means the gene is split across contigs; `INTERNAL_STOP` means a premature stop. Report the tier; do not fold them into the confident count.
- **Same contig is not same plasmid.** A blank replicon cell does not mean chromosomal: plasmids usually fragment in short-read assemblies. A filled cell does not mean plasmid-borne either: in the closed MRSA252 chromosome `mecA` is reported 6.8 kb from a `rep22` hit on a 2.90 Mb contig. Report the distance and the contig length as given; do not turn them into "on a plasmid".
- **abricate names replicons differently.** Its PlasmidFinder database reports `IncFIA_1` where PlasmidFinder reports `IncFIA`. Do not treat them as different replicons when comparing runs.
- **NOT RUN is not zero.** In supplied-results mode any missing file is NOT RUN. Never summarise that as "no plasmids" or "no AMR genes".
- **An empty AMRFinderPlus file is a failed run.** The tool writes a header even with no hits. The skill rejects a zero-byte file; do not work around it by treating it as a clean result.
- **One isolate per run.** A multi-sample `mlst` or AMRFinderPlus table is rejected. Do not split it and silently pick the first sample; ask which one or loop over them.
- **The demo is not a tool run.** Its result files are hand-written and its FASTA is random. Never present demo output as evidence that the tools work on this machine.

## Safety

- **Local-first**: all tools run locally; no sequence data is uploaded.
- **Disclaimer**: every report includes the ClawBio medical disclaimer: *ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions.*
- **Not for treatment decisions**: the output is a genotype, not antimicrobial susceptibility testing.
- **Audit trail**: raw tool outputs, the command, the environment and checksums are saved with every run.
- **No hallucinated science**: calls, classes and subclasses are copied from the tools; the skill invents none.

## Agent Boundary

The agent (LLM) dispatches and explains. The skill (Python) executes.
The agent must NOT add, drop or re-tier a call, infer an organism outside the scheme map, or translate the genotype into a susceptibility or treatment statement.

## Integration with Bio Orchestrator

**Trigger conditions**: the orchestrator routes here when:
- The input is an assembled bacterial genome FASTA and the request mentions AMR, resistance, MLST, sequence type, plasmids or typing
- The input is an AMRFinderPlus, mlst or PlasmidFinder result file

**Chaining partners**: this skill connects with:
- `busco-assessor`: check assembly completeness before trusting a negative AMR result
- `ncbi-datasets`: download a public assembly, then type it here
- `phylogenetics-builder` / `fastreer`: place isolates in a tree; `result.json` gives the ST and genotype to annotate tips
- `claw-metagenomics`: community-level resistome, as opposed to this single-isolate genotype

## Maintenance

- **Review cadence**: every six months, and on any AMRFinderPlus major release.
- **Staleness signals**: AMRFinderPlus renames output columns again (v4 already did); its organism list changes; `mlst` renames a scheme (`ecoli` became `ecoli_achtman_4`, and `salmonella` appeared beside `senterica_achtman_2`); PlasmidFinder changes its output again (v3 dropped `results_tab.tsv` for JSON). The real outputs in `tests/fixtures/` pin the layouts.
- **Deprecation**: archive to `skills/_deprecated/` if a maintained single tool supersedes the three-tool merge.

## Citations

- [AMRFinderPlus](https://doi.org/10.1038/s41598-021-91456-0); Feldgarden et al. 2021, AMR genes and point mutations
- [mlst](https://github.com/tseemann/mlst); Seemann, sequence typing from assemblies
- [PubMLST](https://doi.org/10.12688/wellcomeopenres.14826.1); Jolley et al. 2018, MLST schemes
- [PlasmidFinder](https://doi.org/10.1128/AAC.02412-14); Carattoli et al. 2014, replicon typing and default thresholds
