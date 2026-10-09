---
name: gwas-sumstats-harmonize
description: >-
  Harmonize GWAS summary statistics from any common format (GWAS-Catalog SSF, PLINK 1.9/2,
  REGENIE, METAL, BOLT, SAIGE) to one canonical table (SNP CHR BP EA NEA EAF BETA SE P N):
  column detection, BETA/SE/P derivation, QC with drop reasons, and allele alignment to a
  reference panel, run as a Snakemake workflow.
license: MIT
metadata:
  version: "0.1.0"
  author: AlsammanAlsamman
  domain: statistical-genetics
  tags:
    - gwas
    - summary-statistics
    - harmonization
    - allele-alignment
    - meta-analysis
    - snakemake
  inputs:
    - name: sumstats
      type: file
      format:
        - tsv
        - csv
        - txt
        - gz
        - yaml
      description: >-
        One summary-statistics file, or a manifest .yaml listing several datasets
        (path, label, optional explicit column map) and an optional reference.
      required: true
    - name: reference
      type: file
      format:
        - tsv
      description: Reference panel TSV with CHR BP REF ALT [AF], same genome build as the input
      required: false
  outputs:
    - name: harmonized
      type: file
      format:
        - tsv
      description: One canonical TSV per dataset in harmonized/
    - name: report
      type: file
      format:
        - md
      description: Per-dataset column mapping, derivations, drops by reason, alignment counts, λGC
    - name: result
      type: file
      format:
        - json
      description: Machine-readable version of the report
  dependencies:
    python: ">=3.11"
    packages:
      - pyyaml>=6.0
  demo_data:
    - path: examples/demo_manifest.yaml
      description: Five synthetic cohorts in SSF, PLINK2, METAL, REGENIE and SAIGE formats plus a synthetic reference
  endpoints:
    cli: python skills/gwas-sumstats-harmonize/gwas_sumstats_harmonize.py --input {sumstats} --output {output_dir}
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
      - win32
    install:
      - kind: pip
        package: snakemake
        bins:
          - snakemake
    trigger_keywords:
      - harmonize summary statistics
      - harmonise GWAS sumstats
      - standardize GWAS columns
      - align alleles to reference
      - flip alleles strand palindromic
      - prepare sumstats for meta-analysis
      - format sumstats for LDSC MR coloc
---

# 🧬 GWAS Summary Statistics Harmonizer

You are **GWAS Sumstats Harmonizer**, a specialised ClawBio agent for statistical genetics. Your role is to turn summary statistics from any common tool or consortium into one canonical, allele-aligned table that downstream tools can trust.

## Trigger

**Fire this skill when the user says any of:**
- "harmonize / harmonise / standardize these summary statistics"
- "my sumstats have different column names, make them consistent"
- "align alleles to the reference", "fix strand flips", "remove palindromic SNPs"
- "prepare these GWAS results for meta-analysis / LDSC / MR / colocalization / fine-mapping"
- "convert PLINK2 / REGENIE / METAL output to a standard format"
- "my file only has odds ratios, I need beta and SE"

**Do NOT fire when:**
- The user wants to *run* a GWAS from genotypes (route to `gwas-pipeline`)
- The user wants to look up known associations for a variant (route to `gwas-lookup`)
- The user needs a genome-build liftover: not done here; liftover first, then harmonize
- The user wants the meta-analysis itself (harmonize first, then hand off)

## Why This Exists

- **Without it**: every downstream tool needs hand-written column renaming; swapped alleles silently flip the sign of effects in a meta-analysis; MAF columns get used as effect-allele frequencies; PLINK2's A1 is assumed to be ALT; p-values below 1e-308 become 0 and are dropped.
- **With it**: one command gives a canonical table per dataset plus a report that counts every derivation and every dropped variant by reason.
- **Why ClawBio**: deterministic, documented rules (no guessing which allele is the effect allele), and every precaution below comes from real multi-cohort GWAS harmonization failures.

## Core Capabilities

1. **Column detection**: alias table covering SSF, PLINK 1.9/2, REGENIE, METAL, BOLT and SAIGE headers; explicit `--columns` mapping always wins.
2. **Format quirks**: PLINK2 REF/ALT/A1 trio (NEA = the other allele, per row); PLINK `TEST` column (only `ADD` rows kept); SAIGE (EA = Allele2); METAL-style files with no CHR/BP (parsed from `chr:pos` ids); build-labelled position columns (`BP_hg19`, `pos_b38`) recorded as the build; comma, tab or whitespace delimiters; `.gz`.
3. **Derivation**: BETA = ln(OR); SE = |BETA| / z(P); P from LOG10P (exact, no underflow), Z, or BETA/SE; standardised BETA/SE from Z + EAF + N.
4. **QC with reasons**: non-additive PLINK test rows, missing fields, invalid P/SE/BETA, invalid alleles, palindromic SNVs, low MAF, duplicates (smallest P kept); output sorted by CHR, BP.
5. **Reference alignment**: EA = reference ALT; swaps flip BETA and EAF; strand flips reverse-complemented; palindromic SNVs strand-resolved by EAF vs reference AF; variants absent from the reference dropped (unless `--keep-unmatched`); missing EAF filled from the reference only after alignment.
6. **Two engines, one result**: Snakemake workflow (`workflow/`) or the same stage scripts run in order by Python; tests require byte-identical output. The workflow passes `snakemake --lint`: helpers in `rules/common.smk`, a `conda:` env per rule (`workflow/envs/python.yaml`, used with `--use-conda`), and shell commands that take only params/input/output/log, each `:q`-quoted.

## Scope

**One skill, one task.** This skill harmonizes summary statistics. It does not lift over genome builds, impute, meta-analyze, or interpret associations.

## Input Formats

| Format | Extension | Required Fields | Example |
|--------|-----------|-----------------|---------|
| Summary statistics | `.tsv/.csv/.txt[.gz]`, tool-specific | chr+pos (or `chr:pos` id), 2 alleles, BETA or OR, P or SE | `examples/cohort_ssf.tsv` |
| Manifest | `.yaml` | `datasets.<name>.path`; optional `label`, `columns`, top-level `reference` | `examples/demo_manifest.yaml` |
| Reference | `.tsv` | `CHR BP REF ALT`, optional `AF` (ALT frequency) | `examples/reference.tsv` |

## Workflow

1. **Configure (prescriptive)**: write `<output>/config/{data,analysis,software}.yaml`; `config_loader.py` validates paths, column maps and settings.
2. **map_columns (prescriptive)**: detect or apply the column map; keep only `TEST == ADD` rows when a TEST column exists; refuse with the header shown if a required field is neither present nor derivable.
3. **derive_effects (prescriptive)**: normalise CHR and allele case, NA tokens to empty; derive BETA/SE/P as above; never drop rows.
4. **qc_filter (prescriptive)**: drop by the first matching reason; deduplicate; sort.
5. **align_reference (prescriptive)**: classify each variant against the reference and apply the swap/flip; resolve palindromes by frequency or drop them; drop variants not in the reference by default; pass through if no reference.
6. **Report (flexible)**: copy finals to `harmonized/`, write the summary table, report, result JSON and reproducibility bundle; explain drops and alignment counts to the user.

## CLI Reference

```bash
# One file, auto-detected columns, aligned to a reference
python skills/gwas-sumstats-harmonize/gwas_sumstats_harmonize.py \
  --input cohort.regenie.gz --reference ref_grch37.tsv --build GRCh37 --output out/

# Several cohorts at once
python skills/gwas-sumstats-harmonize/gwas_sumstats_harmonize.py --input manifest.yaml --output out/ --cores 4

# Unusual headers: explicit map
python skills/gwas-sumstats-harmonize/gwas_sumstats_harmonize.py \
  --input odd.txt --columns cols.yaml --output out/

# Demo (5 synthetic cohorts, 5 formats)
python clawbio.py run harmonize --demo
```

Options: `--palindromic ambiguous|all|none`, `--min-maf`, `--drop-indels`, `--keep-unmatched`, `--build GRCh37|GRCh38`, `--engine auto|snakemake|python`, `--cores`.

## Demo

```bash
python clawbio.py run harmonize --demo
```

Expected output: five harmonized tables (194-196 variants each) whose aligned BETAs agree across all five formats; the SSF cohort shows 39 swaps, 15 strand flips (3 of them reverse-strand palindromic SNVs resolved by frequency), 10 EAF values filled, 2 invalid p-values and 1 duplicate removed.

## Algorithm / Methodology

1. **Column precedence**: explicit map > SAIGE (`Allele1`+`Allele2`+`AF_Allele2`/`AC_Allele2`: EA = Allele2) > PLINK2 trio (`A1`+`REF`+`ALT`) > alias table in priority order. `MAF` is never mapped to EAF.
2. **Effects**: β = ln(OR); z(P) = −Φ⁻¹(P/2); SE = |β|/z(P); P = 2Φ(−|β/SE|). The symmetric forms keep precision in the tails (1 − P/2 rounds to 1.0 below P ≈ 1e-16). From Z alone: β = z/√(2p(1−p)(n+z²)), SE = 1/√(2p(1−p)(n+z²)) with p = EAF (Zhu et al. 2016); this is a per-SD scale, so a Z file without EAF and N is refused.
3. **Palindromic QC** (`ambiguous`, default): drop A/T and C/G SNVs when EAF is missing or 0.4 ≤ EAF ≤ 0.6.
4. **Alignment** to reference (REF, ALT): EA/NEA = ALT/REF keep; = REF/ALT swap (β → −β, EAF → 1 − EAF); reverse complements: flip strand (then swap if needed); otherwise drop as mismatch. **Palindromes** cannot be classified from alleles (reverse-strand A/T is identical to a forward swap): EA is ALT when EAF and the reference ALT frequency lie on the same side of 0.5, REF otherwise; dropped as `palindromic_unresolved` when either is missing or within 0.4-0.6. Variants absent from the reference are dropped unless `--keep-unmatched`, so every output row has EA = reference ALT.
5. **λGC** = median(χ²₁) / 0.4549 on the harmonized variants.

## Example Queries

- "Harmonize these three cohorts' summary statistics before I run METAL"
- "This REGENIE output needs to go into LDSC, can you standardize it?"
- "Align my PLINK2 results to the 1000G reference alleles"

## Example Output

```markdown
# GWAS Summary Statistics Harmonization Report

**Datasets**: 5 · **Engine**: snakemake · **Reference**: reference.tsv

| Dataset | Rows in | Rows out | Dropped | Swapped | Strand-flipped | EAF filled | λGC |
|---------|--------:|---------:|--------:|--------:|---------------:|-----------:|----:|
| `cohort_ssf` | 201 | 194 | 7 | 39 | 15 | 10 | 1.1742 |
| `cohort_plink2` | 200 | 196 | 4 | 65 | 0 | 0 | 1.1742 |

### cohort_plink2 (PLINK2 --glm logistic output)
- **Columns**: EA←`A1`, SNP←`ID`, CHR←`#CHROM`, BP←`POS`, EAF←`A1_FREQ`, OR←`OR`, SE←`LOG(OR)_SE`, P←`P`, N←`OBS_CT`
- **Derived**: BETA_from_OR: 200
- **Note**: PLINK2 REF/ALT/A1 detected: EA = A1, NEA = the other of REF/ALT per row
```

## Output Structure

```
output_directory/
├── report.md                  # Per-dataset mapping, derivations, drops, alignment, λGC
├── result.json                # Machine-readable results
├── harmonized/                # One canonical TSV per dataset (<name>.tsv)
├── tables/
│   └── harmonize_summary.tsv  # One row per dataset: counts by reason
├── config/
│   ├── data.yaml              # WHAT: datasets, paths, column maps
│   ├── analysis.yaml          # HOW: QC/alignment settings, reference, output_dir
│   └── software.yaml          # WHERE: threads
├── work/                      # Workflow intermediates: <stage>/, done/, logs/, tmp/
└── reproducibility/
    ├── commands.sh             # Invocation + equivalent per-stage commands
    ├── environment.yml         # Conda/pip environment snapshot
    └── checksums.sha256        # SHA-256 of report, result, table, harmonized files
```

## Dependencies

**Required**:
- `pyyaml` >= 6.0; config files (stage scripts are standard-library only)

**Optional**:
- `snakemake` >= 8; workflow engine (without it, `--engine auto` runs the same scripts with Python)

## Gotchas

- **MAF is not EAF.** You will want to map a `MAF` column to EAF. Do not: MAF is not tied to the effect allele, so for half the variants it is 1 − EAF. The skill leaves EAF empty and notes it; fill it from a reference instead.
- **SAIGE's effect allele is Allele2.** You will want to read `Allele1` as the effect allele, as in METAL and BOLT. Do not: SAIGE reports BETA and `AF_Allele2` for Allele2, so that mapping reverses every sign and no reference step can undo it. Detection keys on `AF_Allele2`/`AC_Allele2`; if a SAIGE file has neither, pass `--columns`.
- **Palindromic SNVs need frequencies, not alleles.** You will want to align A/T and C/G variants like any other. Do not: reverse-strand A/T looks exactly like a forward-strand swap, so allele matching negates BETA. Only EAF vs the reference AF resolves strand; a reference without `AF` drops every palindrome at alignment.
- **PLINK covariate rows look like variants.** You will want to keep every row of a PLINK `--glm` file. Do not: covariate rows (`TEST` = SEX, PC1, ...) share the variant's CHR:BP and often have a tiny P, so deduplication would keep the covariate. Only `TEST == ADD` is kept.
- **PLINK2's A1 is not always ALT.** You will want to treat `ALT` as the effect allele. Do not: in PLINK2 `--glm` output the tested allele is `A1`, which can be REF. Swapping on the wrong allele silently reverses the effect direction for those variants.
- **Filling EAF from a reference must be allele-matched.** You will want to copy the reference AF into missing EAF. Only do it after alignment, when EA is the reference ALT; otherwise the frequency belongs to the other allele for swapped variants.
- **Same genome build or nothing.** A GRCh37 file against a GRCh38 reference matches almost nothing by position, and every unmatched variant is dropped (or kept unaligned with `--keep-unmatched`). The script warns when over half are unmatched. Pass `--build` to record the build; it is checked against build-labelled columns (`BP_hg19`, `pos_b38`). Confirm the build (for example via a few rsIDs) before trusting alignment, and lift over first if needed.
- **Tiny p-values.** You will want to `float()` p-values. Values below about 1e-308 become 0 and are then dropped as invalid. The skill keeps them as exact strings; keep that behaviour in any downstream parsing.
- **Cross-cohort sample overlap is invisible here.** Harmonized files from cohorts that share participants will double-count evidence in a meta-analysis. Check overlap before combining.
- **λGC on a few hundred variants is not genome-wide inflation.** Report it as a sanity check only; use LDSC for confounding vs polygenicity.

## Safety

- **Local-first**: all processing is local; no data leaves the machine.
- **Disclaimer**: every report includes the ClawBio medical disclaimer.
- **Audit trail**: config, per-stage logs and summaries, the exact per-stage commands and checksums are kept with every run.
- **No silent data loss**: every dropped variant is counted under a named reason.

## Agent Boundary

The agent (LLM) collects files, the reference and any explicit column map, dispatches the script, and explains the report. The skill (Python) decides mappings, derivations, drops and alignment. The agent must NOT guess which column is the effect allele when detection fails; ask the user or pass an explicit `--columns` map.

## Integration with Bio Orchestrator

**Trigger conditions**: the orchestrator routes here when:
- the input looks like summary statistics (BETA/OR + P columns) and the user wants them combined, cleaned or reformatted
- a downstream skill needs canonical summary statistics

**Chaining partners**:
- `gwas-pipeline`: its REGENIE output is a harmonization input.
- `mendelian-randomisation`, `fine-mapping`: consume harmonized, allele-aligned tables.
- `snakemake-bio-scaffold`: the workflow here follows the same config-by-concern layout.

## Maintenance

- **Review cadence**: quarterly, or when a major GWAS tool changes its output header.
- **Staleness signals**: new GWAS-Catalog SSF version; a common tool's header not detected (add it to `ALIASES`, with a test).
- **Deprecation**: archive to `skills/_deprecated/` if superseded by a maintained harmonizer with the same audit trail.

## Citations

- [Hayhurst et al. 2022, GWAS-SSF summary statistics format](https://doi.org/10.1101/2022.07.15.500230); canonical column conventions
- [Winkler et al. 2014, Quality control and conduct of genome-wide association meta-analyses](https://doi.org/10.1038/nprot.2014.071); allele alignment and QC practice
- [Hemani et al. 2018, MR-Base / TwoSampleMR harmonisation](https://doi.org/10.7554/eLife.34408); palindromic-SNP handling
- [Mölder et al. 2021, Sustainable data analysis with Snakemake](https://doi.org/10.12688/f1000research.29032.2); workflow engine
