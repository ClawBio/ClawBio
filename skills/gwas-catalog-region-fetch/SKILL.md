---
name: gwas-catalog-region-fetch
description: |
  Fetch a region of GWAS summary statistics from the NHGRI-EBI GWAS Catalog
  harmonised collection via tabix-on-FTP. Use when an agent needs GWAS beta /
  SE / p-value for every variant in a window for one specific study (GCST
  accession). Input: accession, chromosome, start, end. Output: harmonised
  TSV slice in canonical format.
license: MIT
metadata:
  skill-author: Aviv Madar
  version: 0.1.0
  domain: bioinformatics
  tags:
    - gwas
    - gwas-catalog
    - region-fetch
    - tabix
    - summary-statistics
    - harmonised
  inputs:
    - name: accession
      type: string
      description: GWAS Catalog study accession (e.g. GCST90269602 for cholesterol-VLDL).
      required: true
    - name: chromosome
      type: string
      description: Chromosome name without `chr` prefix.
      required: true
    - name: start_bp
      type: integer
      description: Region start, 1-based GRCh38.
      required: true
    - name: end_bp
      type: integer
      description: Region end, 1-based GRCh38 (inclusive).
      required: true
  outputs:
    - name: variants
      type: list
      description: Per-variant rows with variant_id, chromosome, position, ref, alt, beta, se, p_value, allele frequencies.
    - name: release
      type: object
      description: GWASCatalogRelease with accession, harmonised file path, fetched_at_utc.
  dependencies:
    - python>=3.10
    - pysam>=0.22
    - pandas>=2.0
    - requests>=2.28
  demo_data:
    - examples/input.json
  endpoints:
    - https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/    # tabix-on-FTP harmonised
  openclaw:
    requires:
      bins:
        - python3
        - tabix
      env:
      config:
    always: false
    emoji: "🧬"
    homepage: https://github.com/ClawBio/ClawBio
    os:
      - darwin
      - linux
    install: |
      pip install pysam pandas requests
    trigger_keywords:
      - gwas region fetch
      - gwas catalog region
      - gwas sumstats slice
      - GCST harmonised tabix
      - GWAS catalog tabix
---

# 🧬 GWAS Catalog Region Fetch

You are **GWAS Catalog Region Fetch**, a specialised ClawBio agent for pulling per-variant disease/trait GWAS summary statistics from the NHGRI-EBI GWAS Catalog harmonised collection. Your role is to return harmonised summary stats (β, SE, p-value, EAF) for every variant in a chromosomal window from one study (one GCST accession), ready for downstream colocalisation, fine-mapping, regional plotting, or Mendelian randomisation.

## Overview

The NHGRI-EBI GWAS Catalog (Sollis 2023 *NAR*) maintains harmonised summary statistics for ~25,000 published GWAS at `https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/<GCST>/harmonised/<GCST>.h.tsv.gz`. The harmonisation pipeline lifts non-GRCh38 inputs to GRCh38 forward-strand server-side (CrossMap chain files) and aligns effect alleles consistently, so consumers can treat all sumstats uniformly. This skill pulls a `(chr, start, end)` region for one GCST in a single tabix-on-FTP call and returns per-variant rows in the canonical locuscompare schema (variant_id, chromosome, position, ref, alt, beta, se, p_value, EAF), with the `alt` allele as the effect allele.

## Trigger

**Fire when** the user (or upstream agent step) wants:

- A regional slice of GWAS summary statistics (β, SE, p-value, EAF) for variants in a chromosomal window from one GCST study.
- Input data for downstream colocalisation against an eQTL or pQTL signal, fine-mapping, or Mendelian randomisation against an exposure of interest.
- Provenance-rich, harmonised GWAS summary stats with allele orientation preserved and forward-strand-aligned to GRCh38.

**Do NOT fire when** the user wants:

- A **point lookup of one variant in one GWAS** - `database-lookup` or `gwas-lookup` is the right skill for single-variant queries.
- **Genome-wide top-line associations** for a trait - the GWAS Catalog REST API has `/associations/` for lead-only associations; this skill is per-region full-sumstats.
- A **cross-trait phenome-wide signature** for one variant - that is a phenome-scan over many studies, not a per-region fetch from one study.
- **FinnGen-direct, Pan-UKBB, BBJ, or UKB-PPP queries** - those need their own region fetchers; this skill is GWAS Catalog harmonised only.
- **Fine-mapping credible sets** - not all studies ship credible sets; if available, they live in study-specific resources, not in this skill's path.
- **Per-trait genetic correlation** (LDSC, mvLMM) - different upstream tooling.

## Scope

**One skill, one task.** This skill fetches one GCST study's regional summary statistics from the GWAS Catalog harmonised collection and writes them as a harmonised TSV plus a provenance manifest. It does NOT do single-variant lookups, cross-study comparisons, raw-upload fetches, FinnGen-direct fetches, or fine-mapping - see "Do NOT fire when" above for the right skills for those tasks.

## Workflow

When an agent asks for a regional GWAS slice from the GWAS Catalog:

1. **Resolve `accession`**: the canonical `GCST########` identifier. Look up via the GWAS Catalog REST API (`https://www.ebi.ac.uk/gwas/rest/api/studies/<GCST>`) or the web UI at `https://www.ebi.ac.uk/gwas/`. The metadata response includes `hasSummaryStats` (must be `true` to fetch), `pubmedId` (citation), and `ancestries[]` (sample sizes per ancestry bucket).
2. **Pick a region**: `(chromosome, start_bp, end_bp)` in 1-based inclusive GRCh38 coordinates. For LocusCompare-style coloc inspection centre on the lead variant ± 500 kb; for "what does this trait look like in the gene's cis-window" centre on the gene TSS ± 1 Mb.
3. **Tabix range fetch**: the skill performs a single byte-range request against `<GCST>.h.tsv.gz` on the EBI GWAS Catalog FTP. The `harmonised/` subdirectory is the canonical path; do NOT swap to the raw upload (Gotcha #1).
4. **Use `hm_*` columns**: the harmonised TSV emits `hm_chrom`, `hm_pos`, `hm_effect_allele`, `hm_other_allele`, `hm_beta`, `hm_se`, `hm_effect_allele_frequency`. The skill maps these to canonical `(variant_id, chromosome, position, ref, alt, beta, se, p_value, maf)` with `alt` as the effect allele.
5. **Write outputs** to `--output <dir>/`: a flat `variants.tsv` (effect-allele-aligned, GRCh38, ALT-effect β), a `manifest.yaml` with provenance (accession, harmonised file path, harmoniser pipeline version where surfaced, n_variants, source URL, fetched-at UTC timestamp), and a `report.md` human-readable summary.

## CLI Reference

```bash
# Standard usage with a config file
python skills/gwas-catalog-region-fetch/gwas_catalog_region_fetch.py \
    --input <config.json> --output <output_dir>

# Bundled demo (cholesterol-in-medium-VLDL GWAS at the SORT1 locus)
python skills/gwas-catalog-region-fetch/gwas_catalog_region_fetch.py \
    --demo --output /tmp/sort1_vldl_demo

# Via ClawBio runner
python clawbio.py run gwas-region --input <config.json>
python clawbio.py run gwas-region --demo
```

Config schema (JSON or YAML):

```json
{
  "accession": "GCST90269602",
  "chromosome": "1",
  "start_bp": 108774968,
  "end_bp": 109774968
}
```

Bundled biology demos in `examples/`:

- `sort1_cholesterol_vldl.json` - Musunuru 2010 1p13.3 LDL/CHD locus; pairs with the SORT1 ge-eQTL on the exposure side.
- `il6r_crp.yaml` - IL6R × CRP at chr1:154425508; canonical IVW MR demo.
- `tcf7l2_hba1c.json` - TCF7L2 × HbA1c; type-2 diabetes locus.

## Example Output

Running `--demo` (SORT1 × cholesterol-VLDL):

```
info: using bundled demo sort1_cholesterol_vldl.json
gwas-catalog-region-fetch: 2914 variants -> /tmp/sort1_vldl_demo/variants.tsv
  source: GCST90269602 (cholesterol in medium VLDL)
```

`<output_dir>/manifest.yaml`:

```yaml
skill: gwas-catalog-region-fetch
version: 0.1.0
accession: GCST90269602
trait_label: cholesterol in medium VLDL
region:
  chromosome: '1'
  start_bp: 108774968
  end_bp: 109774968
n_variants: 2914
release:
  accession: GCST90269602
  harmonised_path: GCST90269602/harmonised/GCST90269602.h.tsv.gz
  source_url: https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90269602/harmonised/GCST90269602.h.tsv.gz
  fetched_at_utc: '2026-05-09T11:42:06Z'
outputs:
  variants_tsv: variants.tsv
```

`<output_dir>/variants.tsv` (first three rows shown):

```
variant_id              chromosome  position_bp  allele_a  allele_b  beta        se        p          maf       study_id
1_108774974_TCTAC_T     1           108774974    TCTAC     T          0.0123     0.0089     0.165     0.171     GCST90269602
1_108775337_C_T         1           108775337    C         T         -0.0205     0.0095     0.031     0.314     GCST90269602
1_108775606_G_T         1           108775606    G         T          0.0089     0.0188     0.636     0.072     GCST90269602
```

`<output_dir>/report.md`:

```markdown
# gwas-catalog-region-fetch report

- **Accession:** `GCST90269602`
- **Trait:** cholesterol in medium VLDL
- **Region:** chr1:108,774,968-109,774,968
- **Variants returned:** 2914
- **Lead variant:** `1_109274968_G_T` (rs646776), p ≈ 1e-50, β ≈ -0.27 (per Musunuru 2010 inverse SORT1↑→LDL↓ biology)
- **Output TSV:** variants.tsv
```

## Gotchas

1. **Use the `harmonised/` subdirectory, not the raw upload path.** `<GCST>/harmonised/<GCST>.h.tsv.gz` is forward-strand-aligned to GRCh38 with consistent allele orientation. The raw upload (one directory up) can be on GRCh37 with study-specific allele conventions, and is NOT what this skill fetches. Do NOT swap to the raw path.

2. **`hm_*` columns are the canonical fields.** The harmoniser emits `hm_chrom`, `hm_pos`, `hm_effect_allele`, `hm_other_allele`, `hm_beta`, `hm_se`, `hm_effect_allele_frequency`. Use these, not the raw-upload columns. The skill's manifest preserves both for traceability.

3. **Per-study release lag.** GWAS Catalog mirrors a study some time after the upstream release: typically 2-6 months for FinnGen R12 phenotypes; longer for studies that go through deposit-and-curate. If the user references a phenotype that should be in OT, verify presence via the GWAS Catalog API metadata before assuming an arbitrary GCST is fetchable.

4. **Some studies do not have summary statistics deposited at all.** Older or smaller GWAS may have only top-line lead associations but no full sumstats. The GWAS Catalog API exposes `hasSummaryStats` per study; check it before invoking. The fetcher raises `GWASCatalogFetchError` when the harmonised TSV or its `.tbi` is missing on FTP; caller decides whether to fall back (e.g. download the whole file + tabix-index locally - see `references/harmonised_pipeline.md`).

5. **Palindromic SNPs at MAF near 0.5 are dropped by the harmoniser.** A/T and G/C variants with EAF in [0.45, 0.55] cannot be reliably oriented across studies, so the harmoniser excludes them. Expect some variants present in OT credible sets to be missing from the harmonised file. Surface the count to the user when it materially affects the analysis.

6. **β is reported on the ALT allele.** Do NOT compare effect sizes across studies without explicit allele harmonisation. The skill preserves `ref` / `alt` columns; downstream tools (e.g., TwoSampleMR `harmonise_data`) flip signs when alleles are swapped. Cross-study sign-flip risk is real (`references/effect_allele_harmonisation.md`).

## Safety

**Not for clinical decisions.** This skill returns research-grade GWAS summary statistics from public databases. Do not use the output for direct clinical decision-making, diagnosis, or treatment selection without independent validation by a qualified clinician.

**Effect sizes can include winner's-curse bias.** Variants discovered in the same GWAS that produced the summary stats have inflated effect-size estimates. Downstream causal-effect estimation (Mendelian randomisation) should use independent instruments or apply winner's-curse correction.

**Effect estimates may not generalise across ancestries.** The GWAS Catalog records each study's primary ancestry (`ancestries[]` block in the REST metadata); effect sizes from a single-ancestry cohort should not be assumed to apply trans-ancestrally without explicit harmonisation and validation.

## Agent Boundary

The skill returns harmonised GWAS summary statistics (β, SE, p-value, EAF) for variants in a chromosomal window from one GCST study. The agent should:

- **Use the output as input to colocalisation, fine-mapping, or Mendelian randomisation tooling.** These are the appropriate downstream methods for inferring causal effects.
- **NOT make causal claims directly from a single GWAS p-value.** Association is not causation. Causal interpretation requires colocalisation or MR analysis with proper instrumental-variable assumptions.
- **NOT cherry-pick variants by p-value alone.** Statistical inference requires the full window context (credible set, posterior inclusion probabilities, joint conditional analyses).
- **NOT compare effect sizes across studies without harmonising effect alleles.** Cross-study comparison requires a step like TwoSampleMR's `harmonise_data`. Sign-flip risk is real for swapped alleles and palindromic ambiguity.
- **Surface the GCST id, trait label, sample size (cases / controls if binary, total N if continuous), ancestry, and consortium (when present)** alongside any β / p-value the agent quotes. Per the user-friendly enum-expansion rule (`AGENTS.md`), expand each field: `study GCST90269602; trait cholesterol in medium VLDL; ancestry European; N=44,000` (not just `GCST90269602`).
- **NOT report a binary-trait OR (odds ratio) as if it were a continuous-trait β.** The harmonised file's `hm_beta` for binary traits is log(OR); the manifest carries the trait type so the agent can disambiguate.

## Citations

- Sollis et al. (2023). *The NHGRI-EBI GWAS Catalog: knowledgebase and deposition resource.* Nucleic Acids Res 51, D977-D985. doi:10.1093/nar/gkac1010
- Per-study original GWAS publication (cited from the GWAS Catalog metadata `pubmedId` field).
