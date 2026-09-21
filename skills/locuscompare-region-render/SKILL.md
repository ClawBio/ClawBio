---
name: locuscompare-region-render
description: |
  Render a 4-panel regional LocusCompare diagnostic for one (lead variant,
  exposure study, outcome study) tuple - overlays GWAS Manhattan, QTL Manhattan,
  GENCODE gene track, and cross-trait scatter colored by LD r². Use when an
  agent needs visual confirmation that two GWAS / QTL signals share the same
  causal variant (the Liu 2019 LocusCompare convention). Inputs: lead variant +
  two pre-fetched harmonised sumstats slices (or eQTL Catalogue / GWAS Catalog
  identifiers for bundled fetch). Output: PNG + JSON manifest.
license: MIT
metadata:
  skill-author: Aviv Madar
  version: 0.1.0
  domain: bioinformatics
  tags:
    - regional-plot
    - coloc
    - eqtl
    - gwas
    - locuscompare
    - mendelian-randomisation
    - target-validation
    - drug-discovery
  inputs:
    - name: input_file
      type: file
      description: |
        Run config (JSON or YAML) describing the lead variant + region,
        exposure and outcome data sources (either pre-fetched harmonised TSV
        paths OR a bundled-fetcher block referencing eQTL Catalogue / GWAS
        Catalog), and optional LD + gene-track sources. See INPUT_SCHEMA.md
        for the canonical sumstats-slice TSV column spec.
      required: true
  outputs:
    - name: report
      type: file
      description: Markdown report describing the rendered plot + caveats
    - name: figure
      type: file
      description: 4-panel regional LocusCompare PNG (exposure Manhattan + outcome Manhattan + gene track + cross-trait scatter)
    - name: manifest
      type: file
      description: Reproducibility manifest (YAML) with source releases, LD panel id, plink version, n_pairs, palindromic-exclusion count
  dependencies:
    - python>=3.10
    - numpy>=1.24
    - scipy>=1.10
    - pandas>=2.0
    - matplotlib>=3.7
    - pysam>=0.22
    - pyyaml>=6.0
    - pydantic>=2.0
    - requests>=2.28
  demo_data:
    - examples/01_synthetic_demo/config.json
    - examples/02_eqtl_catalogue_x_gwas_catalog/config.yaml
  endpoints:
    - https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/         # eQTL Catalogue tabix-on-FTP
    - https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/     # GWAS Catalog harmonised tabix-on-FTP
    - https://ftp.1000genomes.ebi.ac.uk/                                # 1000G phased VCFs (LD compute)
    - https://rest.ensembl.org/                                         # GENCODE gene-track REST
  openclaw:
    requires:
      bins:
        - python3
        - tabix
        - plink
      env:
      config:
    always: false
    emoji: "🎯"
    homepage: https://github.com/ClawBio/ClawBio
    os:
      - darwin
      - linux
    install: |
      pip install 'numpy>=1.24' 'scipy>=1.10' 'pandas>=2.0' 'matplotlib>=3.7' \
                  'pysam>=0.22' 'pyyaml>=6.0' 'pydantic>=2.0' 'requests>=2.28'
      # plus a system plink 1.9 binary (brew install plink on macOS, apt-get install plink1.9 on Linux)
      # Pins mirror the `dependencies:` block above; keep both in sync.
    trigger_keywords:
      - regional locuscompare
      - regional coloc plot
      - two-trait coloc visualization
      - GWAS eQTL overlay
      - liu 2019 locuscompare convention
      - colocalization plot for a lead variant
      - regional plot for coloc lead
---

# 🎯 LocusCompare Regional Diagnostic

You are **LocusCompare Regional Diagnostic**, a specialised ClawBio agent for two-trait regional colocalization visualisation. Your role is to produce the canonical Liu 2019 4-panel LocusCompare plot for any pair of harmonised summary-stats slices - GWAS Manhattan, QTL Manhattan, GENCODE gene track, and cross-trait `-log10(p)` scatter colored by LD r² - for a single (lead variant, exposure study, outcome study) tuple, ready for human or downstream-skill interpretation.

## Overview

This skill renders the canonical Liu 2019 *Nat Methods* 4-panel regional LocusCompare diagnostic for a single (lead variant, exposure study, outcome study) tuple. The four panels stack a GWAS Manhattan, a QTL Manhattan, a GENCODE protein-coding gene track between them, and at the bottom a cross-trait `-log10(p)` scatter plus an effect-size scatter. Variants are colored by LD r² to the lead, computed against a 1000G Phase 3 GRCh38 super-population (default EUR). Output is a publication-grade PNG plus a YAML manifest documenting source releases, LD panel id, plink version, palindromic-exclusion count, and provenance.

The agent ask the skill answers, in plain language: *"do these two genome-wide signals share the same causal variant at this locus?"* The answer is visual and auditable. The skill does NOT compute coloc.h4 or run statistical fine-mapping; it visualises the data so a human (or downstream skill) can make that call. Pass `coloc.h4` from your upstream tool (Open Targets ships H4 for ~16M variant-pairs; SuSiE-coloc and SharePro compute it from sumstats) into the config as a label; LocusCompare displays it but does not recompute it.

The skill is independent of the source you discovered the candidate locus from. Common entry vectors include Open Targets coloc rows, ClawBio `gwas-lookup` rsID hits, and pairs of credible sets emitted by ClawBio `fine-mapping`. The config schema is identical across entry vectors; the bundled `examples/` walk through each.

## Trigger

**Fire when** the user (or upstream agent step) wants:

- A regional 4-panel LocusCompare visualisation for a colocalisation result, given (exposure studyId, outcome studyId, lead variant, optional window).
- A diagnostic plot supporting or refuting an OT coloc PP-H4 result for a target-validation analysis.
- Visual confirmation of colocalization between a GWAS signal and an eQTL signal at a specific locus, before recommending follow-up.
- The visual companion to a `mendelian-randomisation` IVW run.
- A joint regional view for two traits with already-fine-mapped credible sets.

**Do NOT fire when** the user wants:

- **Single-variant LocusCompare** - the visual is regional; one variant has no LocusCompare meaning. Use a 2-variant association lookup instead.
- **Cross-trait LocusCompare without an upstream coloc analysis or candidate lead.** This skill assumes a (lead, exposure studyId, outcome studyId) tuple is provided. Without a lead, chain to `gwas-lookup` first; do not invent one.
- **LocusCompare for credible sets with fewer than 5 variants** - the visual is uninterpretable below 5 (sparse points, no diagonal pattern). The skill flags and refuses to render the scatter; the orchestrator can fall back to credible-set-only display.
- **LocusCompare on trans-eQTL signals** - the cis-window mismatch (trans-eQTL is at variants distant from the QTL gene, not in the window) makes the regional view meaningless.
- **Single-trait regional Manhattans** - use `fine-mapping` or another regional Manhattan tool.
- **Numeric coloc statistics only** - use `coloc.abf`, SuSiE-coloc, SharePro, or take H4 from Open Targets.
- **Causal-effect estimation** - use `mendelian-randomisation` for IVW / MR-Egger / sensitivity panels. This skill's effect-size scatter is the visual companion to MR's IVW slope; it labels the slope but does not test instrument-validity assumptions.
- **Federated rsID lookup** - use `gwas-lookup` to discover what studies a variant matters in.
- **Run a new GWAS** - use `gwas-pipeline` to run PLINK2 / REGENIE.
- **Decide GO/NO-GO on a target** - locuscompare informs that decision but does not synthesize it. See `target-validation-scorer` for that.

## Scope

**One skill, one task.** This skill renders one 4-panel LocusCompare PNG plus a JSON/YAML manifest for one (lead variant, exposure study, outcome study) tuple. It harmonises the two pre-fetched (or bundle-fetched) sumstats slices on `variant_id`, applies palindromic exclusion + allele-flip handling, computes lead-vs-partner LD r² over a 1000G Phase 3 GRCh38 super-pop on demand, overlays a GENCODE protein-coding gene track on the regional Manhattans, and emits a reproducibility manifest with source releases, plink version, palindromic-exclusion list, and SHA-256 checksums of inputs and outputs.

It does NOT compute coloc.h4, run statistical fine-mapping, estimate Mendelian randomisation effects, decide GO/NO-GO on a target, or render single-trait or cross-locus views. Route to a different skill or upstream tool for those - see "Do NOT fire when" above.

## Workflow

When an agent asks for a regional LocusCompare for a (lead, exposure, outcome) tuple:

1. **Validate**: parse the config (JSON or YAML); validate against the embedded pydantic schema; check input files exist and conform to `INPUT_SCHEMA.md`. The lead variant must be in `chr_pos_ref_alt` GRCh38 form; the window defaults to ±500 kb of lead.
2. **Fetch (optional)**: for each side using a bundled fetcher (`source: eqtl_catalogue` or `source: gwas_catalog`), tabix-fetch the region from the source's FTP via the corresponding sibling skill (`eqtl-catalogue-region-fetch` or `gwas-catalog-region-fetch`). Cache to `~/.clawbio/locuscompare_cache/<source>/`.
3. **LD compute (optional)**: pull the 1000G region VCF (~5-50 MB) via the `ld-1000g-region-compute` sibling skill; run `plink --r2` (plink 1.9) against the lead variant; cache to `~/.clawbio/locuscompare_cache/1000g/`. Falls back to grey coloring when the LD source is unavailable (caveat surfaced in manifest).
4. **Gene track (optional)**: fetch the locus's GENCODE GTF region via Ensembl REST; parse exon / intron / strand for protein-coding genes (default; configurable via `gene_track.biotypes`). Cache to `~/.clawbio/locuscompare_cache/gencode/`.
5. **Harmonise + render + emit**: join exposure × outcome on `variant_id`; flip `beta_outcome` for swapped alleles; exclude palindromic A/T or G/C variants from the scatter (count surfaced in manifest; configurable). Render the 4-panel matplotlib PNG. Emit `manifest.yaml` with full provenance, `report.md` human-readable summary, `tables/pairs.tsv` with all joined pairs, `tables/palindromic_excluded.tsv`, and a `reproducibility/` subdir with `commands.sh`, `environment.yml`, `input_config.yaml` (resolved defaults), and `checksums.sha256`.

## CLI Reference

```bash
# Bundled offline demo (synthetic 200-variant locus)
python skills/locuscompare-region-render/locuscompare_region_render.py \
    --demo --output /tmp/locuscompare_demo

# Real data (config drives everything)
python skills/locuscompare-region-render/locuscompare_region_render.py \
    --input config.yaml --output runs/sort1_vldl/

# Via the ClawBio runner
python clawbio.py run locuscompare-region --demo
python clawbio.py run locuscompare-region --input config.yaml --output runs/
```

Config schema (JSON or YAML); two ways to supply each side:

**Pre-fetched harmonised TSV** (canonical):

```yaml
lead:
  variant_id: "1_109274968_G_T"
  window_bp: 1000000
exposure:
  trait_label: "SORT1 expression - minor salivary gland"
  sumstats_path: "path/to/exposure_slice.tsv"
outcome:
  trait_label: "cholesterol in medium VLDL"
  sumstats_path: "path/to/outcome_slice.tsv"
ld:
  super_pop: "EUR"
```

**Bundled fetcher** (eQTL Catalogue or GWAS Catalog):

```yaml
lead:
  variant_id: "1_109274968_G_T"
  window_bp: 1000000
exposure:
  trait_label: "SORT1 expression - minor salivary gland"
  fetch:
    source: "eqtl_catalogue"
    dataset_id: "QTD000276"
    molecular_trait_id: "ENSG00000134243"
outcome:
  trait_label: "cholesterol in medium VLDL"
  fetch:
    source: "gwas_catalog"
    accession: "GCST90269602"
ld:
  super_pop: "EUR"
```

Mix-and-match (one side pre-fetched, other side via fetcher) is supported. See `INPUT_SCHEMA.md` for the canonical TSV column spec; `examples/recipes/` for harmonisation scripts converting other sources (FinnGen, Pan-UKBB, UKB-PPP, GTEx v10) into this format.

Optional `lead.rs_id: "rsXXXXXX"` (string) propagates into the emitted manifest + report alongside the canonical `variant_id`. Upstream agents (e.g. ai_scientist `coloc_with_mr` workflow) resolve and populate this automatically from OT / dbSNP; manual configs may set it for human readability. The orchestrator joins on `variant_id` (chr_pos_ref_alt), not rs_id; the rs_id is human-readability metadata only.

Bundled examples in `examples/` (each runnable via `--demo <NAME>` once any required upstream data is in place):

- `01_synthetic_demo/` - offline 200-variant synthetic locus; no network required
- `02_eqtl_catalogue_x_gwas_catalog/` - SORT1 × cholesterol-VLDL canonical demo (Musunuru 2010)
- `03_open_targets_followup/` - entry from an OT coloc row
- `04_gwas_lookup_followup/` - entry from a `gwas-lookup` rsID hit
- `05_sqtl_sort1_liver_txrev/` - SORT1 sQTL in GTEx liver (txrev quant method)
- `06_sceqtl_sort1_onek1k_cd14_mono/` - SORT1 sceQTL in OneK1K CD14+ monocytes
- `07_pqtl_sort1_ukbppp_eur/` - SORT1 plasma sortilin via UKB-PPP (pQTL dispatch path)

The `examples/chains/` and `examples/recipes/` subdirectories are documentation patterns, not numbered demos: `chains/finemapping_chain/` documents the `fine-mapping` skill -> locuscompare chain (requires `fine-mapping` to have run first; not in `--list-demos`); `recipes/` ships harmonisation scripts that convert FinnGen / Pan-UKBB / UKB-PPP / GTEx v10 native formats into the canonical INPUT_SCHEMA TSV for use via `sumstats_path:`.

## Example Output

The focal QTL gene (carried on the spec via `exposure.gene_symbol`, or auto-resolved from the OT exposure studyId in OT-followup configs) is rendered bold and tinted on the gene track so it anchors the visual in 10-30-gene windows. Other genes render in the default grey-blue.

A real run on the SORT1 × cholesterol-VLDL canonical demo (`examples/02_eqtl_catalogue_x_gwas_catalog/config.yaml`) produces:

`<output_dir>/report.md`:

```markdown
# locuscompare-region-render report

- **Lead variant:** `1_109274968_G_T` (rs12740374; chr1:109274968, ±1000 kb)
- **Exposure:** SORT1 expression - minor salivary gland (eQTL Catalogue QTD000276)
- **Outcome:** cholesterol in medium VLDL (GWAS Catalog GCST90269602)
- **n_pairs joined:** 2547
- **n_palindromic_excluded:** 333
- **LD reference:** 1000G Phase 3 v5b GRCh38, EUR (n=503 samples)
- **Plot:** 1_109274968_G_T_full_locuscompare.png

## Notes
- fetched 1000G region VCF to ~/.clawbio/locuscompare_cache/1000g/chr1_108774968_109774968.vcf.gz
- effect-size slope = -0.177 (consistent with the inverse SORT1↑→LDL↓ direction reported by Musunuru 2010)
```

`<output_dir>/manifest.yaml` (excerpt):

```yaml
skill: locuscompare-region-render
version: 0.1.0
lead_variant_id: 1_109274968_G_T
lead_rs_id: rs12740374
n_pairs: 2547
n_palindromic_excluded: 333
plot_path: 1_109274968_G_T_full_locuscompare.png
render_block:
  exposure_source: eqtl_catalogue
  exposure_source_release: v7+
  exposure_study_id: QTD000276
  outcome_source: gwas_catalog_harmonised
  outcome_source_release: '2026-05-06'
  outcome_study_id: GCST90269602
  ld_panel: 1000g_phase3_v5b_grch38_basic
  ld_panel_super_pop: EUR
  plink_version: PLINK v1.90b6.27 64-bit (2023-05-09)
  window_bp: 1000000
  ancestry_caveats: []
  fetched_at: '2026-05-09T11:42:06Z'
```

`<output_dir>/figure.png`: four panels stacked - outcome Manhattan (cholesterol-VLDL), exposure Manhattan (SORT1 eQTL), GENCODE gene track (CELSR2, PSRC1, SORT1, SARS1, SYPL2, ATXN7L2, ...), cross-trait `-log10p` scatter (tight diagonal at the top-right where the lead lives), effect-size scatter (negative slope). A representative rendering is bundled at `examples/02_eqtl_catalogue_x_gwas_catalog/expected_render.png`.

Output directory layout:

```
<output_dir>/
├── report.md                      # Markdown summary + caveats
├── figure.png                     # 4-panel LocusCompare PNG
├── manifest.yaml                  # Reproducibility provenance
├── tables/
│   ├── pairs.tsv                  # All joined pairs with both sides' beta/se/p + r2_lead
│   └── palindromic_excluded.tsv   # Variants dropped for palindromic ambiguity
└── reproducibility/
    ├── commands.sh                # Exact command to reproduce
    ├── environment.yml            # Pinned package versions
    ├── input_config.yaml          # Resolved config (defaults filled in)
    └── checksums.sha256           # Hashes of input + output artefacts
```

## Gotchas

1. **Palindromic variants (~10-20 % of common variants).** A/T or G/C variants have ambiguous strand alignment. The skill excludes them from the LocusCompare scatter by default (count surfaced in `manifest.yaml: n_palindromic_excluded`). They still appear in the per-trait Manhattans because the position + p-value is informative even without strand. To override (only safe with allele-frequency cross-checks confirming strand alignment), add `harmonisation: {exclude_palindromic: false}` to the config.

2. **Ancestry mismatch with 1000G EUR (the default LD reference).** If your sumstats are from a non-EUR cohort (FinnGen Finnish-EUR, Pan-UKBB AFR/EAS/CSA, BBJ EAS), the LD reference picks up subtle r² differences that color points slightly off. For Finnish-EUR vs 1000G EUR the divergence is ~0.05 r² on common variants per Locke 2019; for AFR vs EUR it can be substantial. Set `ld.super_pop: AFR | AMR | EAS | SAS` to match. Add an explicit caveat to the rendered plot via `caveats: [...]` so consumers see the assumption.

3. **All four inputs MUST be region-aligned to the same window.** Mismatched windows produce a silent garbage plot. The skill validates that the exposure slice, outcome slice, LD panel, and gene track all overlap the requested `(lead, window_bp)` and refuses to render with mismatched coverage.

4. **Lead variant must be present in BOTH the eQTL and GWAS slices.** If the lead is missing from one side (commonly: low-MAF variant dropped by one harmoniser), the renderer falls back to a proxy variant in LD (r² > 0.8); the manifest records the substitution. If no proxy exists, the skill refuses to render and the orchestrator falls back to a credible-set-only view.

5. **Effect-allele harmonisation across the four inputs is the renderer's input contract, not its job.** The two sumstats slices and the LD panel must arrive in the canonical (chr, pos, ref, alt) GRCh38 ALT-effect form - this is what the bundled fetchers (`eqtl-catalogue-region-fetch`, `gwas-catalog-region-fetch`) emit. User-supplied TSVs not in this form should be normalised upstream (`bcftools norm` for indels). The renderer's harmonisation step is for cross-trait flip / palindromic handling, not for single-trait normalisation.

6. **Render time is bounded by tabix fetches and plink LD compute.** Cold-cache (first render of a region): typically 10-50 s. Warm-cache (region already fetched and LD computed): sub-second. Surface the timing to the user when relevant (multi-render comparisons, demo loops).

7. **Visual diagonal-vs-two-cluster interpretation is informative but not definitive.** A clean diagonal SUPPORTS H4 (single shared causal variant) but does not prove it; two distinct causal variants in tight LD can produce a visually diagonal pattern. The agent must always cite the formal coloc PP-H4 (and PP-H3) alongside the visual interpretation.

8. **Wide windows (> 2 Mb) bloat memory and obscure structure.** The default 1 Mb window is calibrated for typical cis-eQTL × GWAS coloc loci. Going wider rarely helps interpretation and slows the LD compute (plink scales with #variants²). For multi-locus views, render multiple plots and compose them.

9. **`p = 0` on rare extremely-significant variants.** Some sources emit `p = 0` when the actual value is below floating-point precision. The renderer substitutes the underflow floor (`5e-324`) before plotting on `-log10`. The reported `-log10(p)` for such variants is ~323, not their true magnitude.

## Safety

**Not for clinical decisions.** This skill returns a visualisation of summary-statistic colocalisation; it is an interpretation aid, not a clinical or regulatory artifact. Do not use rendered plots for direct clinical decision-making.

**Visual interpretation of LocusCompare patterns is informative but not definitive.** Coloc PP-H4 from a formal Bayesian colocalisation analysis (COLOC, COLOC+SuSiE, SharePro, eCAVIAR) and supporting checks (`n_colocalising_variants` ≥ 5, credible-set width < 200 kb, ancestry match between exposure and outcome) should always be cited alongside the visual.

**Local-first; no data upload.** All computation runs locally. Input sumstats stay on the user's machine. Outbound network calls hit only public-domain databases (eQTL Cat FTP, GWAS Cat FTP, 1000G FTP, Ensembl REST for GENCODE). Cache stored at `~/.clawbio/locuscompare_cache/`; the on-demand GENCODE gene-track cache reads `LOCUSCOMPARE_CACHE_DIR` at call time, so CI / sandboxed environments can redirect with `LOCUSCOMPARE_CACHE_DIR=/path/to/cache` per invocation.

**Reproducibility.** Every run emits a SHA-256 checksum manifest covering inputs and outputs. Failed LD compute or missing gene track is surfaced in the manifest's caveats list, not silently degraded.

## Agent Boundary

The skill renders a 4-panel LocusCompare visualisation for a colocalisation result. The agent should:

- **Use the rendered figure to support or refute the formal coloc PP-H4 result.** Visual confirmation strengthens a high PP-H4 claim; a visual two-cluster pattern weakens it.
- **NOT replace formal coloc analysis with the visual alone.** PP-H4 and the supporting checks are the canonical decision criteria; the visual is an interpretation aid.
- **Surface the H3 vs H4 ambiguity when the visual is two-cluster.** Even with PP-H4 > 0.8, a two-cluster visual pattern indicates likely distinct causal variants in LD; flag for human review.
- **NOT cherry-pick variants outside the rendered window for downstream interpretation.** The visual establishes context for the rendered window only.
- **Cite the rendered window, LD reference + super-pop, OT release (when applicable), exposure / outcome study ids, and lead variant** in the user-facing reply. Per the user-friendly enum-expansion rule (`AGENTS.md`): `STRN × heart failure (FINNGEN_R12_I9_HEARTFAIL); window ±500 kb of lead chr2:36910110:C>T; LD = 1000G Phase 3 EUR; OT release 26.03`.
- **Surface the manifest's caveats list** (palindromic exclusions, missing-lead proxy notes, ancestry mismatches) verbatim in the user-facing reply.
- **NOT decide GO/NO-GO on a target** based on the visual alone. Chain to `target-validation-scorer` for synthesis; this skill is one input among many.
- **NOT silently swap super-populations.** If the upstream cohort's ancestry does not match the requested LD super-pop, surface explicitly and ask the user to confirm before proceeding.

## Citations

- Liu et al. (2019). *Abundant associations with gene expression complicate GWAS follow-up.* Nat Methods 16, 749-755. doi:10.1038/s41592-019-0431-x - LocusCompare convention.
- Hemani et al. (2018). *The MR-Base platform supports systematic causal inference across the human phenome.* eLife 7:e34408. doi:10.7554/eLife.34408 - palindromic-exclusion in MR/coloc.
- Kerimov et al. (2021). *A compendium of uniformly processed human gene expression and splicing quantitative trait loci.* Nat Genet 53, 1290-1299. doi:10.1038/s41588-021-00924-w - eQTL Catalogue.
- Sollis et al. (2023). *The NHGRI-EBI GWAS Catalog: knowledgebase and deposition resource.* Nucleic Acids Res 51, D977-D985. doi:10.1093/nar/gkac1010 - GWAS Catalog.
- 1000 Genomes Project Consortium (2015). *A global reference for human genetic variation.* Nature 526, 68-74. doi:10.1038/nature15393 - 1000G Phase 3 reference.
- Frankish et al. (2021). *GENCODE 2021.* Nucleic Acids Res 49, D916-D923. doi:10.1093/nar/gkaa1087 - GENCODE.
- Wang et al. (2020). *A simple new approach to variable selection in regression, with application to genetic fine mapping.* JRSS-B 82, 1273-1300. doi:10.1111/rssb.12388 - SuSiE (chained-skill: `fine-mapping`).
