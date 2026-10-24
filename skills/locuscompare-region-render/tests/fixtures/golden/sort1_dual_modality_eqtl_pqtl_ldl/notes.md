# SORT1 dual-modality (eQTL + pQTL) golden-parity fixture

**What this fixture locks.** One canonical seed for the v1.3 release: the
SORT1 1p13.3 LDL/CHD locus rendered TWICE against the same outcome
(GWAS Catalog `GCST90269602`, cholesterol in medium VLDL) via the two
exposure backends the orchestrator dispatches between:

1. **eQTL Catalogue** path: `EXPOSURE_KIND_EQTL_CATALOGUE` -> `EQTLCatalogueClient`,
   pulling `QTD000276` (GTEx minor salivary gland ge-eQTL) for ENSG00000134243.
2. **UKB-PPP** path: `EXPOSURE_KIND_UKB_PPP` -> `UKBPPPClient`, pulling
   the bundled SORT1 EUR slice (`syn51469328`, `OID20213` Olink reagent).

Both paths share the canonical SORT1 lead variant `1_109274968_G_T`
(`rs12740374`, Musunuru 2010) and a +/-500 kb window.

## Class of bug this catches

A golden-parity fixture earns its keep when it exercises cross-source /
cross-skill behaviour that toy unit fixtures normalise away. This fixture
covers:

- Cross-backend exposure dispatch: a regression in `dispatch_exposure_kind`
  or the `_render_for_spec` branch on `spec.exposure_kind` flips the
  fetcher and the manifest_block diverges (e.g. `exposure_source` switches
  to the wrong backend, `exposure_protein_label` empties out for pQTL).
- Shared-outcome harmonisation: both renders consume the same GWAS slice
  but with two different exposure variant sets, exercising the
  `wald_ratio.harmonise_regions_for_locuscompare` join across two exposure
  shapes. A regression in palindrome handling, ALT-effect alignment, or
  the join key (`chr_pos_ref_alt`) flips `n_pairs` /
  `n_palindromic_excluded`.
- Manifest-block correctness for pQTL extras: `exposure_protein_label`,
  `exposure_ancestry`, `exposure_ancestry_label` are pQTL-only fields
  added in v1.3. The pQTL block locks them; the eQTL block locks them
  as empty strings.

## What is and is not locked

See `expected.yaml` header. Locked: study identifiers, ancestry,
window, lead variant, n_pairs / n_palindromic_excluded,
scatter_downsample state, ancestry caveats, ld_panel state
(`"none"`, since the dev box lacks plink2 and the test passes
`ld_client=None`). NOT locked: `fetched_at` (render timestamp),
`plot_artifact` (filename basename).

## Refresh

A failing assertion on this fixture is by default a code regression,
not a fixture-staleness signal. Refresh requires an explicit commit
with rationale; the fixture is the source of truth for cross-backend
dispatch.

## Four keys added (0.2.2)

`exposure_modality`, `exposure_modality_label`, `exposure_quant_method` and
`exposure_quant_method_label` are now locked on both blocks (eQTL: `eqtl` /
`ge`; pQTL: `pqtl` / null). New keys, no existing value moved.


## Caveat wording (0.2.2)

The no-LD caveat on both blocks changed from "points show as grey" to "points drawn hollow (no LD value)". The renderer now draws a variant with no r² to the lead as a hollow point under its own legend entry, so "grey" became false: filled grey now means a measured r² of 0.0-0.2. Wording only; no number moved.
