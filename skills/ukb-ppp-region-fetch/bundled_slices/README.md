# Bundled UKB-PPP regional slices

Pre-computed, harmonised regional pQTL summary statistics for the canonical demo cohort. Shipped inside the skill so end users running the bundled demos do NOT need a Synapse PAT.

## Slice file convention

```
<PROTEIN>__<ANCESTRY>__chr<C>__<start>_<end>.json.gz
```

Examples:

- `SORT1__EUR__chr1__108774968_109774968.json.gz`
- `IL6R__EUR__chr1__154404000_155404000.json.gz`

Each file is a gzipped JSON serialisation of a `RegionResult` (see `ukb_ppp_region_fetch.RegionResult.to_dict`); per-variant pQTL rows compress ~8.5x, so a 5,000-variant slice goes from ~3.5 MB on disk to ~430 KB in the PR. The file contains:

- `protein_label_short`, `ancestry`, `chromosome`, `region_start_bp`, `region_end_bp`, `n_variants`
- `variants`: list of harmonised `RegionVariant` rows (`variant_id`, `chromosome`, `position`, `ref`, `alt`, `beta`, `se`, `p_value`, `maf`, `effect_allele_frequency`, raw REGENIE fields)
- `release`: `UKBPPPRelease` block with `release_label`, `protein_hgnc`, `protein_uniprot`, `olink_reagent_id`, `olink_panel`, `ancestry`, `ancestry_label`, `n_samples` (largest per-variant REGENIE `N` over the rows), `n_samples_min` (smallest), `cohort_n` (published cohort size for the ancestry, Sun et al. 2023), `synapse_id`, `source_url`, `fetched_at_utc`
- `notes`: any parser caveats from the original live fetch

## How slices are generated

Maintainer-only path, requires a Synapse PAT one time:

```python
import os
import sys
sys.path.insert(0, "skills/ukb-ppp-region-fetch")
from ukb_ppp_region_fetch import (
    UKBPPPClient, save_region_result_as_bundled_slice,
)
os.environ["SYNAPSE_AUTH_TOKEN"] = "..."  # one-time
client = UKBPPPClient()
result = client.fetch_region(
    protein_label="SORT1", ancestry="EUR",
    chromosome="1", start_bp=108_774_968, end_bp=109_774_968,
)
save_region_result_as_bundled_slice(result)
# -> bundled_slices/SORT1__EUR__chr1__108774968_109774968.json.gz
```

## License + attribution

The pre-computed slices are direct derivatives of UKB-PPP summary statistics (Sun et al. 2023 Nature, PMID 37794186), redistributed under **CC-BY 4.0** with the upstream attribution string carried in each slice's `release` block:

> "UK Biobank Pharma Proteomics Project (Sun et al. 2023 Nature; PMID 37794186); accessed via Synapse syn51364943."

CC-BY 4.0 explicitly permits redistribution with attribution; no separate UK Biobank Application is required for the summary-statistics layer.

## When the slice you need is not here

The bundled inventory covers the **canonical demo cohort only**. For arbitrary protein x ancestry x region queries, the skill falls through to the live Synapse fetcher, which requires a free Synapse Personal Access Token in `SYNAPSE_AUTH_TOKEN`. See the top-level `SKILL.md` "First-time setup" section for the PAT walkthrough.

## Inventory (current)

v0.1.0 ships the SORT1 / EUR slice (chr1:108774968-109774968, the canonical 1p13.3 LDL/CHD locus). Additional slices land as the bundled-demo cohort grows; the slice convention scales by dropping further files into this directory.

- `SORT1__EUR__chr1__108774968_109774968.json.gz`: 5,205 variants; every row carries REGENIE `N` = 33,628.

2026-10-21 (0.1.2): `release.n_samples` 46673 -> 33628, `n_samples_min` 33628 and `cohort_n` 34557 added. The old value was a hardcoded cohort size cited to a table that does not exist in Sun et al. 2023; the new value is the REGENIE `N` carried by all 5,205 rows of this slice. The variant rows, `fetched_at_utc` (2026-05-15T16:36:28Z) and `synapse_id` (syn51469328) are unchanged.
