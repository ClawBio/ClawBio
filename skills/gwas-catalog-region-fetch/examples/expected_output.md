# Expected output — cholesterol-VLDL GWAS slice at the SORT1 locus

Input: `examples/input.json` — fetches a 1 Mb window of GCST90269602
("cholesterol in medium VLDL" GWAS).

## What you get

`RegionResult` dataclass with:

- `n_variants`: ~3000 variants
- `release.accession`: `GCST90269602`
- `release.fetched_at_utc`: ISO timestamp of the fetch
- `variants`: per-variant rows (variant_id, chromosome, position, ref, alt,
  beta, se, p_value, eaf if present, odds_ratio and its CI where reported)

The lead variant `1_109274968_G_T` (rs646776) should appear with p-value
near 1e-50 and a strongly-negative beta (per Musunuru 2010 inverse SORT1↑→LDL↓
biology).

## Reproducing

```bash
bash run_example.sh
```

Network required (~1-5 MB tabix fetch from EBI GWAS Catalog FTP).
