# Expected output of the two demos

Both demos read synthetic tables from `sumstats/`, written by `make_demo_data.py`
from a fixed seed. No network, no real data.

## `sort1_eqtl_x_ldl_gwas`: one shared causal variant

```bash
python skills/coloc-abf-region-run/coloc_abf_region_run.py \
  --demo sort1_eqtl_x_ldl_gwas --output /tmp/coloc_shared
```

```
[coloc-abf-region-run] PP.H4=0.990 PP.H3=0.010 shared_variants=200 \
  lead=1_109300000_A_C sensitivity=robust (the H4 call holds across p12 1e-06 to 1e-04)
```

`tables/posteriors.tsv`:

| hypothesis | posterior |
|---|---|
| PP.H0 | 1.23327e-33 |
| PP.H1 | 6.52889e-17 |
| PP.H2 | 2.01008e-19 |
| PP.H3 | 0.00965095 |
| PP.H4 | 0.990349 |

The variant made causal for both traits sits at position 109,300,000, and it is
the one reported. `excluded_variants.tsv` is empty apart from its header: both
tables carry the same 200 variants.

## `distinct_signals_negative`: two causal variants, 60 kb apart

```bash
python skills/coloc-abf-region-run/coloc_abf_region_run.py \
  --demo distinct_signals_negative --output /tmp/coloc_distinct
```

```
[coloc-abf-region-run] PP.H4=0.000 PP.H3=1.000 shared_variants=200 \
  lead=1_109361000_C_A sensitivity=no shared signal (H4 below the rule at the default p12)
```

Both traits are strongly associated across this window, and a regional p-value
plot would show two overlapping mountains. PP.H3 of 1.00 is the answer: they are
not the same signal. The report says so under "What this does not say", and the
reported lead variant is meaningless when H4 is this low, which is why the
sensitivity verdict leads rather than the lead variant.

The posteriors above are copied from a run of the committed tables. They
reproduce because the tables do; the tests assert the calls rather than every
digit, so a change in the arithmetic would be caught by the frozen R-parity
fixture rather than by this file.
