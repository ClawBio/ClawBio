# Measuring the legacy-EFF gap

Target: `skills/rare-high-impact-variants/rare_high_impact_variants.py`, unmodified.  
Input: `challenge1-b37-segregation.tsv` — the same records in all three runs.

| Metric | raw `EFF=` only | `MC=` unmapped | `MC=` mapped |
|---|---|---|---|
| `variants_processed` | 68 | 68 | 68 |
| `carried_variants` | 68 | 68 | 68 |
| `high_impact_carried` | 0 | 29 | 68 |
| `rare_high_impact_count` | 0 | 2 | 12 |
| `high_impact_common` | 0 | 26 | 54 |
| `high_impact_frequency_unknown` | 0 | 1 | 2 |
| **exit code** | **0** | **0** | **0** |

## Reading the first column

`variants_processed` and `carried_variants` are correct. Every impact metric is 0. The exit code is 0 and the report is well formed. There is no warning on stderr. This is a false negative that reads as a clean run — which is why it is worth measuring rather than describing.

- `rare_high_impact_variants.py:126` — `consequence = info.get("MC", "") or info.get("Consequence", "") or info.get("ANN", "")`  
  EFF is not in this chain, so consequence is the empty string
- `rare_high_impact_variants.py:130` — `gene = (info.get("GENEINFO", "") or info.get("SYMBOL", "")).split(":")[0]...`  
  the gene symbol inside the EFF field is never recovered
- `rare_high_impact_variants.py:160` — `if not _is_high_impact(consequence): continue`  
  the record is dropped here, before any impact metric is incremented

## Reading the second column — the part that matters

Setting `MC=` to the SnpEff effect name is the obvious one-line fix, and it is worse than the bug. The matcher is an unanchored, case-insensitive substring test over eight terms, and SnpEff's names agree with three of them by coincidence:

| SnpEff name | contains a matching substring? |
|---|---|
| `FRAME_SHIFT` | no |
| `SPLICE_SITE_ACCEPTOR` | no |
| `SPLICE_SITE_DONOR` | no |
| `START_LOST` | yes, by accident |
| `STOP_GAINED` | no |
| `STOP_LOST` | yes, by accident |

So the unmapped run reports **29 of 68** rather than 68. A zero is obviously broken. A partial count is not, and nothing in the output distinguishes it from a correct one.

## Frequencies

66 of 68 records received a documented frequency from the cached build-matched layer. The rest carry no frequency key at all, so they land in the skill's own `frequency_unknown` bucket. Supplying a placeholder would manufacture the very claim this project refuses to make.

Unmapped effect names in the mapped run: none.

---

*The brief said this skill could not read the data. It can now, and the interesting result was not the fix but the near-miss beside it.*
