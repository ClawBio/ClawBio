# variant_annotation (demo)

> **Driver:** see [USAGE.md](../../USAGE.md) for the `clawbio` CLI lifecycle and refusals. Then `clawbio --help`. This file only describes *this* skill's DAG, inputs, and outputs.

Two-step demo skill that exercises the `clawbio` CLI spine. Not a port of the real
ClawBio 1.x `variant_annotation` skill — the commands here are placeholder `cp`/`echo`
shell snippets so the pipeline runs without bcftools or VEP installed.

## DAG

```
normalize  →  annotate
```

## Steps

1. **normalize** — stand-in for `bcftools norm`. Emits `normalized_vcf`.
2. **annotate** — stand-in for VEP. Requires `normalize` done and its `normalized_vcf`
   emit present. Emits `annotated_vcf`.

## Inputs

| Name | Required | Description |
|---|---|---|
| `vcf_path` | yes | path to the input VCF |
| `reference` | yes | path to the reference FASTA |

## Outputs

| Emit | Producer | Description |
|---|---|---|
| `normalized_vcf` | `normalize` | placeholder copy of the input VCF |
| `annotated_vcf` | `annotate`  | placeholder copy of `normalized_vcf` |

## Quick run

```bash
clawbio start example/variant_annotation \
  --input vcf_path=/path/to/input.vcf \
  --input reference=/path/to/ref.fa
# then follow the loop in USAGE.md: next → open → done → next → … → bundle
```
