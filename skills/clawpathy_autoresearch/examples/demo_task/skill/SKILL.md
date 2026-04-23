# GWAS sumstats analysis — seed

## Workflow

1. Read the summary statistics file.
2. Identify genome-wide significant lead variants.
3. Compute the genomic inflation factor lambda.
4. Produce a QQ plot under `output/`.
5. Write a final summary JSON to `output/summary.json` with fields for
   lead count, lambda, and plot path.
