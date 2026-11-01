# igv-validator CLI options

Every option of `igv_validator.py`. Run `python skills/igv-validator/igv_validator.py --help` for the authoritative list; the skill's SKILL.md `## CLI Reference` has the everyday commands.

| Option | What it does |
|---|---|
| `--vcf` (`--input`) | VCF/VCF.gz of SNV, indel or SV calls (any caller) |
| `--tumor` | tumor BAM/CRAM, indexed (required except for `--summarize`) |
| `--normal` | matched normal BAM/CRAM; leave out for tumor-only |
| `--tumor-only` | ignore any normal (e.g. with `--demo`) |
| `--reference` | FASTA the BAMs were aligned to (screenshots, indel rescue, CRAM) |
| `--genes` | gene symbols. With `--samplesheet`: the genes to check, found in `--annotation` (default: the curated calls' genes). In a single check: selects VCF records by their gene names (VEP, SnpEff, ANNOVAR, GENE) |
| `--regions` | BED of genes (chrom, start, end, name): selects variants/SVs, sets copy-number genes, overview and interactive windows |
| `--variants` | TSV of chrom, pos, [ref, alt] to check exactly (e.g. a filtered table) |
| `--pass-only` | only FILTER=PASS calls |
| `--max-variants` | cap per run (default 50) |
| `--caller-sample` | VCF column whose caller counts to compare (default: the tumor's sample name) |
| `--tumor-name`, `--normal-name` | labels in reports and images (use your sample names) |
| `--cnv`, `--cnv-sample` | copy-number mode: segment file, and the sample in a multi-sample file |
| `--samplesheet` | one command for a whole run: every check for every sample, into a new dated folder with its summary |
| `--skip-unknown-genes` | with `--samplesheet`: run without gene names the annotation does not know (default: stop with suggestions) |
| `--reports-dir` | with `--samplesheet`: where runs are kept (default `igv_reports/`) |
| `--summarize` | build the summary page over every run in a folder (full: with images, interactive pages, zip) |
| `--index` | write `<folder>/index.html` listing every batch folder with a summary; it then keeps itself updated |
| `--project` | project folder of all runs: its `summary.html` is refreshed after this check (automatic for later checks once the folder has a summary) |
| `--curated-calls` (`--heatmap`) | curated table to compare with IGV (sample, gene, alteration) |
| `--overview` | with `--summarize` and `--regions`: one IGV image per gene and sample |
| `--interactive` | with `--summarize` and `--regions`: zoomable pages per sample and gene (needs igv-reports; contain read data) |
| `--annotation` | GTF/GFF3/BED of genes (e.g. GENCODE) drawn as a genes track in every image |
| `--no-bundle` | with `--summarize`: skip `igv_validation_full.zip` (made by default: the whole report, linked from the page) |
| `--no-igv` | counts and reports only |
| `--igv-path`, `--igv-timeout` | IGV to use (found automatically) and how long to wait for it (default 300 s) |
| `--output` | output folder |
| `--demo`, `--demo-cnv` | run on the bundled synthetic data |
