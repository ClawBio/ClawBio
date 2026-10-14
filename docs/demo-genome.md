# Demo Genome — The Corpasome

ClawBio's demo data comes from a **real human genome**: the 23andMe genotype of Manuel Corpas, published under a CC0 (public domain) license as part of the [Corpasome](https://link.springer.com/article/10.1186/1751-0473-8-13) — one of the first personal genomes to be made fully open for research.

> Corpas, M. (2013). Crowdsourcing the Corpasome. *Source Code for Biology and Medicine*, **8**, 13.
> [doi:10.1186/1751-0473-8-13](https://link.springer.com/article/10.1186/1751-0473-8-13)

## Why a real genome matters

Synthetic test data can only exercise the happy path. A real genome contains the messy reality of consumer genotyping: missing SNPs, heterozygous calls across multiple pharmacogenes, and combinations that produce genuinely actionable clinical alerts.

In Manuel's case, the combination of **VKORC1 TT** (high warfarin sensitivity) and **CYP2C9 \*1/\*2** (intermediate metabolizer) triggers an **AVOID** recommendation for warfarin — the most commonly prescribed anticoagulant worldwide. This is exactly the kind of finding that pharmacogenomic testing is designed to catch.

## Download the full genome

The complete 23andMe SNP chip file (~600,000 variants) is available on Figshare:

**[23andMe SNP chip genotype data — doi:10.6084/m9.figshare.693052](https://figshare.com/articles/dataset/23andMe_SNP_chip_genotype_data/92682)**

ClawBio ships a compressed copy at `skills/genome-compare/data/manuel_corpas_23andme.txt.gz`.

## File format

23andMe files are tab-separated with four columns: rsID, chromosome, position, and genotype. Lines starting with `#` are comments.

```
# rsid  chromosome  position  genotype
rs4244285   10  96541616    AG
rs4986893   10  96540410    GG
rs12248560  10  96522463    CT
rs3892097   22  42524947    CC
rs16947     22  42523943    AG
rs1065852   22  42526694    CC
rs28371725  22  42524175    CT
rs1799853   10  96702047    CT
rs1057910   10  96741053    AA
rs9923231   16  31107689    TT
```

Each row is one **SNP** (single nucleotide polymorphism) — a position in the genome where people commonly differ by a single DNA letter. You inherit one copy from each parent, so the genotype column always shows two letters: `AG` means you got an A from one parent and a G from the other (heterozygous — two different copies), while `TT` means both parents passed on the same T variant (homozygous — two identical copies).

Out of the ~600,000 SNPs in a typical 23andMe file, ClawBio's PharmGx reporter focuses on 21 that sit inside genes encoding drug-metabolising enzymes. It translates each genotype into a **star allele** — a standardised label like CYP2C9 \*2 that geneticists use to describe a known functional variant. It then combines the two star alleles you carry (one per chromosome) into a **diplotype** (e.g. \*1/\*2), which determines how fast or slow your body processes a given drug. Finally, it looks up that diplotype in the [CPIC guidelines](https://cpicpgx.org/) — peer-reviewed, evidence-based rules that map each gene–drug pair to a dosing recommendation: standard dose, use with caution, or avoid.

## Try it

```bash
python clawbio.py run pharmgx --demo
```

This runs the PharmGx reporter against the 21 PGx SNPs extracted from the Corpasome. The output includes a warfarin AVOID alert, gene profile table, and drug summary for 51 medications across 12 pharmacogenes.

## Run PharmGx with your own genotype file

From a source checkout, run `uv sync` first, then use the following Bash/Zsh
example from the repository root. Replace the input filename with your own
23andMe or AncestryDNA file:

```bash
INPUT_FILE="$HOME/genomes/my_genotype.txt"
REPORT_DIR="$PWD/output/my_pharmgx_$(date +%Y%m%d_%H%M%S)"

if [ ! -f "$INPUT_FILE" ]; then
  printf 'Input file not found: %s\n' "$INPUT_FILE"
elif [ -e "$REPORT_DIR" ]; then
  printf 'Output already exists; choose a new directory: %s\n' "$REPORT_DIR"
elif uv run python clawbio.py run pharmgx \
  --input "$INPUT_FILE" \
  --output "$REPORT_DIR"; then
  printf 'Report directory: %s\n' "$REPORT_DIR"
  # Optional on macOS: open the report only after a successful run.
  if [ "$(uname -s)" = "Darwin" ]; then
    if [ -f "$REPORT_DIR/report.html" ]; then
      open "$REPORT_DIR/report.html"
    elif [ -f "$REPORT_DIR/report.md" ]; then
      open "$REPORT_DIR/report.md"
    fi
  fi
fi
```

### Troubleshooting paths and missing reports

- **Use either an absolute path or a path relative to your home directory.**
  `INPUT_FILE="/path/to/my_genotype.txt"` is absolute;
  `INPUT_FILE="$HOME/genomes/my_genotype.txt"` is home-relative. Do not prepend
  `$HOME` to an already absolute path: `"$HOME/Users/example/genomes/file.txt"`
  joins two paths rather than referring to `/Users/example/genomes/file.txt`.
- **Quote paths.** Quoting `"$INPUT_FILE"` and `"$REPORT_DIR"` preserves spaces
  in filenames. Use `$HOME`, not a quoted literal `"~/genomes/file.txt"`,
  because a quoted tilde does not expand to your home directory.
- **Check the analysis status before opening a report.** A failed run may
  leave an empty output directory. An HTML report that does not exist is not
  evidence of a browser problem; resolve the analysis error first.
- **Keep personal outputs private.** `output/` is ignored by Git in this
  repository. Do not force-add genotype files, generated reports, or archives
  to a public commit or pull request.

ClawBio is a research and educational tool. It is not a medical device and
does not provide clinical diagnoses. Consult a healthcare professional before
making any medical decisions.

## 30x Whole-Genome Sequencing

In addition to the 23andMe SNP chip, ClawBio now ships subsets from Manuel Corpas's **30x Illumina whole-genome sequence** (GRCh37, Dante Labs). The WGS data covers ~4 million SNPs, ~600K indels, structural variants (DEL, DUP, INV, BND, INS), and copy number variants, all licensed CC0.

This dataset is provided for research and educational purposes only; it must not be used for clinical decision-making.

| Subset | Contents | Location |
|--------|----------|----------|
| chr20 SNPs + indels | Chromosome 20 variants | `corpas-30x/subsets/chr20_snps_indels.vcf.gz` |
| PGx loci | 31 pharmacogenomic SNPs (WGS) | `corpas-30x/subsets/pgx_loci.vcf.gz` |
| NutriGx loci | 30 nutrigenomics SNPs (WGS) | `corpas-30x/subsets/nutrigx_loci.vcf.gz` |
| SV calls | Structural variants | `corpas-30x/subsets/sv_calls.vcf.gz` |
| CNV calls | Copy number variants | `corpas-30x/subsets/cnv_calls.vcf.gz` |

Full VCFs are available from Zenodo ([doi:10.5281/zenodo.19297389](https://doi.org/10.5281/zenodo.19297389)). To download and prepare:

```bash
python scripts/prepare_corpas_30x.py --all
```

See [reference-genome.md](reference-genome.md) for full documentation, use cases, and citation details.
