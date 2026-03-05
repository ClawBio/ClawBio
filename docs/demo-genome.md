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
