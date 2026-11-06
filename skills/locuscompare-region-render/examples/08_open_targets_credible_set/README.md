# Example 08 - from an Open Targets credible-set page (IL23R plasma protein × inflammatory bowel disease)

You are looking at this Open Targets Platform page:

    https://platform.opentargets.org/credible-set/2f6755df0a7737586c9ff36fa2643c61

It shows an inflammatory bowel disease signal at IL23R and, under its colocalisations, one molecular QTL. You want the regional LocusCompare plot for that pair.

The page itself is a JavaScript application, so its content cannot be read from the URL. The id at the end of the URL is enough: `references/from_open_targets.md` turns it into this config in four Open Targets API steps.

## In a conversation

The skill is meant to be driven by an assistant (Claude, ChatGPT or another model that can run the skill's scripts). The user pastes the page and says what they want; the assistant runs `ot_credible_set.py`, relays what it found, and renders. An example exchange (the assistant's lines restate what the helper prints for this page):

> **User:** I'm looking at https://platform.opentargets.org/credible-set/2f6755df0a7737586c9ff36fa2643c61 . Can you make the LocusCompare plot for the molecular QTL there?
>
> **Assistant:** This is an inflammatory bowel disease signal (GCST90225550, Mize 2022) led by rs11209026, the IL23R R381Q missense variant. Open Targets lists one molecular QTL for it: the IL23R plasma protein QTL from Sun 2018 (colocalisation H4 0.991). The GWAS file can be read by region, so I'm plotting that pair. One caution: rs11209026 changes the protein sequence, which can affect how the assay binds as well as how much protein there is.

With a single QTL and a readable GWAS there is nothing to choose, so the assistant does not ask.

## Resolved pair

| | Open Targets | Config |
|---|---|---|
| GWAS (outcome) | `GCST90225550`, inflammatory bowel disease (Mize 2022) | `outcome.fetch.accession: GCST90225550` |
| Lead variant | `1_67240275_G_A` (rs11209026, IL23R R381Q) | `lead.variant_id` |
| Molecular QTL (exposure) | `sun_2018_aptamer_plasma_il23r_5088_175_3_1`, pQTL, H4 0.991; the only molecular QTL listed | `dataset_id: QTD000584` (from `summarystatsLocation`), `molecular_trait_id: ENSG00000162594` (from `target.id`) |

The GWAS file passes Step 4: `harmonised/GCST90225550.h.tsv.gz` and its `.tbi` are on the GWAS Catalog FTP.

## Run

```bash
python skills/locuscompare-region-render/cli.py \
    --demo 08_open_targets_credible_set \
    --output runs/08_open_targets_credible_set/
```

## What you should see

- Both regional panels peak at rs11209026.
- The effect-size panel puts the lead on the Wald-ratio line (slope about +1.1): the A allele lowers the measured IL23R protein and lowers IBD risk. The study reports odds ratios, so the outcome axis reads "log odds ratio, derived from the reported odds ratio".
- On the LocusCompare panel, other variants sit along each axis rather than on the diagonal: strong for one trait and weak for the other, so not every association in the window is shared with the protein QTL.

Reading it: rs11209026 is a missense variant. The study that produced this protein QTL notes that "genetic variants that change protein structure may result in apparent cis pQTLs owing to altered aptamer binding rather than true quantitative differences in protein levels" (Sun BB, Maranville JC, Peters JE, et al. Genomic atlas of the human plasma proteome. Nature. 2018;558(7708):73-79. doi:10.1038/s41586-018-0175-2. PMID: 29875488). Treat the protein-level reading of this pair with that caveat.

Rendering needs network access to the eQTL Catalogue and GWAS Catalog FTP, the 1000 Genomes FTP and Ensembl REST, plus `plink` 1.9 on PATH.
