# Example 09 - from an Open Targets credible-set page that lists many QTLs (PCSK9 lung expression × hypercholesterolemia)

You are looking at this Open Targets Platform page:

    https://platform.opentargets.org/credible-set/e723d7b2d797c4c3d9e3889ef3e79997

It shows a high-cholesterol signal near PCSK9 (GCST90013932, Mbatchou 2021, lead rs4500361; Open Targets maps the trait to the disease Hypercholesterolemia) and, under its colocalisations, ten molecular QTLs. You ask an assistant:

> I'm on https://platform.opentargets.org/credible-set/e723d7b2d797c4c3d9e3889ef3e79997 ; plot LocusCompare for gtex-v10_ge_lung_ensg00000169174.

locuscompare draws one pair. The QTL is named exactly, so nothing is asked about it; the page's GWAS study turns out not to be readable by region, so a substitute is proposed and confirmed. `ot_credible_set.py` (in the skill directory) makes both steps explicit; `references/from_open_targets.md` describes the rules.

## In a conversation

An example exchange with an assistant running the skill (its lines restate what `ot_credible_set.py` prints for this page):

> **User:** I'm on https://platform.opentargets.org/credible-set/e723d7b2d797c4c3d9e3889ef3e79997 . Plot LocusCompare for the PCSK9 lung QTL.
>
> **Assistant:** This high-cholesterol signal near PCSK9 (GCST90013932, lead rs4500361) colocalises with 10 molecular QTLs. Two match "PCSK9 lung": gene-level expression (H4 0.9992) and exon-level expression (H4 0.9991), both GTEx lung. I'd use gene-level expression. Which one?
>
> **User:** Gene level.
>
> **Assistant:** Two things before I plot. The page's GWAS study has no indexed summary statistics, so it can't be read by region; the closest substitute is GCST90475721, hypercholesterolemia in the Million Veteran Program (Verma 2024), which colocalises with the same lung QTL at H4 0.989. And the lung data the skill can read is GTEx v8 (n=510); Open Targets used GTEx v10 (n=652). Use GCST90475721?
>
> **User:** Yes.
>
> **Assistant:** Here's the plot. Its lead is 1_55023869_C_T, a member of the page's credible set: rs4500361 isn't in the 1000 Genomes LD panel, and the next member up is a C/G variant whose strand can't be resolved. All three substitutions are in the caption.

Had the user pasted the study id (`gtex-v10_ge_lung_ensg00000169174`) the first question would not be asked. The assistant never picks silently: one match is used, several are listed with a recommendation, a substitution is proposed and confirmed.

## Three substitutions, stated

- **The outcome is not the page's study.** The plot shows GCST90475721 (hypercholesterolemia, PheCode 272.11, Verma 2024), not GCST90013932: the page study's GWAS Catalog folder has no indexed harmonised file under the name the fetcher reads. GCST90475721 has the same Open Targets disease (Hypercholesterolemia) and colocalises with the chosen QTL at H4 0.989.
- **The QTL data are an older release.** Open Targets computed the colocalisation on GTEx v10 lung (n=652); the bundled fetcher reads GTEx v8 lung from eQTL Catalogue r7 (dataset QTD000271, n=510). Same gene and tissue, earlier release.
- **The plot's lead is another member of the credible set, chosen at render time.** The config keeps the page's lead rs4500361 (`1_55025188_T_C`) and lists the five credible-set members, PIP high to low, under `lead.candidates`. The render then finds that the page's lead is absent from the 1000 Genomes LD panel (centred on it, no point would have LD colouring), that the next member, `1_55021673_C_G` (PIP 0.153), is a C/G variant whose strand cannot be aligned, and uses `1_55023869_C_T` (PIP 0.117), present in the eQTL data, the GWAS and the panel. The x-range and title stay on the window fetched around the page's lead, so the diamond sits 1.3 kb left of centre rather than implying data the window never fetched. `1_55026242_C_T`, the lead of GCST90475721's own credible set for this signal, is absent from GTEx v8 lung.

The first two are written into `config.yaml` under `provenance.substitutions` and `caveats`; the third is decided by the render, which prints it with the other caveats and records it in the manifest (`lead_fallback`, `lead_candidates`). Delete `candidates` from the config to keep the plot on the page's lead.

## Step 1: the QTLs

```
$ python ot_credible_set.py options e723d7b2d797c4c3d9e3889ef3e79997
Open Targets credible set e723d7b2d797c4c3d9e3889ef3e79997
  GWAS: GCST90013932, High cholesterol (SPA correction) (Mbatchou J 2021), n=407746
  lead variant: 1_55025188_T_C (rs4500361)

Molecular QTLs that colocalise with this signal: 10 listed by Open Targets, 10 shown.
#   gene   kind         tissue                         study     H4      direction  n    readable by the bundled fetchers
1   PCSK9  eQTL (gene)  upper lobe of left lung        GTEx-v10  0.9992  same       652  yes: eQTL Catalogue QTD000271 (older release: GTEx v8, n=510)
2   PCSK9  exon eQTL    upper lobe of left lung        GTEx-v10  0.9991  same       652  yes: eQTL Catalogue QTD000272 (older release: GTEx v8, n=510)
3   PCSK9  exon eQTL    tibial nerve                   GTEx-v10  0.9988  same       708  yes: eQTL Catalogue QTD000287 (older release: GTEx v8, n=532)
4   PCSK9  eQTL (gene)  tibial nerve                   GTEx-v10  0.9988  same       708  yes: eQTL Catalogue QTD000286 (older release: GTEx v8, n=532)
5   BSND   exon eQTL    blood                          GTEx-v10  0.9976  same       853  yes: eQTL Catalogue QTD000357 (older release: GTEx v8, n=670)
6   PCSK9  exon eQTL    blood                          GTEx-v10  0.9974  same       853  yes: eQTL Catalogue QTD000357 (older release: GTEx v8, n=670)
7   PCSK9  eQTL (gene)  blood                          GTEx-v10  0.9969  same       853  yes: eQTL Catalogue QTD000356 (older release: GTEx v8, n=670)
8   PCSK9  eQTL (gene)  adipose tissue                 GTEx-v10  0.9943  same       749  yes: eQTL Catalogue QTD000116 (older release: GTEx v8, n=581)
9   BSND   eQTL (gene)  blood                          GTEx-v10  0.9928  same       853  yes: eQTL Catalogue QTD000356 (older release: GTEx v8, n=670)
10  PCSK9  exon eQTL    induced pluripotent stem cell  HipSci    0.7463  opposite   322  yes: eQTL Catalogue QTD000362

Recommended QTL: row 1, PCSK9 eQTL (gene) in upper lobe of left lung (gtex-v10_ge_lung_ensg00000169174)
  - highest colocalisation probability (H4 0.9992) among the 10 QTLs shown that the bundled fetchers can read
  - Open Targets computed this colocalisation on GTEx v10 (n=652); the bundled fetcher reads GTEx v8 from eQTL Catalogue r7 (n=510), so the plot shows the same gene and tissue from an older release

Next: choose a QTL, then run
  python ot_credible_set.py options e723d7b2d797c4c3d9e3889ef3e79997 --qtl <row number or study id>
```

With the exact study id the request names one row, so the assistant goes straight to Step 2. Had the user said only "the PCSK9 lung QTL", `--qtl "PCSK9 lung"` lists the two lung rows (gene-level and exon-level eQTL) with a recommendation among them, and the assistant asks which one.

Eight of the ten rows are PCSK9, in five tissues. Gene-level and exon-level rows are listed separately because they are different measurements, so nothing is folded here. The highest H4 is not evidence that lung is where PCSK9 acts; see "Which tissue" in `references/from_open_targets.md`.

## Step 2: the GWAS side, for the chosen QTL

```
$ python ot_credible_set.py options e723d7b2d797c4c3d9e3889ef3e79997 --qtl gtex-v10_ge_lung_ensg00000169174
...
Chosen QTL: PCSK9 eQTL (gene) in upper lobe of left lung (gtex-v10_ge_lung_ensg00000169174)
Page GWAS file check (GCST90013932): no: no indexed harmonised file under the name the fetcher reads
The page study GCST90013932 cannot be read by region, so these GWAS studies that colocalise with the chosen QTL could stand in for it:
123 listed by Open Targets, 15 shown.
#   study         trait                                          first author              vs page          H4      direction  n       readable by region
1   GCST90013882  High cholesterol (Firth correction)            Mbatchou J 2021           same disease     0.9993  same       407746  no: no indexed harmonised file under the name the fetcher reads
2   GCST90475721  Hypercholesterolemia (PheCode 272.11)          Verma A 2024              same disease     0.9886  same       417271  yes: harmonised file with a tabix index
3   GCST90435755  Hypercholesterolemia (PheCode 272.11)          Zhou W 2018               same disease     0.9644  same       406276  yes: harmonised file with a tabix index
4   GCST90019501  Total cholesterol levels                       Sinnott-Armstrong N 2021  different trait  0.9999  same       355858  yes: harmonised file with a tabix index
5   GCST90019512  Low density lipoprotein cholesterol levels     Sinnott-Armstrong N 2021  different trait  0.9999  same       355197  yes: harmonised file with a tabix index
6   GCST90019496  Apolipoprotein B levels                        Sinnott-Armstrong N 2021  different trait  0.9999  same       354097  yes: harmonised file with a tabix index
7   GCST90475415  low density lipoprotein cholesterol (LDLC, me  Verma A 2024              different trait  0.9999  same       113213  yes: harmonised file with a tabix index
8   GCST90002412  Low density lipoprotein cholesterol levels     Klimentidis YC 2020       different trait  0.9998  same       431167  no: no indexed harmonised file under the name the fetcher reads
9   GCST90501203  Phospholipids to Total Lipids in Medium LDL p  Zoodsma M 2025            different trait  0.9995  same       450015  yes: harmonised file with a tabix index
10  GCST90498842  Phospholipids to Total Lipids in Medium LDL p  Zoodsma M 2025            different trait  0.9995  same       434646  yes: harmonised file with a tabix index
11  GCST90449583  Phospholipids to total lipids ratio in very l  Tambets R 2026            different trait  0.9993  same       185352  yes: harmonised file with a tabix index
12  GCST90092963  Concentration of small LDL particles           Richardson TG 2022        different trait  0.9980  same       115082  no: no indexed harmonised file under the name the fetcher reads
13  GCST90269623  Cholesterol in IDL (UKB data field 23526)      Davyson E 2023            different trait  0.9975  same       88329   yes: harmonised file with a tabix index
14  GCST90269624  Cholesteryl esters in IDL (UKB data field 235  Davyson E 2023            different trait  0.9975  same       88329   yes: harmonised file with a tabix index
15  GCST90269621  Total lipids in IDL (UKB data field 23524)     Davyson E 2023            different trait  0.9975  same       88329   yes: harmonised file with a tabix index
1 row(s) are the page's own study and are not offered.
107 more not shown (lower in the order, not checked); --json writes every row.
Recommended GWAS: GCST90475721, Hypercholesterolemia (PheCode 272.11) (Verma A 2024)
  - same Open Targets disease as the page study
  - highest H4 (0.9886) among the same-disease studies whose file is indexed
  - substitution: the plot would show GCST90475721, not the page study GCST90013932

Next: confirm the pair with the user, then run
  python ot_credible_set.py config e723d7b2d797c4c3d9e3889ef3e79997 --qtl 227c96a4630d2dbb1cb3cf1e637ebc49 --gwas GCST90475721 --out config.yaml
```

GCST90013882 colocalises more strongly but has no indexed file either, so it is not recommended. This output was recorded live on 2026-10-13; a file check that fails (a timeout, say) shows as "unknown" and is never recommended.

## Step 3: config and render

```bash
python ot_credible_set.py config e723d7b2d797c4c3d9e3889ef3e79997 \
    --qtl gtex-v10_ge_lung_ensg00000169174 --gwas GCST90475721 --out config.yaml
python skills/locuscompare-region-render/cli.py \
    --demo 09_open_targets_credible_set_many_qtls \
    --output runs/09_open_targets_credible_set_many_qtls/
```

## What you should see

Rendered on 2026-10-13: 4,146 variants joined in the 1 Mb window. The lead diamond is on `1_55023869_C_T` in both regional panels and both scatters, points are coloured by LD with it, and the effect-size panel's Wald-ratio slope at that lead is +0.209 (outcome axis in log odds ratio, derived from the reported odds ratio). The plot file keeps the requested lead's name, `1_55025188_T_C_full_locuscompare.png`. The three substitution sentences are printed as warnings and in the figure caption; the manifest's `lead_variant_id` is `1_55023869_C_T` and its `lead_requested_variant_id` is `1_55025188_T_C`.

Rendering needs network access to the eQTL Catalogue and GWAS Catalog FTP, the 1000 Genomes FTP and Ensembl REST, plus `plink` 1.9 on PATH.
