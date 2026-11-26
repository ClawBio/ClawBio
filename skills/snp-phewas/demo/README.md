# snp-phewas demo data

Real source rows, not synthetic values. Used by `--demo` and by the offline tests.

| File | Source | Content |
|---|---|---|
| `pan_ukbb/biomarkers-30780-both_sexes-irnt.tsv.bgz` (+ `.tbi`) | Pan-UK Biobank, release 2024-06-26 (`sumstats_flat_files/`) | GRCh37 1:109,807,590-109,827,590, every column of the source file |
| `pan_ukbb/icd10-E78-both_sexes.tsv.bgz` (+ `.tbi`) | same | same window |
| `pan_ukbb/phecode-411-both_sexes.tsv.bgz` (+ `.tbi`) | same | same window |
| `finngen/finngen_R12_I9_CHD.gz` (+ `.tbi`) | FinnGen R12 (`summary_stats/release/`) | GRCh38 1:109,264,968-109,284,968, every column of the source file |
| `finngen/finngen_R12_E4_LIPOPROT.gz` (+ `.tbi`) | same | same window |
| `finngen/finngen_R12_I9_IHD.gz` (+ `.tbi`) | same | same window |
| `mvp/GCST90475001-GCST90476000/GCST90475929/harmonised/GCST90475929.h.tsv.gz` (+ `.tbi`) | Million Veteran Program via the GWAS Catalog (`summary_statistics/`, harmonized), ischemic heart disease (phecode 411), European ancestry | GRCh38 1:109,264,968-109,284,968, every column, the source's CRLF line ends kept |
| `mvp/GCST90475001-GCST90476000/GCST90475416/harmonised/GCST90475416.h.tsv.gz` (+ `.tbi`) | same, LDL cholesterol mean, European ancestry | same window |
| `mvp/GCST90475001-GCST90476000/GCST90475416/harmonised/GCST90475416.running.log` | the GWAS Catalog's harmonization log for that file, `https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90475001-GCST90476000/GCST90475416/harmonised/GCST90475416.running.log`, downloaded 2026-10-24 | the whole file (2,244 bytes): per-`hm_code` counts of the sites harmonized and of those that failed |
| `mvp/GCST90475001-GCST90476000/GCST90475929/harmonised/GCST90475929.running.log` | same, `.../GCST90475929/harmonised/GCST90475929.running.log`, downloaded 2026-10-24 | the whole file (2,224 bytes) |
| `finngen/api_variant_1-109274968-G-T.json` | `https://r12.finngen.fi/api/variant/1:109274968-G-T` | the API answer trimmed to the three endpoints above plus one censored and one uncatalogued endpoint |

Each slice is +/-10 kb around rs12740374, cut from the public files (Pan-UK Biobank and
FinnGen on 2026-10-17, MVP on 2026-10-19), keeping
each file's own header line, and re-compressed and re-indexed with htslib `bgzip` / `tabix`
(the slices are tiny, so their indexes are correct; see the skill's Gotchas for why the
full Pan-UK Biobank files are read with this skill's own index walker). The window also
holds the strand-ambiguous A/T variant rs41279716 (GRCh38 1:109,270,888), which the
palindromic tests use.

## Terms and acknowledgments

- **Pan-UK Biobank** results are released under CC BY 4.0
  (https://pan.ukbb.broadinstitute.org/downloads/). Cite: Karczewski KJ, Gupta R, Kanai M,
  et al. Pan-UK Biobank genome-wide association analyses enhance discovery and resolution
  of ancestry-enriched effects. Nat Genet. 2025;57(10):2408-2417.
  doi:10.1038/s41588-025-02335-7. PMID: 40968291. This research has been conducted using the
  UK Biobank Resource.
- **Million Veteran Program** summary statistics and the GWAS Catalog's harmonization logs
  beside them are served by the NHGRI-EBI GWAS Catalog
  under the EMBL-EBI terms of use (https://www.ebi.ac.uk/about/terms-of-use/), which place no
  restrictions beyond the data owner's; MVP deposited no license text. Cite: Verma A, Huffman
  JE, Rodriguez A, et al. Diversity and scale: genetic architecture of 2068 traits in the VA
  Million Veteran Program. Science. 2024;385(6706):eadj1182. doi:10.1126/science.adj1182.
  PMID: 39024449.
- **FinnGen** R12 public results are free to download and use with citation and
  acknowledgment (https://www.finngen.fi/en/access_results). Cite: Kurki MI, Karjalainen J,
  Palta P, et al. FinnGen provides genetic insights from a well-phenotyped isolated
  population. Nature. 2023;613(7944):508-518. doi:10.1038/s41586-022-05473-8.
  PMID: 36653562. We want to acknowledge the participants and investigators of the FinnGen
  study.
