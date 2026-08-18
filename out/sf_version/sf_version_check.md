# Is the secondary-findings list we screen against current?

**Verdict: `CONFIRMED`**

evidence names ACMG SF v3.3 with 84 genes; the bundled list is v3.2 with 81. Genes named as added: ABCD1, CYP27A1, PLN.

- Bundled list: ACMG SF **v3.2**, 81 genes
- Provenance: ACMG_SF_V32_GENES imported from skills/clinical-variant-reporter/acmg_engine.py (81 genes, ACMG SF v3.2)

## Evidence retrieved

- [The ACMG releases 2025 update to secondary findings gene list; SF v3.3 | EurekAlert!](https://www.eurekalert.org/news-releases/1090415)
- [ACMG SF v3.3 list for reporting of secondary findings in clinical exome and genome sequencing: A policy statement of the American College of Medical Genetics and Genomics (ACMG) - PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC12318660)
- [ACMG SF v3.3 list for reporting of secondary findings in clinical exome and genome sequencing: A policy statement of the American College of Medical Genetics and Genomics (ACMG)](https://www.gimjournal.org/article/S1098-3600(25)00101-7/fulltext)
- [ACMG SF Genes](https://search.clinicalgenome.org/kb/genes/acmgsf)
- [ACMG Secondary Findings v3.3 Released](https://3billion.io/blog/acmg-secondary-findings-v3-3-released)
- [ACMG SF v3.3 list for reporting of secondary findings in clinical exome and genome sequencing: A policy statement of the American College of Medical Genetics and Genomics (ACMG) - Icahn School of Medicine at Mount Sinai](https://scholars.mssm.edu/en/publications/acmg-sf-v33-list-for-reporting-of-secondary-findings-in-clinical-)

## Why this matters here rather than in general

This skill withholds records for screening against stale evidence. Screening them against a stale gene list would be the same error one level up. The version therefore travels with every result, and where it could not be established the report says so instead of assuming.
