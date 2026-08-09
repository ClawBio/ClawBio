---
name: clinical-variant-prioritizer
description: Screen a genotype set (array or WGS-derived) against OMIM-morbid, ACMG-SF and Hereditary-Cancer gene panels and prioritise carried variants by ClinVar significance, gnomAD frequency, inheritance model and zygosity, following the pathogenicity-screening method of Corpas et al. 2021 (Whole Genome Interpretation for a Family of Five).
license: MIT
metadata:
  version: 0.1.0
  author: Manuel Corpas
  domain: genomics
  reference: "Corpas M, Megy K, Mistry V, Metastasio A, Lehmann E. Whole Genome Interpretation for a Family of Five. Front Genet. 2021;12:535123. doi:10.3389/fgene.2021.535123"
  tags:
  - clinical-genomics
  - variant-prioritisation
  - clinvar
  - acmg
  - pathogenicity
  - carrier-screening
  openclaw:
    emoji: "🩺"
    os: [darwin, linux]
    trigger_keywords:
    - variant prioritisation
    - clinical variants
    - pathogenic variant
    - ClinVar
    - carrier status
    - disease risk variants
---

# clinical-variant-prioritizer

Turn a genotype set into a prioritised list of clinically relevant variants, the
way a clinical genome analyst would: screen catalogued disease-gene panels,
then rank what is carried by **how much it matters**, not by how loud the raw
ClinVar label is.

This skill implements the pathogenicity-screening stage of *Whole Genome
Interpretation for a Family of Five* (Corpas et al., Front Genet 2021): variants
are filtered through **OMIM-morbid**, **ACMG-SF** and **Hereditary-Cancer**
panels, intersected with **ClinVar** significance and **gnomAD** population
frequency, and classified by **inheritance model** and **zygosity**.

## Why it is not a raw ClinVar lookup

A raw lookup reports a label. This skill reports *actionability*. The same
"pathogenic" allele means very different things depending on context:

| Context | Category |
|---|---|
| Dominant / risk gene, allele carried | `actionable` |
| Recessive gene, homozygous | `affected` |
| Recessive gene, heterozygous | `carrier` (reproductive-risk only) |
| Uncertain / conflicting ClinVar | `uncertain` (flagged, not acted on) |
| Benign allele carried | `benign` |
| Variant not carried | `reference` |

A heterozygous carrier of a common, recessive, benign-spectrum allele is *not*
an actionable finding, even when ClinVar shows "pathogenic" submissions. Saying
so plainly is the point.

## Interface

```python
from api import run

result = run(
    {"rs28941785": "CT", "rs1800562": "GG"},   # rsid -> genotype
    options={"panel_path": "..."},              # optional custom panel
)
```

`run()` returns:

- `summary`: `panel_size`, `loci_tested`, `loci_carried`, `reference`,
  `not_tested`, and per-category counts (`actionable`, `affected`, `carriers`,
  `uncertain`, `benign`).
- `findings`: ranked list (highest priority first); each carries gene, HGVS,
  consequence, genotype, zygosity, ClinVar significance + review status, gnomAD
  frequency, condition, inheritance, panel membership, category and a
  plain-language `rationale`.
- `headline`, `method`, `disclaimer`.

## Panel

`data/clinical_panel.json` is a curated set of catalogued clinical loci, each
shipping its ClinVar significance, ClinVar review status, gnomAD frequency,
consequence, condition and inheritance model, so the screen is deterministic and
offline-reproducible (no per-call ClinVar/gnomAD/VEP network round-trips). Extend
it by adding entries; keys may be rsids or stable variant ids for WGS-only
variants not present on arrays.

## Limitations

Array-based input covers only catalogued loci and misses most rare variants; a
clean screen is not a clean genome. Heterozygous calls do not establish phase.
Confirm any finding with an accredited clinical assay. Research and educational
use only; not a clinical diagnosis.

## Test

```bash
python -m pytest tests/ -q
```
