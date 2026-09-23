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
    os:
    - darwin
    - linux
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

## Optional disease-severity evidence

Use this offline research annotation with independently curated HPO tier
assignments. It does not classify HPO terms with an LLM or establish clinical
validity. Pass an evidence dictionary as `options["severity_evidence"]`.

Run the fully synthetic example from the repository root:

```bash
python - <<'PY'
import importlib.util
import json
from pathlib import Path
skill = Path("skills/clinical-variant-prioritizer").resolve()
spec = importlib.util.spec_from_file_location("clinical_prioritizer", skill / "api.py")
api = importlib.util.module_from_spec(spec)
spec.loader.exec_module(api)
result = api.run(
    {"synthetic-variant": "AG"},
    {"panel_path": skill / "examples/synthetic_panel.json",
     "severity_evidence": json.loads((skill / "examples/severity_evidence.json").read_text())},
)
print(json.dumps(result, indent=2))
PY
```

Expected: carrier category and Pathogenic ClinVar label unchanged; separate
Severe disease-severity annotation and severity_disagreement warning against
the synthetic panel's mild label. Example identifiers and assertions are
synthetic placeholders, not evidence about a real disease.

### Input contract (schema 1.0)

- Envelope: schema_version "1.0"; hpo_release as an ISO date;
  upstream_commit "03bc5f6be6456a3ca2c5206e399f0a38f879ce57"; nonempty records list.
- Each record: gene, disease_id (MONDO: plus seven digits), disease name,
  hpo_id (HP: plus seven digits), integer tier 1 through 4, frequency,
  qualifier (empty or NOT), citation (PMID:digits, doi:DOI, or HTTP(S) URL),
  source_span, source_version, and hpo_release matching the envelope.
- Supply a custom panel with a curated disease_id for each relevant entry.
  Match exact gene and MONDO identity; never infer disease from gene alone.
  Duplicate panel IDs and duplicate HPO terms within a gene-disease pair
  cause abstention. Resolve conflicting sources through curation.
- Frequencies: proportions in [0,1], fractions ("3/10"), explicit percentages
  ("30%"), or HP:0040280 through HP:0040285. Bare numbers above 1, missing
  frequencies, Excel dates and unknown codes are rejected.

Identifier syntax, declared versions and provenance presence are checked.
HPO/MONDO membership, actual release availability, source-span entailment,
tier accuracy and completeness are NOT verified. HPO dates are caller-declared
provenance, not loaded ontology releases. The SHA-256 records the canonical
input JSON; it does not authenticate evidence. Source spans are retained in
output and must be handled according to their sensitivity.

### Aggregation and output

Exclude NOT annotations and frequencies below 30%. Count unique eligible HPO
terms: more than one Tier 1 gives Profound; one gives Severe. At least one
Tier 2 gives Severe when Tier 2 + Tier 3 total at least four, otherwise
Moderate. Tier 3 alone gives Moderate. Eligible Tier 4-only evidence gives
Mild. Empty or entirely excluded evidence gives severity_unknown.
Mild describes only supplied eligible evidence; omitted severe phenotypes can
change the result. Do not interpret it as proof of a mild disease.

Rules follow the MIT-licensed
[upstream code](https://github.com/T0hid/hpo-classification-agent/blob/03bc5f6be6456a3ca2c5206e399f0a38f879ce57/severity_classification.py),
associated with [arXiv:2609.19569v1](https://arxiv.org/abs/2609.19569v1).
The exact pinned function and licence are retained in test fixtures.
This adapter deliberately abstains where upstream retains missing frequencies
or attempts to recover Excel-corrupted dates. No upstream datasets or
restricted production prompts are copied.

Findings gain disease_severity: label, reasons, warnings, tier counts,
included/excluded evidence, HPO release, upstream commit and canonical JSON
SHA-256. Top-level severity_evidence reports classified/unknown counts.
Malformed evidence invalidates the entire bundle, avoiding selective omission.
Incompatible schema/algorithm versions or missing identity yield unknown.
Without the option, the original output is unchanged. Severity never changes
ranking, category, ClinVar significance, summary or headline. Existing legacy
free-text severity logic is unchanged and is not validated by this feature.

### Verification and limitations

Run `python -m pytest skills/clinical-variant-prioritizer/tests/test_severity_evidence.py -q`.
Synthetic tests compare 216 tier-count combinations against the pinned
upstream function. They establish software behaviour only, not biological
correctness, clinical safety, population fairness or adoption.

Known baseline defect: main does not track data/clinical_panel.json, required
by the original default API. Use an explicit custom panel as above. Nine
existing default-panel tests fail independently of this extension.

ClawBio is a research and educational tool. It is not a medical device and
does not provide clinical diagnoses. Consult a healthcare professional before
making any medical decisions.
